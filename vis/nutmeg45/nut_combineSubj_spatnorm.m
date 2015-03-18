function nut_combineSubj_spatnorm(files,subjuse,do_aal_crop,outfilename,aalpath)
% a spin-off of nut_average_spatnorm, but saves mean, stdev, and each
% subject value for further stats
% thanks to Mich for aal speedup

if ~exist('subjuse')
    subjuse=1; %which subject's normalised MRI to use for final display.
end
if ~exist('do_aal_crop')
    do_aal_crop=0; %crop out voxels outside of AAL labels
end

MNIlivevox=[];
for ii=1:length(files)
    beam=load(files(ii).name);
    if isfield(beam,'beam')
        beam=beam.beam;
    end
    length(beam.voxels)
    if ii==1
        MNIlivevox=beam.voxels;
    else
        MNIlivevox=intersect(MNIlivevox,beam.voxels,'rows');
    end
end

beam=nut_beam_legacy_compatibility(beam);
[~,nums]=size(beam.s);
for ss=1:nums
    ssum{ss}=zeros(size(MNIlivevox,1),size(beam.s{ss},2),size(beam.s{ss},3));
    for ii=1:length(files)
        seach{ii,ss}=zeros(size(MNIlivevox,1),size(beam.s{ss},2),size(beam.s{ss},3));
    end
end

clear beam
for ii=1:length(files)
    tmp{ii}.beam=load(files(ii).name);
    voxuse{ii}=ismember(tmp{ii}.beam.voxels,MNIlivevox,'rows');
    
    for ss=1:nums
        ssum{ss}=ssum{ss}+tmp{ii}.beam.s{ss}(voxuse{ii},:,:);
        seach{ii,ss}=tmp{ii}.beam.s{ss}(voxuse{ii},:,:);
    end
end
beam=tmp{subjuse}.beam;
voxels=tmp{subjuse}.beam.voxels(voxuse{subjuse},:);
for ss=1:nums
    s{ss}=ssum{ss}/length(files);
end

% compute standard deviation
for ss=1:nums
    vssum{ss}=zeros(size(MNIlivevox,1),size(beam.s{ss},2),size(beam.s{ss},3));
    for ii=1:length(files)
        vssum{ss}=vssum{ss}+(tmp{ii}.beam.s{ss}(voxuse{ii},:,:)-s{ss}).^2;
    end
    vs{ss}=vssum{ss}/length(files);
end

beam=rmfield(beam,'s');

for jj=1:size(s{1},3),
    for ss=1:nums
        beam.s{ss}(:,:,jj)=s{ss}(find(s{ss}(:,1,1)),:,jj);
        beam.vs{ss}(:,:,jj)=vs{ss}(find(vs{ss}(:,1,1)),:,jj);
        for ii=1:length(files)
            beam.seach{ii,ss}(:,:,jj)=seach{ii,ss}(find(s{ss}(:,1,1)),:,jj);
        end
    end
    beam.voxels=voxels(find(s{1}(:,1,1)),:);
end

[pathstr,filnam,fext]=fileparts(files(subjuse).name);
% save([pathstr(1:(end-2)) filnam(1:8) filnam(12:end) '_avg' fext],'-v7.3','-struct','beam');
save(outfilename,'-struct','beam');

if do_aal_crop
    aal=load([aalpath '/aal_info.mat']);
    aalxyzflip = aal.xyz';       %saves tons and tons of flipping operations
    voxaal = zeros(3,size(beam.voxels,1)); %pre-allocates memory
    for ii=1:size(beam.voxels,1)
        mni8_2aal=dsearchn(aalxyzflip,beam.voxels(ii,:));
        voxaal(:,ii)=nut_coordtfm(aal.xyz(:,mni8_2aal)',aal.mm2vox);
%         if mod(ii,100) == 0
%             disp (['AALSTEP 1 vox: ' num2str(ii) ' of ' num2str(size(beam.voxels,1))]);
%         end;
    end
    clear aalxyzflip; %clears memory
    
    for ii=1:size(beam.voxels,1)
%         disp (['AALSTEP 2 vox: ' num2str(ii) ' of ' num2str(size(beam.voxels,1))]);
        aalreg(ii)=aal.y(voxaal(1,ii),voxaal(2,ii),voxaal(3,ii));
        if aalreg(ii)==0
            for ss=1:nums
            beam.s{ss}(ii,:,:)=0;
            end
        end
    end
    clear s voxnew
    for ss=1:nums
    s{ss}=beam.s{ss}(find(beam.s{ss}(:,1,1)),:,:);
    vs{ss}=beam.vs{ss}(find(beam.s{ss}(:,1,1)),:,:);
    end
    voxnew=beam.voxels(find(beam.s{1}(:,1,1)),:);
    clear seach
    for ii=1:length(files)
        for ss=1:nums
        seach{ii,ss}=beam.seach{ii,ss}(find(beam.s{ss}(:,1,1)),:,:);
%         disp (['AALSTEP 3 ii: ' num2str(ii) ' of ' num2str(length(files))]);
        end
    end
    beam=rmfield(beam,'s');
    beam=rmfield(beam,'vs');
    beam=rmfield(beam,'seach');
    beam=rmfield(beam,'voxels');
    beam.s=s;
    beam.vs=vs;
    beam.voxels=voxnew;
    beam.seach=seach;
    
    beam.coreg.brainrender_path='/net/cador/data_local/Anatomical_Library/JZumer/render_c1wzumer_P3T_N.mat';
    save(outfilename,'-struct','beam');
end





