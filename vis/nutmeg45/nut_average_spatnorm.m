function nut_average_spatnorm(stats,subjuse)
% stats are from nut_stat_settings*.mat

if ~exist('subjuse')
    subjuse=1; %which subject's normalised MRI to use for final display.
end

MNIlivevox=[];
for ii=1:length(stats.files)
    beam=load(stats.files{ii});
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

if ~isfield(beam,'sa')
    beam=nut_oldbeam2newbeam(beam);
end
if isfield(beam,'s')
    beam=rmfield(beam,'s');
end
sasum{1}=zeros(size(MNIlivevox,1),size(beam.sa{1},2),size(beam.sa{1},3));
scsum{1}=zeros(size(MNIlivevox,1),size(beam.sa{1},2),size(beam.sa{1},3));
nsum{1}=zeros(size(MNIlivevox,1),size(beam.sa{1},2),size(beam.sa{1},3));
zsum{1}=zeros(size(MNIlivevox,1),size(beam.sa{1},2),size(beam.sa{1},3));
sasdsum{1}=zeros(size(MNIlivevox,1),size(beam.sa{1},2),size(beam.sa{1},3));
scsdsum{1}=zeros(size(MNIlivevox,1),size(beam.sa{1},2),size(beam.sa{1},3));

clear beam
for ii=1:length(stats.files)
    tmp{ii}.beam=load(stats.files{ii});
    voxuse{ii}=ismember(tmp{ii}.beam.voxels,MNIlivevox,'rows');
    sasum{1}=sasum{1}+tmp{ii}.beam.sa{1}(voxuse{ii},:,:);
    scsum{1}=scsum{1}+tmp{ii}.beam.sc{1}(voxuse{ii},:,:);
    if isfield(tmp{ii}.beam,'n')
    nsum{1}=nsum{1}+tmp{ii}.beam.n{1}(voxuse{ii},:,:);
    end
    if isfield(tmp{ii}.beam,'sasd')
    sasdsum{1}=sasdsum{1}+tmp{ii}.beam.sasd{1}(voxuse{ii},:,:);
    scsdsum{1}=scsdsum{1}+tmp{ii}.beam.scsd{1}(voxuse{ii},:,:);
    zsum{1}=zsum{1}+tmp{ii}.beam.z{1}(voxuse{ii},:,:);
    end
end
beam=tmp{subjuse}.beam;
voxels=tmp{1}.beam.voxels(voxuse{1},:);
sa{1}=1e30*sasum{1}/length(stats.files);
sc{1}=1e30*scsum{1}/length(stats.files);
n{1}=1e30*nsum{1}/length(stats.files);
sasd{1}=1e30*sasdsum{1}/length(stats.files);
scsd{1}=1e30*scsdsum{1}/length(stats.files);
z{1}=zsum{1}/length(stats.files);

% compute standard deviation
vsasum{1}=zeros(size(MNIlivevox,1),size(beam.sa{1},2),size(beam.sa{1},3));
vscsum{1}=zeros(size(MNIlivevox,1),size(beam.sa{1},2),size(beam.sa{1},3));
vnsum{1}=zeros(size(MNIlivevox,1),size(beam.sa{1},2),size(beam.sa{1},3));
vsasdsum{1}=zeros(size(MNIlivevox,1),size(beam.sa{1},2),size(beam.sa{1},3));
vscsdsum{1}=zeros(size(MNIlivevox,1),size(beam.sa{1},2),size(beam.sa{1},3));
vzsum{1}=zeros(size(MNIlivevox,1),size(beam.sa{1},2),size(beam.sa{1},3));
for ii=1:length(stats.files)
%     tmp{ii}=load(stats.files{ii});
%     voxuse{ii}=ismember(tmp{ii}.beam.voxels,MNIlivevox,'rows');
    vsasum{1}=vsasum{1}+(1e30*tmp{ii}.beam.sa{1}(voxuse{ii},:,:)-sa{1}).^2;
    vscsum{1}=vscsum{1}+(1e30*tmp{ii}.beam.sc{1}(voxuse{ii},:,:)-sc{1}).^2;
    if isfield(tmp{ii}.beam,'n')
    vnsum{1}=vnsum{1}+(1e30*tmp{ii}.beam.n{1}(voxuse{ii},:,:)-n{1}).^2;
    vsasdsum{1}=vsasdsum{1}+(1e30*tmp{ii}.beam.sasd{1}(voxuse{ii},:,:)-sasd{1}).^2;
    vscsdsum{1}=vscsdsum{1}+(1e30*tmp{ii}.beam.scsd{1}(voxuse{ii},:,:)-scsd{1}).^2;
    vzsum{1}=vzsum{1}+(tmp{ii}.beam.z{1}(voxuse{ii},:,:)-z{1}).^2;
    end
end
vsa{1}=vsasum{1}/length(stats.files);
vsc{1}=vscsum{1}/length(stats.files);
vn{1}=vnsum{1}/length(stats.files);
vsasd{1}=vsasdsum{1}/length(stats.files);
vscsd{1}=vscsdsum{1}/length(stats.files);
vz{1}=vzsum{1}/length(stats.files);


beam=rmfield(beam,'sa');
beam=rmfield(beam,'sc');
if isfield(beam,'n')
beam=rmfield(beam,'n');
beam=rmfield(beam,'sasd');
beam=rmfield(beam,'scsd');
beam=rmfield(beam,'z');
end

for jj=1:size(sa{1},3),
    beam.sa{1}(:,:,jj)=sa{1}(find(sa{1}(:,:,1)),:,jj);
    beam.sc{1}(:,:,jj)=sc{1}(find(sc{1}(:,:,1)),:,jj);
    beam.n{1}(:,:,jj)=n{1}(find(n{1}(:,:,1)),:,jj);
    beam.vsa{1}(:,:,jj)=vsa{1}(find(vsa{1}(:,:,1)),:,jj);
    beam.vsc{1}(:,:,jj)=vsc{1}(find(vsc{1}(:,:,1)),:,jj);
    beam.vn{1}(:,:,jj)=vn{1}(find(vn{1}(:,:,1)),:,jj);
    beam.sasd{1}(:,:,jj)=sasd{1}(find(sasd{1}(:,:,1)),:,jj);
    beam.scsd{1}(:,:,jj)=scsd{1}(find(scsd{1}(:,:,1)),:,jj);
    beam.z{1}(:,:,jj)=z{1}(find(z{1}(:,:,1)),:,jj);
    beam.vsasd{1}(:,:,jj)=vsasd{1}(find(vsasd{1}(:,:,1)),:,jj);
    beam.vscsd{1}(:,:,jj)=vscsd{1}(find(vscsd{1}(:,:,1)),:,jj);
    beam.vz{1}(:,:,jj)=vz{1}(find(vz{1}(:,:,1)),:,jj);
    beam.voxels=voxels(find(sa{1}(:,:,1)),:);
end

[pathstr,filnam,fext]=fileparts(stats.files{1});
% save([pathstr(1:(end-2)) filnam(1:8) filnam(12:end) '_avg' fext],'-v7.3','-struct','beam');
save([filnam(1:8) filnam(12:end) '_avg' fext],'-struct','beam');

% ttest
if ~isempty(beam.sasd{1})
beam.sa{1}=(beam.sa{1}-beam.sc{1})./(((beam.sasd{1}+beam.scsd{1})/2)/sqrt(length(stats.files)));
beam.sc{1}=[];beam.n{1}=[];beam.sav{1}=[];beam.scv{1}=[];beam.nv{1}=[];
save([filnam(1:8) filnam(12:end) '_avg_ttest' fext],'-struct','beam');
end

% z
if isfield(beam,'z')
beam.sa=beam.z;
save([filnam(1:8) filnam(12:end) '_avg_z' fext],'-struct','beam');
end