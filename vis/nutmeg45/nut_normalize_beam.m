function iserr=nut_normalize_beam(beamname,normedMRI,newvoxsize)
%function nut_normalize_beam(beamname,normedMRI,newvoxsize)
% beamname is saved s_beam_*.mat file
% normedMRI (optional input): path to spatially normalized MRI
% newvoxsize (optional input): voxel size after spatial normalization (e.g., [5 5 5])

ONLYS1=false;

global defaults
if isempty(defaults)
    spm('Defaults','fMRI');
end
if (strcmp(spm('ver'),'SPM2') && defaults.analyze.flip) || (strncmp(spm('ver'),'SPM8',4) && spm_flip_analyze_images)
    fprintf('WARNING IF USING ANALYZE-FORMAT MRIs (*.hdr/*.img):\nYour SPM default settings assume that MRIs have left and right sides flipped. This setting may be\nincompatible with NUTMEG!\n')
end
iserr=false;

% functional data
%----------------
beam=load(beamname);
beam=nut_beam_legacy_compatibility(beam);

if(isfield(beam,'orientation'))
    if beam.orientation==2
        errordlg('throw me a frickin'' bone here! you can''t spatially normalize if you aren''t using a neurological MRI'); iserr=true; return
    end
end
if nargin>2 && ~isempty(newvoxsize)
    beam.voxelsize = newvoxsize;
end

if isfield(beam,'rois')  % do not normalize voxel data when using ROIs, therefore replace here with dummy 0's
    beam.s={zeros(size(beam.s{1},1),1)};
end

% numor=length(beam.s);
[numvox,numtim,numfrq,numor]=size(beam.s{1});
if ~ONLYS1, [dum,nums]=size(beam.s); else nums=1; end
if isfield(beam,'rois')
    noiseflag=false; zflag=false; varflag=false;
else
    noiseflag= isfield(beam,'n');
    zflag= isfield(beam,'z');
    varflag= isfield(beam,'sasd');
end

[voxelsblob,blob2mri_tfm]=nut_voxels2blob(beam);  % same as rivets.voxelsblob, no need to have global rivets just for this one line
%voxelsblob_round=round(voxelsblob);  %this does seem necessary..once in awhile, something is not perfect integer

if( max(max(abs( round(voxelsblob) - voxelsblob ))) <= 1e-3 )
    dimxyz = max(round(voxelsblob));
else    % this happens when working with some EEG headmodels
    disp('nut_normalize_beam: the blob coordinates are not close to integer. Using a linear interpolation.')
    dimxyz = max(ceil(voxelsblob));
end

% structural data
%----------------
if nargin==1 || isempty(normedMRI)
    if isfield(beam.coreg,'norm_mripath')
        normedMRI=beam.coreg.norm_mripath;
    else
        [dum,mrifile,dum] = fileparts(beam.coreg.mripath);
        if ( mrifile(1)=='w' || strcmp(mrifile,'T1') || strncmp(mrifile,'avg',3) )  % if struct MRI is already normalized.     
            normedMRI=beam.coreg.mripath;                                           % This may happen in studies where no individual MRI is available.
        else
            [mrifile, norm_mripath]=uigetfile('w*.img','Please Select Normed MRI...');
            if isequal(mrifile,0) || isequal(norm_mripath,0)
                iserr=true;
                return
            end
            normedMRI = fullfile(norm_mripath,mrifile);
        end
    end
end
if ~exist(normedMRI,'file')
    errordlg('Normalized structural MRI not found.')
    iserr=true;
    return
end
if strcmp(normedMRI,beam.coreg.mripath)         % If struct MRI is already normalized, we do not need to transform, only interpolate voxel coordinates to common grid.
    normmats.Affine=eye(4);
    normmats.VF=spm_vol(beam.coreg.mripath);
    if strncmp(spm('ver'),'SPM8',4)
        normmats.VG=spm_vol(which('T1.nii'));
    elseif strcmp(spm('ver'),'SPM2')
        normmats.VG=spm_vol(which('T1.mnc'));
    else
        error('This version of SPM is not supported.')
    end
    normmats.Tr=[];
else
    [path,name]=fileparts(normedMRI);
    snmat=fullfile(path,[name(2:end) '_sn.mat']);
    normmats=load(snmat);
end

% settings
%---------
[bdir, bname]=fileparts(beamname);
SN.fname=fullfile(bdir, [bname '_tmp.img']);
% SN.mat=st.vols{1}.blobs{1}.mat;
SN.mat=blob2mri_tfm;  %so that we don't need global st for batch-running
SN.pinfo=[1;0;0];
SN.dim=dimxyz;
% datatype handled differently between SPM2 and later versions
% 64 means double precision
if(strcmp(spm('ver'),'SPM2'))
    SN.dim(4) = 64;
else
    SN.dt = [64 0];
end

flags.vox  = beam.voxelsize .* [-1 1 1]; % don't ask about that first -1. that's just the way it is. (some things will never change.)
flags.bb   = [min(-80,defaults.normalise.write.bb(1,1)) min(-115,defaults.normalise.write.bb(1,2)) min(-65,defaults.normalise.write.bb(1,3));
    max( 80,defaults.normalise.write.bb(2,1)) max(  80,defaults.normalise.write.bb(2,2)) max(100,defaults.normalise.write.bb(2,3))];
% the boundary box set by SPM tends to be too small, therefore we set here the minimal box size.

% loop through data dimensions
%-----------------------------
for tt=1:numtim  % tt num of time points
    for ff=1:numfrq  % ff num of frequency points
        for oo=1:numor
            %if isround
            for ss=1:nums
                grids{ss} = nut_vector2vol(beam.s{ss}(:,tt,ff,oo),voxelsblob,dimxyz);
            end
            if noiseflag
                gridn = nut_vector2vol(beam.n{1}(:,tt,ff,oo),voxelsblob,dimxyz);
            end
            if zflag
                gridz = nut_vector2vol(beam.z{1}(:,tt,ff,oo),voxelsblob,dimxyz);
            end
            if varflag
                gridsasd = nut_vector2vol(beam.sasd{1}(:,tt,ff,oo),voxelsblob,dimxyz);
                gridscsd = nut_vector2vol(beam.scsd{1}(:,tt,ff,oo),voxelsblob,dimxyz);
            end
%             else
%                 % Create an integer image grid, that matches the non-integer coords without leaving
%                 % blanks
%                 if(tt==1 && ff==1)
%                     gridlow=floor(voxelsblob);
%                     gridhigh=ceil(voxelsblob);
%                     gridmid=round(voxelsblob);
%                     gridall=union(gridlow,gridhigh,'rows');
%                     gridall=union(gridall,gridmid,'rows');
%                     gridall=gridall( find(all(gridall,2)) , : );        % Remove 0 coordinates to avoid error below
%                     clear gridlow gridhigh gridmid
%                 end
%                 
%                 % Nearest neighbor interpolation
%                 nn=dsearchn(gridall,voxelsblob);
%                 for ss=1:nums
%                     grids{ss}=nan(dimxyz);
%                     for vv=1:size(voxelsblob_round,1)
%                         grids{ss}(gridall(nn(vv),1),gridall(nn(vv),2),gridall(nn(vv),3))=beam.s{ss}(vv,tt,ff,oo);
%                     end
%                 end
%             end

            for ss=1:nums
                SNn{ss} = spm_write_vol(SN,grids{ss});
                SNout=spm_write_sn(SNn{ss},normmats,flags);
                Z=SNout.dat(:);
                if ss==1 && tt==1 && ff==1 && oo==1
                    keep = find(isfinite(Z));
                    [s_norm{1:nums}] = deal( zeros( length(keep),numtim,numfrq,numor) );
                    [crap,voxels]=spm_read_vols(SNout);
                    voxels = voxels';
                    beam.voxels = voxels(keep,:);
                    clear crap voxels
                end
                s_norm{ss}(:,tt,ff,oo) = Z(keep);
            end
            
            if noiseflag
                SN = spm_write_vol(SN,gridn);
                SNout=spm_write_sn(SN,normmats,flags);
                Z=SNout.dat(:);
                n_norm{1}(:,tt,ff,oo) = Z(keep);
            end
            if zflag
                SN = spm_write_vol(SN,gridz);
                SNout=spm_write_sn(SN,normmats,flags);
                Z=SNout.dat(:);
                z_norm{1}(:,tt,ff,oo) = Z(keep);
            end
            if varflag
                SN = spm_write_vol(SN,gridsasd);
                SNout=spm_write_sn(SN,normmats,flags);
                Z=SNout.dat(:);
                sasd_norm{1}(:,tt,ff,oo) = Z(keep);
                SN = spm_write_vol(SN,gridscsd);
                SNout=spm_write_sn(SN,normmats,flags);
                Z=SNout.dat(:);
                scsd_norm{1}(:,tt,ff,oo) = Z(keep);
            end
        end
    end
end

% delete temp files
%------------------
bimg=fullfile(bdir, [bname '_tmp.img']);
bhdr=fullfile(bdir, [bname '_tmp.hdr']);
bmat=fullfile(bdir, [bname '_tmp.mat']);
warning('off','MATLAB:DELETE:FileNotFound')     % the _tmp.mat file does not seem to be created in SPM8.
delete(bimg,bhdr,bmat);
warning('on','MATLAB:DELETE:FileNotFound')

% finish and save spat norm file
%-------------------------------
% beam.s=s_norm; clear s_norm
beam.s=s_norm;
if noiseflag
    beam.n=n_norm;
end
if zflag
    beam.z=z_norm;
end
if varflag
    beam.sasd=sasd_norm;
    beam.scsd=scsd_norm;
end

beam.coreg.mripath=normedMRI;
beam.coreg.meg2mri_tfm=eye(4);  %hack cuz we can't go back to meg coords anymore with an affine tfm

[bpath,bname,bext]=fileparts(beamname);
norm_sbeamname=fullfile(bpath,[bname '_spatnorm.mat']);
save(norm_sbeamname,'-struct','beam');
