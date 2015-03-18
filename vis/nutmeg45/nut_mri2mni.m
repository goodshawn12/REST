function [xyz_o]=nut_mri2mni(xyz_i,coreg,doaffine)
% [xyz_o]=nut_mri2mni(xyz_i,coreg,doaffine)
% Takes MRI coords (mm) and converts them to MNI coordinates (mm)
%
% XYZ_I is Nx3 list of coords
%
% coreg:  optional input (from either beam.coreg or nuts.coreg)
%         else we will find it from global beam or global nuts
%         Either way, it must contain coreg.norm_mripath and coreg.mripath
%
% doaffine: (optional) applies affine transform for area outside of
%           bounding box specified by SPM warping. Ok to =1 during
%           nut_results_viewer, but should be =0 when calling from
%           nut_makeMNIvoi

if nargin<2 || isempty(coreg)
    global beam
    if(isempty(beam))
        global nuts
        if(isfield(nuts,'coreg'))
            coreg = nuts.coreg;
        else
            error('no coregistration info found in global beam or nuts');
        end
    else
        coreg = beam.coreg;
    end
end
    
if isstruct(coreg)
    if ~isfield(coreg,'norm_mripath')
        [dum,mrifile,dum] = fileparts(coreg.mripath);
        if ( mrifile(1)=='w' || strcmp(mrifile,'T1') || strncmp(mrifile,'avg',3) )  % if struct MRI is already normalized.     
            norm_mripath=coreg.mripath;       
        else
            error('must have a normed MRI loaded');
        end
    else
        norm_mripath=coreg.norm_mripath;
    end
    mripath=coreg.mripath;
else
    error('second input, if specified, must be either blank or coreg structure')
end


if size(xyz_i,2)~=3, error('coord must be an N x 3 matrix'); end;

if strcmp(mripath,norm_mripath)
    % This may happen if the displayed data is already spatially
    % normalized.
    xyz_o = xyz_i;
    return
end

[path,name,ext]=fileparts(norm_mripath);
Vorig=spm_vol(mripath);
norminput=fullfile(path,[name(2:end) ext]);
Vnorminput=spm_vol(norminput);
xyz=nut_coordtfm(xyz_i,Vnorminput.mat*inv(Vorig.mat));

if(strcmp(spm('ver'),'SPM8b') || strcmp(spm('ver'),'SPM8'))
    %     iyimg=fullfile(path,['y_i' name(2:end) ext]);
    iyimg=fullfile(path,['y_i' name(2:end) '.nii']); % ext might be '.img', but y_i* file will always be *.nii if dealing with SPM8
    if(exist(iyimg,'file'))
        iyimg = [repmat(iyimg,3,1) [',1,1';',1,2';',1,3']];
    else
        iyimg=fullfile(path,['iy_' name(2:end) ext]);
        iyimg = [repmat(iyimg,3,1) [',1';',2';',3']];
    end
elseif(strcmp(spm('ver'),'SPM2'))
    iyimg=fullfile(path,['iy_' name(2:end) ext]);
    iyimg = [repmat(iyimg,3,1) [',1';',2';',3']];
else
    errordlg(['NUTMEG is not compatible with ' spm('ver') '. Please place SPM2 or SPM8 in your path, which may be downloaded from http://www.fil.ion.ucl.ac.uk/spm'])
    return
end

% this is thanks to john's gem #5  http://www.sph.umich.edu/~nichols/JohnsGems2.html
Viy = spm_vol(iyimg);
vx = double(nut_coordtfm(xyz,inv(Viy(1).mat)));  % The voxel in the deformation to sample
xyz_o = [ spm_sample_vol(Viy(1),vx(:,1),vx(:,2),vx(:,3),1) spm_sample_vol(Viy(2),vx(:,1),vx(:,2),vx(:,3),1) spm_sample_vol(Viy(3),vx(:,1),vx(:,2),vx(:,3),1)];

if ~exist('doaffine') | isempty(doaffine)
    doaffine=0;
end

if doaffine==1
    for ii=1:size(xyz_o,1)
        if (sum(xyz_o(ii,:)==0)==3 || any(isnan(xyz_o(ii,:))))
            % this will happen if voxel outside of region of Viy, and spm_sample_vol returns 'zeros' for every element
            % then we do the next best thing of applying Affine normalisation transform, rather than warped transform
            load([path '/' name(2:end) '_sn.mat']);
            vxorig=nut_coordtfm(xyz_i(ii,:),inv(Vorig.mat));
            xyz_o(ii,:)=nut_coordtfm(nut_coordtfm(vxorig,inv(Affine)),VG.mat);
            disp(['warning: the ' num2str(ii) ' coordinate was out of bounds and so transformed using Affine only']);
        end
    end
end

