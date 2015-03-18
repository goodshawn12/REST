function [I,C] = fcm_roiidx(roi)
% fcm_roiidx(roi)

global nuts fuse

if iscell(roi) || ischar(roi)  % if roi label given
    nuts.selvox.ipsilab = roi;
    return
end

if isfield(roi,'MEGvoxels')
    ivox = roi.MEGvoxels;
elseif isfield(roi,'VOIvoxels')
    ivox = roi.VOIvoxels;
elseif isfield(roi,'MRIvoxels')
    ivox = nut_coordtfm(roi.MRIvoxels,inv(nuts.coreg.meg2mri_tfm));
elseif isfield(roi,'voxels')
    ivox = nut_mni2mri(roi.voxels,nuts.coreg);
    ivox = nut_coordtfm(ivox,inv(nuts.coreg.meg2mri_tfm));
elseif isfield(roi,'MNIvoxels')
    ivox = nut_mni2mri(roi.MNIvoxels,nuts.coreg);
    ivox = nut_coordtfm(ivox,inv(nuts.coreg.meg2mri_tfm));    
else
    error('Invalid roi structure')
end

[idx,dist]=dsearchn(nuts.voxels,ivox);
if length(idx)>1
    idx=unique(idx(dist<=min(nuts.voxelsize)));
end
nuts.selvox.ipsi = idx;
if isfield(roi,'label')
    nuts.selvox.ipsilab = roi.label;
end
     
if nargout>0
    I=idx;
end
if ( isfield(fuse,'seed') && ~isempty(strfind(fuse.seed,'Contralateral')) )  % homologous contralateral voxel indices
    if isequal(nuts.coreg.meg2mri_tfm,eye(4))  % if voxels in MRI space
        cvox = [-ivox(:,1) ivox(:,2) ivox(:,3)];
    else                                       % if voxels in MEG space
        cvox = [ivox(:,1) -ivox(:,2) ivox(:,3)];
    end
    [idxc,dist]=dsearchn(nuts.voxels,cvox);
    idxc=unique(idxc(dist<=min(nuts.voxelsize)));
    nuts.selvox.contra = idxc;
    if nargout>1
        C=idxc;
    end
end