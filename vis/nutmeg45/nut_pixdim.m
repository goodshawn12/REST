function [pixdim] = nut_readdim(vol)
% pixdim = nut_readdim(vol)
% returns voxelsize from either SPM2 or SPM8/SPM8b

if(strcmp(spm('ver'),'SPM2'))
%    ndim=vol.private.hdr.dime.dim(2:4);
    pixdim=vol.private.hdr.dime.pixdim(2:4);
elseif(strcmp(spm('ver'),'SPM8') || strcmp(spm('ver'),'SPM8b'))
%    ndim=vol.private.hdr.dim(2:4);
    pixdim=vol.private.hdr.pixdim(2:4);
end
