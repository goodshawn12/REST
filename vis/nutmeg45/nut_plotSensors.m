function nut_plotSensors
% NUT_PLOTSENSORS
% Plots MEG sensors on SPM volume.
% Often used as a sanity check for meg2mri transformation.

global nuts st;

% sensor coordinates remain in MEG coords, so we should convert to MRI
% coords for plotting
coilcoord = nut_meg2mri(nuts.meg.sensorCoord(:,:,1));

% following jazz is needed because SPM's blob freaks out when given
% negative coordinates, and needs a dilation to know these are big voxels

% (chose 4 mm^3 voxel size)
voxelsize = [4 4 4];
translation_tfm = [ voxelsize(1)            0                 0 min(coilcoord(:,1,1))-voxelsize(1)
                               0 voxelsize(2)                 0 min(coilcoord(:,2,1))-voxelsize(2)
                               0                 0 voxelsize(3) min(coilcoord(:,3,1))-voxelsize(3)
                               0                 0            0                                             1 ];

coil_blobs = nut_coordtfm(coilcoord(:,:,1),inv(translation_tfm));

% keep=find(prod(double(coil_coord > 0.5),2)); % discard coords with nonpositive voxels
% coil_blobs  = coil_coord(keep,:);
spm_orthviews('rmblobs',1);
spm_orthviews('addblobs',1,coil_blobs',zeros(size(coil_blobs,1),1),translation_tfm);
if strcmp(spm('ver'),'SPM2')
    delete(st.vols{1}.blobs{1}.cbar);
end
spm_orthviews('redraw',1);
