function [voxelsblob,blob2mri_tfm]=nut_voxels2blob(beam)
% [voxelsblob,blob2mri_tfm]=nut_voxels2blob(beam)


beam.voxels = double(beam.voxels);
voxelsMRI = double(nut_coordtfm(beam.voxels,beam.coreg.meg2mri_tfm));

% if(round(voxelsMRI)-voxelsMRI<1e-12)
if (round(voxelsMRI)-voxelsMRI<1e-4)        % if voxels in MRI space
    blob2mri_tfm = [ beam.voxelsize(1)        0                 0             min(voxelsMRI(:,1))-beam.voxelsize(1)
                            0          beam.voxelsize(2)        0             min(voxelsMRI(:,2))-beam.voxelsize(2)
                            0                 0         beam.voxelsize(3)     min(voxelsMRI(:,3))-beam.voxelsize(3)
                            0                 0                 0                                1 ];
    voxelsblob = nut_coordtfm(voxelsMRI,inv(blob2mri_tfm));
else  %if voxels in MEG space
    blob2meg_tfm = [ beam.voxelsize(1)        0                 0           min(beam.voxels(:,1))-beam.voxelsize(1)
                            0           beam.voxelsize(2)       0           min(beam.voxels(:,2))-beam.voxelsize(2)
                            0                 0         beam.voxelsize(3)   min(beam.voxels(:,3))-beam.voxelsize(3)
                            0                 0                 0                                1 ];
    voxelsblob = nut_coordtfm(beam.voxels,inv(blob2meg_tfm));
    blob2mri_tfm=beam.coreg.meg2mri_tfm*blob2meg_tfm;
% else
%     warning('all your voxels are not in either MRI or MEG space??')
%     return
end

blob2mri_tfm=double(blob2mri_tfm);

