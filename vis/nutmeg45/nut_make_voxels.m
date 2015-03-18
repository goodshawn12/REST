function voxels=nut_make_voxels(voxelsize);

global nuts ndefaults

voxelsize = [voxelsize voxelsize voxelsize];

if( isfield(nuts,'voxels') && all(nuts.voxelsize==voxelsize) )
       disp('Using existing nuts.voxels and nuts.voxelsize...');
     voxels = nuts.voxels;
    voxelsize = nuts.voxelsize;
else
    if(isfield(nuts,'VOIvoxels'))
       disp('Using manually selected VOI (nuts.VOIvoxels)...');
         VOIvoxels = nuts.VOIvoxels;

        switch(ndefaults.lf.voxelMRIalign)
            case 1
                VOIvoxels = nut_meg2mri(VOIvoxels);
            case 2
                warning('MNI voxel alignment not supported for manual VOI selection...')
        end

        % fit to desired voxel grid
        for ii = 1:3
            voxels(:,ii) = voxelsize(ii)*round(VOIvoxels(:,ii)/voxelsize(ii));
        end
        clear VOIvoxels;

        % remove duplicate voxels from fit
        voxels = unique(voxels,'rows');
        switch(ndefaults.lf.voxelMRIalign)
            case 1
                voxels = nut_mri2meg(voxels);
            case 2
                
        end
        nuts.voxelsize = voxelsize;
    elseif(isfield(nuts.coreg,'norm_mripath'))
        disp('Generating VOI from MNI template...');
        tic;
        switch(ndefaults.lf.voxelMRIalign)
            case 0
                voxels = nut_makeMNIvoi(voxelsize,'wholehead');
            case 1
                warning('MRI voxel alignment not yet supported for MNI-based VOI selection.')
            case 2
                load MNIvoxels
                MNIvoxels = double(MNIvoxels);
                voxelsMRI = nut_mni2mri(unique(round(MNIvoxels/voxelsize(1))*voxelsize(1),'rows'));
                voxels = nut_mri2meg(voxelsMRI);
        end
        disp('done generating VOI.'); toc
        nuts.voxelsize = voxelsize;
    else
        errordlg('You must load a spatially normalized MRI or manually define a VOI to generate a lead field.')
    end
end
