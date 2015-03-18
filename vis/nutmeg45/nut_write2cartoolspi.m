function vox_trans = nut_write2cartoolspi(spifilename,beam,mriconvention)
% nut_write2cartoolspi  exports MEG voxel coordinates to Cartool spi files.
%
% vox_trans = nut_write2cartoolspi(spifilename,beam,mriconvention)
%
%   spifilename     name of file to be written
%   beam            NUTMEG results structure
%   mriconvention   'spm' (default): Structural MRIs ara in RAS convention
%                   used in SPM and NUTMEG.
%                   'smac': MRIs are in SMAC convention.
%                   'origin': MRIs are saved with same origin as calculated
%                   by SPM and have 1 mm voxel size.

if nargin<3, mriconvention='SMAC'; end

[pIMG,fIMG,eIMG]=fileparts(beam.coreg.mripath);
if( fIMG(1)=='w' && exist(fullfile(pIMG,[fIMG(2:end) eIMG]),'file') )    % Bring to non-normalized if possible
    beam.coreg.norm_mripath=beam.coreg.mripath;
    beam.coreg.mripath = fullfile(pIMG,[fIMG(2:end) eIMG]);
    voxels = nut_mni2mri(beam.voxels,beam.coreg);
else                % Convert from MEG to MRI space
    voxels = nut_coordtfm(beam.voxels,beam.coreg.meg2mri_tfm);
end

if strcmpi(mriconvention,'origin')
    vox_trans = voxels;
        
else
    V = spm_vol(beam.coreg.mripath);

    % Bring mm to voxel pixels
    % This corresponds to: voxels = nut_mm2voxels(voxels);
    voxels = [voxels ones(size(voxels,1),1)] * inv(V.mat)';
    voxels(:,4)=[];   % destroy last column of ones

    % Rotate and transform functional data to SMAC space if necessary
    if strcmpi(mriconvention,'smac')
        vox_trans=zeros(size(voxels));
        vox_trans(:,3) = voxels(:,1);
        vox_trans(:,1) = V.dim(2) - voxels(:,2);
        vox_trans(:,2) = V.dim(3) - voxels(:,3);
    else
        vox_trans=voxels;
    end
end

% Write
dlmwrite(spifilename,vox_trans,' ');

