function nuts = nut_cartoollf2nuts
% Converts Cartool solution points (=voxels) and leadfield to NUTMEG format

mult=2; % Newer versions of Cartool produce voxel coordinates in downsampled space, so we have to correct this.

global nuts

% Coregistered structural MRI
if ( ~isfield(nuts,'coreg') || ~isfield(nuts.coreg,'mripath') || strcmp(nuts.coreg.mripath,which('blank.img')) )
    [fIMG, pIMG] = uigetfile('*.img;*.nii','Select structural MRI');
    if isequal(fIMG,0), return, end
    [dum,fIMG,ext]=fileparts(fIMG);
    if strcmpi(ext,'.hdr'), ext='.img'; end
    nuts.coreg.mripath=fullfile(pIMG,[fIMG ext]);
    if exist(fullfile(pIMG,['w' fIMG ext]),'file')
        nuts.coreg.norm_mripath = fullfile(pIMG,['w' fIMG ext]);
    end
    nuts.coreg.meg2mri_tfm=eye(4);
    nuts.coreg.orientation=1;
else
    [pIMG,fIMG,ext]=fileparts(nuts.coreg.mripath);
    if strcmpi(ext,'.hdr'), ext='.img'; end
    if ( ~isfield(nuts.coreg,'norm_mripath') && exist(fullfile(pIMG,['w' fIMG ext]),'file') )
        nuts.coreg.norm_mripath = fullfile(pIMG,['w' fIMG ext]);
    end
    nuts.coreg.meg2mri_tfm=eye(4);              % overwrite previous settings!
    if isfield(nuts.coreg,'fiducials_mri_mm')   % no fiducials may be specified!
        nuts.coreg = rmfield(nuts.coreg,'fiducials_mri_mm');
    end
end
% answer = questdlg('Is your MRI a template or the MRI of the subject?', 'NUTMEG question', 'Template', 'Individual', 'Individual');
% domult2 = strcmp(answer,'Individual')+1;

mult=2; % This seems to be required in newer versions of Cartool, but not in older...?

% Get all the filenames
d=dir([pIMG filesep '*.spi']);
if length(d)~=1
    [fSPI, pSPI] = uigetfile([pIMG filesep '*.spi'],'Select SolutionPoints file created by Cartool');
    if isequal(fSPI,0), return, end
else
    fSPI=d(1).name;
    pSPI=pIMG;
end

d=dir([pIMG filesep '*.ris']);
if length(d)~=1
    [fLF, pLF] = uigetfile([pIMG filesep '*.ris'],'Select leadfield file created by Cartool');
    if isequal(fLF,0), return, end
else
    fLF=d(1).name;
    pLF = pIMG;
end


% Import voxels and leadfield
[x y z]         = textread(fullfile(pSPI,fSPI),'%f %f %f %*s');
voxelsMRI       = [x y z]; clear x y z
voxelsizeMRI    = [mean(diff(unique(voxelsMRI(:,1)))) mean(diff(unique(voxelsMRI(:,2)))) mean(diff(unique(voxelsMRI(:,3))))];
nuts.voxels     = voxelsMRI .* mult;
nuts.voxelsize  = voxelsizeMRI .* mult;
nuts.Lp         = nut_import_cartoolris(fullfile(pLF,fLF));

% Make integer
% nuts.voxelsize  = round(voxelsizeMRI);
% vsr = voxelsizeMRI ./ nuts.voxelsize;
% tfm = diag([vsr 1]);
% nuts.voxels = nut_coordtfm(voxelsMRI,inv(tfm));
% nuts.coreg.meg2mri_tfm=tfm;

% Update NUTMEG GUI
if ~isempty(findobj('tag','nutmegfig'))
    nut_refresh_image;
    nut_enabler;
end

