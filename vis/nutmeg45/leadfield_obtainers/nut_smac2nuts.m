function nuts = nut_smac2nuts(voxelsize)
% Converts SMAC solution points (=voxels) and leadfield to NUTMEG format
%   ¦nuts¦ = nut_smac2nuts(¦voxelsize¦)
%      (parameters in ¦¦ are optional)
%
% voxelsize     in mm, e.g. 5. Corresponds to the "sampling rate" for
%               solution points in SMAC if you resampled your MRI to a
%               voxel size of 1x1x1 mm.

global nuts

if nargin<1   
    global temp
    if ( ~isfield(nuts,'fig') || ~ishandle(nuts.fig) )
        error('No voxelsize specified.') 
    end
    vsfig = figure('units','normalized','position',[0.45 0.45 0.1 0.1], ...
        'toolbar','none','menubar','none','numbertitle','off');
    uicontrol(vsfig,'Style','text','string','Enter voxel size in mm:', ...
        'units','normalized','position',[.1 .85 .8 .1]);
    uicontrol(vsfig,'Style','edit','units','normalized','tag','vsedit', ...
        'position',[.3 .55 .4 .2]);
    uicontrol(vsfig,'Style','pushbutton','string','Ok', ...
        'callback','global temp; handles=guihandles(gcf); temp=get(handles.vsedit,''string''); delete(gcf);', ...
        'units','normalized','position',[.3 .25 .4 .2]);
    uiwait(vsfig);
    voxelsize = str2num(temp);
    if isempty(voxelsize), return, end
    clear global temp
end
if length(voxelsize)==1, voxelsize=voxelsize*ones(1,3); end

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

% Get all the filenames
d=dir([fullfile(pIMG,fIMG) '*.spi']);
if length(d)~=1
    [fSPI, pSPI] = uigetfile('*.spi','Select SolutionPoints file created by SMAC');
    if isequal(fSPI,0), return, end
else
    fSPI=d(1).name;
    pSPI=pIMG;
end

d=dir([fullfile(pIMG,fIMG) '*.lf']);
if length(d)~=1
    [fLF, pLF] = uigetfile('*.lf','Select leadfield file created by SMAC');
    if isequal(fLF,0), return, end
else
    fLF=d(1).name;
    pLF = pIMG;
end

% Import voxels and leadfield
if strcmp(ext,'.img'), ext='.hdr'; end
nuts.voxels     = nut_import_smacsolutionpoints(fullfile(pSPI,fSPI),fullfile(pIMG,[fIMG ext]));
nuts.voxelsize  = voxelsize;
nuts.Lp         = nut_import_smacleadfield(fullfile(pLF,fLF));

if ~isempty(findobj('tag','nutmegfig'))
    nut_refresh_image;
    nut_enabler;
end
