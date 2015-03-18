function nut_create_nuts(datafile,varargin)
% datafile: full path name of where data file is (*.meg4, *.fif, etc)
% mrifile: full path name of where individual subject MRI
% lsc_file: full path name of 

nutmeg;
global nuts;

mrifiletmp=keyval('mrifile',varargin);
lscfiletmp=keyval('lscfile',varargin);
orientationtmp=keyval('orient',varargin);
normmritmp=keyval('mnifile',varargin);
voxelstmp=keyval('voxelsfile',varargin);
lfcomptmp=keyval('lfcompmode',varargin);
voxelsizetmp=keyval('voxsize',varargin);
sessfiletmp=keyval('sessionfile',varargin);

if ~isempty(mrifiletmp)
    nuts.coreg.mripath = mrifiletmp;
    nut_refresh_image
end

if ~isempty(orientationtmp)
    nuts.coreg.orientation=orientationtmp;
else
    nuts.coreg.orientation=1;
end

global coreg
coreg=nuts.coreg;
if isempty(lscfiletmp)
    nut_importmeg(datafile);
else
    nut_importmeg(datafile,lsc_file);
end

if ~isempty(normmritmp)
    nuts.coreg.norm_mripath=normmritmp;
elseif ~isempty(voxelstmp)
    load(voxelstmp);
    nuts.VOIvoxels=VOIvoxels;
end

if isempty(voxelsizetmp)
    voxelsize=5;
else
    voxelsize=voxelsizetmp;
end

if ~isempty(lfcomptmp)
    lfcomp=lfcompttmp;
    [nuts,nuts.Lp,nuts.voxels]=nut_obtain_lead_field(nuts,voxelsize,lfcomp);
else
    [nuts,nuts.Lp,nuts.voxels]=nut_obtain_lead_field(nuts,voxelsize);
end


if ~isempty(sessfiletmp)
    save(sessfiletmp,'nuts');
end


    




