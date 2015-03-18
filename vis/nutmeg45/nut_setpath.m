function nut_setpath(dosave)
% NUT_SETPATH   sets correct NUTMEG path.
% Assumes that base NUTMEG directory is already in path or current working directory.
%
% Usage:
%   nut_setpath(dosave)
%
% dosave      (optional) save new path to Matlab's pathdef.m file.
%             Default is true.
% 
% NOTE: This does not set the SPM path which is also required for NUTMEG!


if nargin<1, dosave=true; end

NUT_BASE = fileparts(which('nutmeg.m'));

% ndirname = fliplr(strtok(fliplr(NUT_BASE),filesep));

% Remove old (and possibly wrong) entries, if existing
P = path;
P = textscan(P,'%s','delimiter',pathsep);
P = P{1};
idx = strmatch(NUT_BASE,P);
if ~isempty(idx)
    rmpath(P{idx});
end

% Add correct path
newpath = { 
     NUT_BASE
    [NUT_BASE filesep 'beamformers'];
    [NUT_BASE filesep 'data_importers'];
    [NUT_BASE filesep 'denoisers'];
    [NUT_BASE filesep 'examples'];
    [NUT_BASE filesep 'external'];
    [NUT_BASE filesep 'external' filesep 'ctf'];
    [NUT_BASE filesep 'external' filesep 'openmeeg'];
    [NUT_BASE filesep 'fcm'];
    [NUT_BASE filesep 'fcm' filesep 'params'];
    [NUT_BASE filesep 'fcm' filesep 'src'];
    [NUT_BASE filesep 'leadfield_obtainers'];
    [NUT_BASE filesep 'mesh'];
    [NUT_BASE filesep 'nuteeg'];
    [NUT_BASE filesep 'nuteeg' filesep 'util'];
    [NUT_BASE filesep 'nuteeg' filesep 'util' filesep 'Sphere_Tessellation'];
    [NUT_BASE filesep 'openmeeglink'];
    [NUT_BASE filesep 'templates'];
    [NUT_BASE filesep 'tfbf'];
    [NUT_BASE filesep 'tfbf' filesep 'params'];
    [NUT_BASE filesep 'tfbf' filesep 'src'];
    };
addpath(newpath{:});

if exist('filt')
    warning('off');rmpath(fileparts(which('filt')));warning('on');
    warning('we have removed control toolbox from your path, as there is a variable - file name conflict');
end
% if strfind(P,[NUT_BASE '/deprecated'])
%     rmpath([NUT_BASE '/deprecated'])
% end

% save path
if dosave
    savepath;
end

