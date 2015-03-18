function nut_reset
% NUT_RESET
% Reset NUTMEG

clear global nuts beam bolts;  % flush any previous nutmeg leftovers

global defaults

global nuts screws st defaults
% nuts will contain all data loaded for analysis
% screws will contain any GUI incidentals
% st is from SPM and contains MRI info
% defaults is from SPM

spm('Defaults','fMRI');

if(strcmp(spm('ver'),'SPM2'))
    if defaults.analyze.flip
        warndlg('WARNING: Your SPM default settings assume that MRIs have left and right sides flipped. This setting may be incompatible with NUTMEG! Please type "edit spm_defaults" in the MATLAB command window, and change the line "defaults.analyze.flip=1" to "defaults.analyze.flip=0".')
    end
    defaults.analyze.flip=0; % needed to placate SPM2's trigger-happy flipping
elseif(strcmp(spm('ver'),'SPM8') || strcmp(spm('ver'),'SPM8b'))
    if spm_flip_analyze_images
        fprintf('WARNING IF USING ANALYZE-FORMAT MRIs (*.hdr/*.img):\nYour SPM default settings assume that MRIs have left and right sides flipped. This setting may be\nincompatible with NUTMEG! Please type "edit spm_flip_analyze_images.m" in the MATLAB command window,\nand change the line "flip=1" to "flip=0".\n')
    end
    defaults.analyze.flip=0; % needed to placate SPM2's trigger-happy flipping
else
    warning('not sure what to do about defaults.analyze.flip for this version of SPM');
end


nuts.coreg.mripath = which('blank.img');
% spm_image('init',nuts.coreg.mripath);
nuts.coreg.meg2mri_tfm = eye(4);

% if(isempty(gcbf))  % this happens when opening nutmeg
%     nuts.fig = gcf;
%     handles = guihandles(nuts.fig);
% else
%     nuts.fig = gcbf;
%     handles = guihandles(nuts.fig);
% end
nuts.fig = findobj('tag','nutmegfig');
handles = guihandles(nuts.fig);

set(handles.nut_megfile,'String','(none loaded)');
nut_enabler(handles);

nut_refresh_image;

