function nut_refresh_image
% NUT_REFRESH_IMAGE

global nuts;
if(isempty(nuts))
    msgbox('Go to Load Image','Image was never loaded in workspace','warn')
    return;
end       
if(isfield(nuts.coreg,'orientation'))
    set(findobj('Tag','nut_orientation_menu'),'Value',nuts.coreg.orientation);
else
    set(findobj('Tag','nut_orientation_menu'),'Value',1);
end
if(~exist(nuts.coreg.mripath,'file'))
    warndlg(['Cannot find MRI:' nuts.coreg.mripath]);
    nuts.coreg.mripath = which('blank.img');
end
V = spm_vol(nuts.coreg.mripath);

if (strcmp(spm('ver'),'SPM2'))
if(V.private.hdr.dime.dim(1) < 1 | V.private.hdr.dime.dim(1) > 15)
    warning('looks like you''s guys got you''s endians all flipped around n'' whatnot');
end;
end
% DIM = V.dim(1:3);
% VOX = V.private.hdr.dime.pixdim(2:4);
% ORIGIN = V.private.hdr.hist.origin(1:3);

%display image
spm_image('init',nuts.coreg.mripath);
disp('SPM done initializing...');
nut_spmfig_setup;
