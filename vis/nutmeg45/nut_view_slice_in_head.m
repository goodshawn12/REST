function nut_view_slice_in_head
% NUT_VIEW_SLICE_IN_HEAD
%
% This function will allow you to see the portions inside the head surface.

global nuts;

%get the volume from the image
V=spm_vol(nuts.coreg.mripath);
vol = spm_read_vols(V);
dim=size(vol);
% generate PROMPT for the user
def1=int2str([1 dim(1)]);
def2=int2str([1 dim(2)]);
def3=int2str([1 dim(3)]);
prompt   = {'Select range for sagittal view, x','Select range for coronal view, y',...
        'Select range for axial view, z'};
title    = 'It slices! It dices!';
lines = 1;
def      = {def1 def2 def3};   
answer   = inputdlg(prompt,title,lines,def);
if (isempty(answer)) 
    return;
end;
% read the values from PROMPT 
sag=['sagittal=[' answer{1} '];'];
eval(sag);
cor=['coronal=[' answer{2} '];'];
eval(cor);
ax=['axial=[' answer{3} '];'];
eval(ax);
% smooth the 3D image
v=smooth3(vol(sagittal(1):sagittal(2),coronal(1):coronal(2),axial(1):axial(2)));
figure
% create an isosurface and its patch
p=patch(isosurface(v,'noshare'), 'FaceColor', 'yellow', 'EdgeColor', 'none');
% get isocaps 
patch(isocaps(vol(sagittal(1):sagittal(2),coronal(1):coronal(2),axial(1):axial(2)),'noshare'),...
    'FaceColor', 'interp', 'EdgeColor', 'none');
% find normals to the isosurface
isonormals(v,p);
% set display property
view(3); 
axis tight;  
daspect([1 1 1])
colormap(gray)
camlight; lighting phong 
return