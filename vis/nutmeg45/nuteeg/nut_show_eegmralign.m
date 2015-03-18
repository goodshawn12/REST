function nut_show_eegmralign
% Displays sensor positions over head.

global nuts 

if ~isfield(nuts,'meg') || ~isfield(nuts.meg,'sensorCoord')
    return
end

handles.fig = findobj('tag','nut_show_eegmralign');
if isempty(handles.fig)
    handles.fig = figure('units','normalized','tag','nut_show_eegmralign','CloseRequestFcn',@closeit);
else
    figure(handles.fig);
    cla
end

handles.hires = uicontrol('Style','checkbox','Units','normalized',...
    'Position',[.1 .05 .1 .05], 'string', 'Hi Res', ...
    'Callback',@plotdahead);

handles.toggle = uicontrol('Style','pushbutton','Units','normalized',...
    'Position',[.25 .05 .2 .05], 'string', 'Toggle Labels', ...
   'Callback',@toggle);

handles.project = uicontrol('Style','pushbutton','Units','normalized',...
    'Position',[.5 .05 .2 .05], 'string', 'Project Sensors', ...
   'Callback',@project);

handles.warp = uicontrol('Style','pushbutton','Units','normalized',...
    'Position',[.75 .05 .2 .05], 'string', 'Warp MR', ...
   'Callback',@warp);

guidata(handles.fig,handles);

plotdahead([],[],handles);

function plotdahead(hObject,dum2,handles)

global st mesh nuts
if nargin<3
    handles = guidata(hObject);
end

sensCoord = nuts.meg.sensorCoord .* 1.05; % multiply to keep them visible on head

if get(handles.hires,'Value')
    sensCoord = nut_meg2mri(sensCoord);
    
    Y = spm_read_vols(st.vols{1});
    thresh = nut_find_thresh(Y);

    position=nut_vol2cloud6(thresh,Y);
    position = nut_voxels2mm(position); 

    mesh = nut_cloud2mesh(position);
    mesh = nut_mesh_smoothing(mesh);
    mesh = nut_optmesh(mesh);
      
	nut_show_head(mesh);
    
else
      
    if ~isfield(nuts.coreg,'volfile')
        [file,path] = uigetfile('*_vol.mat','Open Mesh File');
        if isequal(file,0)
            return
        else
            nuts.coreg.volfile = fullfile(path,file);
        end
    end
    load(nuts.coreg.volfile);
    if isfield(vol.bnd(1),'tri'); vol.bnd(1).faces = vol.bnd(1).tri; end;   % Handle fieldtrip vols
    if isfield(vol.bnd(1),'pnt'); vol.bnd(1).vertices = vol.bnd(1).pnt; end;
    mesh = PrepareTriangleMesh(nut_mri2meg(vol.bnd(1).vertices),vol.bnd(1).faces(:,[1 3 2]));

    figure(handles.fig);
    trisurf(mesh.e,mesh.p(:,1),mesh.p(:,2),mesh.p(:,3),max(mesh.p,[],2),'EdgeColor','none');
    colormap copper
end

hold on
plot3(sensCoord(:,1),sensCoord(:,2),sensCoord(:,3),'wo','MarkerFaceColor','w','MarkerSize',10);
hold off
set(gca,'visible','off')
rotate3d on
%TODO: how to prevent rotation from changing size of head?

handles.label=zeros(1,size(sensCoord,1));
for k=1:size(sensCoord,1)
    handles.label(k)=text(sensCoord(k,1),sensCoord(k,2),sensCoord(k,3),nuts.meg.sensor_labels{k}, ...
        'visible','off','horizontalalignment','center');
end
guidata(handles.fig,handles);

function project(hObject,dum2)

nut_eegmralign('Project Sensors');
plotdahead(hObject);

function warp(hObject,dum2)

clear global mesh
nut_eegmralign('Warp MR');
plotdahead(hObject);

function toggle(hObject,dum2)

handles = guidata(hObject);
setting=setdiff(['on ';'off'],get(handles.label(1),'visible'),'rows');
set(handles.label,'visible',setting)

function closeit(dum1,dum2)

clear global mesh
closereq

%------------
function [Y,XYZ] = spm_read_vols(V,mask)
% Read in entire image volumes
% FORMAT [Y,XYZ] = spm_read_vols(V,mask)
% V    - vector of mapped image volumes to read in (from spm_vol)
% mask - implicit zero mask?
% Y    - 4D matrix of image data, fourth dimension indexes images
% XYZ  - 3xn matrix of XYZ locations returned
%_______________________________________________________________________
%
% For image data types without a representation of NaN (see spm_type),
% implicit zero masking can be used. If mask is set, then zeros are
% treated as masked, and returned as NaN.
%_______________________________________________________________________
% @(#)spm_read_vols.m	2.5 Andrew Holmes 00/03/21
% 
% Modified for NUTEEG - Daniel D.E. Wong


%-Argument checks
%-----------------------------------------------------------------------
if nargin<2, mask = 0; end
if nargin<1, error('insufficient arguments'), end

%-Image dimension, orientation and voxel size checks
%-----------------------------------------------------------------------
if length(V)>1 & any(any(diff(cat(1,V.dim),1,1),1)&[1,1,1,0])
	error('images don''t all have the same dimensions'), end
if any(any(any(diff(cat(3,V.mat),1,3),3)))
	error('images don''t all have same orientation & voxel size'), end

%-Read in image data
%-----------------------------------------------------------------------
n  = prod(size(V));			%-#images
Y = zeros([V(1).dim(1:3),n]);		%-image data matrix

for i=1:n, for p=1:V(1).dim(3)
	Y(:,:,p,i) = spm_slice_vol(V(i),spm_matrix([0 0 p]),V(i).dim(1:2),0);
end, end

%-Apply implicit zero mask for image datatypes without a NaNrep
%-----------------------------------------------------------------------
if mask
	%-Work out images without NaNrep
	im = logical(zeros(n,1));
	for i=1:n, im(i)=~spm_type(V(i).dim(4),'NaNrep'); end
	%-Mask
	Y(Y(:,:,:,im)==0)=NaN;
end

%-Return as 3D matrix if single image
%-----------------------------------------------------------------------
% if n==1; Y=Y(:,:,:,1); end

%-Compute XYZ co-ordinates (if required)
%-----------------------------------------------------------------------
if nargout>1
	[R,C,P]=ndgrid(1:V(1).dim(1),1:V(1).dim(2),1:V(1).dim(3));
    R = single(R); C = single(C); P = single(P);
	RCP = [R(:)';C(:)';P(:)'];
	RCP(4,:)=1;
	XYZ = V(1).mat*RCP;
	XYZ=XYZ(1:3,:);
end
