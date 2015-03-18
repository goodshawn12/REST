function varargout = nut_GenMesh_tool(varargin)
% NUT_GENMESH_TOOL M-file for nut_GenMesh_tool.fig
%   A standard .fig-associated M-file that can launch the GUI, return the
%   figure handle and process some callbacks.

warning('off')  % hack so warning for "created by matlab 7.0" is turned off
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @nut_GenMesh_tool_OpeningFcn, ...
                   'gui_OutputFcn',  @nut_GenMesh_tool_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
warning('on')  % we want to see all other warnings.


%%-----------------------------------------------
function nut_GenMesh_tool_OpeningFcn(hObject, eventdata, handles, varargin)
% --- Executes just before nut_GenMesh_tool is made visible.
handles.output = hObject;
guidata(hObject, handles);
global nuts
if isempty(nuts)
   msgbox('Load the image','Volume not found','warn')
   return;
end

%%-----------------------------------------------
function varargout = nut_GenMesh_tool_OutputFcn(hObject, eventdata, handles)
% --- Outputs from this function are returned to the command line.
varargout{1} = handles.output;



function nut_renderhead_button_Callback(hObject, eventdata, handles)
% --- [Generate Mesh] button callback.
global st coreg
threshold = str2num(get(findobj('Tag','nut_vol2cloud_text'),'String'));
if (isempty(threshold))
    msgbox('Not a valid threshold value',...
   'THRESHOLD VALUE','warn');
    return;
else            
    position=nut_vol2cloud6(threshold,spm_read_vols(st.vols{1}));
end
 
if(isfield(coreg,'mesh'))
    coreg = rmfield(coreg,'mesh');
end
coreg.mesh = nut_cloud2mesh(position);
coreg.mesh=nut_mesh_smoothing(coreg.mesh);
coreg.mesh=nut_optmesh(coreg.mesh);
position = nut_voxels2mm(position);  % makes more sense to plot mm
figure('Units','normalized','Position',[.05 .2 .6 .6]);
subplot(2,3,1)
plot(position(:,1),position(:,2),'.'); axis equal off;
title('axial view')
subplot(2,3,2);
plot(position(:,1),position(:,3),'.'); axis equal off;
title('coronal view')
subplot(2,3,3)
plot(position(:,2),position(:,3),'.'); axis equal off;
title('sagittal view')
subplot('position',[.25 .05 .5 .5]);
[m1,m2,m3]=size(coreg.mesh);
coreg.mesh=reshape(coreg.mesh,m1*m2,m3);
coreg.mesh=nut_voxels2mm(coreg.mesh);  %we'll keep mesh in mri mm coords from now on
coreg.mesh=reshape(coreg.mesh,m1,m2,m3);
nut_show_head(coreg.mesh);
nut_coreg_enabler


function nut_vol2cloud_text_Callback(hObject, eventdata, handles)
nut_renderhead_button_Callback(hObject);


function nut_thresholdselect_Callback(hObject, eventdata, handles)
% -- Autoselect threshold button.
global st
thresh = nut_find_thresh(spm_read_vols(st.vols{1}));
threshold=num2str(thresh);
set(findobj('Tag','nut_vol2cloud_text'),'String',threshold);
