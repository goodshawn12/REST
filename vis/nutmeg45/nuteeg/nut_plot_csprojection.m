function varargout = nut_plot_csprojection(varargin)
% NUT_PLOT_CSPROJECTION M-file for nut_plot_csprojection.fig
%      NUT_PLOT_CSPROJECTION, by itself, creates a new NUT_PLOT_CSPROJECTION or raises the existing
%      singleton*.
%
%      H = NUT_PLOT_CSPROJECTION returns the handle to a new NUT_PLOT_CSPROJECTION or the handle to
%      the existing singleton*.
%
%      NUT_PLOT_CSPROJECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NUT_PLOT_CSPROJECTION.M with the given input arguments.
%
%      NUT_PLOT_CSPROJECTION('Property','Value',...) creates a new NUT_PLOT_CSPROJECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before nut_plot_csprojection_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to nut_plot_csprojection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% @author Daniel D.E. Wong
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help nut_plot_csprojection

% Last Modified by GUIDE v2.5 10-Nov-2009 16:47:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @nut_plot_csprojection_OpeningFcn, ...
                   'gui_OutputFcn',  @nut_plot_csprojection_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before nut_plot_csprojection is made visible.
function nut_plot_csprojection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to nut_plot_csprojection (see VARARGIN)

% Choose default command line output for nut_plot_csprojection
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes nut_plot_csprojection wait for user response (see UIRESUME)
% uiwait(handles.csprojectionfig);
global csprojection
% Set to initial values from gui
csprojection.beam = [];
csprojection.mesh = [];
csprojection.ori = [];
csprojection.t1 = str2num(get(handles.t1,'String'));
csprojection.t2 = str2num(get(handles.t2,'String'));
csprojection.scale1 = str2num(get(handles.scale1,'String'));
csprojection.scale2 = str2num(get(handles.scale2,'String'));
csprojection.depth = str2num(get(handles.depth,'String'));
csprojection.alpha = str2num(get(handles.alpha,'String'));

set(handles.csprojectionfig,'CloseRequestFcn',@closeGUI);


% --- Outputs from this function are returned to the command line.
function varargout = nut_plot_csprojection_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function cs_Callback(hObject, eventdata, handles)
% hObject    handle to cs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cs as text
%        str2double(get(hObject,'String')) returns contents of cs as a double
global csprojection
filename = get(hObject,'String');
try; load(filename); catch; msgbox('No such file.', 'CS Projection','error'); mesh = []; end;
if ~exist('cs')
    mesh = [];
else
    mesh = PrepareTriangleMesh(nut_mri2meg(cs.vertices),cs.faces(:,[1 3 2]));
end
csprojection.mesh = mesh;
if isempty(mesh); set(hObject,'String',''); guidata(hObject, handles); end;


% --- Executes during object creation, after setting all properties.
function cs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in open_cs.
function open_cs_Callback(hObject, eventdata, handles)
% hObject    handle to open_cs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global csprojection
[file,path] = uigetfile('*_cs.mat','Open Cortical Surface Mesh File');
if isequal(file,0) | isequal(path,0); return; end;
filename = strcat(path,file);
set(handles.cs,'String',filename);
load(filename);
if ~exist('cs'); csprojection.mesh = []; return; end;
mesh = PrepareTriangleMesh(nut_mri2meg(cs.vertices),cs.faces(:,[1 3 2]));
csprojection.mesh = mesh;
guidata(hObject, handles);


function sbeam_Callback(hObject, eventdata, handles)
% hObject    handle to sbeam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sbeam as text
%        str2double(get(hObject,'String')) returns contents of sbeam as a double
global csprojection
filename = get(hObject,'String');
try
    beam = load(filename);
    beam = nut_beam_legacy_compatibility(beam);
catch
    msgbox('No such file.', 'CS Projection','error');
    set(handles.sbeam,'String','');
end
if ~exist('beam')
    csprojection.beam = [];
    return;
else
    if isfield(beam,'s')
        beam.s{1} = sqrt(sum(beam.s{1}(:,:,1,:).^2,4));
    else
        errordlg('Invalid s_beam file.'); return
    end
end
csprojection.beam = beam;
set(handles.t1,'String','0'); csprojection.t1 = 0;
set(handles.t2,'String',num2str(beam.timepts(end))); csprojection.t2 = beam.timepts(end);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function sbeam_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sbeam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in open_sbeam.
function open_sbeam_Callback(hObject, eventdata, handles)
% hObject    handle to open_sbeam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global csprojection
[file,path] = uigetfile('s_beam*.mat','Open Cortical Surface Mesh File');
if isequal(file,0) | isequal(path,0); return; end;
filename = strcat(path,file);
set(handles.sbeam,'String',filename);
beam = load(filename);
beam = nut_beam_legacy_compatibility(beam);
if ~exist('beam')
    beam = [];
    return;
end
if isfield(beam,'s')
    beam.s{1} = sqrt(sum(beam.s{1}(:,:,1,:).^2,4));
else
    errordlg('Invalid s_beam file.'); return
end

csprojection.beam = beam;
set(handles.t1,'String','0'); csprojection.t1 = 0;
set(handles.t2,'String',num2str(beam.timepts(end))); csprojection.t2 = beam.timepts(end);
guidata(hObject, handles);



function depth_Callback(hObject, eventdata, handles)
% hObject    handle to depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of depth as text
%        str2double(get(hObject,'String')) returns contents of depth as a double
global csprojection
csprojection.depth = str2num(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function depth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to depth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in display.
function display_Callback(hObject, eventdata, handles)
% hObject    handle to display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global csprojection

if isempty(csprojection.beam) | isempty(csprojection.mesh)
    msgbox('Beamformer output or cortical surface not loaded.','CS Projection','error');
    return;
end

beam = csprojection.beam;
mesh = csprojection.mesh;
ori = csprojection.ori;

t1 = dsearchn(beam.timepts,csprojection.t1);
t2 = dsearchn(beam.timepts,csprojection.t2);

samin = min(mean(csprojection.beam.s{1}(:,t1:t2),2));
samax = max(mean(csprojection.beam.s{1}(:,t1:t2),2));

scale1 = samin + csprojection.scale1/1000 * (samax-samin);
scale2 = samin + csprojection.scale2/1000 * (samax-samin);

voxval = mean(beam.s{1}(:,t1:t2),2);
inclvox = find(voxval >= scale1)';  % Voxels with values above threshold

% Map voxels to faces
fvalfaces = zeros(size(mesh.e,1),1);
fvalvertices = zeros(size(mesh.p,1),1);
for v = inclvox
    vec = mesh.mp-repmat(beam.voxels(v,:),size(mesh.mp,1),1);
    dist = nut_rownorm(vec);
    if ~isempty(ori)
        dp = abs(dot(mesh.un,repmat(ori(v,:),size(mesh.un,1),1),2));
    else
        dp = ones(size(mesh.un,1),1);
    end
    closefaces = find((dist < csprojection.depth) & (voxval(v)*dp > fvalfaces));   % CS faces within projection depth of voxel with fvals less than voxels
    uvec = vec(closefaces,:)./repmat(dist(closefaces),1,3);
    dp2 = abs(dot(uvec,mesh.un(closefaces,:),2));
    losfaces = find(dp2.^2./dist(closefaces).^2 < 1);    % CS faces with LOS to voxel (as a subset of closefaces)
    fvalfaces(closefaces(losfaces)) = voxval(v)*dp(closefaces(losfaces));
    
    %fvalvertices(reshape(mesh.e(closefaces(losfaces),:),length(losfaces)*3,1)) = voxval(v);
    fvv = voxval(v)*dp(closefaces(losfaces));
    fvalvertices(mesh.e(closefaces(losfaces),1)) = fvv;
    fvalvertices(mesh.e(closefaces(losfaces),2)) = fvv;
    fvalvertices(mesh.e(closefaces(losfaces),3)) = fvv;
end

incl = find(fvalfaces);
excl = find(~fvalfaces);

fvalvertices(find(fvalvertices > scale2)) = scale2;

h = figure('DoubleBuffer','on');

colormap(hot(50));


if ~isempty(ori)
    %caxis([scale1 scale2])
    caxis([min(fvalvertices) max(fvalvertices)])
else
    caxis([scale1 scale2])
end
cameratoolbar('SetCoordSys','y')
cameratoolbar('SetMode','orbit') 
camzoom(1);
light
axis equal;
axis off;

p = patch('faces',mesh.e(incl,:),'vertices',mesh.p,'facevertexcdata',fvalvertices,'facecolor','interp','edgecolor','none','facelighting','phong');
colorbar('ycolor',[1 1 1])
hold on
patch('faces',mesh.e(excl,:),'vertices',mesh.p,'facecolor',[0.75 0.75 0.75],'edgecolor','none','facelighting','phong','facealpha',csprojection.alpha);
hold off
set(h,'color',[0 0 0])

dcm_obj = datacursormode(h);
set(dcm_obj,'UpdateFcn',@myupdatefcn);


function t1_Callback(hObject, eventdata, handles)
% hObject    handle to t1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t1 as text
%        str2double(get(hObject,'String')) returns contents of t1 as a double
global csprojection
csprojection.t1 = str2num(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function t1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function t2_Callback(hObject, eventdata, handles)
% hObject    handle to t2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t2 as text
%        str2double(get(hObject,'String')) returns contents of t2 as a double
global csprojection
csprojection.t2 = str2num(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function t2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function scale1_Callback(hObject, eventdata, handles)
% hObject    handle to scale1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scale1 as text
%        str2double(get(hObject,'String')) returns contents of scale1 as a double
global csprojection
csprojection.scale1 = str2num(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function scale1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scale1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function scale2_Callback(hObject, eventdata, handles)
% hObject    handle to scale2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scale2 as text
%        str2double(get(hObject,'String')) returns contents of scale2 as a double
global csprojection
csprojection.scale2 = str2num(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function scale2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scale2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function alpha_Callback(hObject, eventdata, handles)
% hObject    handle to alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha as text
%        str2double(get(hObject,'String')) returns contents of alpha as a double
global csprojection
csprojection.alpha = str2num(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function ori_Callback(hObject, eventdata, handles)
% hObject    handle to ori (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ori as text
%        str2double(get(hObject,'String')) returns contents of ori as a double
global csprojection
filename = get(hObject,'String');
try; load(filename); catch; msgbox('No such file.', 'CS Projection','error'); ori = []; end;
if ~exist('ori'); ori = []; end;
csprojection.ori = ori;
if isempty(ori); set(hObject,'String',''); guidata(hObject, handles); end;


% --- Executes during object creation, after setting all properties.
function ori_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ori (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in open_ori.
function open_ori_Callback(hObject, eventdata, handles)
% hObject    handle to open_ori (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global csprojection
[file,path] = uigetfile('*_ori.mat','Open Dipole Orientation File');
if isequal(file,0) | isequal(path,0); return; end;
filename = strcat(path,file);
load(filename);
if ~exist('ori'); ori = []; filename=''; end;
csprojection.ori = ori;
set(handles.ori,'String',filename);
guidata(hObject, handles);


function closeGUI(hObject,eventdata)
clear global csprojection
delete(gcf)