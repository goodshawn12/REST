function varargout = fcm_roiselect_gui(varargin)
% FCM_ROISELECT_GUI 

% Last Modified by GUIDE v2.5 16-Nov-2011 18:25:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fcm_roiselect_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @fcm_roiselect_gui_OutputFcn, ...
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


% --- Executes just before fcm_roiselect_gui is made visible.
function fcm_roiselect_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fcm_roiselect_gui (see VARARGIN)

% Choose default command line output for fcm_roiselect_gui
%handles.output = '';

% Update handles structure
%guidata(hObject, handles);

% UIWAIT makes fcm_roiselect_gui wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = fcm_roiselect_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if isstruct(handles) && isfield(handles,'list_sel')
    varargout{1} = get(handles.list_sel,'string');
    delete(handles.figure1);
else
    varargout{1} = '';
end


% --- Executes on selection change in list_unsel.
function list_unsel_Callback(hObject, eventdata, handles)
% hObject    handle to list_unsel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns list_unsel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_unsel


% --- Executes during object creation, after setting all properties.
function list_unsel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_unsel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
L = fcm_anatroi('-list');
set(hObject,'string',L(:,end))

% --- Executes on selection change in list_sel.
function list_sel_Callback(hObject, eventdata, handles)
% hObject    handle to list_sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns list_sel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_sel


% --- Executes during object creation, after setting all properties.
function list_sel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_select.
function push_select_Callback(hObject, eventdata, handles)
% hObject    handle to push_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global st

S = get(handles.list_sel,'string');
U = get(handles.list_unsel,'string');
u = get(handles.list_unsel,'value');
S = [S; U(u)];
set(handles.list_sel,'string',S)
U(u) = [];
set(handles.list_unsel,'string',U,'Value',1)

roi=fcm_anatroi(S);

if get(handles.view,'value')
    vmm = nut_mni2mri(double(roi.MNIvoxels));
    XYZ= nut_mm2voxels(vmm);
    %vvox=nut_mm2voxels(vmm);
    %XYZ = unique([floor(vvox);ceil(vvox)],'rows');
    %XYZ = unique([XYZ;round(vvox)],'rows'); clear vvox
    Z=ones(size(XYZ,1),1);

    figure(st.fig)
    spm_orthviews('rmblobs',1)
    spm_orthviews('addcolouredblobs',1,XYZ.',Z,st.vols{1}.mat,[1 0 0])
    spm_orthviews('reposition',mean(vmm))
end

% --- Executes on button press in push_unselect.
function push_unselect_Callback(hObject, eventdata, handles)
% hObject    handle to push_unselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global st

S = get(handles.list_sel,'string');
s = get(handles.list_sel,'value');
S(s) = [];
set(handles.list_sel,'string',S,'Value',1)

U = fcm_anatroi('-list');
U = U(:,end);
U(ismember(U,S))=[];
set(handles.list_unsel,'string',U)

if get(handles.view,'value')
    spm_orthviews('rmblobs',1)
    if ~isempty(S)
        roi=fcm_anatroi(S);
        vmm = nut_mni2mri(double(roi.MNIvoxels));
        XYZ= nut_mm2voxels(vmm);
        Z=ones(size(XYZ,1),1);

        figure(st.fig)
        spm_orthviews('addcolouredblobs',1,XYZ.',Z,st.vols{1}.mat,[1 0 0])
        spm_orthviews('reposition',mean(vmm))
    end
end

% --- Executes on button press in view.
function view_Callback(hObject, eventdata, handles)
% hObject    handle to view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global st nuts

spm_orthviews('rmblobs',1)
if get(hObject,'Value')
    if isempty(findobj('tag','nutmegfig'))
        if isempty(nuts)
            error('You must load session file first.')
        end
        nutmeg(nuts);
    end

    S = get(handles.list_sel,'string');
    if ~isempty(S)
        roi=fcm_anatroi(S);
        vmm = nut_mni2mri(double(roi.MNIvoxels));
        XYZ= nut_mm2voxels(vmm);
        Z=ones(size(XYZ,1),1);

        figure(st.fig)
        spm_orthviews('addcolouredblobs',1,XYZ.',Z,st.vols{1}.mat,[1 0 0])
        spm_orthviews('reposition',mean(vmm))
    end
end

% Hint: get(hObject,'Value') returns toggle state of view
% --- Executes on button press in push_ok.
function push_ok_Callback(hObject, eventdata, handles)
% hObject    handle to push_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);

% --- Executes on button press in push_cancel.
function push_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to push_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%set(handles.list_sel,'string','')
delete(handles.figure1);
