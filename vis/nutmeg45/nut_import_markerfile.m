function varargout = nut_import_markerfile(varargin)
% NUT_IMPORT_MARKERFILE imports CTF MarkerFiles via GUI interface.
%    Is called by nutmeg GUI.

% Last Modified by GUIDE v2.5 14-May-2007 15:29:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @nut_import_markerfile_OpeningFcn, ...
                   'gui_OutputFcn',  @nut_import_markerfile_OutputFcn, ...
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


% --- Executes just before nut_import_markerfile is made visible.
function nut_import_markerfile_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to nut_import_markerfile (see VARARGIN)

% Choose default command line output for nut_import_markerfile
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global nuts markers

if ~isfield(nuts,'meg')
    delete(gcf)
    errordlg('You need to load MEG data first.')
    return
end

if ~isfield(nuts.meg,'markers')
    nut_loadmarkerfile(handles);
else
    [mkpath,mkfile]=fileparts(nuts.meg.markers.file);
    set(handles.text_markerfilename,'string',mkfile)
    if isfield(nuts.meg.markers,'activemarker')
        set(handles.check_active,'Enable','on','Value',1)
        set(handles.pop_active,'Enable','on','string',nuts.meg.markers.activemarker)
    end
    if isfield(nuts.meg.markers,'controlmarker')
        set(handles.check_control,'Enable','on','Value',1)
        set(handles.pop_control,'Enable','on','string',nuts.meg.markers.controlmarker)
    end
end
uiwait

% UIWAIT makes nut_import_markerfile wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = nut_import_markerfile_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;


%-------------
function nut_loadmarkerfile(handles) 

global nuts markers

dum=pwd;
if exist([nuts.meg.filename '.ds'],'dir')
    cd([nuts.meg.filename '.ds'])
end

[mkfile, mkpath] = uigetfile('*.mrk;*.cls','Load CTF MarkerFile/ClassFile...');
if ~ischar(mkfile)
    return
elseif strcmpi(mkfile(end-3:end),'.mrk')
    markers=ctf_read_markerfile([],fullfile(mkpath,mkfile));
    set([handles.text_classtime handles.text3 handles.text4],'Enable','off')
elseif strcmpi(mkfile(end-3:end),'.cls')
    class=ctf_read_classfile([],fullfile(mkpath,mkfile));
    set([handles.text_classtime handles.text3 handles.text4],'Enable','on')
    tim=str2num(get(handles.text_classtime','string'));
    for k=1:class.number
        markers(k).file=[class.path filesep mkfile];
        markers(k).marker_names=class.data(k).name;
        markers(k).trial_times=[class.data(k).trials' tim*ones(class.data(k).Ntrials,1)];
    end
else
    warndlg('File is invalid.')
    return
end

set(handles.text_markerfilename,'string',mkfile)
set([handles.check_active handles.check_control handles.pop_active handles.pop_control],'Enable','on')
set([handles.pop_active handles.pop_control],'string',{markers.marker_names}')

cd(dum)

% --- Executes on button press in push_ok.
function push_ok_Callback(hObject, eventdata, handles)
% hObject    handle to push_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global nuts markers

if ( ~isempty(markers) && (get(handles.check_active,'Value') || get(handles.check_control,'Value')) )
    nuts.meg.markers.file=markers(1).file;
    
    % clear previous settings
    nuts.meg.markers.active=zeros(size(nuts.meg.data,3),1);
    nuts.meg.markers.control=zeros(size(nuts.meg.data,3),1);
    if isfield(nuts.meg.markers,'activemarker'), nuts.meg.markers=rmfield(nuts.meg.markers,'activemarker'); end
    if isfield(nuts.meg.markers,'controlmarker'), nuts.meg.markers=rmfield(nuts.meg.markers,'controlmarker'); end
    
    if get(handles.check_active,'Value')
        selmrk=get(handles.pop_active,'Value');
        t_t=markers(selmrk).trial_times;
        goodact=[1;1+find(diff(t_t(:,1)))];    % only use first occurrance of marker in each trial
        nuts.meg.markers.active=nan(size(nuts.meg.data,3),1);
        nuts.meg.markers.active(t_t(goodact,1))=t_t(goodact,2)*1000;  % set values to milliseconds
        nuts.meg.markers.activemarker=markers(selmrk).marker_names;
    end
    if get(handles.check_control,'Value')
        selmrk=get(handles.pop_control,'Value');
        t_t=markers(selmrk).trial_times;
        goodcon=[1;1+find(diff(t_t(:,1)))];    % only use first occurrance of marker in each trial
        nuts.meg.markers.control=nan(size(nuts.meg.data,3),1);
        nuts.meg.markers.control(t_t(goodcon,1))=t_t(goodcon,2)*1000;  % set values to milliseconds
        nuts.meg.markers.controlmarker=markers(selmrk).marker_names;
    end
    if abs(length(goodact)-length(goodcon))>min([length(goodact) length(goodcon)])/10
        warndlg('WARNING: The number of valid trials for active and control windows differs by more than 10%.')
    end
elseif ( isfield(nuts.meg,'markers') && (~get(handles.check_active,'Value') && ~get(handles.check_control,'Value')) )    % clear markers if user deselected checkboxes
    nuts.meg=rmfield(nuts.meg,'markers');
end
clear global markers
delete(handles.output)

% --- Executes on button press in push_cancel.
function push_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to push_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear global markers
delete(gcf)

% --- Executes on button press in check_active.
function check_active_Callback(hObject, eventdata, handles)
% hObject    handle to check_active (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_active
global nuts

if ~get(hObject,'Value')
    if isfield(nuts.meg.markers,'active')
        nuts.meg.markers.active=zeros(size(nuts.meg.data,3),1);
    end
    if isfield(nuts.meg.markers,'activemarker')
        nuts.meg.markers=rmfield(nuts.meg.markers,'activemarker');
    end
end

% --- Executes on button press in check_control.
function check_control_Callback(hObject, eventdata, handles)
% hObject    handle to check_control (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_control
global nuts

if ~get(hObject,'Value')
    if isfield(nuts.meg.markers,'control')
        nuts.meg.markers.control=zeros(size(nuts.meg.data,3),1);
    end
    if isfield(nuts.meg.markers,'controlmarker')
        nuts.meg.markers=rmfield(nuts.meg.markers,'controlmarker');
    end
end
    

% --- Executes on selection change in pop_active.
function pop_active_Callback(hObject, eventdata, handles)
% hObject    handle to pop_active (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns pop_active contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_active


% --- Executes during object creation, after setting all properties.
function pop_active_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_active (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_control.
function pop_control_Callback(hObject, eventdata, handles)
% hObject    handle to pop_control (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns pop_control contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_control


% --- Executes during object creation, after setting all properties.
function pop_control_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_control (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_load_markerfile.
function push_load_markerfile_Callback(hObject, eventdata, handles)
% hObject    handle to push_load_markerfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nut_loadmarkerfile(handles);


function text_classtime_Callback(hObject, eventdata, handles)
% hObject    handle to text_classtime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_classtime as text
%        str2double(get(hObject,'String')) returns contents of text_classtime as a double
global markers
tim=str2num(get(handles.text_classtime','string'));
for k=1:length(markers)
    markers(k).trial_times(:,2)=tim;
end

% --- Executes during object creation, after setting all properties.
function text_classtime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_classtime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


