function varargout = visORICA(varargin)
% VISORICA MATLAB code for visORICA.fig
%      VISORICA, by itself, creates a new VISORICA or raises the existing
%      singleton*.
%
%      H = VISORICA returns the handle to a new VISORICA or the handle to
%      the existing singleton*.
%
%      VISORICA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VISORICA.M with the given input arguments.
%
%      VISORICA('Property','Value',...) creates a new VISORICA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before visORICA_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to visORICA_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help visORICA

% Last Modified by GUIDE v2.5 05-Nov-2014 17:08:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @visORICA_OpeningFcn, ...
                   'gui_OutputFcn',  @visORICA_OutputFcn, ...
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


% --- Executes just before visORICA is made visible.
function visORICA_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to visORICA (see VARARGIN)

inlet = vis_stream_ORICA('StreamName','EEGLAB','figurehandles',handles.figure1,'axishandles',handles.axisEEG);
1;
eegTimer = timerfind;
eegTimer = eegTimer(strcmp(eegTimer.Name,'eeg_timer'));
% if switch to components, change eegTimer.TimerFcn and add that fcn to
% vis_stream_ORICA




% Choose default command line output for visORICA
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes visORICA wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = visORICA_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenuEEG.
function popupmenuEEG_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuEEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuEEG contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuEEG


% --- Executes during object creation, after setting all properties.
function popupmenuEEG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuEEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonIC.
function pushbuttonIC_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonIC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in togglebuttonSaveIC1.
function togglebuttonSaveIC1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonSaveIC1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonSaveIC1


% --- Executes on button press in togglebuttonExcludeIC1.
function togglebuttonExcludeIC1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonExcludeIC1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonExcludeIC1


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in togglebuttonExcludeIC8.
function togglebuttonExcludeIC8_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonExcludeIC8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonExcludeIC8


% --- Executes on button press in togglebuttonSaveIC8.
function togglebuttonSaveIC8_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonSaveIC8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonSaveIC8


% --- Executes on button press in togglebuttonExcludeIC7.
function togglebuttonExcludeIC7_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonExcludeIC7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonExcludeIC7


% --- Executes on button press in togglebuttonSaveIC7.
function togglebuttonSaveIC7_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonSaveIC7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonSaveIC7


% --- Executes on button press in togglebuttonExcludeIC6.
function togglebuttonExcludeIC6_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonExcludeIC6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonExcludeIC6


% --- Executes on button press in togglebuttonSaveIC6.
function togglebuttonSaveIC6_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonSaveIC6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonSaveIC6


% --- Executes on button press in togglebuttonExcludeIC5.
function togglebuttonExcludeIC5_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonExcludeIC5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonExcludeIC5


% --- Executes on button press in togglebuttonSaveIC5.
function togglebuttonSaveIC5_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonSaveIC5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonSaveIC5


% --- Executes on button press in togglebuttonExcludeIC4.
function togglebuttonExcludeIC4_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonExcludeIC4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonExcludeIC4


% --- Executes on button press in togglebuttonSaveIC4.
function togglebuttonSaveIC4_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonSaveIC4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonSaveIC4


% --- Executes on button press in togglebuttonExcludeIC3.
function togglebuttonExcludeIC3_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonExcludeIC3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonExcludeIC3


% --- Executes on button press in togglebuttonSaveIC3.
function togglebuttonSaveIC3_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonSaveIC3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonSaveIC3


% --- Executes on button press in togglebuttonExcludeIC2.
function togglebuttonExcludeIC2_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonExcludeIC2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonExcludeIC2


% --- Executes on button press in togglebuttonSaveIC2.
function togglebuttonSaveIC2_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonSaveIC2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonSaveIC2
