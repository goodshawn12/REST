function varargout = nut_maketrial_gui(varargin)
% NUT_MAKETRIAL_GUI M-file for nut_maketrial_gui.fig
%      NUT_MAKETRIAL_GUI, by itself, creates a new NUT_MAKETRIAL_GUI or raises the existing
%      singleton*.
%
%      H = NUT_MAKETRIAL_GUI returns the handle to a new NUT_MAKETRIAL_GUI or the handle to
%      the existing singleton*.
%
%      NUT_MAKETRIAL_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NUT_MAKETRIAL_GUI.M with the given input arguments.
%
%      NUT_MAKETRIAL_GUI('Property','Value',...) creates a new NUT_MAKETRIAL_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before nut_maketrial_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to nut_maketrial_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help nut_maketrial_gui

% Last Modified by GUIDE v2.5 03-Aug-2011 15:21:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @nut_maketrial_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @nut_maketrial_gui_OutputFcn, ...
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


% --- Executes just before nut_maketrial_gui is made visible.
function nut_maketrial_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to nut_maketrial_gui (see VARARGIN)

% Choose default command line output for nut_maketrial_gui
%handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes nut_maketrial_gui wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function out = nut_maketrial_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
global MT
out=MT;
clear global MT

% --- Executes on button press in push_ok.
function push_ok_Callback(hObject, eventdata, handles)
% hObject    handle to push_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global MT

MT=[];  % re-init
MT.ok = 1;  % In case no options were selected
if get(handles.check_maketrial,'Value')
    if get(handles.radio_fixed,'Value')
        MT.trialdef=struct('prestim',str2num(get(handles.edit_pre,'string')),'poststim',str2num(get(handles.edit_post,'string')));
    elseif get(handles.radio_trialfun,'Value')
        MT.trialfun=get(handles.edit_callback,'string');
    end
elseif get(handles.check_rema,'Value')
    MT.rema.secpertrial = str2num(get(handles.edit_secpertrial,'string'));
end

delete(gcf)

% --- Executes on button press in push_cancel.
function push_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to push_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global MT
MT=[];
delete(gcf)

% --- Executes on button press in check_maketrial.
function check_maketrial_Callback(hObject, eventdata, handles)
% hObject    handle to check_maketrial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_maketrial
set(handles.check_rema,'Value',0)

% --- Executes on button press in check_rema.
function check_rema_Callback(hObject, eventdata, handles)
% hObject    handle to check_rema (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_rema
set(handles.check_maketrial,'Value',0)

% --- Executes when selected object is changed in uipanel1.
function uipanel1_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel1 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

function edit_pre_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pre (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pre as text
%        str2double(get(hObject,'String')) returns contents of edit_pre as a double


% --- Executes during object creation, after setting all properties.
function edit_pre_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pre (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_post_Callback(hObject, eventdata, handles)
% hObject    handle to edit_post (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_post as text
%        str2double(get(hObject,'String')) returns contents of edit_post as a double


% --- Executes during object creation, after setting all properties.
function edit_post_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_post (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_callback_Callback(hObject, eventdata, handles)
% hObject    handle to edit_callback (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_callback as text
%        str2double(get(hObject,'String')) returns contents of edit_callback as a double

% --- Executes during object creation, after setting all properties.
function edit_callback_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_callback (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_secpertrial_Callback(hObject, eventdata, handles)
% hObject    handle to edit_secpertrial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_secpertrial as text
%        str2double(get(hObject,'String')) returns contents of edit_secpertrial as a double

% --- Executes during object creation, after setting all properties.
function edit_secpertrial_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_secpertrial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


