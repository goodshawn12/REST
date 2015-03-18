function varargout = nut_normsettings_gui(varargin)
% NUT_NORMSETTINGS_GUI M-file for nut_normsettings_gui.fig
%      NUT_NORMSETTINGS_GUI, by itself, creates a new NUT_NORMSETTINGS_GUI or raises the existing
%      singleton*.
%
%      H = NUT_NORMSETTINGS_GUI returns the handle to a new NUT_NORMSETTINGS_GUI or the handle to
%      the existing singleton*.
%
%      NUT_NORMSETTINGS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NUT_NORMSETTINGS_GUI.M with the given input arguments.
%
%      NUT_NORMSETTINGS_GUI('Property','Value',...) creates a new NUT_NORMSETTINGS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before nut_normsettings_gui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to nut_normsettings_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help nut_normsettings_gui

% Last Modified by GUIDE v2.5 08-Jan-2007 17:10:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @nut_normsettings_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @nut_normsettings_gui_OutputFcn, ...
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


% --- Executes just before nut_normsettings_gui is made visible.
function nut_normsettings_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to nut_normsettings_gui (see VARARGIN)

% Choose default command line output for nut_normsettings_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
global settings
settings=struct('all',logical(0),'newvoxsize',[]);

% UIWAIT makes nut_normsettings_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = nut_normsettings_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

uiwait
global settings
varargout{1} = settings;
clear global settings


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ok_push.
function ok_push_Callback(hObject, eventdata, handles)
% hObject    handle to ok_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global settings

if get(handles.current_radio,'Value')
    settings.all=logical(0);
elseif get(handles.all_radio,'Value')
    settings.all=logical(1);
else
    return
end

if get(handles.nonnorm_radio,'Value')
    settings.newvoxsize=[];
elseif get(handles.struc_radio,'Value')
    global defaults
    settings.newvoxsize=defaults.normalise.write.vox;
elseif get(handles.custom_radio,'Value')
    settings.newvoxsize=str2num(get(handles.edit1,'string'));
else
    return
end
delete(gcf)

% --- Executes on button press in cancel_push.
function cancel_push_Callback(hObject, eventdata, handles)
% hObject    handle to cancel_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global settings
settings=[];
delete(gcf)

% --- Executes on button press in current_radio.
function current_radio_Callback(hObject, eventdata, handles)
% hObject    handle to current_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of current_radio

% --- Executes on button press in all_radio.
function all_radio_Callback(hObject, eventdata, handles)
% hObject    handle to all_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of all_radio

% --- Executes on button press in nonnorm_radio.
function nonnorm_radio_Callback(hObject, eventdata, handles)
% hObject    handle to nonnorm_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nonnorm_radio


% --- Executes on button press in struc_radio.
function struc_radio_Callback(hObject, eventdata, handles)
% hObject    handle to struc_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of struc_radio


% --- Executes on button press in custom_radio.
function custom_radio_Callback(hObject, eventdata, handles)
% hObject    handle to custom_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of custom_radio


