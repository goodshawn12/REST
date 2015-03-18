function varargout = nut_select_beamratio(varargin)
% NUT_SELECT_BEAMRATIO M-file for nut_select_beamratio.fig
%      NUT_SELECT_BEAMRATIO, by itself, creates a new NUT_SELECT_BEAMRATIO or raises the existing
%      singleton*.
%
%      H = NUT_SELECT_BEAMRATIO returns the handle to a new NUT_SELECT_BEAMRATIO or the handle to
%      the existing singleton*.
%
%      NUT_SELECT_BEAMRATIO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NUT_SELECT_BEAMRATIO.M with the given input arguments.
%
%      NUT_SELECT_BEAMRATIO('Property','Value',...) creates a new NUT_SELECT_BEAMRATIO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before nut_select_ratio4stats_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to nut_select_beamratio_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help nut_select_beamratio

% Last Modified by GUIDE v2.5 10-Nov-2010 13:31:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @nut_select_beamratio_OpeningFcn, ...
                   'gui_OutputFcn',  @nut_select_beamratio_OutputFcn, ...
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


% --- Executes just before nut_select_beamratio is made visible.
function nut_select_beamratio_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to nut_select_beamratio (see VARARGIN)

% Choose default command line output for nut_select_beamratio
handles.output = hObject;
%handles.output2= []; 
%set(handles.ok_pushbutton,'UserData',0);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes nut_select_beamratio wait for user response (see UIRESUME)
% uiwait(handles.nut_select_beamratio);


% --- Outputs from this function are returned to the command line.
function varargout = nut_select_beamratio_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on selection change in ratio_menu.
function ratio_menu_Callback(hObject, eventdata, handles)
% hObject    handle to ratio_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns ratio_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ratio_menu

% --- Executes during object creation, after setting all properties.
function ratio_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ratio_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in noise_checkbox.
function noise_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to noise_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of noise_checkbox

% --- Executes on button press in ok_pushbutton.
function ok_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to ok_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global stats
stats.ratio=get(handles.ratio_menu,'Value');
stats.noise=get(handles.noise_checkbox,'Value');
delete(handles.output)

% --- Executes when user attempts to close nut_select_beamratio.
function nut_select_ratio4stats_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to nut_select_beamratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);

% --- Executes during object creation, after setting all properties.
function noise_checkbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noise_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function ok_pushbutton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ok_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
