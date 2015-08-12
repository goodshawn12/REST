function varargout = selStream(varargin)
% SELSTREAM MATLAB code for selStream.fig
%      SELSTREAM, by itself, creates a new SELSTREAM or raises the existing
%      singleton*.
%
%      H = SELSTREAM returns the handle to a new SELSTREAM or the handle to
%      the existing singleton*.
%
%      SELSTREAM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELSTREAM.M with the given input arguments.
%
%      SELSTREAM('Property','Value',...) creates a new SELSTREAM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before selStream_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to selStream_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help selStream

% Last Modified by GUIDE v2.5 28-Jan-2015 13:19:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @selStream_OpeningFcn, ...
                   'gui_OutputFcn',  @selStream_OutputFcn, ...
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


% --- Executes just before selStream is made visible.
function selStream_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to selStream (see VARARGIN)

% Choose default command line output for selStream
handles.output = hObject;
if length(varargin) == 1
    handles.streamnames = varargin{1};
else
    errordlg('Input correct streamnames sturcture.');
end

nStream = length(handles.streamnames);
menuOptions = [];
for it=1:nStream
   menuOptions = [menuOptions;handles.streamnames(it)];
end
set(handles.popupmenu1,'String',menuOptions);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes selStream wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = selStream_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

% items = get(hObject,'String');
% idx_sel = get(hObject,'Value');
% item_sel = items(idx_sel);
% 
% assignin('base','streamname',item_sel);


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

items = get(handles.popupmenu1,'String');
idx_sel = get(handles.popupmenu1,'Value');
item_sel = items(idx_sel);
assignin('base','streamname',item_sel);

% Hint: delete(hObject) closes the figure
if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
%     guidata(hObject,handles);
end
delete(hObject);

