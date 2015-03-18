function varargout = nut_connectivity_gui(varargin)
% NUT_CONNECTIVITY_GUI M-file for nut_connectivity_gui.fig
%      NUT_CONNECTIVITY_GUI, by itself, creates a new NUT_CONNECTIVITY_GUI or raises the existing
%      singleton*.
%
%      H = NUT_CONNECTIVITY_GUI returns the handle to a new NUT_CONNECTIVITY_GUI or the handle to
%      the existing singleton*.
%
%      NUT_CONNECTIVITY_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NUT_CONNECTIVITY_GUI.M with the given input arguments.
%
%      NUT_CONNECTIVITY_GUI('Property','Value',...) creates a new NUT_CONNECTIVITY_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before nut_connectivity_gui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to nut_connectivity_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help nut_connectivity_gui

% Last Modified by GUIDE v2.5 17-Nov-2010 20:27:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @nut_connectivity_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @nut_connectivity_gui_OutputFcn, ...
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


% --- Executes just before nut_connectivity_gui is made visible.
function nut_connectivity_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to nut_connectivity_gui (see VARARGIN)

% Choose default command line output for nut_connectivity_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

if size(varargin)>0
    global conset
    conset=varargin{1};
    set(handles.do_hec,'Value',conset.do.hec);
    set(handles.do_shec,'Value',conset.do.shec);
    set(handles.do_aCoh,'Value',conset.do.coh);
    set(handles.do_iCoh,'Value',conset.do.icoh);
    
    
    if conset.seed.type==1
        set(handles.seed_voxel,'Value',1);
        set(handles.edit_voxel,'String',num2str(conset.seed.index));
    elseif conset.seed.type==2
        set(handles.seed_mni,'Value',1)
        set(handles.edit_mni,'String',num2str(conset.seed.index));
    elseif conset.seed.type==3
        set(handles.seed_all,'Value',1)
    else
        warning('conset.seed.type not valid')
        set(handles.seed_voxel,'Value',1);
    end
        
    set(handles.twindows,'String',[ '[' num2str(conset.windows(1,:)) ']' ]);

end


% UIWAIT makes nut_connectivity_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = nut_connectivity_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

uiwait
global conset
varargout{1} = conset;
if ~isempty(conset)
    disp('running connectivity');
    nut_sourceConnectivity(conset);
end
clear global conset;


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

global conset

conset.do.hec=get(handles.do_hec,'Value');
conset.do.shec=get(handles.do_shec,'Value');
conset.do.coh=get(handles.do_aCoh,'Value');
conset.do.icoh=get(handles.do_iCoh,'Value');

tmp=get(handles.edit_inputfile,'String')
if ~strcmp(tmp(end-3:end),'.mat')
    push_input_Callback(hObject, eventdata, handles)
end
conset.inputfile=get(handles.edit_inputfile,'String');

if get(handles.radio_dofilter,'Value')
    conset.filter=1;
elseif get(handles.radio_isfiltered,'Value')
    conset.filter=0;
else
    return
end
conset.filterband=str2num(get(handles.edit_filter,'String'));


if get(handles.seed_voxel,'Value')
    conset.seed.type=1;
    conset.seed.index=eval(get(handles.edit_voxel,'String'));
elseif get(handles.seed_mni,'Value')
    conset.seed.type=2;
    conset.seed.index=eval(get(handles.edit_mni,'String'));
elseif get(handles.seed_all,'Value')
    conset.seed.type=3;
else
    return
end

if get(handles.edit_window_input,'Value')
    conset.windows=eval(get(handles.twindows,'String'));
elseif get(handles.input_window_file)
    try
        conset.windows=load(get(handles.edit_window_file,'String'));
    catch
        windowfile=uigetfile('*.mat','Select Window file');
        conset.windows=load(windowfile);
    end
else
    return
end

conset.chunktime=str2num(get(handles.edit_chunk,'String'));
conset.schunktime=str2num(get(handles.edit_schunk,'String'));

if strcmp(get(handles.edit_trials,'String'),'all')
    conset.trials=0;
elseif isnumeric(str2num(get(handles.edit_trials,'String')))
    conset.trials=str2num(get(handles.edit_trials,'String'));
else
    warning('input for trials to use is wrong format');
end


delete(gcf)

% --- Executes on button press in cancel_push.
function cancel_push_Callback(hObject, eventdata, handles)
% hObject    handle to cancel_push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global conset
conset=[];
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


% --- Executes on button press in do_hec.
function do_hec_Callback(hObject, eventdata, handles)
% hObject    handle to do_hec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of do_hec


% --- Executes on button press in do_shec.
function do_shec_Callback(hObject, eventdata, handles)
% hObject    handle to do_shec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of do_shec


% --- Executes on button press in do_aCoh.
function do_aCoh_Callback(hObject, eventdata, handles)
% hObject    handle to do_aCoh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of do_aCoh


% --- Executes on button press in do_iCoh.
function do_iCoh_Callback(hObject, eventdata, handles)
% hObject    handle to do_iCoh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of do_iCoh



function twindows_Callback(hObject, eventdata, handles)
% hObject    handle to twindows (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of twindows as text
%        str2double(get(hObject,'String')) returns contents of twindows as a double


% --- Executes during object creation, after setting all properties.
function twindows_CreateFcn(hObject, eventdata, handles)
% hObject    handle to twindows (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_voxel_Callback(hObject, eventdata, handles)
% hObject    handle to edit_voxel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_voxel as text
%        str2double(get(hObject,'String')) returns contents of edit_voxel as a double


% --- Executes during object creation, after setting all properties.
function edit_voxel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_voxel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_mni_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mni (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mni as text
%        str2double(get(hObject,'String')) returns contents of edit_mni as a double


% --- Executes during object creation, after setting all properties.
function edit_mni_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mni (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_window_file_Callback(hObject, eventdata, handles)
% hObject    handle to edit_window_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_window_file as text
%        str2double(get(hObject,'String')) returns contents of edit_window_file as a double


% --- Executes during object creation, after setting all properties.
function edit_window_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_window_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_window_file.
function push_window_file_Callback(hObject, eventdata, handles)
% hObject    handle to push_window_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[windowfile,windowpath]=uigetfile('*.mat','Select Window file');
set(handles.edit_window_file,'String',[windowpath windowfile])


function edit_chunk_Callback(hObject, eventdata, handles)
% hObject    handle to edit_chunk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_chunk as text
%        str2double(get(hObject,'String')) returns contents of edit_chunk as a double


% --- Executes during object creation, after setting all properties.
function edit_chunk_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_chunk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_schunk_Callback(hObject, eventdata, handles)
% hObject    handle to edit_schunk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_schunk as text
%        str2double(get(hObject,'String')) returns contents of edit_schunk as a double


% --- Executes during object creation, after setting all properties.
function edit_schunk_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_schunk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_input.
function push_input_Callback(hObject, eventdata, handles)
% hObject    handle to push_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global conset
conset.freq=str2num(get(handles.edit_filter,'String'));
freqstring=[num2str(conset.freq(1)) 'to' num2str(conset.freq(2))];

[inputfile,inputpath]=uigetfile(['weights_*' freqstring '*.mat'],'Select Window file');
conset.inputfile=[inputpath inputfile];
set(handles.edit_inputfile,'String',conset.inputfile)

if ~strcmp(freqstring,conset.inputfile(strfind(conset.inputfile,'Hz')-length(freqstring):strfind(conset.inputfile,'Hz')-1))
    error('you must choose a weights file corresponding to the filter settings')
end
load(conset.inputfile);
if ~exist('W','var')
    error('wrong input file type. enter weights file')
end
conset.W=W;


function edit_inputfile_Callback(hObject, eventdata, handles)
% hObject    handle to edit_inputfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_inputfile as text
%        str2double(get(hObject,'String')) returns contents of edit_inputfile as a double


% --- Executes during object creation, after setting all properties.
function edit_inputfile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_inputfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_trials_Callback(hObject, eventdata, handles)
% hObject    handle to edit_trials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_trials as text
%        str2double(get(hObject,'String')) returns contents of edit_trials as a double


% --- Executes during object creation, after setting all properties.
function edit_trials_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_trials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_filter_Callback(hObject, eventdata, handles)
% hObject    handle to edit_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_filter as text
%        str2double(get(hObject,'String')) returns contents of edit_filter as a double


% --- Executes during object creation, after setting all properties.
function edit_filter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
