function varargout = nut_simulation(varargin)
% NUT_SIMULATION M-file for nut_simulation.fig
%      NUT_SIMULATION, by itself, creates a new NUT_SIMULATION or raises the existing
%      singleton*.
%
%      H = NUT_SIMULATION returns the handle to a new NUT_SIMULATION or the handle to
%      the existing singleton*.
%
%      NUT_SIMULATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NUT_SIMULATION.M with the given input arguments.
%
%      NUT_SIMULATION('Property','Value',...) creates a new NUT_SIMULATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before nut_simulation_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to nut_simulation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help nut_simulation

% Last Modified by GUIDE v2.5 12-Mar-2006 21:25:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @nut_simulation_OpeningFcn, ...
                   'gui_OutputFcn',  @nut_simulation_OutputFcn, ...
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


% --- Executes just before nut_simulation is made visible.
function nut_simulation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to nut_simulation (see VARARGIN)

% Choose default command line output for nut_simulation
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes nut_simulation wait for user response (see UIRESUME)
% uiwait(handles.figure1);

global sim nuts
if(isempty(varargin))
sim.num_gauss_sources=str2num(get(handles.nut_num_gauss_sources_edit,'String'));
sim.num_gaussdamp_sin=str2num(get(handles.nut_num_gaussdamp_sin_edit,'String'));
sim.num_sin_source=str2num(get(handles.nut_num_sin_source_edit,'String'));
sim.num_timepoints=str2num(get(handles.nut_num_timepoints_edit,'String'));
sim.percent_stim=str2num(get(handles.nut_percent_stim_edit,'String'));
sim.freqrange=str2num(get(handles.nut_freqrange_edit,'String'));
sim.sinphase=str2num(get(handles.nut_sinphase_edit,'String'));
sim.sindamp=str2num(get(handles.nut_sindamp_edit,'String'));
sim.sourceloc=str2num(get(handles.nut_sourceloc_edit,'String'));
sim.source_separation=str2num(get(handles.nut_source_separation_edit,'String'));
sim.source_orientation=str2num(get(handles.nut_source_orientation_edit,'String'));
sim.source_corr=str2num(get(handles.nut_source_corr_edit,'String'));
sim.num_inter=str2num(get(handles.nut_num_inter_edit,'String'));
sim.num_sensors=str2num(get(handles.nut_num_sensors_edit,'String'));
sim.SNR=str2num(get(handles.nut_SNR_edit,'String'));
sim.SNIR=str2num(get(handles.nut_SNIR_edit,'String'));
sim.seed=str2num(get(handles.nut_seed_edit,'String'));

sim.gauss_rand_source=get(handles.nut_gauss_rand_source_popupmenu,'Value');
sim.randfreq=get(handles.nut_randfreq_popupmenu,'Value');
sim.sinphase=get(handles.nut_sinphase_popupmenu,'Value');
sim.damploc=get(handles.nut_damploc_popupmenu,'Value');
sim.specify_sourceloc=get(handles.nut_specify_sourceloc_popupmenu,'Value');
sim.randoriented=get(handles.nut_randoriented_popupmenu,'Value');
sim.intertype=get(handles.nut_intertype_popupmenu,'Value');
sim.randinterloc=get(handles.nut_randinterloc_popupmenu,'Value');
sim.randinterorient=get(handles.nut_randinterorient_popupmenu,'Value');
sim.gauss_rand_inter=get(handles.nut_gauss_rand_inter_popupmenu,'Value');
sim.gauss_rand_noise=get(handles.nut_gauss_rand_noise_popupmenu,'Value');
sim.sensor_select=get(handles.nut_sensor_select_popupmenu,'Value');

sim.inter_rotate=get(handles.nut_inter_rotate_checkbox,'Value');
sim.source_rotate=get(handles.nut_source_rotate_checkbox,'Value');
else
    simpath=varargin{1};
    load(simpath);
end

sim.srate=nuts.meg.srate;
sim.handles=handles;


% --- Outputs from this function are returned to the command line.
function varargout = nut_simulation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

global sim nuts
% [sources,t]=nut_create_sim_sources(sim);
% nuts.meg.latency=t;
% [data]=nut_sources2sensors(sources,nuts.Lp,sim);
% nuts.meg.data=data;

function nut_num_gauss_sources_edit_Callback(hObject, eventdata, handles)
% hObject    handle to nut_num_gauss_sources_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_num_gauss_sources_edit as text
%        str2double(get(hObject,'String')) returns contents of nut_num_gauss_sources_edit as a double


% --- Executes during object creation, after setting all properties.
function nut_num_gauss_sources_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_num_gauss_sources_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nut_num_timepoints_edit_Callback(hObject, eventdata, handles)
% hObject    handle to nut_num_timepoints_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_num_timepoints_edit as text
%        str2double(get(hObject,'String')) returns contents of nut_num_timepoints_edit as a double


% --- Executes during object creation, after setting all properties.
function nut_num_timepoints_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_num_timepoints_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nut_num_sensors_edit_Callback(hObject, eventdata, handles)
% hObject    handle to nut_num_sensors_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_num_sensors_edit as text
%        str2double(get(hObject,'String')) returns contents of nut_num_sensors_edit as a double


% --- Executes during object creation, after setting all properties.
function nut_num_sensors_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_num_sensors_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in nut_specify_sourceloc_togglebutton.
function nut_specify_sourceloc_togglebutton_Callback(hObject, eventdata, handles)
% hObject    handle to nut_specify_sourceloc_togglebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nut_specify_sourceloc_togglebutton
aa=get(hObject,'Value')



function nut_SNR_edit_Callback(hObject, eventdata, handles)
% hObject    handle to nut_SNR_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_SNR_edit as text
%        str2double(get(hObject,'String')) returns contents of nut_SNR_edit as a double


% --- Executes during object creation, after setting all properties.
function nut_SNR_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_SNR_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nut_sourceloc_edit_Callback(hObject, eventdata, handles)
% hObject    handle to nut_sourceloc_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_sourceloc_edit as text
%        str2double(get(hObject,'String')) returns contents of nut_sourceloc_edit as a double


% --- Executes during object creation, after setting all properties.
function nut_sourceloc_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_sourceloc_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in nut_sourcetype_popupmenu.
function nut_sourcetype_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to nut_sourcetype_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns nut_sourcetype_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nut_sourcetype_popupmenu


% --- Executes during object creation, after setting all properties.
function nut_sourcetype_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_sourcetype_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in nut_intertype_popupmenu.
function nut_intertype_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to nut_intertype_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns nut_intertype_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nut_intertype_popupmenu


% --- Executes during object creation, after setting all properties.
function nut_intertype_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_intertype_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nut_source_separation_edit_Callback(hObject, eventdata, handles)
% hObject    handle to nut_source_separation_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_source_separation_edit as text
%        str2double(get(hObject,'String')) returns contents of nut_source_separation_edit as a double


% --- Executes during object creation, after setting all properties.
function nut_source_separation_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_source_separation_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nut_source_corr_edit_Callback(hObject, eventdata, handles)
% hObject    handle to nut_source_corr_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_source_corr_edit as text
%        str2double(get(hObject,'String')) returns contents of nut_source_corr_edit as a double


% --- Executes during object creation, after setting all properties.
function nut_source_corr_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_source_corr_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nut_SNIR_edit_Callback(hObject, eventdata, handles)
% hObject    handle to nut_SNIR_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_SNIR_edit as text
%        str2double(get(hObject,'String')) returns contents of nut_SNIR_edit as a double


% --- Executes during object creation, after setting all properties.
function nut_SNIR_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_SNIR_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nut_num_inter_edit_Callback(hObject, eventdata, handles)
% hObject    handle to nut_num_inter_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_num_inter_edit as text
%        str2double(get(hObject,'String')) returns contents of nut_num_inter_edit as a double


% --- Executes during object creation, after setting all properties.
function nut_num_inter_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_num_inter_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in nut_source_rotate_checkbox.
function nut_source_rotate_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to nut_source_rotate_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nut_source_rotate_checkbox


% --- Executes on button press in nut_random_sensors_checkbox.
function nut_random_sensors_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to nut_random_sensors_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nut_random_sensors_checkbox



function nut_source_orientation_edit_Callback(hObject, eventdata, handles)
% hObject    handle to nut_source_orientation_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_source_orientation_edit as text
%        str2double(get(hObject,'String')) returns contents of nut_source_orientation_edit as a double


% --- Executes during object creation, after setting all properties.
function nut_source_orientation_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_source_orientation_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in nut_sensor_select_popupmenu.
function nut_sensor_select_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to nut_sensor_select_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns nut_sensor_select_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nut_sensor_select_popupmenu


% --- Executes during object creation, after setting all properties.
function nut_sensor_select_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_sensor_select_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in nut_specify_sourceloc_popupmenu.
function nut_specify_sourceloc_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to nut_specify_sourceloc_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns nut_specify_sourceloc_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nut_specify_sourceloc_popupmenu


% --- Executes during object creation, after setting all properties.
function nut_specify_sourceloc_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_specify_sourceloc_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nut_seed_edit_Callback(hObject, eventdata, handles)
% hObject    handle to nut_seed_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_seed_edit as text
%        str2double(get(hObject,'String')) returns contents of nut_seed_edit as a double


% --- Executes during object creation, after setting all properties.
function nut_seed_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_seed_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in nut_randoriented_popupmenu.
function nut_randoriented_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to nut_randoriented_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns nut_randoriented_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nut_randoriented_popupmenu


% --- Executes during object creation, after setting all properties.
function nut_randoriented_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_randoriented_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in nut_randinterloc_popupmenu.
function nut_randinterloc_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to nut_randinterloc_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns nut_randinterloc_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nut_randinterloc_popupmenu


% --- Executes during object creation, after setting all properties.
function nut_randinterloc_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_randinterloc_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in nut_randinterorient_popupmenu.
function nut_randinterorient_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to nut_randinterorient_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns nut_randinterorient_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nut_randinterorient_popupmenu


% --- Executes during object creation, after setting all properties.
function nut_randinterorient_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_randinterorient_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nut_num_gaussdamp_sin_edit_Callback(hObject, eventdata, handles)
% hObject    handle to nut_num_gaussdamp_sin_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_num_gaussdamp_sin_edit as text
%        str2double(get(hObject,'String')) returns contents of nut_num_gaussdamp_sin_edit as a double


% --- Executes during object creation, after setting all properties.
function nut_num_gaussdamp_sin_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_num_gaussdamp_sin_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nut_num_sin_source_edit_Callback(hObject, eventdata, handles)
% hObject    handle to nut_num_sin_source_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_num_sin_source_edit as text
%        str2double(get(hObject,'String')) returns contents of nut_num_sin_source_edit as a double


% --- Executes during object creation, after setting all properties.
function nut_num_sin_source_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_num_sin_source_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in nut_inter_rotate_checkbox.
function nut_inter_rotate_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to nut_inter_rotate_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nut_inter_rotate_checkbox


% --- Executes on selection change in nut_gauss_rand_source_popupmenu.
function nut_gauss_rand_source_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to nut_gauss_rand_source_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns nut_gauss_rand_source_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nut_gauss_rand_source_popupmenu


% --- Executes during object creation, after setting all properties.
function nut_gauss_rand_source_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_gauss_rand_source_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in nut_gauss_rand_inter_popupmenu.
function nut_gauss_rand_inter_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to nut_gauss_rand_inter_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns nut_gauss_rand_inter_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nut_gauss_rand_inter_popupmenu


% --- Executes during object creation, after setting all properties.
function nut_gauss_rand_inter_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_gauss_rand_inter_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in nut_gauss_rand_noise_popupmenu.
function nut_gauss_rand_noise_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to nut_gauss_rand_noise_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns nut_gauss_rand_noise_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nut_gauss_rand_noise_popupmenu


% --- Executes during object creation, after setting all properties.
function nut_gauss_rand_noise_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_gauss_rand_noise_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nut_percent_stim_edit_Callback(hObject, eventdata, handles)
% hObject    handle to nut_percent_stim_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_percent_stim_edit as text
%        str2double(get(hObject,'String')) returns contents of nut_percent_stim_edit as a double


% --- Executes during object creation, after setting all properties.
function nut_percent_stim_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_percent_stim_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nut_freqrange_edit_Callback(hObject, eventdata, handles)
% hObject    handle to nut_freqrange_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_freqrange_edit as text
%        str2double(get(hObject,'String')) returns contents of nut_freqrange_edit as a double


% --- Executes during object creation, after setting all properties.
function nut_freqrange_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_freqrange_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nut_sinphase_edit_Callback(hObject, eventdata, handles)
% hObject    handle to nut_sinphase_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_sinphase_edit as text
%        str2double(get(hObject,'String')) returns contents of nut_sinphase_edit as a double


% --- Executes during object creation, after setting all properties.
function nut_sinphase_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_sinphase_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in nut_randfreq_popupmenu.
function nut_randfreq_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to nut_randfreq_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns nut_randfreq_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nut_randfreq_popupmenu


% --- Executes during object creation, after setting all properties.
function nut_randfreq_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_randfreq_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in nut_sinphase_popupmenu.
function nut_sinphase_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to nut_sinphase_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns nut_sinphase_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nut_sinphase_popupmenu


% --- Executes during object creation, after setting all properties.
function nut_sinphase_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_sinphase_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nut_sindamp_edit_Callback(hObject, eventdata, handles)
% hObject    handle to nut_sindamp_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_sindamp_edit as text
%        str2double(get(hObject,'String')) returns contents of nut_sindamp_edit as a double


% --- Executes during object creation, after setting all properties.
function nut_sindamp_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_sindamp_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in nut_damploc_popupmenu.
function nut_damploc_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to nut_damploc_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns nut_damploc_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nut_damploc_popupmenu


% --- Executes during object creation, after setting all properties.
function nut_damploc_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_damploc_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in nut_done_button.
function nut_done_button_Callback(hObject, eventdata, handles)
% hObject    handle to nut_done_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in nut_save_button.
function nut_save_button_Callback(hObject, eventdata, handles)
% hObject    handle to nut_save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in nut_load_button.
function nut_load_button_Callback(hObject, eventdata, handles)
% hObject    handle to nut_load_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


