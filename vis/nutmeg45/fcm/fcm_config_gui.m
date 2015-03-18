function varargout = fcm_config_gui(varargin)
% FCM_CONFIG_GUI GUI for defining and saving FCM configuration.
% Usage: fcm_config_gui
%
% See also  FCM_START

% Last Modified by GUIDE v2.5 05-Sep-2011 16:08:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fcm_config_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @fcm_config_gui_OutputFcn, ...
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


% --- Executes just before fcm_config_gui is made visible.
function fcm_config_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fcm_config_gui (see VARARGIN)

% Choose default command line output for fcm_config_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fcm_config_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);

global fuse

if ( ~isstruct(fuse) || ~isfield(fuse,'seed') || ~isfield(fuse,'connection') )
    % fuse = struct('seed','All','connection','All','connectionavg','All','output',{{'Mean' 'Z'}});
    defaultfile=[fileparts(which('fcm_config_gui.m')) filesep 'params' filesep 'fcmConfig.ini'];
    fn=textread(defaultfile,'%c','headerlines',1)';
    load([fileparts(which('fcm_config_gui.m')) filesep 'params' filesep fn]);
end

update_gui(handles);
    

% --- Outputs from this function are returned to the command line.
function varargout = fcm_config_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in pop_seed.
function pop_seed_Callback(hObject, eventdata, handles)
% hObject    handle to pop_seed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_seed contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_seed
global fuse
seedC=get(hObject,'String');
fuse.seed=seedC{get(hObject,'Value')};

if strcmp(fuse.seed,'Extracerebral')
    set(handles.pop_conntype,'value',6,'enable','off')
    fuse.connectionavg='None';
else
    set(handles.pop_conntype,'value',1,'enable','on')
    fuse.connectionavg='All';
end
if strcmp(fuse.seed,'All')
    set([handles.check_roi handles.edit_roi handles.browse_roidef],'Enable','on')
else
    set([handles.check_roi handles.edit_roi handles.browse_roidef],'Enable','off')
end

% --- Executes during object creation, after setting all properties.
function pop_seed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_seed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in pop_connection.
function pop_connection_Callback(hObject, eventdata, handles)
% hObject    handle to pop_connection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_connection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_connection
global fuse

connC=get(hObject,'String');
fuse.connection=connC{get(hObject,'Value')};

if ~strcmp(fuse.connection,'All')
    set(handles.pop_conntype,'Value',1,'enable','off')
    fuse.connectionavg='All';
    set([handles.text_grid1 handles.edit_grid handles.text_grid2],'enable','on')    
else
    set(handles.pop_conntype,'enable','on')
    set([handles.text_grid1 handles.edit_grid handles.text_grid2],'enable','off')
end

% --- Executes during object creation, after setting all properties.
function pop_connection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_connection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in pop_conntype.
function pop_conntype_Callback(hObject, eventdata, handles)
% hObject    handle to pop_conntype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_conntype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_conntype
global fuse
ctC={'all' 'interhomo' 'interhetero' 'intra' 'seed' 'none'};
fuse.connectionavg=ctC{get(hObject,'Value')};

% --- Executes during object creation, after setting all properties.
function pop_conntype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_conntype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in pop_datatype.
function pop_datatype_Callback(hObject, eventdata, handles)
% hObject    handle to pop_datatype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_datatype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_datatype
global fuse
switch get(hObject,'Value')
    case 1
        fuse.datatype='1';
    case 2
        fuse.datatype='31';
    case 3
        fuse.datatype='3';
end

% --- Executes during object creation, after setting all properties.
function pop_datatype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_datatype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in pop_funconn.
function pop_funconn_Callback(hObject, eventdata, handles)
% hObject    handle to pop_funconn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_funconn contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_funconn
global fuse
fuco = {'ccohere' 'nccohere' 'glcohere' 'pli' 'ampcorr'};
fuse.funconn = fuco{get(hObject,'Value')};

% --- Executes during object creation, after setting all properties.
function pop_funconn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_funconn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in list_output.
function list_output_Callback(hObject, eventdata, handles)
% hObject    handle to list_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns list_output contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_output
global fuse

outC = {'Mean' 'Z' 'L'}; %get(hObject,'String');
outS = get(hObject,'value')';

if any(outS==2)
    outS = unique([1;outS]);
    set(hObject,'Value',outS)
end
if any(outS==3)
    if isempty(strmatch(fuse.seed,{'All' 'Selected+Contralateral'}))
        errordlg('L-images can only be computed with seed voxels set to ''All'' or ''Selected+Contralateral''');
        outS(outS==3)=[];
        set(hObject,'Value',outS)
        return
    end
end
            
fuse.output = outC(outS);

% --- Executes during object creation, after setting all properties.
function list_output_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in push_close.
function push_close_Callback(hObject, eventdata, handles)
% hObject    handle to push_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(gcf)

% --- Executes on button press in push_cancel.
function push_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to push_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear global fuse
%delete(gcf)

% --- Executes on button press in push_load.
function push_load_Callback(hObject, eventdata, handles)
% hObject    handle to push_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global fuse

temp=pwd;
cd([fileparts(which('fcm_config_gui')) filesep 'params']); 
[ff,fp]=uigetfile('*.mat','Load FCM configuration file...');
if isequal(ff,0), cd(temp), return, end
load(fullfile(fp,ff));
cd(temp)

update_gui(handles);

% --- Executes on button press in push_save.
function push_save_Callback(hObject, eventdata, handles)
% hObject    handle to push_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global fuse

temp=pwd;
cd([fileparts(which('fcm_config_gui')) filesep 'params']); 
[ff,fp]=uiputfile('*.mat','Save FCM configuration as...');
if isequal(ff,0), cd(temp), return, end
save(fullfile(fp,ff),'fuse');
cd(temp);

answer = questdlg('Would you like to make this file the default?','FCM','Yes','No','No');
if strcmp(answer,'Yes')
    defaultfile=[fileparts(which('fcm_config_gui.m')) filesep 'params' filesep 'fcmConfig.ini'];
    delete(defaultfile);
    fid = fopen(defaultfile,'wt');
    fprintf(fid,'[DefaultConfigFilename]\n');
    fprintf(fid,'%s\n',ff);
    fclose(fid);
    fprintf('New default file is %s.\n',ff)
end

%-----------------------------------
function update_gui(handles)

global fuse

try
    set(handles.pop_seed,'Value',strmatch(fuse.seed,get(handles.pop_seed,'string'),'exact'));
    set(handles.pop_connection,'Value',strmatch(fuse.connection,get(handles.pop_connection,'string')));
    if strcmp(fuse.connection,'Grid')
        set(handles.pop_conntype,'Value',1,'enable','off')
        set([handles.text_grid1 handles.edit_grid handles.text_grid2],'enable','on')
        if isfield(fuse,'gridspacing'), set(handles.edit_grid,'string',int2str(fuse.gridspacing)), end
        set(handles.edit_grid,'String',int2str(fuse.gridspacing));
    elseif strcmp(fuse.seed,'Extracerebral')
        set(handles.pop_conntype,'Value',6,'enable','off')
        set([handles.text_grid1 handles.edit_grid handles.text_grid2],'enable','off')
    else
        set(handles.pop_conntype,'Value',strmatch(lower(fuse.connectionavg),{'all' 'interhomo' 'interhetero' 'intra' 'seed' 'none'}),'enable','on');
        set([handles.text_grid1 handles.edit_grid handles.text_grid2],'enable','off')
    end
    if isfield(fuse,'funconn')
        fucoC= {'ccohere' 'nccohere' 'glcohere' 'pli' 'ampcorr'};
        fucoS= find(ismember(fucoC,fuse.funconn));
        set(handles.pop_funconn,'Value',fucoS);
    end
    if isfield(fuse,'datatype')
        dtS = strmatch(fuse.datatype,{'1' '31' '3'},'exact');
        set(handles.pop_datatype,'Value',dtS)
    end
    if strcmp(fuse.seed,'All')
        set(handles.check_roi,'Enable','on')
        if isfield(fuse,'roi') && fuse.roi>0
            set(handles.check_roi,'Value',1);
            set([handles.edit_roi handles.browse_roidef],'Enable','on')
            [dum,fi]=fileparts(fuse.roidef);
            set(handles.edit_roi,'string',fi)
        else
            set(handles.check_roi,'Value',0);
            set([handles.edit_roi handles.browse_roidef],'Enable','off')
        end
    else
        set(handles.check_roi,'Value',0);
        set([handles.check_roi handles.edit_roi handles.browse_roidef],'Enable','off')
    end      
        
    outC = {'Mean' 'Z' 'L'};
    outS = find(ismember(outC,fuse.output));
    set(handles.list_output,'Value',outS);
        
catch
    errordlg('Invalid FCM configuration file.')
end

function edit_grid_Callback(hObject, eventdata, handles)
% hObject    handle to edit_grid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_grid as text
%        str2double(get(hObject,'String')) returns contents of edit_grid as a double

global fuse
fuse.gridspacing = str2double(get(hObject,'string'));

% --- Executes during object creation, after setting all properties.
function edit_grid_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_grid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_roi_Callback(hObject, eventdata, handles)
% hObject    handle to edit_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_roi as text
%        str2double(get(hObject,'String')) returns contents of edit_roi as a double
global fuse

ipa= get(hObject,'string');
[pa,fi]=fileparts(ipa);
if isempty(pa) 
    if isfield(fuse,'roidef') && ~isempty(strfind(fuse.roidef,fi))
        pa = fileparts(fuse.roidef);
    else
        pa=[fileparts(which('fcm_gui')) filesep 'templates'];
    end
end
fu = fullfile(pa,fi);
load(fu);
if ~exist('ROI','var') || ~isfield(ROI,'voxels')
    errordlg('Invalid ROI definition file.'), return
end

fuse.roi = length(ROI.label);
fuse.roidef = fu;

% --- Executes during object creation, after setting all properties.
function edit_roi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in check_roi.
function check_roi_Callback(hObject, eventdata, handles)
% hObject    handle to check_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global fuse

if get(hObject,'Value')
    ipa= get(handles.edit_roi,'string');
    [pa,fi]=fileparts(ipa);
    if isempty(pa) 
        if isfield(fuse,'roidef') && ~isempty(strfind(fuse.roidef,fi))
            pa = fileparts(fuse.roidef);
        else
            pa=[fileparts(which('fcm_gui')) filesep 'templates'];
        end
    end
    fu = fullfile(pa,fi);
    load(fu);
    if ~exist('ROI','var') || ~isfield(ROI,'voxels')
        errordlg('Invalid ROI definition file.'), return
    end

    fuse.roi = length(ROI.label);
    fuse.roidef = fu;
    
    set([handles.edit_roi handles.browse_roidef],'enable','on')
else
    fuse.roi=0;
    fuse.roidef='';
    
    set([handles.edit_roi handles.browse_roidef],'enable','off')
end

% Hint: get(hObject,'Value') returns toggle state of check_roi
% --- Executes on button press in browse_roidef.
function browse_roidef_Callback(hObject, eventdata, handles)
% hObject    handle to browse_roidef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global fuse

depa = [fileparts(which('fcm_gui')) filesep 'templates' filesep '*.mat'];
[fi,pa]=uigetfile(depa,'Select ROI definition file...');
if isequal(fi,0), return, end
fu = fullfile(pa,fi);
load(fu);
if ~exist('ROI','var') || ~isfield(ROI,'voxels')
    errordlg('Invalid ROI definition file.'), return
end

fuse.roi = length(ROI.label);
fuse.roidef = fu;
set(handles.edit_roi,'string',fi)
