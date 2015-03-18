function varargout = fcm_across_subject_gui(varargin)
% FCM_ACROSS_SUBJECT_GUI  GUI for FCM statistical analyses across subjects.
%      FCM_ACROSS_SUBJECT_GUI, by itself, creates a new FCM_ACROSS_SUBJECT_GUI or raises the existing
%      singleton*.
%
%      H = FCM_ACROSS_SUBJECT_GUI returns the handle to a new FCM_ACROSS_SUBJECT_GUI or the handle to
%      the existing singleton*.
%
%      FCM_ACROSS_SUBJECT_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FCM_ACROSS_SUBJECT_GUI.M with the given input arguments.
%
%      FCM_ACROSS_SUBJECT_GUI('Property','Value',...) creates a new FCM_ACROSS_SUBJECT_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fcm_across_subject_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fcm_across_subject_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fcm_across_subject_gui

% Last Modified by GUIDE v2.5 10-Mar-2011 18:10:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fcm_across_subject_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @fcm_across_subject_gui_OutputFcn, ...
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


% --- Executes just before fcm_across_subject_gui is made visible.
function fcm_across_subject_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fcm_across_subject_gui (see VARARGIN)

% Choose default command line output for fcm_across_subject_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global conns
if isfield(conns,'files') && isfield(conns.files,'beam')
    update_beamlist(handles,fullpath2file(conns.files.beam));
else
    beam.files.beam = [];
    update_beamlist(handles);
end
conns.files.behavioral = 'workspace';

% UIWAIT makes fcm_across_subject_gui wait for user response (see UIRESUME)
% uiwait(handles.fcm_across_subject_gui);

% --- Outputs from this function are returned to the command line.
function varargout = fcm_across_subject_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function update_beamlist(handles,preselect)

dir_struct = dir('s_beam*.mat');
sorted_names = sortrows({dir_struct.name}');
set(handles.list_beam,'String',sorted_names)
if nargin>1
    selidx = strmatch(preselect,sorted_names);
    if ~isempty(selidx)
        set(handles.list_beam,'Value',selidx(1))
    else
        set(handles.list_beam,'Value',1)
    end
elseif isempty(get(handles.list_beam,'Value')) && ~isempty(sorted_names)
    set(handles.list_beam,'Value',1)    
elseif get(handles.list_beam,'Value')>length(sorted_names)
    set(handles.list_beam,'Value',1)
end

% --- Executes on button press in push_refresh.
function push_refresh_Callback(hObject, eventdata, handles)
% hObject    handle to push_refresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_beamlist(handles);

% --- Executes on button press in run_pop.
function run_pop_Callback(hObject, eventdata, handles)
% hObject    handle to run_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[voifile,voipath]=uigetfile('*.mat','Load Pointer File...');
if isequal(voifile,0), return, end
load([voipath voifile]);

% Prepare stats structure, for other settings you must you nut_timef_stats
stats.docomp=false;
stats.subjsel=unique(voi.subjnr);
stats.numsubj=length(stats.subjsel);
stats.condsel=unique(voi.condnr);
stats.numtests=length(stats.condsel);
stats.comps=stats.condsel';
stats.subjnr=repmat(stats.subjsel,[stats.numtests 1]);
stats.condnr=repmat(stats.comps,[1 stats.numsubj]);
stats.files=cell(max(stats.subjsel),max(stats.condsel));
for c=stats.condsel
    for s=stats.subjsel
        if any(any(stats.subjnr==s & stats.condnr==c))
            v=find(voi.subjnr==s & voi.condnr==c);
            if ~isempty(v)
                [fpath,ffile,fext]=fileparts(voi.pathnames{v});
                if ~isempty(cell2mat(regexp(ffile,{'_spatnorm','_subj','_fcmri'})))     % if files are already spatnorm
                    suffix='.mat';
                else                            % usually, _spatnorm has to be added
                    suffix='_spatnorm.mat';
                end
                stats.files{s,c}=fullfile(fpath,[ffile suffix]);
                if ~exist(stats.files{s,c},'file'), errordlg(sprintf('Could not find file %s.',stats.files{s,c})), return, end
            else
                errordlg(sprintf('Could not find condition %d of subject %d.',c,s)), return
            end
        end
    end
end
stats.frqband=[];
stats.timepts=[];
stats.fsel='one';
stats.tsel='one';
stats.frqstats=false;
stats.tsstats=false;

nut_beampopulation(stats);

% ------------
function getbeamfile(handles)
global conns
beamlist = get(handles.list_beam,'string');
if ~isempty(beamlist)
    conns.files.beam = beamlist{get(handles.list_beam,'value')};
else
    conns.files.beam = [];
end

% --- Executes on button press in run_viewer.
function run_viewer_Callback(hObject, eventdata, handles)
% hObject    handle to run_viewer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global conns
getbeamfile(handles);
if ~isempty(conns.files.beam)
    tv(conns.files.beam)
end

% --- Executes on button press in browse_controlpop.
function browse_controlpop_Callback(hObject, eventdata, handles)
% hObject    handle to browse_controlpop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global conns
[fi,pa]=uigetfile('*.mat','Select control population file...');
if isequal(fi,0), return, end
fu=fullfile(pa,fi);
conns.files.controlpop = fu;
set(handles.edit_controlpop,'string',fi)

% --- Executes on button press in menu_controlpop.
function menu_controlpop_Callback(hObject, eventdata, handles)
% hObject    handle to menu_controlpop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global conns

load([fileparts(which('fcm.m')) filesep 'params' filesep 'controlpopulations.mat']);   
if isempty(controls)
    errordlg('No control populations were added to the menu yet. Click on "Get population" first.')
    return
end
menutext=cat(2, '''Which control population would you like to use?'',' , ...
                sprintf('''%s'', ',controls.description) , ...
                '''Cancel''' );
eval(['sel=menu(' menutext ');']);
if sel<=length(controls)
    conns.files.controlpop = controls(sel).path;
    set(handles.edit_controlpop,'string',fullpath2file(conns.files.controlpop))
end
     
% --- Executes on button press in run_pimage.
function run_pimage_Callback(hObject, eventdata, handles)
% hObject    handle to run_pimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global conns
getbeamfile(handles);
if isempty(conns.files.beam)
    errordlg('No s_beam file selected.'), return
end
if ~isfield(conns.files,'controlpop')
    errordlg('No control population specified.'), return
end

beam = load(conns.files.beam);
if isfield(beam,'beam'), beam=beam.beam; end
load(conns.files.controlpop);

beam = fcm_beam2Pimage(beam,pop);

[pa,fi]=fileparts(conns.files.beam);
newfile = strrep(fi,'_spatnorm','P');
conns.files.Pimage = fullfile(pa,newfile);
save(conns.files.Pimage,'beam')

msgbox(['P-image saved as ' newfile]);
update_beamlist(handles,newfile);

% --- Executes on button press in browse_clinical.
function browse_clinical_Callback(hObject, eventdata, handles)
% hObject    handle to browse_clinical (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global conns
[fi,pa]=uigetfile('*.mat','Select behavioral data file...');
if isequal(fi,0), return, end
fu=fullfile(pa,fi);
conns.files.behavioral = fu;
set(handles.edit_clinicalsource,'string',fi)

% --- Executes on button press in push_workspace.
function push_workspace_Callback(hObject, eventdata, handles)
% hObject    handle to push_workspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global conns
conns.files.behavioral = 'workspace';
set(handles.edit_clinicalsource,'string','MATLAB Workspace')

% --- Executes on button press in browse_activepop.
function browse_activepop_Callback(hObject, eventdata, handles)
% hObject    handle to browse_activepop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global conns
[fi,pa]=uigetfile('*.mat','Select active/patient population file...');
if isequal(fi,0), return, end
fu=fullfile(pa,fi);
conns.files.activepop = fu;
set(handles.edit_activepop,'string',fi)

% --- Executes on button press in menu_activepop.
function menu_activepop_Callback(hObject, eventdata, handles)
% hObject    handle to menu_activepop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global conns

load([fileparts(which('fcm.m')) filesep 'params' filesep 'activepopulations.mat']);   
if isempty(actives)
    errordlg('No active populations were added to the menu yet. Click on "Get population" first.')
    return
end
menutext=cat(2, '''Which active/patient population would you like to use?'',' , ...
                sprintf('''%s'', ',actives.description) , ...
                '''Cancel''' );
eval(['sel=menu(' menutext ');']);
if sel<=length(actives)
    conns.files.activepop = actives(sel).path;
    set(handles.edit_activepop,'string',fullpath2file(conns.files.activepop))
end
                
% --- Executes on button press in run_corr.
function run_corr_Callback(hObject, eventdata, handles)
% hObject    handle to run_corr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global conns

if ~isfield(conns.files,'activepop')
    errordlg('No active/patient population specified.'), return
end

try
    if strcmp(conns.files.behavioral,'workspace')
        B = evalin('base',get(handles.edit_clinicalvar,'string'));
    else
        load(conns.files.behavioral)        
        B = eval(get(handles.edit_clinicalvar,'string'));
    end
catch
    errordlg('Invalid behavioral variable.'), return
end

covarstring = get(handles.edit_clinicalcovar,'string');
if ~isempty(covarstring) && ~strcmp(covarstring,'(none)')
    try    
        if strcmp(conns.files.behavioral,'workspace')
            C = evalin('base',covarstring);
        else
            C = eval(covarstring);
        end
    catch
        errordlg('Invalid behavioral covariable(s).'), return
    end
else
    C = [];
end

beam = fcm_corr(conns.files.activepop,B,C);

[fi, pa] = uiputfile('s_beam*.mat', 'Save correlation as...');
if isequal(fi,0), return, end

conns.files.correlation = fullfile(pa,fi);
save(conns.files.correlation,'beam');
update_beamlist(handles,fi);

% --- Executes on button press in push_snpm.
function push_snpm_Callback(hObject, eventdata, handles)
% hObject    handle to push_snpm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global conns

if ~isfield(conns.files,'activepop')
    errordlg('No active/patient population specified.'), return
end
try
    if strcmp(conns.files.behavioral,'workspace')
        B = evalin('base',get(handles.edit_clinicalvar,'string'));
    else
        load(conns.files.behavioral)        
        B = eval(get(handles.edit_clinicalvar,'string'));
    end
catch
    errordlg('Invalid behavioral variable.'), return
end

load(conns.files.activepop);

if length(B)~=length(pop.subjnr)
    errordlg('Behavioral variable does not match number of subjects'), return
end
covarstring = get(handles.edit_clinicalcovar,'string');
if ~isempty(covarstring) && isempty(strfind(covarstring,'none'))
    errordlg('Covariates are not supported for SnPM correlations.'), return
end
if ~exist(conns.files.activepop,'file')
    errordlg(sprintf('Population file %s does not exist',conns.files.activepop)), return
end

stats = struct('type','corr_snpm','docomp',false,'fsel','one','frqstats',false, ...
    'corr4multfreq',false,'tsel','one','tsstats',false,'corr4multtime',false, ...
    'covar',[],'whichtail','both','markmode',2,'doavg',false,'dostat',true, ...
    'usespatnormed',true,'subjsel',unique(pop.subjnr),'numsubj',length(unique(pop.subjnr)), ...
    'condsel',unique(pop.condnr),'numtests',1,'comps',1,'subjnr',pop.subjnr, ...
    'condnr',pop.condnr,'files',{{conns.files.activepop}},'behavdata',B(:), ...
    'frqband',[],'frqtxt','','timepts',[],'timetxt','');

if ~all(isfinite(stats.behavdata))
    bad = ~isfinite(stats.behavdata);
    stats.behavdata(bad)= [];
    stats.subjsel(bad)  = [];
    stats.subj2remove=find(bad);
    stats.numsubj       = length(stats.subjsel);
    stats.subjnr        = stats.subjsel;
end
stats.numperm =  2^length(stats.behavdata);
stats.isaproxperm = (stats.numperm>10000);
if stats.isaproxperm, stats.numperm=10000; end

nut_beamstats(stats);

% --- Executes on selection change in list_beam.
function list_beam_Callback(hObject, eventdata, handles)
% hObject    handle to list_beam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns list_beam contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_beam
update_beamlist(handles)

% --- Executes during object creation, after setting all properties.
function list_beam_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_beam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_clinicalsource_Callback(hObject, eventdata, handles)
% hObject    handle to edit_clinicalsource (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_clinicalsource as text
%        str2double(get(hObject,'String')) returns contents of edit_clinicalsource as a double

% --- Executes during object creation, after setting all properties.
function edit_clinicalsource_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_clinicalsource (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_clinicalvar_Callback(hObject, eventdata, handles)
% hObject    handle to edit_clinicalvar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_clinicalvar as text
%        str2double(get(hObject,'String')) returns contents of edit_clinicalvar as a double

% --- Executes during object creation, after setting all properties.
function edit_clinicalvar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_clinicalvar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_clinicalcovar_Callback(hObject, eventdata, handles)
% hObject    handle to edit_clinicalcovar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_clinicalcovar as text
%        str2double(get(hObject,'String')) returns contents of edit_clinicalcovar as a double
cvstr= get(hObject,'string');
if ~isempty(cvstr) && isempty(strfind(cvstr,'none'))
    set(handles.run_corr,'string','Partial Correlation')
else
    set(handles.run_corr,'string','Pearson Correlation')
end


% --- Executes during object creation, after setting all properties.
function edit_clinicalcovar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_clinicalcovar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_activepop_Callback(hObject, eventdata, handles)
% hObject    handle to edit_activepop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_activepop as text
%        str2double(get(hObject,'String')) returns contents of edit_activepop as a double
global conns
conns.files.activepop = get(hObject,'string');

% --- Executes during object creation, after setting all properties.
function edit_activepop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_activepop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_controlpop_Callback(hObject, eventdata, handles)
% hObject    handle to edit_controlpop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_controlpop as text
%        str2double(get(hObject,'String')) returns contents of edit_controlpop as a double
global conns
conns.files.controlpop = get(hObject,'string');

% --- Executes during object creation, after setting all properties.
function edit_controlpop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_controlpop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%------------------
function fi = fullpath2file(fu)
[pa,fi] = fileparts(fu);
