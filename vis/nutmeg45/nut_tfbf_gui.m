function varargout = nut_tfbf_gui(varargin)
% NUT_TFBF_GUI  Graphical user interface for time-frequency source analysis.

% Last Modified by GUIDE v2.5 01-Nov-2011 18:23:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @nut_tfbf_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @nut_tfbf_gui_OutputFcn, ...
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


% --- Executes just before nut_tfbf_gui is made visible.
function nut_tfbf_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to nut_tfbf_gui (see VARARGIN)

% Choose default command line output for nut_tfbf_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes nut_tfbf_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);

global confi nuts
confi=[];

if ( isstruct(nuts) && isfield(nuts,'sessionfile') )
    confi.sessionfile = nuts.sessionfile;
    set(handles.text_sessionfilename,'string',confi.sessionfile);
end

update_listsbeams(pwd,handles);

% --- Outputs from this function are returned to the command line.
function varargout = nut_tfbf_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in load_sessionfile.
function load_sessionfile_Callback(hObject, eventdata, handles)
% hObject    handle to load_sessionfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global confi

[fMAT, pMAT] = uigetfile('*.mat','Select session file...');
if isequal(fMAT,0), return, end
set(handles.text_sessionfilename,'string',fMAT);
confi.sessionfile=fullfile(pMAT,fMAT);

update_listsbeams(pMAT,handles);

% if ~exist([confi.sessionfile(1:end-4) 'Lp.mat'],'file')
%     nut_liposession(confi.sessionfile);     % create Lp file here to avoid errors later
% end

% --- Executes on button press in load_timewindowfile.
function load_timewindowfile_Callback(hObject, eventdata, handles)
% hObject    handle to load_timewindowfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global confi

[fMAT, pMAT] = uigetfile('*.mat','Select timewindow file...');
if isequal(fMAT,0), return, end
timewindowfile=fullfile(pMAT,fMAT);
load(timewindowfile);

if ~exist('active','var') || ~exist('control','var')
    errordlg('Invalid timewindow file.')
    return
end

set(handles.edit_activestart,'string',num2str(active(1,1)));
set(handles.edit_activeend,'string',num2str(active(end,2)));
set([handles.edit_controlstart handles.edit_controlend],'string','')
if ~isempty(control)
    set(handles.edit_controlstart,'string',num2str(control(1,1)));
    set(handles.edit_controlend,'string',num2str(control(end,2)));
end
set(handles.edit_windowsize,'string',num2str(active(1,2)-active(1,1)));
if size(active,1)>1
    set(handles.edit_stepsize,'string',num2str(active(2,1)-active(1,1)));
else
    set(handles.edit_stepsize,'string',num2str(active(1,2)-active(1,1)));
end

set(handles.text_timewindowfilename,'string',fMAT);
confi.timewindowfile=timewindowfile;

% --- Executes on button press in save_timewindowfile.
function save_timewindowfile_Callback(hObject, eventdata, handles)
% hObject    handle to save_timewindowfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global confi

[fMAT, pMAT] = uiputfile('*.mat', 'Save timewindow settings as...');
if isequal(fMAT,0), return, end
timewindowfile=fullfile(pMAT,fMAT);

activestart = str2num(get(handles.edit_activestart,'string'));
activeend   = str2num(get(handles.edit_activeend,'string'));
controlstart= str2num(get(handles.edit_controlstart,'string'));
controlend  = str2num(get(handles.edit_controlend,'string'));
windowsize  = str2num(get(handles.edit_windowsize,'string'));
stepsize    = str2num(get(handles.edit_stepsize,'string'));
if windowsize>activeend-activestart
    errordlg('Your window length is larger than the time difference between Active End and Active Start.')
    return
end

nut_constructTFBFparams(activestart,activeend,windowsize,stepsize,controlstart,controlend,timewindowfile);
set(handles.text_timewindowfilename,'string',fMAT);
confi.timewindowfile=timewindowfile;

% --- Executes on button press in load_filterconfigfile.
function load_filterconfigfile_Callback(hObject, eventdata, handles)
% hObject    handle to load_filterconfigfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global confi

[fMAT, pMAT] = uigetfile([fileparts(which('nutmeg.m')) filesep 'tfbf' filesep 'params' filesep '*.mat'], 'Select filter configuration file...');
if isequal(fMAT,0), return, end
filterconfigfile=fullfile(pMAT,fMAT);
load(filterconfigfile);

if ~exist('filt','var')
    errordlg('Invalid filter configuration file.')
    return
end

if isfield(filt,'class')
    set(handles.radio_nutmegfilter,'value',1)
    val = strmatch(lower(filt.class),{'butter' 'cheby1' 'cheby2' 'ellip' 'firls' 'firpm'});
    if ~isscalar(val), errordlg('Do not know this filter class: %s',filt.class), return, end
    set(handles.pop_filterclass, 'value', val)
    val = strcmpi(filt.type,'notch')+1;
    set(handles.pop_filtertype, 'value', val)
    if ~isfield(filt,'order') || isempty(filt.order), filt.order='auto';
    elseif isnumeric(filt.order), filt.order=int2str(filt.order); 
    end
    set(handles.edit_filterorder,'string',filt.order)
elseif isfield(filt,'callback')
    set(handles.radio_ownfilter,'value',1)
    set(handles.edit_ownfiltercallback,'string',filt.callback)
else
    errordlg('Invalid filter configuration file.')
    return
end

set(handles.text_filterconfigfile,'string',fMAT);
confi.filterconfigfile=filterconfigfile;


% --- Executes on button press in save_filterconfigfile.
function save_filterconfigfile_Callback(hObject, eventdata, handles)
% hObject    handle to save_filterconfigfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global confi

[fMAT, pMAT] = uiputfile([fileparts(which('nutmeg.m')) filesep 'tfbf' filesep 'params' filesep '*.mat'], 'Save filter configuration as...');
if isequal(fMAT,0), return, end
filterconfigfile=fullfile(pMAT,fMAT);

filt.class = {'butter' 'cheby1' 'cheby2' 'ellip' 'firls' 'firpm'};
filt.class = filt.class{get(handles.pop_filterclass,'value')};
filt.type  = {'band/high/low' 'notch'};
filt.type  = filt.type{get(handles.pop_filtertype, 'value')};
filt.order = get(handles.edit_filterorder,'string');
if ~strcmpi(filt.order,'auto'), filt.order=str2double(filt.order); end

save(filterconfigfile,'filt');
set(handles.text_filterconfigfile,'string',fMAT);
confi.filterconfigfile=filterconfigfile;

% --- Executes on button press in load_beamformerconfigfile.
function load_beamformerconfigfile_Callback(hObject, eventdata, handles)
% hObject    handle to load_beamformerconfigfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global confi

[fMAT, pMAT] = uigetfile([fileparts(which('nutmeg.m')) filesep 'tfbf' filesep 'params' filesep '*.mat'], 'Select beamformer configuration file...');
if isequal(fMAT,0), return, end
beamconfigfile=fullfile(pMAT,fMAT);
load(beamconfigfile);

if ~exist('params','var')
    errordlg('Invalid beamformer configuration file.')
    return
end
if ~isfield(params,'cn'), params.cn=true; end
if ~isfield(params,'wn'), params.wn=~params.cn; end

selalgo=[];
for k=1:length(params.algo)
    selalgo(end+1)=strmatch(lower(params.algo{k}),lower(get(handles.list_algo,'string')),'exact');
end
set(handles.list_algo,'value',selalgo)

seloutput=[];
if params.savepower,        seloutput(end+1)=1; end
if params.saveweights,      seloutput(end+1)=2; end
if params.savefiltereddata, seloutput(end+1)=3; end
set(handles.list_output,'value',seloutput)

set(handles.pop_cn,'value',bin2dec(int2str([params.wn params.cn]))+1)

if isfield(params,'regularization')
    rl = get(handles.pop_regul,'string');
    val = find(strcmpi(params.regularization,rl));
    if isempty(val), val=1; end
    set(handles.pop_regul,'Value',val)
else
    set(handles.pop_regul,'Value',1)
end
if get(handles.pop_regul,'Value')==1
    set(handles.pop_regulthres,'Enable','off')
else
    set(handles.pop_regulthres,'Enable','on')
end
if isfield(params,'regulthres') && strcmpi(params.regulthres,'auto')
    set(handles.pop_regulthres,'Value',2)
else
    set(handles.pop_regulthres,'Value',1)
end
        
set(handles.check_qsub,'value',params.qsub)

set(handles.text_beamformerconfigfile,'string',fMAT);
confi.beamconfigfile=beamconfigfile;


% --- Executes on button press in save_beamformerconfigfile.
function save_beamformerconfigfile_Callback(hObject, eventdata, handles)
% hObject    handle to save_beamformerconfigfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global confi

[fMAT, pMAT] = uiputfile([fileparts(which('nutmeg.m')) filesep 'tfbf' filesep 'params' filesep '*.mat'], 'Save beamformer configuration as...');
if isequal(fMAT,0), return, end
beamconfigfile=fullfile(pMAT,fMAT);

params.algo = get(handles.list_algo,'String');
params.algo = params.algo(get(handles.list_algo,'value'));
seloutput = get(handles.list_output,'value');
params.savepower        = any( seloutput == 1 );
params.saveweights      = any( seloutput == 2 );
params.savefiltereddata = any( seloutput == 3 );
params.cn = (get(handles.pop_cn,'value')==2);
params.wn = (get(handles.pop_cn,'value')==3);
lr = get(handles.pop_regul,'string');
params.regularization = lr{get(handles.pop_regul,'Value')};
if ~strcmpi(params.regularization,'none')
    lrt = {'always' 'auto'};
    params.regulthres = lrt{get(handles.pop_regulthres,'Value')};
end
params.qsub = logical(get(handles.check_qsub,'value'));

save(beamconfigfile,'params');
set(handles.text_beamformerconfigfile,'string',fMAT);
confi.beamconfigfile=beamconfigfile;

%-----------------------------
function iserr = check_files(reqfiles)

if nargin<1, reqfiles=1:4; end
global confi

test = ~isfield(confi,{'sessionfile' 'timewindowfile' 'filterconfigfile' 'beamconfigfile'});
test(setdiff(1:4,reqfiles))=false;
iserr = any(test);

if iserr
    needed={'session' 'timewindow' 'filter configuration' 'beamformer configuration'};
    errortxt=['The following file(s) was/were not specified: ' sprintf('%s file, ',needed{find(test)})];
    errortxt(end-1)='.';
    errordlg(errortxt);
end

%-----------------------------
function update_listsbeams(dir_path,handles)
cd(dir_path);
dir_struct = dir('s_beam*.mat');
sorted_names = sortrows({dir_struct.name}');
set(handles.list_sbeams,'String',sorted_names,'Value',1)

% --- Executes on button press in push_run.
function push_run_Callback(hObject, eventdata, handles)
% hObject    handle to push_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global confi
iserr=check_files;
if iserr, return, end

if isempty(which('tfbf'))
    addpath([fileparts(which('nutmeg')) filesep 'tfbf' filesep 'src']);
end
tfbf(confi.sessionfile,confi.timewindowfile, ...
    get(handles.edit_lowcutoff,'string'),get(handles.edit_highcutoff,'string'), ...
    confi.filterconfigfile,confi.beamconfigfile);

update_listsbeams(pwd,handles);

% --- Executes on button press in push_newbatch.
function push_newbatch_Callback(hObject, eventdata, handles)
% hObject    handle to push_newbatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global confi

[fM, pM] = uiputfile('*.m', 'Create new batch file...');
if isequal(fM,0), return, end
confi.batchfile=fullfile(pM,fM);

% --- Executes on button press in push_add2batch.
function push_add2batch_Callback(hObject, eventdata, handles)
% hObject    handle to push_add2batch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global confi
iserr=check_files;
if iserr, return, end

currcall=sprintf('tfbf(''%s'',''%s'',%s,%s,''%s'',''%s'');',confi.sessionfile, ...
    confi.timewindowfile,get(handles.edit_lowcutoff,'string'), ...
    get(handles.edit_highcutoff,'string'), ...
    confi.filterconfigfile,confi.beamconfigfile);

if ~isfield(confi,'batchfile')
    %[fM, pM] = uiputfile('*.m', 'Create new batch file...');
    %if isequal(fM,0), return, end
    %confi.batchfile=fullfile(pM,fM);
    push_newbatch_Callback(gcbo,[],guidata(gcbo));
end
fid=fopen(confi.batchfile,'a+t');
fseek(fid,0,'eof');
fprintf(fid,'%s\n',currcall);
fclose(fid);

% --- Executes on button press in push_runbatch.
function push_runbatch_Callback(hObject, eventdata, handles)
% hObject    handle to push_runbatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global confi
if ~isfield(confi,'batchfile')
    [fM, pM] = uigetfile('*.m', 'Select batch file...');
    if isequal(fM,0), return, end
    confi.batchfile = fullfile(pM,fM);
end

if isempty(which('tfbf'))
    addpath([fileparts(which('nutmeg')) filesep 'tfbf' filesep 'src']);
end

[pM,fM,ext]=fileparts(confi.batchfile);
notinpath = isempty(which(fM));
if notinpath
    addpath(pM);    % need to add path here in case batch file is stored elsewhere
end

eval(fM);

if notinpath
    rmpath(pM);
end

update_listsbeams(pwd,handles);
msgbox('Batch done.')

% --- Executes on button press in push_assemble.
function push_assemble_Callback(hObject, eventdata, handles)
% hObject    handle to push_assemble (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global confi
iserr=check_files([1 3 4]);
if iserr, return, end

[pa,sf,ext]=fileparts(confi.sessionfile);
[pa,ff,ext]=fileparts(confi.filterconfigfile);
[pa,bf,ext]=fileparts(confi.beamconfigfile);
load(confi.beamconfigfile);
for k=1:length(params.algo)
    nut_tfbf2timef([sf '_' ff '_' bf],params.algo{k});
end

answer = questdlg('Would you like to delete the assembled single output files?','NUTMEG question','Yes','No','Yes');
if strcmp(answer,'Yes')
    for k=1:length(params.algo)
        delete(['s_beamtf_' sf '_' ff '_' bf '_*Hz_*ms_' params.algo{k} '.mat']);
    end
end

update_listsbeams(pwd,handles);
outputfile = ['s_beamtf_' sf '_' ff '_' bf '_' params.algo{1} '_all.mat'];
pos = strmatch(outputfile,get(handles.list_sbeams,'string'));
set(handles.list_sbeams,'value',pos);

% --- Executes on button press in push_view.
function push_view_Callback(hObject, eventdata, handles)
% hObject    handle to push_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

allfiles = get(handles.list_sbeams,'string');
if ~isempty(allfiles)
    selfile  = allfiles{get(handles.list_sbeams,'value')};
    if exist(selfile,'file')
        nut_results_viewer(selfile);
    end
end

% --- Executes on selection change in pop_filterclass.
function pop_filterclass_Callback(hObject, eventdata, handles)
% hObject    handle to pop_filterclass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns pop_filterclass contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_filterclass


% --- Executes during object creation, after setting all properties.
function pop_filterclass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_filterclass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radio_nutmegfilter.
function radio_nutmegfilter_Callback(hObject, eventdata, handles)
% hObject    handle to radio_nutmegfilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_nutmegfilter


% --- Executes on button press in radio_ownfilter.
function radio_ownfilter_Callback(hObject, eventdata, handles)
% hObject    handle to radio_ownfilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_ownfilter


% --- Executes on selection change in pop_filtertype.
function pop_filtertype_Callback(hObject, eventdata, handles)
% hObject    handle to pop_filtertype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns pop_filtertype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_filtertype


% --- Executes during object creation, after setting all properties.
function pop_filtertype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_filtertype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_filterorder_Callback(hObject, eventdata, handles)
% hObject    handle to edit_filterorder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_filterorder as text
%        str2double(get(hObject,'String')) returns contents of edit_filterorder as a double


% --- Executes during object creation, after setting all properties.
function edit_filterorder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_filterorder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ownfiltercallback_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ownfiltercallback (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ownfiltercallback as text
%        str2double(get(hObject,'String')) returns contents of edit_ownfiltercallback as a double


% --- Executes during object creation, after setting all properties.
function edit_ownfiltercallback_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ownfiltercallback (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_activestart_Callback(hObject, eventdata, handles)
% hObject    handle to edit_activestart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_activestart as text
%        str2double(get(hObject,'String')) returns contents of edit_activestart as a double


% --- Executes during object creation, after setting all properties.
function edit_activestart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_activestart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_activeend_Callback(hObject, eventdata, handles)
% hObject    handle to edit_activeend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_activeend as text
%        str2double(get(hObject,'String')) returns contents of edit_activeend as a double


% --- Executes during object creation, after setting all properties.
function edit_activeend_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_activeend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_controlstart_Callback(hObject, eventdata, handles)
% hObject    handle to edit_controlstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_controlstart as text
%        str2double(get(hObject,'String')) returns contents of edit_controlstart as a double


% --- Executes during object creation, after setting all properties.
function edit_controlstart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_controlstart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_controlend_Callback(hObject, eventdata, handles)
% hObject    handle to edit_controlend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_controlend as text
%        str2double(get(hObject,'String')) returns contents of edit_controlend as a double


% --- Executes during object creation, after setting all properties.
function edit_controlend_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_controlend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_windowsize_Callback(hObject, eventdata, handles)
% hObject    handle to edit_windowsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_windowsize as text
%        str2double(get(hObject,'String')) returns contents of edit_windowsize as a double


% --- Executes during object creation, after setting all properties.
function edit_windowsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_windowsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_stepsize_Callback(hObject, eventdata, handles)
% hObject    handle to edit_stepsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_stepsize as text
%        str2double(get(hObject,'String')) returns contents of edit_stepsize as a double

% --- Executes during object creation, after setting all properties.
function edit_stepsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_stepsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in list_algo.
function list_algo_Callback(hObject, eventdata, handles)
% hObject    handle to list_algo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns list_algo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_algo


% --- Executes during object creation, after setting all properties.
function list_algo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_algo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in list_output.
function list_output_Callback(hObject, eventdata, handles)
% hObject    handle to list_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns list_output contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_output


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


% --- Executes on selection change in pop_cn.
function pop_cn_Callback(hObject, eventdata, handles)
% hObject    handle to pop_cn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_cn contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_cn


% --- Executes during object creation, after setting all properties.
function pop_cn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_cn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_regul.
function pop_regul_Callback(hObject, eventdata, handles)
% hObject    handle to pop_regul (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_regul contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_regul
if get(hObject,'Value')==1
    set(handles.pop_regulthres,'Enable','off')
else
    set(handles.pop_regulthres,'Enable','on')
end

% --- Executes during object creation, after setting all properties.
function pop_regul_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_regul (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in pop_regulthres.
function pop_regulthres_Callback(hObject, eventdata, handles)
% hObject    handle to pop_regulthres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_regulthres contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_regulthres


% --- Executes during object creation, after setting all properties.
function pop_regulthres_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_regulthres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in check_qsub.
function check_qsub_Callback(hObject, eventdata, handles)
% hObject    handle to check_qsub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_qsub

function edit_lowcutoff_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lowcutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lowcutoff as text
%        str2double(get(hObject,'String')) returns contents of edit_lowcutoff as a double


% --- Executes during object creation, after setting all properties.
function edit_lowcutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lowcutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_highcutoff_Callback(hObject, eventdata, handles)
% hObject    handle to edit_highcutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_highcutoff as text
%        str2double(get(hObject,'String')) returns contents of edit_highcutoff as a double


% --- Executes during object creation, after setting all properties.
function edit_highcutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_highcutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
global confi
clear global confi
delete(hObject);


% --- Executes on selection change in list_sbeams.
function list_sbeams_Callback(hObject, eventdata, handles)
% hObject    handle to list_sbeams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns list_sbeams contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_sbeams


% --- Executes during object creation, after setting all properties.
function list_sbeams_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_sbeams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


