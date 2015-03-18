function varargout = fcm_gui(varargin)
% FCM_GUI  GUI for the functional connectivity map (FCM) toolbox.
%      FCM_GUI, by itself, creates a new FCM_GUI or raises the existing
%      singleton*.
%
%      H = FCM_GUI returns the handle to a new FCM_GUI or the handle to
%      the existing singleton*.
%
%      FCM_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FCM_GUI.M with the given input arguments.
%
%      FCM_GUI('Property','Value',...) creates a new FCM_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fcm_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fcm_gui_OpeningFcn via varargin.



% Edit the above text to modify the response to help fcm_gui

% Last Modified by GUIDE v2.5 19-Oct-2011 14:05:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fcm_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @fcm_gui_OutputFcn, ...
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


% --- Executes just before fcm_gui is made visible.
function fcm_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fcm_gui (see VARARGIN)

% Choose default command line output for fcm_gui
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

global fuse nuts conns

if isempty(fuse)
    fcm_start
end
if isfield(nuts,'sessionfile')
    conns.files.session=nuts.sessionfile;
end
update_gui(handles);
% UIWAIT makes fcm_gui wait for user response (see UIRESUME)
% uiwait(handles.fcm_gui);


% --- Outputs from this function are returned to the command line.
function varargout = fcm_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in push_setconfig.
%function push_setconfig_Callback(hObject, eventdata, handles)
% hObject    handle to push_setconfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%fcm_start gui


% --- Executes on button press in browse_session.
function browse_session_Callback(hObject, eventdata, handles)
% hObject    handle to browse_session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global conns nuts
[fi,pa]=uigetfile('*.mat','Select NUTMEG session file...');
if isequal(fi,0), return, end
fu=fullfile(pa,fi);
conns.files.session = fu;
nuts = load(conns.files.session,'voxels','voxelsize','coreg');
set(handles.edit_session,'string',fi)

% --- Executes on button press in browse_filtdata.
function browse_filtdata_Callback(hObject, eventdata, handles)
% hObject    handle to browse_filtdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global conns
[fi,pa]=uigetfile('filtdata_*.mat','Select filtered data file...');
if isequal(fi,0), return, end
fu=fullfile(pa,fi);
conns.files.filtdata = fu;
set(handles.edit_filtdata,'string',fi)
    
% --- Executes on button press in browse_W.
function browse_W_Callback(hObject, eventdata, handles)
% hObject    handle to browse_W (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global conns
[fi,pa]=uigetfile('W_*.mat;weight*.mat;*.is','Select weight matrix file...');
if isequal(fi,0), return, end
fu=fullfile(pa,fi);
conns.files.W = fu;
set(handles.edit_W,'string',fi)

% --- Executes on button press in browse_ext.
function browse_ext_Callback(hObject, eventdata, handles)
% hObject    handle to browse_ext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global conns
[fi,pa]=uigetfile('*.mat','Select extracerebral data file...');
if isequal(fi,0), return, end
fu=fullfile(pa,fi);
test=who('-file',fu);
if ~any(strcmp(test,'ext'))
    errordlg('File does not contain variable "ext".')
    return
end
conns.files.ext = fu;
set(handles.edit_ext,'string',fi)

% --- Executes on button press in push_psd.
function push_psd_Callback(hObject, eventdata, handles)
% hObject    handle to push_psd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global nuts conns

wbar = waitbar(0,'Calculating sensor power spectrogram...');
if isfield(nuts,'meg') && isfield(nuts.meg,'data') && ~isempty(nuts.meg.data)
    data=nuts.meg.data; fs=nuts.meg.srate;
    waitbar(.1,wbar);
else
    if isfield(conns,'files') && isfield(conns.files,'session')
        tmp=load(conns.files.session,'meg');
        data=tmp.meg.data; fs=tmp.meg.srate; clear tmp
        waitbar(.4,wbar);
    else
        close(wbar)
        errordlg('Please select NUTMEG session file first.')
        return
    end
end

%periodogram3(data,fs,45);
periodogram3(data(:,1:2:end,1:2:end),fs,45);  % Change here if you want to test every channel and trial at expense of speed.
waitbar(1,wbar);
close(wbar)


% --- Executes on button press in push_selvoxidx.
function push_selvoxidx_Callback(hObject, eventdata, handles)
% hObject    handle to push_selvoxidx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global fuse

if ( isfield(fuse,'roi') && fuse.roi>0 );
    roi = fcm_roiselect_gui;
else
    sel=menu('How do you want to define the ROI?','Paint polygons','Select anatomical ROI','Load MATLAB file with ROI voxel coordinates','Cancel');
    switch sel
        case 1
            roi = fcm_paintroi;      
        case 2
            roi = fcm_anatroi('gui'); 
        case 3
            [fi,pa]=uigefile('*.mat','Load ROI voxel coordinates...');
            if isequal(fi,0), return, end
            roi=load(fullfile(pa,fi));
            if ~any(isfield(roi,{'voxels' 'MEGvoxels' 'MNIvoxels' 'MRIvoxels'}))
                errordlg('File must contain one of the variables "MRIvoxels", "MNIvoxels" or "MEGvoxels".')
                return
            end
        case 4
            return
    end
end
if isempty(roi), return, end

fcm_roiidx(roi);

figure(handles.fcm_gui)
set(handles.edit_seedcoord,'string','ROI');

% --- Executes on button press in push_seedcurr.
function push_seedcurr_Callback(hObject, eventdata, handles)
% hObject    handle to push_seedcurr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global st nuts conns

if isempty(findobj('tag','nutmegfig'))
    errordlg('Open NUTMEG main GUI first.'), return
end
if ~isfield(st,'megp')
    errordlg('NUTMEG has not been properly initialized. Start it over and try again.'), return
end
if ~isfield(nuts,'voxels')
    if issubfield(conns,'files.session')
        nuts=load(conns.files.session,'voxels');
        if ~isfield(nuts,'voxels')
            errordlg('Your session file seems to be invalid or saved in an old format.'), return
        end
    else
        errordlg('Please select NUTMEG session file first.')
        return
    end
end
coord = str2num(get(st.megp,'string'));
nuts.selvox.ipsi = dsearchn(nuts.voxels,coord);
set(handles.edit_seedcoord,'string',sprintf('%1.2f %1.2f %1.2f',nuts.voxels(nuts.selvox.ipsi,:)));

% -------------
function iserr=getparams(handles)

global conns nuts fuse
iserr=false;


conns.freq.ftype = {'s' 'b'};
conns.freq.ftype = conns.freq.ftype{get(handles.pop_ftype,'Value')};
conns.freq.low = get(handles.edit_lofrq,'string'); 
conns.freq.high = get(handles.edit_hifrq,'string');
if any(strcmp(fuse.funconn,{'ampcorr' 'pli'})) && strcmp(fuse.datatype,'3')
    conns.freq.nfft = '';
else
    conns.freq.nfft = get(handles.pop_nfft,'string');
    conns.freq.nfft = conns.freq.nfft{get(handles.pop_nfft,'value')};
end

conns.taper.do = get(handles.check_mtaper,'value');
if conns.taper.do
    conns.taper.width = get(handles.edit_taper,'string');
end

conns.time.len = get(handles.edit_len,'string');
conns.time.step= get(handles.edit_step,'string');
conns.baseline.start = get(handles.edit_basestart,'string');
conns.baseline.stop  = get(handles.edit_basestop,'string');

if strcmpi(get(handles.push_selvoxidx,'enable'),'on')
    coord = get(handles.edit_seedcoord,'string');
    if ~strcmpi(coord,'ROI')
        coord = str2num(coord);
        if isempty(coord)
            errordlg('Invalid seed voxel coordinate(s).')
            iserr=true; return
        end
        nuts.selvox.ipsi = dsearchn(nuts.voxels,coord);
    elseif ~isfield(nuts,'selvox') || ~isfield(nuts.selvox,'ipsi') || isempty(nuts.selvox.ipsi)
        errordlg('ROI voxel definition seems to have failed.')
        iserr=true; return
    end
end

% --- Executes on button press in push_runlocal.
function push_runlocal_Callback(hObject, eventdata, handles)
% hObject    handle to push_runlocal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global conns fuse nuts

if ~isfield(nuts,'voxels')
    if ( isfield(conns,'files') && isfield(conns.files,'session') && exist(conns.files.session,'file') )
        nut_opensession(conns.files.session,true);
    else
        errordlg('You must open a session file first.'), return
    end
end

% get parameters
iserr=getparams(handles);
if iserr, return, end
if isfield(nuts,'meg')
    nuts.meg.data=[]; %save memory
end

% ROI labelling for current subject
useroi=(isfield(fuse,'roi') && fuse.roi>0);
if useroi
    [fi,pa]=uigetfile('*.mat','ROI transformation matrices available? (Cancel otherwise)');
    if isequal(fi,0)
        hw=waitbar(0,'Labelling ROIs... Please wait...');
        R=fcm_voxel2roi(nuts.voxels,nuts.coreg);
        waitbar(1,hw);
        delete(hw);
        [fi,pa]=uiputfile('*.mat','Save ROI transformation matrices as...','roi.mat');
        if isequal(fi,0), return, end
        roifile=fullfile(pa,fi);
        save(roifile,'R')
    else
        roifile=fullfile(pa,fi);
        fcm_useroi(roifile);
    end
end

% define connections between voxel pairs
if strcmpi(fuse.connection,'grid') && ~isfield(fuse,'gridspacing')
    errordlg('You must define grid spacing first (click "Set configuration").')
    return
end
if strcmpi(fuse.seed,'extracerebral')
    load(conns.files.ext)
    comps = fcm_conndef(size(ext,2));
else
    comps = fcm_conndef;
end

% Calculate functional connectivity  
varargin={conns.files.filtdata conns.files.W comps};
if ~strcmp(fuse.datatype,'31')
    varargin = [varargin {conns.time.len conns.time.step}];
end
varargin = [varargin conns.freq.ftype str2num(conns.freq.low) str2num(conns.freq.high) conns.freq.nfft];

if strcmp(fuse.datatype,'3') && ~isempty(conns.baseline.start) && ~isempty(conns.baseline.stop)
    varargin = [varargin {'baseline' conns.baseline.start conns.baseline.stop}];
end
if strcmpi(fuse.seed,'extracerebral')
    varargin = [varargin {'ext' ext}];
end
if conns.taper.do
    varargin = [varargin {'taper' conns.taper.width}];
end
if useroi
    varargin = [varargin {'roi' roifile}];
end
% ttt=clock;
% varargin = [varargin {'savevirt' sprintf('rloc_%d%d%d_%dh%d',ttt(1:5))}];  clear ttt

% output filename
isok=false;
while ~isok
    [fi, pa] = uiputfile([fuse.funconn '*.mat'], 'Save functional connectivity as...',[fuse.funconn '_.mat']);
    if isequal(fi,0), return, end
    isok=strncmpi(fi,fuse.funconn,length(fuse.funconn));
end
conns.files.comcoh = fullfile(pa,fi);

% Calc
fcm_sourceconnectivity(conns.files.comcoh,varargin{:});

update_listcomcoh(handles,fi);

% --- Executes on button press in push_prepqsub.
function push_prepqsub_Callback(hObject, eventdata, handles)
% hObject    handle to push_prepqsub (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global conns fuse

% get parameters
iserr=getparams(handles);
if iserr, return, end
if ~isscalar(str2num(conns.freq.low)) || ~isscalar(str2num(conns.freq.high))
    errordlg('Multiple frequency bands are not possible when using the compiled function.'), return
end
numjobs = str2num(get(handles.edit_numjobs,'string'));
if isempty(numjobs), errordlg('Invalid number of cluster jobs.'), return, end
if numjobs>1 && strcmpi(fuse.seed,'all')
    errordlg('Use one job per subject when calculating all to all connections. Multiple jobs will be slower because each computer runs the same calculations.')
    return
end

% define connections between voxel pairs
if strcmpi(fuse.connection,'grid') && ~isfield(fuse,'gridspacing')
    errordlg('You must define grid spacing first (click "Set configuration").')
    return
end
if strcmpi(fuse.seed,'extracerebral')
    load(conns.files.ext)
    comps = fcm_conndef(size(ext,2));
else
    comps = fcm_conndef;
end
fcm_preparejobs(comps,numjobs);
cd ..

% ROI labelling for current subject
useroi=(isfield(fuse,'roi') && fuse.roi>0);
if useroi
    [fi,pa]=uigetfile('*.mat','ROI transformation matrices available? (Cancel otherwise)');
    if isequal(fi,0)
        hw=waitbar(0,'Labelling ROIs... Please wait...');
        R=fcm_voxel2roi(nuts.voxels,nuts.coreg,fuse.roidef);
        waitbar(1,hw);
        delete(hw);
        [fi,pa]=uiputfile('*.mat','Save ROI transformation matrices as...','roi.mat');
        if isequal(fi,0), return, end
        roifile=fullfile(pa,fi);        
        save(roifile,'R');
        rofime = ['\n' fi];
    else
        roifile=fullfile(pa,fi);
        rofime = ['\n' fi];
    end
else
    rofime = '';
end

% Inform user
fun = ['s' fuse.funconn fuse.datatype];
varargin={conns.files.filtdata conns.files.W comps};
if ~strcmp(fuse.datatype,'31')
    varargin = [varargin {conns.time.len conns.time.step}];
end
varargin = [varargin conns.freq.ftype conns.freq.low conns.freq.high conns.freq.nfft];

if strcmp(fuse.datatype,'3') && ~isempty(conns.baseline.start) && ~isempty(conns.baseline.stop)
    varargin = [varargin {'baseline' conns.baseline.start conns.baseline.stop}];
end
if strcmpi(fuse.seed,'extracerebral')
    varargin = [varargin {'ext' fullpath2file(conns.files.ext)}];
end
if conns.taper.do
    varargin = [varargin {'taper' conns.taper.width}];
end
if isfield(fuse,'roi') && fuse.roi>0
    varargin = [varargin {'roi' fullpath2file(roifile)}];
end
msgbox(sprintf(['Your data is now ready to be processed on a linux cluster. You may need to transfer the following files/directories to the cluster:\n\n' ...
    fullpath2file(conns.files.filtdata) '\n' fullpath2file(conns.files.W) '\ncomps' rofime '\n\n\nEach job can then be called as follows (all on 1 line):\n\n' ...
    fun ' ' sprintf('%s ',varargin{:}) ...
    '\n\n\nSee qCC' fuse.datatype ' in the "scripts" subfolder for an example.\nOnce done, transfer the results in the folder "' fuse.funconn '" back to the local machine and click "Assemble results from linux cluster"']),'FCM message')


% --- Executes on button press in push_getcomcoh.
function push_getcomcoh_Callback(hObject, eventdata, handles)
% hObject    handle to push_getcomcoh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global conns fuse

isok=false;
while ~isok
    [fi, pa] = uiputfile([fuse.funconn '*.mat'], 'Save functional connectivity as...',[fuse.funconn '_.mat']);
    if isequal(fi,0), return, end
    isok=strncmpi(fi,fuse.funconn,length(fuse.funconn));
end
conns.files.comcoh = fullfile(pa,fi);

fcm_assemble(conns.files.comcoh);
save(conns.files.comcoh,'CC');
      
update_listcomcoh(handles,conns.files.comcoh);
set(handles.push_comcoh2beam,'enable','on')


% --- Executes on button press in push_comcoh2beam.
function push_comcoh2beam_Callback(hObject, eventdata, handles)
% hObject    handle to push_comcoh2beam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global conns fuse
wbar = waitbar(0,'Calculating output volume...');

comcohlist = get(handles.list_comcoh,'string');
if isempty(comcohlist), return,  end
conns.files.comcoh = comcohlist{get(handles.list_comcoh,'value')};
waitbar(.1,wbar);
if isfield(conns.files,'W') && ( ~isfield(fuse,'roi') || fuse.roi<1 )
    fcm_comcoh2beam(conns.files.comcoh,conns.files.W);
else
    fcm_comcoh2beam(conns.files.comcoh);
end
waitbar(1,wbar)
close(wbar)

[pa,fi] = fileparts(conns.files.comcoh);
fi = ['s_beamtf_' fi];

conns.files.beam = fullfile(pa,fi);
update_listbeam(handles,fi);

% --- Executes on button press in push_refresh.
function push_refresh_Callback(hObject, eventdata, handles)
% hObject    handle to push_refresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

update_listcomcoh(handles);
update_listbeam(handles);


% --- Executes on button press in push_view.
function push_view_Callback(hObject, eventdata, handles)
% hObject    handle to push_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global conns
beamlist = get(handles.list_beam,'string');
if ~isempty(beamlist)
    conns.files.beam = beamlist{get(handles.list_beam,'value')};
    tv(conns.files.beam)
end


% --- Executes on button press in push_acrosssubj_gui.
function push_acrosssubj_gui_Callback(hObject, eventdata, handles)
% hObject    handle to push_acrosssubj_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fcm_across_subject_gui;
delete(handles.fcm_gui);

%----------------------------------
function update_gui(handles)

global fuse conns

if ~isfield(fuse,'datatype'), fuse.datatype='31'; end
if ~isfield(fuse,'funconn'),  fuse.funconn='ccohere'; end

if isfield(conns,'files') && isfield(conns.files,'session')
    set(handles.edit_session,'string',conns.files.session)
end

set(handles.text_seed,'string',['Seed: ' fuse.seed])
set(handles.text_connections,'string',['Connections: ' fuse.connection])
set(handles.text_funconn,'string',['Measure: ' fuse.funconn])
set(handles.text_datatype,'string',['Datatype: ' fuse.datatype])
if strcmpi(fuse.seed,'Extracerebral')
    set([handles.text22 handles.edit_ext handles.browse_ext],'enable','on')
else
    set([handles.text22 handles.edit_ext handles.browse_ext],'enable','off')
end
if strncmpi(fuse.seed,'Selected',8) || strcmpi(fuse.connectionavg,'Seed')
    set([handles.push_selvoxidx handles.push_seedcurr handles.edit_seedcoord handles.text14],'Enable','on')
else
    set([handles.push_selvoxidx handles.push_seedcurr handles.edit_seedcoord handles.text14],'Enable','off')  
end
if ~strcmp(fuse.datatype,'31')
    set([handles.text38 handles.text39 handles.text40 handles.text41 handles.edit_len handles.edit_step],'enable','on')
else
    set([handles.text38 handles.text39 handles.text40 handles.text41 handles.edit_len handles.edit_step],'enable','off')
end
if strcmp(fuse.datatype,'3')
    set([handles.text45 handles.text46 handles.text47 handles.edit_basestart handles.edit_basestop],'enable','on')
else
    set([handles.text45 handles.text46 handles.text47 handles.edit_basestart handles.edit_basestop],'enable','off')
end
if get(handles.pop_ftype,'Value')==2
    set([handles.text37 handles.check_mtaper handles.edit_taper],'enable','off')
    set(handles.check_mtaper,'Value',0)
else
    set([handles.text37 handles.check_mtaper handles.edit_taper],'enable','on')
end

update_listcomcoh(handles);
update_listbeam(handles);

%-------------------------
function update_listcomcoh(handles,preselect)

tmp = dir('ccohere*.mat'); names = {tmp.name}';
tmp = dir('nccohere*.mat'); names = [names;{tmp.name}'];
tmp = dir('glcohere*.mat'); names = [names;{tmp.name}'];
tmp = dir('comcoh*.mat'); names = [names;{tmp.name}'];
tmp = dir('ampcorr*.mat'); names = [names;{tmp.name}'];
tmp = dir('pli*.mat'); names = [names;{tmp.name}']; clear tmp
sorted_names = sortrows(names);
set(handles.list_comcoh,'String',sorted_names)
if nargin>1
    selidx = strmatch(preselect,sorted_names);
    if ~isempty(selidx)
        set(handles.list_comcoh,'Value',selidx(1))
    else
        set(handles.list_comcoh,'Value',1)
    end
elseif isempty(get(handles.list_comcoh,'Value')) && ~isempty(sorted_names)
    set(handles.list_comcoh,'Value',1)
elseif get(handles.list_comcoh,'Value')>length(sorted_names)
    set(handles.list_comcoh,'Value',1)
end

%-----------------
function update_listbeam(handles,preselect)

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

% --- Executes on selection change in pop_ftype.
function pop_ftype_Callback(hObject, eventdata, handles)
% hObject    handle to pop_ftype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns pop_ftype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_ftype
global fuse
if get(hObject,'Value')==2 
    set([handles.text37 handles.check_mtaper handles.edit_taper],'enable','off')
    set(handles.check_mtaper,'Value',0)
else
    set([handles.text37 handles.check_mtaper handles.edit_taper],'enable','on')
end

% --- Executes during object creation, after setting all properties.
function pop_ftype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_ftype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_lofrq_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lofrq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lofrq as text
%        str2double(get(hObject,'String')) returns contents of edit_lofrq as a double


% --- Executes during object creation, after setting all properties.
function edit_lofrq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lofrq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_hifrq_Callback(hObject, eventdata, handles)
% hObject    handle to edit_hifrq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_hifrq as text
%        str2double(get(hObject,'String')) returns contents of edit_hifrq as a double


% --- Executes during object creation, after setting all properties.
function edit_hifrq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_hifrq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_nfft.
function pop_nfft_Callback(hObject, eventdata, handles)
% hObject    handle to pop_nfft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns pop_nfft contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_nfft


% --- Executes during object creation, after setting all properties.
function pop_nfft_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_nfft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in check_mtaper.
function check_mtaper_Callback(hObject, eventdata, handles)
% hObject    handle to check_mtaper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_mtaper

function edit_taper_Callback(hObject, eventdata, handles)
% hObject    handle to edit_taper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_taper as text
%        str2double(get(hObject,'String')) returns contents of edit_taper as a double


% --- Executes during object creation, after setting all properties.
function edit_taper_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_taper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_session_Callback(hObject, eventdata, handles)
% hObject    handle to edit_session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_session as text
%        str2double(get(hObject,'String')) returns contents of edit_session as a double


% --- Executes during object creation, after setting all properties.
function edit_session_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_filtdata_Callback(hObject, eventdata, handles)
% hObject    handle to edit_filtdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_filtdata as text
%        str2double(get(hObject,'String')) returns contents of edit_filtdata as a double


% --- Executes during object creation, after setting all properties.
function edit_filtdata_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_filtdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_W_Callback(hObject, eventdata, handles)
% hObject    handle to edit_W (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_W as text
%        str2double(get(hObject,'String')) returns contents of edit_W as a double


% --- Executes during object creation, after setting all properties.
function edit_W_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_W (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_ext_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ext as text
%        str2double(get(hObject,'String')) returns contents of edit_ext as a double


% --- Executes during object creation, after setting all properties.
function edit_ext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_seedcoord_Callback(hObject, eventdata, handles)
% hObject    handle to edit_seedcoord (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_seedcoord as text
%        str2double(get(hObject,'String')) returns contents of edit_seedcoord as a double


% --- Executes during object creation, after setting all properties.
function edit_seedcoord_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_seedcoord (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_len_Callback(hObject, eventdata, handles)
% hObject    handle to edit_len (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_len as text
%        str2double(get(hObject,'String')) returns contents of edit_len as a double


% --- Executes during object creation, after setting all properties.
function edit_len_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_len (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_step_Callback(hObject, eventdata, handles)
% hObject    handle to edit_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_step as text
%        str2double(get(hObject,'String')) returns contents of edit_step as a double


% --- Executes during object creation, after setting all properties.
function edit_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_basestart_Callback(hObject, eventdata, handles)
% hObject    handle to edit_basestart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_basestart as text
%        str2double(get(hObject,'String')) returns contents of edit_basestart as a double

% --- Executes during object creation, after setting all properties.
function edit_basestart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_basestart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_basestop_Callback(hObject, eventdata, handles)
% hObject    handle to edit_basestop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_basestop as text
%        str2double(get(hObject,'String')) returns contents of edit_basestop as a double

% --- Executes during object creation, after setting all properties.
function edit_basestop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_basestop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_numjobs_Callback(hObject, eventdata, handles)
% hObject    handle to edit_numjobs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_numjobs as text
%        str2double(get(hObject,'String')) returns contents of edit_numjobs as a double


% --- Executes during object creation, after setting all properties.
function edit_numjobs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_numjobs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in list_comcoh.
function list_comcoh_Callback(hObject, eventdata, handles)
% hObject    handle to list_comcoh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns list_comcoh contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_comcoh
update_listcomcoh(handles);

% --- Executes during object creation, after setting all properties.
function list_comcoh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_comcoh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in list_beam.
function list_beam_Callback(hObject, eventdata, handles)
% hObject    handle to list_beam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns list_beam contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_beam
update_listbeam(handles);

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

function periodogram3(data,fs,maxfrq)
% plots average periodgram across all sensors and trials

[len,numchan,numtrial]=size(data);
nfft=2^nextpow2(2*fs);

P = zeros(nfft/2+1,numchan,numtrial);
for k=1:numtrial
    for c=1:numchan
        data(:,c,k)=data(:,c,k)-mean(data(:,c,k));
        [P(:,c,k),Pf]=periodogram(data(:,c,k),hamming(len),nfft,fs);      % [Pxx,F] = PWELCH(X,WINDOW,NOVERLAP,NFFT,Fs)
    end
end
f = find(Pf<=maxfrq);
P = P(f,:,:);
Pf= Pf(f);

PP = mean(mean(10.*log10(P),3),2);

figure;
plot(Pf,PP);
xlabel('Frequency [Hz]','fontweight','bold')
ylabel('Power [dB]','fontweight','bold')
title('Average power spectrogram of all sensors and trials','fontweight','bold')
set(gca,'xlim',[min(Pf) max(Pf)])

%------------------
function fi = fullpath2file(fu)
[dum,fi] = fileparts(fu);



