function varargout = visORICA(varargin)
% VISORICA MATLAB code for visORICA.fig
%      VISORICA, by itself, creates a new VISORICA or raises the existing
%      singleton*.
%
%      H = VISORICA returns the handle to a new VISORICA or the handle to
%      the existing singleton*.
%
%      VISORICA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VISORICA.M with the given input arguments.
%
%      VISORICA('Property','Value',...) creates a new VISORICA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before visORICA_OpeningFcn gets called. An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to visORICA_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help visORICA

% Last Modified by GUIDE v2.5 20-Nov-2014 13:50:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @visORICA_OpeningFcn, ...
    'gui_OutputFcn',  @visORICA_OutputFcn, ...
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


% --- Executes just before visORICA is made visible.
function visORICA_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to visORICA (see VARARGIN)
%   1: channel locations


% Parse varargsin
handles.chanlocs = evalin('base','chanlocs');
handles.ntopo = 8;
handles.nic = length(handles.chanlocs);
handles.ics = 1:handles.ntopo;
handles.streamName = 'visORICAst'; %'EEGLAB'; %
handles.curIC = 1;
handles.lock = [];
handles.exclude = [];
calibData = varargin{1};

% Start EEG stream
vis_stream_ORICA('StreamName',handles.streamName,'figurehandles',handles.figure1,'axishandles',handles.axisEEG);
% vis_stream_ORICA('figurehandles',handles.figure1,'axishandles',handles.axisEEG); % assume only one stream existed (FIXME: multiple streams case)
eegTimer = timerfind('Name','eegTimer');

% Intialize ORICA
initializeORICA(handles,calibData);

% Create ORICA timer
oricaTimer = timer('Period',.1,'ExecutionMode','fixedSpacing','TimerFcn',{@onl_filtered_ORICA,1},'StartDelay',0.1,'Tag','oricaTimer','Name','oricaTimer');

% Populate scalp maps
for it = 1:handles.ntopo
    set(handles.figure1, 'CurrentAxes', handles.(['axesIC' int2str(it)]))
    topoplotFast(rand(size(handles.chanlocs)), handles.chanlocs);
end

% Create scalp map timer
topoTimer = timer('Period',.5/handles.ntopo,'ExecutionMode','fixedRate','TimerFcn',{@vis_topo,hObject},'StartDelay',0.2,'Tag','topoTimer','Name','topoTimer');

% Create data timer (starts as power spectrum)
infoTimer = timer('Period',1,'ExecutionMode','fixedRate','TimerFcn',{@icPS,hObject},'StartDelay',0.2,'Tag','infoTimer','Name','infoTimer');

% Set panel and button colors
names = fieldnames(handles);
ind = find(any([strncmpi(names,'panel',5),strncmpi(names,'toggle',6),strncmpi(names,'push',4),strncmpi(names,'popup',5)],2));
for it = 1:length(ind)
set(handles.(names{ind(it)}),'BackgroundColor',get(handles.figure1,'Color'))
end

% Save timers
handles.pauseTimers = [eegTimer,topoTimer,infoTimer];

% Choose default command line output for visORICA
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Start timers
start(oricaTimer);
start(eegTimer);
start(topoTimer);
start(infoTimer);

% UIWAIT makes visORICA wait for user response (see UIRESUME)
% uiwait(handles.figure1);


function initializeORICA(handles,calibData)

run_readlsl_ORICA('MatlabStream',handles.streamName,'MarkerStreamQuery', []);

opts.lsl.StreamName = handles.streamName;
opts.BCILAB_PipelineConfigFile = 'ORICA_pipeline_config_realtime.mat'; % make sure this file doesn't have 'signal' entry

% define the pipeline configuration
try    fltPipCfg = exp_eval(io_load(opts.BCILAB_PipelineConfigFile));
catch, disp('-- no existing pipeline --'); fltPipCfg = {}; end
fltPipCfg = arg_guipanel('Function',@flt_pipeline, ...
    'Parameters',[{'signal',calibData} fltPipCfg], ...
    'PanelOnly',false);

% run pipline on calibration data
cleaned_data = exp_eval(flt_pipeline('signal',calibData,fltPipCfg));

% initialize the pipeline for streaming data
pipeline     = onl_newpipeline(cleaned_data,opts.lsl.StreamName);
assignin('base','pipeline',pipeline);


function icPS(varargin)

% plot PSD of selected IC
try
secs2samp = 5; % seconds

W = evalin('base','W');
% if isempty(W), W = evalin('base','Wn'); end
sphere = evalin('base','sphere');
handles = guidata(varargin{3});

set(handles.panelInfo,'Title',['Power spectral density of IC' int2str(handles.curIC)])

srate = evalin('base',['lsl_' handles.streamName '_stream.srate']);
data = evalin('base',['lsl_' handles.streamName '_stream.data']);
if all(data(:,end)==0)
    mark=1;
    while true
        ind = find(data(1,mark:end)==0,1);
        mark = mark+ind;
        if all(data(:,mark-1)==0)
            break; end
    end
    mark = mark-2;
    data = data(:,max(1,mark-srate*secs2samp+1):mark);
else
    data = data(:,max(1,end-srate*secs2samp+1):end);
end

data = bsxfun(@minus,data,mean(data,2));
data = W*sphere*data;
data = data(handles.curIC,:);

    [data,f] = pwelch(data,[],[],[],srate);
    
    plot(handles.axisInfo,f,db(data))
    grid(handles.axisInfo,'on');
    xlabel(handles.axisInfo,'Frequency (Hz)')
    ylabel(handles.axisInfo,'Power/Frequency (dB/Hz)')
    % axis(handles.axisInfo,[0 srate/2 -60 40])
    axis(handles.axisInfo,'tight')
    set(handles.axisInfo,'XTick',[0 10:10:f(end)])
end


function vis_topo(varargin)
% get the updated stream buffer
W = evalin('base','W');

% if isempty(W)
%     W = evalin('base','Wn');
% end
handles = guidata(varargin{3});

% update topo plot
it = mod(get(varargin{1},'TasksExecuted')-1,handles.ntopo)+1;
hstr = ['axesIC' int2str(it)];
hand = get(handles.(hstr),'children');

try
    sphere = evalin('base','sphere');
    Winv = inv(W*sphere);
    
    [map, cmin, cmax] = topoplotUpdate(Winv(:,handles.ics(it)), handles.chanlocs,'electrodes','off','gridscale',32);
    set(hand(end),'CData',map);
    set(handles.(hstr),'CLim',[cmin cmax]);
end

% update name and buttons
% old = get(handles.(['panelIC' int2str(it)]),'Title');
% old = str2num(old(regexp(old,'[\d]+'):end));
set(handles.(['panelIC' int2str(it)]),'Title',['IC' int2str(handles.ics(it))])
set(handles.(['togglebuttonLock' int2str(it)]),'Value',any(handles.lock==handles.ics(it)))
set(handles.(['togglebuttonExclude' int2str(it)]),'Value',any(handles.exclude==handles.ics(it)))


% --- Outputs from this function are returned to the command line.
function varargout = visORICA_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenuEEG.
function popupmenuEEG_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuEEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
eegTimer = timerfind('Name','eegTimer');
stop(eegTimer)
if get(hObject,'value')==1
    set(eegTimer,'UserData',0);
else
    set(eegTimer,'UserData',1);
end
start(eegTimer)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuEEG contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuEEG


% --- Executes during object creation, after setting all properties.
function popupmenuEEG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuEEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonIC.
function pushbuttonIC_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonIC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'figIC')
    figure(handles.figIC.handle)
    return
end

if strcmpi(get(handles.pushbuttonPause,'string'),'Resume')
    genICSelectGUI(hObject,handles)
else
    stop(handles.pauseTimers)
    genICSelectGUI(hObject,handles)
    start(handles.pauseTimers)
end


function genICSelectGUI(hObject,handles)
temp = get(handles.figure1,'Position');
fhandle = figure('toolbar','none','Menubar','none','Name','IC Select','position',[1 1 temp(3:4)],'Resize','on','DeleteFcn',{@closeFigIC,hObject});
W = evalin('base','W');
sphere = evalin('base','sphere');
Winv = inv(W*sphere);

rowcols(2) = ceil(sqrt(handles.nic));
rowcols(1) = ceil(handles.nic/rowcols(2));
scaleMatTopo = [1 0 0 0;0 1 0 0;0 0 1 0;0 .2 0 .8];
buttonGap = .1;
scaleMatExclude = [1 0 0 0;0 1 0 0;0 0 .5-buttonGap/2 0;.5+buttonGap/2 0 0 .2];
for it = 1:handles.nic
    h(it) = subaxis(rowcols(1),rowcols(2),it,'MR',.025,'ML',.025,'MT',.025,'MB',.025,'SH',0,'SV',0.02);
    tempPos = get(h(it),'Position');
    set(h(it),'position',get(h(it),'position')*scaleMatTopo)
    topoplotFast(Winv(:,it),handles.chanlocs);
    title(['IC' int2str(it)])
    
    buttonLock(it) = uicontrol('Style', 'togglebutton', 'String', 'Lock','Units','normalize','Position', tempPos.*[1 1 .5-buttonGap/2 .2],'Callback', {@lockIC,hObject,it},'Value',any(handles.lock==it));
    buttonExclude(it) = uicontrol('Style', 'togglebutton', 'String', 'Exclude','Units','normalize','Position', tempPos*scaleMatExclude,'Callback', {@excludeIC,hObject,it},'Value',any(handles.exclude==it));
end
handles.figIC.buttonLock = buttonLock;
handles.figIC.buttonExclude = buttonExclude;
handles.figIC.handle = fhandle;
guidata(hObject,handles);


function closeFigIC(varargin)
hObject = varargin{3};
% load handles
handles = guidata(hObject);
if isfield(handles,'figIC')
    handles = rmfield(handles,'figIC'); end
guidata(hObject,handles);


function lockIC(varargin)
hObject = varargin{3};
ic = varargin{4};
button = varargin{1};
% load handles
handles = guidata(hObject);
if get(button,'Value') % turned lock on
    handles.lock = sort([handles.lock ic]);
    if any(handles.exclude==ic)
        handles.exclude(handles.exclude==ic) = [];
        % update fig
        if isfield(handles,'figIC')
            set(handles.figIC.buttonExclude(ic),'Value',0); end
    end
else % turned lock off
    handles.lock(handles.lock==ic) = [];
end
% save handles
guidata(hObject,handles);
% update ics to plot
updateICs(hObject)


function excludeIC(varargin)
hObject = varargin{3};
ic = varargin{4};
button = varargin{1};
% load handles
handles = guidata(hObject);
if get(button,'Value') % turned exclude on
    handles.exclude = sort([handles.exclude ic]);
    if any(handles.lock==ic)
        handles.lock(handles.lock==ic) = [];
        % update fig
        if isfield(handles,'figIC')
            set(handles.figIC.buttonLock(ic),'Value',0); end
    end
else % turned exclude off
    handles.exclude(handles.exclude==ic) = [];
end
% save handles
guidata(hObject,handles);
% update ics to plot
updateICs(hObject)


function updateICs(hObject)
handles = guidata(hObject);
temp = [handles.lock setdiff(1:handles.nic,[handles.lock,handles.exclude]) handles.exclude];
handles.ics = temp(1:handles.ntopo);
guidata(hObject,handles);


% function updateButtons(handles,ic,onFlag,lockFlag)
% if lockFlag
%     % update gui
%     ind = find(handles.ics==ic);
%     if ind
%         set(handles.(['togglebuttonLock' int2str(ind)]),'Value',onFlag); end
%     % update fig
%     if isfield(handles,'figIC')
%         set(handles.figIC.buttonLock(ic),'Value',onFlag); end
% else
%     % update gui
%     ind = find(handles.ics==ic);
%     if ind
%         set(handles.(['togglebuttonExclude' int2str(ind)]),'Value',onFlag); end
%     % update fig
%     if isfield(handles,'figIC')
%         set(handles.figIC.buttonExclude(ic),'Value',onFlag); end
% end


% --- Executes on selection change in popupmenuInfo.
function popupmenuInfo_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuInfo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuInfo


% --- Executes during object creation, after setting all properties.
function popupmenuInfo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in togglebuttonExclude8.
function togglebuttonExclude8_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonExclude8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
excludeIC(hObject,[],hObject,handles.ics(8))


% --- Executes on button press in togglebuttonLock8.
function togglebuttonLock8_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lockIC(hObject,[],hObject,handles.ics(8))


% --- Executes on button press in togglebuttonExclude7.
function togglebuttonExclude7_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonExclude7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
excludeIC(hObject,[],hObject,handles.ics(7))


% --- Executes on button press in togglebuttonLock7.
function togglebuttonLock7_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lockIC(hObject,[],hObject,handles.ics(7))


% --- Executes on button press in togglebuttonExclude6.
function togglebuttonExclude6_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonExclude6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
excludeIC(hObject,[],hObject,handles.ics(6))


% --- Executes on button press in togglebuttonLock6.
function togglebuttonLock6_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lockIC(hObject,[],hObject,handles.ics(6))


% --- Executes on button press in togglebuttonExclude5.
function togglebuttonExclude5_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonExclude5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
excludeIC(hObject,[],hObject,handles.ics(5))


% --- Executes on button press in togglebuttonLock5.
function togglebuttonLock5_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lockIC(hObject,[],hObject,handles.ics(5))


% --- Executes on button press in togglebuttonExclude4.
function togglebuttonExclude4_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonExclude4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
excludeIC(hObject,[],hObject,handles.ics(4))


% --- Executes on button press in togglebuttonLock4.
function togglebuttonLock4_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lockIC(hObject,[],hObject,handles.ics(4))


% --- Executes on button press in togglebuttonExclude3.
function togglebuttonExclude3_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonExclude3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
excludeIC(hObject,[],hObject,handles.ics(3))


% --- Executes on button press in togglebuttonLock3.
function togglebuttonLock3_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lockIC(hObject,[],hObject,handles.ics(3))


% --- Executes on button press in togglebuttonExclude2.
function togglebuttonExclude2_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonExclude2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
excludeIC(hObject,[],hObject,handles.ics(2))


% --- Executes on button press in togglebuttonLock2.
function togglebuttonLock2_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lockIC(hObject,[],hObject,handles.ics(2))


% --- Executes on button press in togglebuttonExclude1.
function togglebuttonExclude1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonExclude1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
excludeIC(hObject,[],hObject,handles.ics(1))


% --- Executes on button press in togglebuttonLock1.
function togglebuttonLock1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lockIC(hObject,[],hObject,handles.ics(1))

% --- Executes on mouse press over axes background.
function axesIC1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesIC2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.curIC = handles.ics(1);
% Update handles structure
guidata(hObject, handles);

% --- Executes on mouse press over axes background.
function axesIC2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesIC2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.curIC = handles.ics(2);
% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axesIC3_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesIC3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.curIC = handles.ics(3);
% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axesIC4_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesIC4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.curIC = handles.ics(4);
% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axesIC5_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesIC5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.curIC = handles.ics(5);
% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axesIC6_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesIC6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.curIC = handles.ics(6);
% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axesIC7_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesIC7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.curIC = handles.ics(7);
% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function axesIC8_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesIC8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.curIC = handles.ics(8);
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in pushbuttonPause.
function pushbuttonPause_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonPause (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmpi(get(handles.pushbuttonPause,'string'),'Pause');
    set(handles.pushbuttonPause,'string','Resume');
    stop(handles.pauseTimers)
else
    set(handles.pushbuttonPause,'string','Pause');
    start(handles.pauseTimers);
end


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
temp = timerfindall;
warning off MATLAB:timer:deleterunning
delete(temp)

if isfield(handles,'figIC')
    close(handles.figIC.handle); end
