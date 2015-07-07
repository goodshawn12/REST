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

% Last Modified by GUIDE v2.5 22-Jan-2015 21:14:45

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
end


% --- Executes just before visORICA is made visible.
function visORICA_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to visORICA (see VARARGIN)
%   1: channel locations



% Parse varargsin
handles.chanlocs = evalin('base','chanlocs'); % change !!!
handles.ntopo = 8;
handles.nic = length(handles.chanlocs);
handles.ics = 1:handles.ntopo;
% handles.streamName = 'visORICAst';
handles.curIC = 1;
handles.lock = [];
handles.color_lock = [0.5 1 0.5];
handles.exclude = [];
handles.color_exclude = [1 0.5 0.5];
calibData = varargin{1};

% Check for orica.m and arica.m in bcilab path
startup_check_flt_files();

% Check if localization is possible and adjust GUI accordingly
% handles = startup_check_localization(handles,calibData); % !!! turn back on after testing

% Intialize ORICA
handles = initializeORICA(handles,calibData);

% Create ORICA timer 
oricaTimer = timer('Period',.1,'ExecutionMode','fixedSpacing','TimerFcn',{@onl_filtered_ORICA,handles.streamName},'StartDelay',0.1,'Tag','oricaTimer','Name','oricaTimer');%,'BusyMode','queue');

% Start EEG stream
[~, handles.bufferName] = vis_stream_ORICA('figurehandles',handles.figure1,'axishandles',handles.axisEEG,'streamname',handles.streamName);
eegTimer = timerfind('Name','eegTimer');
set(eegTimer,'UserData',{hObject,1})

% Populate scalp maps
for it = 1:handles.ntopo
    set(handles.figure1, 'CurrentAxes', handles.(['axesIC' int2str(it)]))
    topoplotFast(rand(size(handles.chanlocs)), handles.chanlocs);
end

% Create scalp map timer
topoTimer = timer('Period',round(1/handles.ntopo*1000)/1000,'ExecutionMode','fixedRate','TimerFcn',{@vis_topo,hObject},'StartDelay',0.2,'Tag','topoTimer','Name','topoTimer');

% Create data timer (starts as power spectrum)
infoTimer = timer('Period',1,'ExecutionMode','fixedRate','TimerFcn',{@infoPSD,hObject},'StartDelay',0.2,'Tag','infoTimer','Name','infoTimer');

% Set panel and button colors
handles.color_bg = get(handles.figure1,'Color');
names = fieldnames(handles);
ind = find(any([strncmpi(names,'panel',5),strncmpi(names,'toggle',6),strncmpi(names,'push',4),strncmpi(names,'popup',5)],2));
for it = 1:length(ind)
    set(handles.(names{ind(it)}),'BackgroundColor',handles.color_bg)
end

% Gather and disperse pipeline function names
funs = get_pipeline_functions();
set(handles.popupmenuEEG,'String',[funs; 'ICA cleaned'])
buffer = evalin('base',handles.bufferName);
buffer.funs = funs;
assignin('base',handles.bufferName,buffer);

% Save timers
handles.pauseTimers = [eegTimer,topoTimer,infoTimer];

% Choose default command line output for visORICA
handles.output = hObject;

% Update handles structure
handles.intialized = true;
guidata(hObject, handles);

% Start timers
start(oricaTimer);
start(eegTimer);
start(topoTimer);
start(infoTimer);
end

% UIWAIT makes visORICA wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function [funs] = get_pipeline_functions(p)
if ~exist('p','var')
    p = evalin('base','pipeline'); end
funs = {};
if p.subnodes
    % recursive call for deeper parts of pipeline
    for k=p.subnodes
        [funs] = get_pipeline_functions(p.parts{k}); end
end
% build the outputs
funs = [funs;func2str(p.head)];
end


function startup_check_flt_files
contents = dir(env_translatepath('functions:/filters/flt_*.m'));
flag_refreshBCILAB = false;
if ~any(strcmp({contents.name},'flt_orica.m'))
    display('visORICA: flt_orica.m not present in BCILAB path. Attempting to create copy. This will only happen the first time visEEG is run.')
    destination = env_translatepath('functions:/filters/flt_orica.m');
    source = [fileparts(fileparts(which('visORICA'))) '/flt_orica.m'];
    if ~exist(source,'file')
        warning('visORICA cannot find flt_orica.m. visORICA will likely error soon.')
    else
        [status, message] = copyfile(source,destination);
        if ~status
            warning('visORICA cannot copy flt_orica.m into BCILAB. visORICA will likely error soon.\n%s',message)
        else
            flag_refreshBCILAB = true;
        end
    end
end
if ~any(strcmp({contents.name},'flt_arica.m'))
    display('visORICA: flt_arica.m not present in BCILAB path. Attempting to create copy. This will only happen the first time visEEG is run.')
    destination = env_translatepath('functions:/filters/flt_arica.m');
    source = [fileparts(fileparts(which('visORICA'))) '/flt_arica.m'];
    if ~exist(source,'file')
        warning('visORICA cannot find flt_arica.m. visORICA will likely error soon.')
    else
        [status, message] = copyfile(source,destination);
        if ~status
            warning('visORICA cannot copy flt_arica.m into BCILAB. visORICA will likely error soon.\n%s',message)
        else
            flag_refreshBCILAB = true;
        end
    end
end

if flag_refreshBCILAB
    flt_pipeline('update');
end
end

function handles = startup_check_localization(handles,calibData)

if ~isfield(calibData,'headModel') || isempty(calibData.headModel)
    set(handles.pushbuttonLocalize,'HitTest','off','visible','off')
elseif ischar(calibData.headModel)
    handles.headModel = headModel.loadFromFile(calibData.headModel);
    temp = load(handles.headModel.surfaces);
    handles.nVertices = size(temp.surfData(3).vertices,1);
    [~,handles.K,handles.L,rmIndices] = getSourceSpace4PEB(handles.headModel);
    handles.hmInd = setdiff(1:handles.nVertices,rmIndices);
%     temp = load(handles.headModel.leadFieldFile);
%     handles.K = temp.K;
else
    handles.headModel = calibData.headModel;
    % if calibData does not contain field 'localization' with lead field
    % matrix, laplacian matrix, number of vertices in cortex model, and
    % valid vertex indeces, then calulate them. (takes a while).
    if ~isfield(calibData.localization,'nVertices')
        if ischar(handles.headModel.surfaces)
            temp = load(handles.headModel.surfaces);
            handles.nVertices = size(temp.surfData(3).vertices,1);
        else
            handles.nVertices = size(handles.headModel.surfaces(3).vertices,1);
        end
    else
        handles.nVertices = calibData.localization.nVertices;
    end
    
    if ~isfield(calibData,'localization') || ...
            ~isfield(calibData.localization,'L')
            
        [~,handles.K,handles.L,rmIndices] = ...
            getSourceSpace4PEB(handles.headModel);
        handles.hmInd = setdiff(1:handles.nVertices,rmIndices);
    else
        if ~isfield(calibData.localization,'K')
            temp = load(handles.headModel.leadFieldFile);
            handles.K = temp.K;
        else
            handles.K = calibData.localization.K;
        end
        handles.L = calibData.localization.L;
        if isfield(calibData.localization,'ind') || isfield(calibData.localization,'rmIndices')
            try
                handles.hmInd = calibData.localization.ind;
            catch
                handles.hmInd = setdiff(1:handles.nVertices,calibData.localization.rmIndices);
            end
        else
            handles.hmInd = 1:handles.nVertices;
        end
    end
end
end

function handles = initializeORICA(handles,calibData)

% create/refresh convergence buffers
bufflen = 60; % seconds
assignin('base','conv_statIdx',zeros(1,bufflen*calibData.srate)); % !!! clean output
assignin('base','conv_mir',zeros(1,bufflen*calibData.srate));
assignin('base','learning_rate',zeros(1,bufflen*calibData.srate));

% load LSL
if ~exist('env_translatepath','file')
    % standalone case
    lib = lsl_loadlib();
else
    % if we're within BCILAB we want to make sure that the library is also found if the toolbox is compiled
    lib = lsl_loadlib(env_translatepath('dependencies:/liblsl-Matlab/bin'));
end

% find streams
streams = lsl_resolve_all(lib,3);
streamnames = cellfun(@(s)s.name(),streams ,'UniformOutput',false);
if isempty(streamnames)
    error('There is no stream visible on the network.');
elseif length(streamnames) == 1
    assignin('base','streamname',streamnames); % save stream name in base workspace
else
    % if more than 2 (EEG) streams, pop up a GUI to select one.
    selStream(streamnames); % result 'streamname' is saved in base workspace
    streamnames = evalin('base','streamname');
    streamnames = streamnames{1};
end
run_readlsl('MatlabStream',streamnames,'DataStreamQuery', ['name=''' streamnames '''']);
handles.streamName = streamnames;

opts.lsl.StreamName = handles.streamName;
opts.BCILAB_PipelineConfigFile = 'data/ORICA_pipeline_config_realtime.mat'; % make sure this file doesn't have 'signal' entry

% define the pipeline configuration
try    fltPipCfg = exp_eval(io_load(opts.BCILAB_PipelineConfigFile));
catch, disp('-- no existing pipeline --'); fltPipCfg = {}; end
fltPipCfg = arg_guipanel('Function',@flt_pipeline, ...
    'Parameters',[{'signal',calibData} fltPipCfg], ...
    'PanelOnly',false);

% save the configuration
if ~isempty(fltPipCfg)
    if isfield(fltPipCfg,'signal'); fltPipCfg = rmfield(fltPipCfg,'signal'); end
    save(env_translatepath(opts.BCILAB_PipelineConfigFile),...
        '-struct','fltPipCfg');
end


% grab calib data from online stream
disp('Collecting calibration data from online stream... please wait 20 seconds...');
pause(10); % uh oh!
calibData = onl_peek(opts.lsl.StreamName,10,'seconds');
calibData = warmStartWithBadChRemoved(calibData);
assignin('base','badChIndex',calibData.etc.badChIndex);
assignin('base','badChLabels',calibData.etc.badChLabels);

% run pipline on calibration data
cleaned_data = exp_eval(flt_pipeline('signal',calibData,fltPipCfg));

% initialize the pipeline for streaming data
pipeline     = onl_newpipeline(cleaned_data,opts.lsl.StreamName);
assignin('base','pipeline',pipeline);
end


function infoPSD(varargin)

% plot PSD of selected IC
try
    secs2samp = 5; % seconds
    
    W = evalin('base','pipeline.state.icaweights');
    % if isempty(W), W = evalin('base','Wn'); end
    sphere = evalin('base','pipeline.state.icasphere');
    handles = guidata(varargin{3});
    
    set(handles.panelInfo,'Title',['Power spectral density of IC' int2str(handles.curIC)])
    
    srate = evalin('base',[handles.bufferName '.srate']);
    data = evalin('base',[handles.bufferName '.data{end}']); % !!! make this more robust
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
%     [data,f,conf] = pwelch(data,[],[],[],srate,'ConfidenceLevel',0.95);!!!
    
    plot(handles.axisInfo,f,db(data)) % !!!
    grid(handles.axisInfo,'on');
    xlabel(handles.axisInfo,'Frequency (Hz)')
    ylabel(handles.axisInfo,'Power/Frequency (dB/Hz)')
    axis(handles.axisInfo,'tight')
    set(handles.axisInfo,'XTick',[0 10:10:f(end)])
end
end


function infoConverge(varargin)

% plot convergence statistics

% parse inputs
% handle_statIdx = varargin{3};
handle_learning_rate = varargin{3};
% handle_mir = varargin{4};

% load convergence statistics
% conv_statIdx = evalin('base','conv_statIdx');
conv_mir = evalin('base','conv_mir');
learning_rate = evalin('base','learning_rate');

ylim_mir = minmax(conv_mir);
ylim_lr = minmax(learning_rate);
% set(get(handle_mir,'parent'),'YLim',ylim_mir,'YTick',linspace(ylim_mir(1),ylim_mir(2),5))
set(get(handle_learning_rate,'parent'),'YLim',ylim_lr,'YTick',linspace(ylim_lr(1),ylim_lr(2),5))

% update plots
% set(handle_statIdx,'YData',conv_statIdx)
% set(handle_mir,'YData',conv_mir)
% set(handle_learning_rate,'YData',learning_rate)
set(handle_learning_rate,'YData',learning_rate)



end


function vis_topo(varargin)
% get the updated stream buffer
W = evalin('base','pipeline.state.icaweights');

% get handles
handles = guidata(varargin{3});

% update topo plot
it = mod(get(varargin{1},'TasksExecuted')-1,handles.ntopo)+1;
hstr = ['axesIC' int2str(it)];
hand = get(handles.(hstr),'children');

try
    sphere = evalin('base','pipeline.state.icasphere');
    Winv = inv(W*sphere);
    
    [map, cmin, cmax] = topoplotUpdate(Winv(:,handles.ics(it)), handles.chanlocs,'electrodes','off','gridscale',32);
    surfInd = strcmp(get(hand,'type'),'surface');
    set(hand(surfInd),'CData',map);
    set(handles.(hstr),'CLim',[cmin cmax]);
end

% update name and buttons
lock = any(handles.lock==handles.ics(it));
exclude = any(handles.exclude==handles.ics(it));
set(handles.(['panelIC' int2str(it)]),'Title',['IC' int2str(handles.ics(it))])
set(handles.(['togglebuttonLock' int2str(it)]),'Value',lock,...
    'BackgroundColor',handles.color_lock*lock + handles.color_bg*(1-lock))
set(handles.(['togglebuttonExclude' int2str(it)]),'Value',exclude,...
    'BackgroundColor',handles.color_exclude*exclude + handles.color_bg*(1-exclude))
end


% --- Outputs from this function are returned to the command line.
function varargout = visORICA_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end


% --- Executes on selection change in popupmenuEEG.
function popupmenuEEG_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuEEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
eegTimer = timerfind('Name','eegTimer');
stop(eegTimer)
temp = get(eegTimer,'UserData');
set(eegTimer,'UserData',[temp(1) get(hObject,'value')]);
start(eegTimer)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuEEG contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuEEG
end


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
end


% --- Executes on button press in pushbuttonLocalize.
function pushbuttonLocalize_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonLocalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% if it already exists then just focus on it
if isfield(handles,'figLoc')
    figure(handles.figLoc.handle.hFigure)
    return
end

% otherwise pause if necessary and then create the window
method =  questdlg('Choose method for source Localization','Localization', ...
    'LORETA','Dipole Fit','LORETA');

flag_resume = false;
if strcmpi(get(handles.pushbuttonPause,'string'),'Pause')
    stop(handles.pauseTimers);
    flag_resume = true;
end

switch method
    case 'LORETA'
        figLoc_gen_LORETA(hObject,handles)
    case 'Dipole Fit'
        figLoc_gen_dipolefit(hObject,handles)
end

if flag_resume
    handles = guidata(hObject);
    start(handles.pauseTimers); end

end


function figLoc_gen_LORETA(hObject,handles)
% get ica decomposition
W = evalin('base','pipeline.state.icaweights');
sphere = evalin('base','pipeline.state.icasphere');
Winv = inv(W*sphere);

% run dynamicLoreta once to generate state and initial localization
[U,S,V] = svd(handles.K/handles.L,'econ');
Ut = U';
s2 = diag(S).^2;
iLV = handles.L\V;
options = struct('maxTol',1e-3,'maxIter',100,'gridSize',100,'verbose',false,'history',true,'useGPU',false,'initNoiseFactor',0.001);
[J,sigma2,tau2] = dynamicLoreta(Winv(:,handles.curIC), Ut, s2,iLV,[],[], options);

% create headModel plot
if length(handles.hmInd) == length(J)
    Jest = zeros(handles.nVertices,1);
    Jest(handles.hmInd,:) = J;
else
    Jest = zeros(handles.nVertices,3);
    Jest(handles.hmInd,:) = reshape(J,[],3);
end
fhandle = handles.headModel.plotOnModel(Jest(:),Winv(:,handles.curIC),sprintf('IC %d Localization (LORETA)',handles.curIC));
set(fhandle.hFigure,'DeleteFcn',{@closeFigLoc,hObject},'name',['IC' num2str(handles.curIC)]);


% create timer
locTimer = timer('Period',3,'StartDelay',3,'ExecutionMode','fixedRate','TimerFcn',{@figLoc_update_LORETA,hObject},'Tag','locTimer','Name','locTimer');

% save headModel plot, timer, and dynamicLoreta parameters to handles
handles.pauseTimers = [handles.pauseTimers,locTimer];
handles.figLoc.handle = fhandle;
handles.figLoc.IC = handles.curIC;
handles.figLoc.Ut = Ut;
handles.figLoc.s2 = s2;
handles.figLoc.iLV = iLV;
handles.figLoc.tau2 = tau2;
handles.figLoc.sigma2 = sigma2;
handles.figLoc.options = options;
temp = load(handles.headModel.surfaces);
handles.figLoc.scalp = geometricTools.localGaussianInterpolator(handles.headModel.channelSpace,temp.surfData(1).vertices,32);
if range(handles.figLoc.scalp)<.5
%     warning('visORICA - Localization: Unsure about channel space units. Trying meters instead of millimeters.')
    handles.figLoc.scalp = geometricTools.localGaussianInterpolator(handles.headModel.channelSpace,temp.surfData(1).vertices,32/1000);
end
guidata(hObject,handles);

end


function figLoc_update_LORETA(varargin)
% get ica decomposition
W = evalin('base','pipeline.state.icaweights');
sphere = evalin('base','pipeline.state.icasphere');
Winv = inv(W*sphere);

% parse inputs
handles = guidata(varargin{3});

% run dynamicLoreta
Ut = handles.figLoc.Ut;
s2 = handles.figLoc.s2;
iLV = handles.figLoc.iLV;
tau2 = handles.figLoc.tau2;
sigma2 = handles.figLoc.sigma2;
options = handles.figLoc.options;
[J,sigma2,tau2] = dynamicLoreta(Winv(:,handles.figLoc.IC), Ut, s2,iLV,sigma2,tau2, options);

% update figure and related object values
if handles.nVertices == length(J)
    Jest = J;
    handles.figLoc.handle.sourceMagnitud = Jest;
    set(handles.figLoc.handle.hCortex,'FaceVertexCData',handles.figLoc.handle.sourceMagnitud)
else
    Jest = zeros(handles.nVertices,3);
    Jest(handles.hmInd,:) = reshape(J,[],3);
    handles.figLoc.handle.sourceOrientation = Jest;
    handles.figLoc.handle.sourceMagnitud = squeeze(sqrt(sum(Jest.^2,2)));
    set(handles.figLoc.handle.hVector,'udata',Jest(:,1),'vdata',Jest(:,2),'wdata',Jest(:,3))
    set(handles.figLoc.handle.hCortex,'FaceVertexCData',handles.figLoc.handle.sourceMagnitud)
end
scalp_val = handles.figLoc.scalp*Winv(:,handles.figLoc.IC);
set(handles.figLoc.handle.hScalp,'FaceVertexCData',scalp_val)

% upade color limits
maxabs = max(abs(scalp_val));
handles.figLoc.handle.clim.scalp = [-maxabs maxabs];
maxabs = max(abs(handles.figLoc.handle.sourceMagnitud));
handles.figLoc.handle.clim.source = [-maxabs maxabs];


% save dynamicLoreta parameters to handles
handles.figLoc.tau2 = tau2;
handles.figLoc.sigma2 = sigma2;
guidata(varargin{3},handles);
end

function figLoc_gen_dipolefit(hObject,handles)
% get ica decomposition
W = evalin('base','pipeline.state.icaweights');
sphere = evalin('base','pipeline.state.icasphere');
Winv = inv(W*sphere);

% load surfaces
sd = load(handles.headModel.surfaces);
vertices = sd.surfData(3).vertices(handles.hmInd,:);

% calculate scalp potential transform and guess units
scalp = geometricTools.localGaussianInterpolator(handles.headModel.channelSpace,sd.surfData(1).vertices,32);
if range(scalp)<.5
%     warning('visORICA - Localization: Unsure about channel space units. Trying meters instead of millimeters.')
    scalp = geometricTools.localGaussianInterpolator(handles.headModel.channelSpace,sd.surfData(1).vertices,32/1000);
    Q_location = .01*eye(3)/1000;
else
    Q_location = .01*eye(3);
end

% create figure and dipole plot
fhandle = handles.headModel.plotDipoles([0 0 0],[0 0 0]);
set(fhandle.hFigure,'DeleteFcn',{@closeFigLoc,hObject},'name',['IC' num2str(handles.curIC)]);
hold(fhandle.hAxes,'on');

% initialize bfpf with all vertices active
nParticles = 500;
[dipoles, L, moments, weights, rv, state] = ...
    bfpf(Winv(:,handles.curIC),handles.K,vertices,1,nParticles,Q_location,[],[],1);

% adjust state to only contain the nParticles best particles
[~,ind] = sort(state.weights,2,'descend');
state.weights = state.weights(ind(1:nParticles));
state.L = state.L(ind(1:nParticles));

% set scalp potentials
set(fhandle.hScalp,'FaceVertexCData',scalp*Winv(:,handles.curIC),'facecolor','interp')

if length(moments)/3==length(vertices)
    handles.figLoc.fixed_dip = false;
    moments = reshape(moments,[],3);
    moments = bsxfun(@rdivide,moments,row_pnorm(moments));
    dmnorm = norm(dipoles.moment)/50;
    arrows_dip = quiver3(dipoles.location(1),dipoles.location(2),dipoles.location(3), ...
        dipoles.moment(1)/dmnorm,dipoles.moment(2)/dmnorm,dipoles.moment(3)/dmnorm,'ko','filled');
%     arrows = quiver3(vertices(L,1),vertices(L,2),vertices(L,3), ...
%         weights'.*moments(:,1),weights'.*moments(:,2),weights'.*moments(:,3));
else
    handles.figLoc.fixed_dip = true;
    normals = geometricTools.getSurfaceNormals(sd.surfData(end).vertices,sd.surfData(end).faces,false);
    arrows_dip = quiver3(dipoles.location(1),dipoles.location(2),dipoles.location(3), ...
        dipoles.moment*normals(dipoles.L,1)/50,dipoles.moment*normals(dipoles.L,2)/50,dipoles.moment*normals(dipoles.L,3)/50,'ko','filled');
%     arrows = quiver3(vertices(L,1),vertices(L,2),vertices(L,3), ...
%         weights'.*moments.*normals(L,1),weights'.*moments.*normals(L,2),weights'.*moments.*normals(L,3)); 
end
colormap(bipolar(512, 0.99))
maxabs = max(abs(vec(scalp*Winv(:,handles.curIC))));
caxis(fhandle.hAxes,[-maxabs maxabs]);

% create Residual Variance text
handles.figLoc.axisRV = axes('parent',fhandle.hFigure,'position',[.05 .9 .3 .1],'hittest','off');
axis(handles.figLoc.axisRV,'off')
handles.figLoc.textRV = text(0,0.5,sprintf('Residual Variance: %04.1f%%',rv*100), ...
    'parent',handles.figLoc.axisRV,'fontweight','bold','fontsize',16,'hittest','off');

% create timer
locTimer = timer('Period',3,'StartDelay',3,'ExecutionMode','fixedRate', ...
    'TimerFcn',{@figLoc_update_dipolefit,hObject},'Tag','locTimer','Name','locTimer');

% save headModel plot, timer, and bfpf state to handles
handles.figLoc.state = state;
handles.pauseTimers = [handles.pauseTimers,locTimer];
handles.figLoc.handle = fhandle;
handles.figLoc.IC = handles.curIC;
handles.figLoc.arrows_dip = arrows_dip;
% handles.figLoc.arrows = arrows;
handles.figLoc.scalp = scalp;
handles.figLoc.Q_location = Q_location;
if handles.figLoc.fixed_dip
    handles.figLoc.normals = normals; end
handles.figLoc.nParticles = nParticles;
guidata(hObject,handles);

end


function figLoc_update_dipolefit(varargin)
% get ica decomposition
W = evalin('base','pipeline.state.icaweights');
sphere = evalin('base','pipeline.state.icasphere');
Winv = inv(W*sphere);

% parse inputs
handles = guidata(varargin{3});
Q_location = handles.figLoc.Q_location;
if handles.figLoc.fixed_dip
    normals = handles.figLoc.normals; end
nParticles = handles.figLoc.nParticles;

% load surfaces
temp = load(handles.headModel.surfaces);
vertices = temp.surfData(3).vertices(handles.hmInd,:);

% update bfpf
[dipoles, L, moments, weights, rv, handles.figLoc.state] = ...
    bfpf(Winv(:,handles.figLoc.IC),handles.K,vertices,1,handles.figLoc.state.nParticles,Q_location,[],handles.figLoc.state,1);
if ~handles.figLoc.fixed_dip
    moments = reshape(moments,[],3);
    moments = bsxfun(@rdivide,moments,row_pnorm(moments));
    dmnorm = norm(dipoles.moment)/50;
    set(handles.figLoc.arrows_dip,'XData',dipoles.location(1), ...
                                  'YData',dipoles.location(2), ...
                                  'ZData',dipoles.location(3), ...
                                  'UData',dipoles.moment(1)/dmnorm, ...
                                  'VData',dipoles.moment(2)/dmnorm, ...
                                  'WData',dipoles.moment(3)/dmnorm);
%     set(handles.figLoc.arrows,'XData',vertices(L,1), ...
%                               'YData',vertices(L,2), ...
%                               'ZData',vertices(L,3), ...
%                               'UData',weights'.*moments(:,1), ...
%                               'VData',weights'.*moments(:,2), ...
%                               'WData',weights'.*moments(:,3));
else
    set(handles.figLoc.arrows_dip,'XData',dipoles.location(1), ...
                                  'YData',dipoles.location(2), ...
                                  'ZData',dipoles.location(3), ...
                                  'UData',dipoles.moment*normals(dipoles.L,1)/50, ...
                                  'VData',dipoles.moment*normals(dipoles.L,2)/50, ...
                                  'WData',dipoles.moment*normals(dipoles.L,3)/50);
%     set(handles.figLoc.arrows,'XData',vertices(L,1), ...
%                               'YData',vertices(L,2), ...
%                               'ZData',vertices(L,3), ...
%                               'UData',weights'.*moments.*normals(L,1), ...
%                               'VData',weights'.*moments.*normals(L,2), ...
%                               'WData',weights'.*moments.*normals(L,3));
    
end
set(handles.figLoc.handle.hScalp,'FaceVertexCData', ...
    handles.figLoc.scalp*Winv(:,handles.figLoc.IC))
maxabs = max(abs(vec(handles.figLoc.scalp*Winv(:,handles.figLoc.IC))));
caxis(handles.figLoc.handle.hAxes,[-maxabs maxabs]);

% update RV text
set(handles.figLoc.textRV,'string',sprintf('Residual Variance: %04.1f%%',rv*100));

end


function closeFigLoc(varargin)
hObject = varargin{3};
% load handles
handles = guidata(hObject);
% delete figure handle from handles
if isfield(handles,'figLoc')
    handles = rmfield(handles,'figLoc'); end
% delete timer and remove from pauseTimers
locTimerInd = strcmp(get(handles.pauseTimers,'Name'),'locTimer');
delete(handles.pauseTimers(locTimerInd));
handles.pauseTimers(locTimerInd) = [];
% save handles
guidata(hObject,handles);
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
end


function genICSelectGUI(hObject,handles)
temp = get(handles.figure1,'Position');
fhandle = figure('toolbar','none','Menubar','none','Name','IC Select','position',[1 1 temp(3:4)],'Resize','on','DeleteFcn',{@closeFigIC,hObject});
W = evalin('base','pipeline.state.icaweights');
sphere = evalin('base','pipeline.state.icasphere');
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
    
    lock = any(handles.lock==it);
    exclude = any(handles.exclude==it);
    buttonLock(it) = uicontrol('Style', 'togglebutton', 'String', 'Lock',...
        'Units','normalize','Position', tempPos.*[1 1 .5-buttonGap/2 .2],...
        'Callback', {@lockIC,it,hObject},'Value',lock,...
        'BackgroundColor',handles.color_lock*lock + handles.color_bg*(1-lock));
    buttonExclude(it) = uicontrol('Style', 'togglebutton', 'String', 'Exclude',...
        'Units','normalize','Position', tempPos*scaleMatExclude,...
        'Callback', {@excludeIC,it,hObject},'Value',exclude,...
        'BackgroundColor',handles.color_exclude*exclude + handles.color_bg*(1-exclude));
end
handles.figIC.buttonLock = buttonLock;
handles.figIC.buttonExclude = buttonExclude;
handles.figIC.handle = fhandle;
guidata(hObject,handles);
end


function closeFigIC(varargin)
hObject = varargin{3};
% load handles
handles = guidata(hObject);
if isfield(handles,'figIC')
    handles = rmfield(handles,'figIC'); end
guidata(hObject,handles);
end


function lockIC(varargin)
ic = varargin{3};
button = varargin{1};
% load handles
if numel(varargin)>3
    hObject = varargin{4};
else
    hObject = get(button,'parent');
end
handles = guidata(hObject);
if get(button,'Value') % turned lock on
    handles.lock = sort([handles.lock ic]);
%     set(button,'BackgroundColor',[0.5 1 0.5])
    if isfield(handles,'figIC')
        set(handles.figIC.buttonLock(ic),'Value',1,'BackgroundColor',handles.color_lock); end
%     if any(handles.exclude==ic)
%         handles.exclude(handles.exclude==ic) = [];
%         % update fig
%         if isfield(handles,'figIC')
%             set(handles.figIC.buttonExclude(ic),'Value',0,...
%                 'BackgroundColor',handles.color_bg);
%         end
%     end
else % turned lock off
    handles.lock(handles.lock==ic) = [];
    if isfield(handles,'figIC')
        set(handles.figIC.buttonLock(ic),'value',0,'BackgroundColor',handles.color_bg); end
end
% save handles
guidata(hObject,handles);
% update ics to plot
updateICs(hObject)
end


function excludeIC(varargin)
ic = varargin{3};
button = varargin{1};
% load handles
if numel(varargin)>3
    hObject = varargin{4};
else
    hObject = get(button,'parent');
end
handles = guidata(hObject);
if get(button,'Value') % turned exclude on
    handles.exclude = sort([handles.exclude ic]);
%     set(button,'BackgroundColor',[1 0.5 0.5])
    if isfield(handles,'figIC')
        set(handles.figIC.buttonExclude(ic),'Value',1,'BackgroundColor',handles.color_exclude); end
%     if any(handles.lock==ic)
%         handles.lock(handles.lock==ic) = [];
%         % update fig
%         if isfield(handles,'figIC')
%             set(handles.figIC.buttonLock(ic),'Value',0,...
%                 'BackgroundColor',handles.color_bg);
%         end
%     end
else % turned exclude off
    handles.exclude(handles.exclude==ic) = [];
    if isfield(handles,'figIC')
        set(handles.figIC.buttonExclude(ic),'value',0,'BackgroundColor',handles.color_bg); end
end
% save handles
guidata(hObject,handles);
% update ics to plot
updateICs(hObject)
end


function updateICs(hObject)
handles = guidata(hObject);
temp = [handles.lock setdiff(1:handles.nic,[handles.lock])];
% temp = [handles.lock setdiff(1:handles.nic,[handles.lock,handles.exclude]) handles.exclude];
handles.ics = temp(1:handles.ntopo);
guidata(hObject,handles);
end



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

infoTimer = timerfind('name','infoTimer');
timerFcn = subsref(get(infoTimer,'TimerFcn'), substruct('{}',{1}));

contents = get(hObject,'String');
switch contents{get(handles.popupmenuInfo,'Value')}
    case 'Power Spectrum'
        % changed?
        if isequal(timerFcn,@infoPSD)
            return
        end
        % if so...
        stop(infoTimer)
        set(infoTimer,'Period',1,'ExecutionMode','fixedRate','TimerFcn',{@infoPSD,handles.figure1},'StartDelay',0);
        handles.axisInfo = handles.axisInfo(1);
        start(infoTimer)
    case 'Convergence'
        % changed?
        if isequal(timerFcn,@infoConverge)
            return
        end
        % if so...
        stop(infoTimer)
%         conv_statIdx = evalin('base','pipeline.state.statIdx');
%         conv_statIdx = evalin('base','conv_statIdx');
        conv_mir = evalin('base','conv_mir');
        learning_rate = evalin('base','learning_rate');
        srate = evalin('base',[handles.streamName '.srate']);
        x = -(length(conv_mir)-1)/srate:1/srate:0;
        axes(handles.axisInfo)
        [handles.axisInfo] = plot(x,learning_rate);
        line1 = get(handles.axisInfo,'children');
%         [handles.axisInfo,line1,line2] = plotyy(x,learning_rate,x,conv_mir);
%         [handles.axisInfo,line1,line2] = plotyy(x,conv_statIdx,x,conv_mir);
        set(get(get(handles.axisInfo(1),'parent'),'XLabel'),'String','Time (seconds)')
        set(get(get(handles.axisInfo(1),'parent'),'YLabel'),'String','Learning Rate (dB)')
%         set(get(handles.axisInfo(1),'YLabel'),'String','Convergence Index')
%         set(get(handles.axisInfo(2),'YLabel'),'String','Mutual Information Reduction')
        set(infoTimer,'Period',1,'ExecutionMode','fixedRate','TimerFcn',{@infoConverge,handles.axisInfo},'StartDelay',0);
%         set(infoTimer,'Period',1,'ExecutionMode','fixedRate','TimerFcn',{@infoConverge,line1,line2},'StartDelay',0);
        axis(get(handles.axisInfo,'parent'),'tight')
        start(infoTimer)
    otherwise
        warning('visORICA: popupmenuInfo recieved a strange input')
end


% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuInfo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuInfo
end


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
end


% --- Executes on button press in togglebuttonExclude8.
function togglebuttonExclude8_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonExclude8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
excludeIC(hObject,[],handles.ics(8))
end


% --- Executes on button press in togglebuttonLock8.
function togglebuttonLock8_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lockIC(hObject,[],handles.ics(8))
end


% --- Executes on button press in togglebuttonExclude7.
function togglebuttonExclude7_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonExclude7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
excludeIC(hObject,[],handles.ics(7))
end


% --- Executes on button press in togglebuttonLock7.
function togglebuttonLock7_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lockIC(hObject,[],handles.ics(7))
end


% --- Executes on button press in togglebuttonExclude6.
function togglebuttonExclude6_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonExclude6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
excludeIC(hObject,[],handles.ics(6))
end


% --- Executes on button press in togglebuttonLock6.
function togglebuttonLock6_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lockIC(hObject,[],handles.ics(6))
end


% --- Executes on button press in togglebuttonExclude5.
function togglebuttonExclude5_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonExclude5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
excludeIC(hObject,[],handles.ics(5))
end


% --- Executes on button press in togglebuttonLock5.
function togglebuttonLock5_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lockIC(hObject,[],handles.ics(5))
end


% --- Executes on button press in togglebuttonExclude4.
function togglebuttonExclude4_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonExclude4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
excludeIC(hObject,[],handles.ics(4))
end


% --- Executes on button press in togglebuttonLock4.
function togglebuttonLock4_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lockIC(hObject,[],handles.ics(4))
end


% --- Executes on button press in togglebuttonExclude3.
function togglebuttonExclude3_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonExclude3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
excludeIC(hObject,[],handles.ics(3))
end


% --- Executes on button press in togglebuttonLock3.
function togglebuttonLock3_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lockIC(hObject,[],handles.ics(3))
end


% --- Executes on button press in togglebuttonExclude2.
function togglebuttonExclude2_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonExclude2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
excludeIC(hObject,[],handles.ics(2))
end


% --- Executes on button press in togglebuttonLock2.
function togglebuttonLock2_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lockIC(hObject,[],handles.ics(2))
end


% --- Executes on button press in togglebuttonExclude1.
function togglebuttonExclude1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonExclude1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
excludeIC(hObject,[],handles.ics(1))
end


% --- Executes on button press in togglebuttonLock1.
function togglebuttonLock1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lockIC(hObject,[],handles.ics(1))
end

% --- Executes on mouse press over axes background.
function axesIC1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesIC2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.curIC = handles.ics(1);
% Update handles structure
guidata(hObject, handles);
end

% --- Executes on mouse press over axes background.
function axesIC2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesIC2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.curIC = handles.ics(2);
% Update handles structure
guidata(hObject, handles);
end


% --- Executes on mouse press over axes background.
function axesIC3_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesIC3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.curIC = handles.ics(3);
% Update handles structure
guidata(hObject, handles);
end


% --- Executes on mouse press over axes background.
function axesIC4_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesIC4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.curIC = handles.ics(4);
% Update handles structure
guidata(hObject, handles);
end


% --- Executes on mouse press over axes background.
function axesIC5_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesIC5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.curIC = handles.ics(5);
% Update handles structure
guidata(hObject, handles);
end


% --- Executes on mouse press over axes background.
function axesIC6_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesIC6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.curIC = handles.ics(6);
% Update handles structure
guidata(hObject, handles);
end


% --- Executes on mouse press over axes background.
function axesIC7_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesIC7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.curIC = handles.ics(7);
% Update handles structure
guidata(hObject, handles);
end


% --- Executes on mouse press over axes background.
function axesIC8_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesIC8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.curIC = handles.ics(8);
% Update handles structure
guidata(hObject, handles);
end


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
end


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'intialized') && handles.intialized
    if isfield(handles,'figIC')
        try
            close(handles.figIC.handle); end, end
    if isfield(handles,'figLoc')
        try
            close(handles.figIC.handle); end, end

    timerNames = {'eegTimer','oricaTimer','topoTimer','infoTimer','locTimer',[handles.streamName '_timer']};
    % warning off MATLAB:timer:deleterunning
    for it = 1:length(timerNames)
        temp = timerfind('name',timerNames{it});
        if isempty(temp)
            continue; end
        stop(temp)
        delete(temp)
    end
end
end

