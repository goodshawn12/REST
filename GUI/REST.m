function varargout = REST(varargin)
% REST MATLAB code for REST.fig
%
%      REST(opts) starts the Real-time EEG Source-mapping Toolbox.
%       - OPTS is a structure with 2 possible fields:
%           - chanlocs (possibly required): EEGALAB channel locations
%             structure. This information could also be present in
%             calibration_data or in headModel, in which case chanlocs is
%             not required.
%           - calibration_data (optional): EEGLAB EEG structure with some
%             initaial data or a string containing the path to the *.set
%             file containing the calibration data (optional).
%           - headModel (optional): MoBILAB headModel object or path to 
%             saved headModel file. headModels be made using the included 
%             function make_headModel. Required to do source localization.
%             - K (optional): if headModel is an object, then it is 
%               advisable to attach the cropped lead-field matrix
%             - L (optional): if headModel is an object, then it is 
%               advisable to attach the cropped laplacian matrix
%             - hmInd (optional): if headModel is an object, then it is 
%               advisable to attach the remaining indices used in the 
%               headModel
%           - customize_pipeline (optional): a value of "true" will bring
%             up a dialog for customizing the pipeline. Be warned: the
%             pipeline must end with flt_orica or flt_arica and may be
%             brittle to other changes as well. This is largely untested
%             but provided as an option for the brave who wish to
%             experiment (default: false)
%           - playback (optional): a value of "true" will cause REST to
%             stream the calibration data (looping) and treat it as a
%             real-time data stream for testing purposes. (default: false)
%           - calibration_window (optional): if you do not want to use all
%             the provided calibration data to initialze the online ICA
%             pipeline, you may use this option to shorten it. A scalar
%             will take that many seconds from the recording and throw away
%             the rest. A vector with two elements defines a start and stop
%             time (in seconds) to take from the data. This is usefule for
%             playback as you can have a large playback dataset and a small
%             calibration dataset.
%
%      REST('Property','Value',...) creates a new REST or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before REST_OpeningFcn gets called. An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to REST_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @REST_OpeningFcn, ...
    'gui_OutputFcn',  @REST_OutputFcn, ...
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

%% Setup

% --- Executes just before REST is made visible.
function REST_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to REST (see VARARGIN)

% Set parameters
handles.ntopo = 8;
handles.curIC = 1;
handles.lock = [];
handles.color_lock = [0.5 1 0.5];
handles.reject = [];
handles.color_reject = [1 0.5 0.5];

% Parse varagin
[handles,calibData,config] = startup_check_inputs(handles,varargin);
    
% Intialize ORICA
handles = startup_initializeORICA(handles,calibData,config);

% Create ORICA timer 
oricaTimer = timer('Period',.1,'ExecutionMode','fixedSpacing','TimerFcn',{@onl_filtered_ORICA,parseStreamName(handles.streamName)},'StartDelay',0.1,'Tag','oricaTimer','Name','oricaTimer');%,'BusyMode','queue');

% Start EEG stream
[~, handles.bufferName] = vis_stream_ORICA('figurehandles',handles.figure1,'axishandles',handles.axisEEG,'streamname',handles.streamName);
eegTimer = timerfind('Name','eegTimer');
set(eegTimer,'UserData',{hObject,1})

% Gather and disperse pipeline function names
funs = get_pipeline_functions();
funnames = cell(length(funs)-2,1);
for it = 2:length(funs)
    temp = arg_report('properties',funs{it});
    funnames{it-1} = temp.name;
    if iscell(funnames{it-1})
        funnames{it-1} = funnames{it-1}{1}; end
end
set(handles.popupmenuEEG,'String',['Raw Data'; funnames; 'ICA Cleaned'])
buffer = evalin('base',handles.bufferName);
buffer.funs = funs;
assignin('base',handles.bufferName,buffer);

% Find if channels have been removed
funsstr = cellfun(@func2str,funs,'uniformoutput',false');
if any(strcmp(funsstr,'flt_selchans'))
    pipeline = evalin('base','pipeline'); %#ok<NASGU>
    numparts = find(flipud(strcmp(funsstr,'flt_selchans')))-1;
    eval(['removed = pipeline' repmat('.parts{2}',1,numparts) '.parts{4};']);
    handles.rmchan_index = ismember({handles.chanlocs.labels},removed);
    % adjust chanlocs
    handles.urchanlocs = handles.chanlocs;
    handles.chanlocs(handles.rmchan_index) = [];
    handles.nic = length(handles.chanlocs);
    handles.ics = 1:handles.nic;
    
    if isfield(handles,'headModel')
        % adjust headModel
        handles.urheadModel = handles.headModel;
        handles.headModel.dropChannels(handles.rmchan_index); % !!! had to change the headModel contructor
        handles.K(handles.rmchan_index,:) = [];
        %     LFM = load(handles.headModel.leadFieldFile);
        %     LFM.K(handles.rmchan_index,:) = [];
        %     save(handles.headModel.leadFieldFile,'-struct','LFM')
    end
end

% Populate scalp maps
for it = 1:handles.ntopo
    set(handles.figure1, 'CurrentAxes', handles.(['axesIC' int2str(it)]))
    [~,Zi,~,Xi,Yi,intx,inty] = topoplotFast_LRBF(rand(size(handles.chanlocs)), handles.chanlocs);
end

% Generate scalp map interpolation matrix (jerry rigged)
nChan = length(handles.chanlocs);
in = eye(nChan);
out = zeros(32^2,nChan);
for it = 1:nChan
    op = rbfcreate(double([inty;intx]),in(:,it)','RBFFunction', 'linear');
    out(:,it) = rbfinterp(double([Xi(:),Yi(:)]'), op);
end
handles.topoMat = out/in;
handles.topoNaNMask = isnan(Zi);
handles.topoNPixel = size(out,1);
handles.topoMat(handles.topoNaNMask,:) = [];

% Create scalp map timer
topoTimer = timer('Period',round(1/handles.ntopo*1000)/1000,'ExecutionMode','fixedRate','TimerFcn',{@vis_topo,hObject},'StartDelay',0.2,'Tag','topoTimer','Name','topoTimer');

% Create data timer (starts as power spectrum)
infoTimer = timer('Period',0.25,'ExecutionMode','fixedRate','TimerFcn',{@infoPSD,hObject}, ...
    'StartDelay',0.2,'Tag','infoTimer','Name','infoTimer');

% PSD settings
sec2samp = 5;
sec4fft = 0.5;
overlap = 0.5;
fmax = 45;
nwindows = floor((sec2samp / sec4fft - 1) / (1 - overlap)) + 1;
handles.psd = struct('curIC', 1, 'ffts', zeros(fmax, 0), 'inds', [], ...
    'sec2samp', sec2samp, 'sec4fft', sec4fft, 'overlap', overlap, 'fmax', fmax, ...
    'nwindows', nwindows, 'line', []);

% Set panel and button colors
handles.color_bg = get(handles.figure1,'Color');
names = fieldnames(handles);
ind = find(any([strncmpi(names,'panel',5),strncmpi(names,'toggle',6),strncmpi(names,'push',4),strncmpi(names,'popup',5)],2));
for it = 1:length(ind)
    set(handles.(names{ind(it)}),'BackgroundColor',handles.color_bg)
end

% Save timers
handles.pauseTimers = [eegTimer,topoTimer,infoTimer];

% Choose default command line output for REST
handles.output = hObject;

% Update handles structure
handles.intialized = true;
guidata(hObject, handles);

% % for testing
% set(handles.figure1, 'visible', 'on')
% set(handles.pushbuttonSortICs, 'callback',  ...
%     @(hObject,eventdata) REST('pushbuttonLSLOut_Callback', ...
%                               hObject, eventdata, guidata(hObject)))

% Start timers
start(oricaTimer);
start(eegTimer);
start(topoTimer);
start(infoTimer);
end


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
funs = [funs;{p.head}];
end


function [handles, calibData, config] = startup_check_inputs(handles,in)
% !!! need to add appropriate errors and explanations

% check channel locations
if isfield(in{1},'chanlocs')
    handles.chanlocs = in{1}.chanlocs; end

% check calibration data
if isfield(in{1},'calibration_data')
    if isstruct(in{1}.calibration_data)
        calibData = in{1}.calibration_data;
    elseif ischar(in{1}.calibration_data)
        calibData = pop_loadset(in{1}.calibration_data);
    end
    if ~isfield(handles,'chanlocs') && isfield(calibData,'chanlocs') && ~isempty(calibData.chanlocs)
        handles.chanlocs = calibData.chanlocs; end
else
    calibData = [];
end

% check if playback is requested
if isfield(in{1},'playback') && in{1}.playback
    OMIT_MARKERS = true;
    % find existing streams
    lib = lsl_loadlib(env_translatepath('dependencies:/liblsl-Matlab/bin'));
    streams = lsl_resolve_all(lib,0.1);
    streamnames = cellfun(@(s)s.name(),streams ,'UniformOutput',false)';
    ind = 1;
    while any(strcmp(['REST_playback_data' num2str(ind)], streamnames)) ...
            || any(strcmp(['REST_playback_markers' num2str(ind)], streamnames))
        ind = ind + 1; end
    % create playback stream
    handles.playbackStream = play_eegset_lsl('Dataset', calibData, ...
        'DataStreamName', ['REST_playback_data' num2str(ind)], ...
        'EventStreamName', ['REST_playback_markers' num2str(ind)], ...
        'Background', true, 'NoMarkers', OMIT_MARKERS);
    handles.streamName = ['REST_playback_data' num2str(ind)];
end

% shorten calibration data if requested
if isfield(in{1},'calibration_window')
    if isscalar(in{1}.calibration_window)
        calibData = pop_select(calibData,'time',[0 in{1}.calibration_window]);
    else
        calibData = pop_select(calibData,'time',in{1}.calibration_window);
    end
end

% check if localization is possible and adjust GUI accordingly
handles = startup_check_localization(handles,in{1});

if ~isfield(handles,'chanlocs') && isfield(handles,'headModel')
    for it = 1:size(handles.headModel.channelSpace,1)
        handles.chanlocs(it).labels = handles.headModel.channelLabel{it};
        handles.chanlocs(it).X = handles.headModel.channelSpace(it,1);
        handles.chanlocs(it).Y = handles.headModel.channelSpace(it,2);
        handles.chanlocs(it).Z = handles.headModel.channelSpace(it,3);
    end
    handles.chanlocs = convertlocs(handles.chanlocs);
end

% if still no chanlocs, error
if ~isfield(handles,'chanlocs')
    error('REST: No channel location information provided!')
end

% check whether to open pipeline arg_guipanel
if isfield(in{1},'customize_pipeline')
	config.customize_pipeline = in{1}.customize_pipeline;
else
	config.customize_pipelin = false;
end

% check whether to save the pipeline configuration
if isfield(in{1},'save_config')
	config.save_config = in{1}.save_config;
else
	config.save_config = false;
end

% check if config file is defined
if isfield(in{1},'config')
    handles.config = in{1}.config;
else
    handles.config = 'Config_ORICA';
end
end


function handles = startup_check_localization(handles,in) % !!! combine headmodel in localization
% if no headModel provided, remove localization button
if ~isfield(in,'headModel') || isempty(in.headModel)
    set(handles.pushbuttonLocalize,'HitTest','off','visible','off')
    
% if the provided headModel is a string, load the headModel
elseif isa(in.headModel,'char')
    handles.headModel = headModel.loadFromFile(in.headModel);
    temp = load(handles.headModel.surfaces);
    handles.nVertices = size(temp.surfData(3).vertices,1);
    
    % if there is not an accompanying *_SSPEB.mat file containing the
    % cropped lead-field matrix, laplacian matrix, and valid vertex indices,
    % then calulate them. (takes a while).
    if ~exist([in.headModel '_SSPEB.mat'],'file')
        [~,K,L,rmIndices] = getSourceSpace4PEB(handles.headModel);
        hmInd = setdiff(1:handles.nVertices,rmIndices);
        save([in.headModel '_SSPEB.mat'],'K','L','hmInd')

        handles.K = K;
        handles.L = L;
        handles.hmInd = hmInd;
    else
        temp = load([in.headModel '_SSPEB.mat']);
        handles.K = temp.K;
        handles.L = temp.L;
        handles.hmInd = temp.hmInd;
    end
    
% if the provided headModel is an object, use it
elseif isa(in.headModel,'headModel')
    handles.headModel = in.headModel;
    
    % find the number of vertices in the model
    if ischar(handles.headModel.surfaces)
        temp = load(handles.headModel.surfaces);
        handles.nVertices = size(temp.surfData(3).vertices,1);
    else
        handles.nVertices = size(handles.headModel.surfaces(3).vertices,1);
    end
    
    % if K, L, and hmInd are not provided in opts, then calulate them.
    if ~isfield(in,'K') || ~isfield(in,'L') || ~isfield(in,'hmInd')
        [~,handles.K,handles.L,rmIndices] = ...
            getSourceSpace4PEB(handles.headModel);
        handles.hmInd = setdiff(1:handles.nVertices,rmIndices);
    else
        handles.K = in.K;
        handles.L = in.L;
        handles.hmInd = in.hmInd;
    end
end
end


function handles = startup_initializeORICA(handles,calibData,config)

% load LSL
if ~exist('env_translatepath','file')
    % standalone case
    lib = lsl_loadlib();
else
    % if we're within BCILAB we want to make sure that the library is also found if the toolbox is compiled
    lib = lsl_loadlib(env_translatepath('dependencies:/liblsl-Matlab/bin'));
end

% find streams
if ~isfield(handles, 'streamName')
    streams = lsl_resolve_all(lib,0.1);
    streamnames = cellfun(@(s)s.name(),streams ,'UniformOutput',false);
    if isempty(streamnames)
        error('There is no stream visible on the network.');
    elseif length(streamnames) == 1
        assignin('base','streamname',streamnames); % save stream name in base workspace
    else
        % if more than 2 (EEG) streams, pop up a GUI to select one.
    %     selStream(streamnames); % result 'streamname' is saved in base workspace
    %     streamnames = evalin('base','streamname');
    %     streamnames = streamnames{1};
        [streamind, ok] = listdlg('ListString', streamnames, ...
            'SelectionMode', 'Single', 'PromptString', 'Select which LSL stream to use.');
        assert(ok && ~isempty(streamind), 'No LSL stream selection was made.')
        streamnames = streamnames{streamind};
    end
    handles.streamName = streamnames;
end
lslin(handles)
% run_readlsl('MatlabStream', handles.streamName, ...
%     'DataStreamQuery', ['name=''' handles.streamName ''''], ...
%     'MarkerStreamQuery', '');
opts.lsl.StreamName = parseStreamName(handles.streamName);

% create learning rate buffer
bufflen = 60; % seconds
handles.srate = getfield(onl_peek(opts.lsl.StreamName,1,'samples'),'srate');
assignin('base','learning_rate',nan(1,bufflen*handles.srate));

% look for pre-existing config file for pipeline
REST_path = fileparts(fileparts(which('REST')));
opts.BCILAB_PipelineConfigFile = ...
    [REST_path filesep 'data' filesep 'config' filesep handles.config '.mat']; % make sure this file doesn't have 'signal' entry

% define the pipeline configuration
tic
try
    fltPipCfg = exp_eval(io_load(opts.BCILAB_PipelineConfigFile));
    if config.customize_pipeline
        fltPipCfg = arg_guipanel('Function',@flt_pipeline, ...
            'Parameters',[{'signal',onl_peek(opts.lsl.StreamName,1,'samples')} fltPipCfg], ...
            'PanelOnly',false);
    end
catch
    disp('-- no existing pipeline or fail loading pipeline--'); 
    fltPipCfg = {};
end

% open pipeline configuration gui if no settings found or if user requested
if isempty(fltPipCfg)
    fltPipCfg = arg_guipanel('Function',@flt_pipeline, ...
        'Parameters',[{'signal',onl_peek(opts.lsl.StreamName,1,'samples')} fltPipCfg], ...
        'PanelOnly',false);
end

if isfield(fltPipCfg,'pselchans')
    if isfield(calibData.etc,'badChLabels')
        fltPipCfg.pselchans.channels = calibData.etc.badChLabels;
    end
end

% save the configuration %!!! maybe disable this?
if config.save_config
    if ~isempty(fltPipCfg)
        if isfield(fltPipCfg,'signal')
            fltPipCfg = rmfield(fltPipCfg,'signal'); end
        save(env_translatepath(opts.BCILAB_PipelineConfigFile), ...
            '-struct','fltPipCfg');
    end
end

% grab calib data from online stream if there is none
if isempty(calibData)
disp('Collecting calibration data from online stream... please wait 10 seconds...');
pause(10-toc); % uh oh!
calibData = onl_peek(opts.lsl.StreamName,10,'seconds');
end

% check for bad channels
calibData = warmStartWithBadChRemoved(calibData);

% run pipline on calibration data
cleaned_data = exp_eval(flt_pipeline('signal',calibData,fltPipCfg));

% initialize the pipeline for streaming data
pipeline     = onl_newpipeline(cleaned_data,opts.lsl.StreamName);
assignin('base','pipeline',pipeline);
end

%% timer functions

% plot PSD of selected IC
function infoPSD(varargin)
try %#ok<*TRYNC>
% load handles
handles = guidata(varargin{3});

% load buffer
buffer = evalin('base',handles.bufferName);

% clear old or invalid fft estimates
if handles.curIC ~= handles.psd.curIC
    handles.psd.ffts = [];
    handles.psd.inds = [];
    handles.psd.curIC = handles.curIC;
    
    % set string to correct IC number
    set(handles.panelInfo,'Title',['Power spectral density of IC' int2str(handles.curIC)])
else
    delete_ind = handles.psd.inds < buffer.smax - handles.psd.sec2samp * buffer.srate;
    handles.psd.ffts(:, delete_ind) = [];
    handles.psd.inds(:, delete_ind) = [];
end

% get indices for windows
winlen = round(buffer.srate * handles.psd.sec4fft);
winspacing = round(buffer.srate * handles.psd.sec4fft * (1 - handles.psd.overlap));
if isempty(handles.psd.inds)
    % calculate all fft windows starting at current time
    ind =  max(1, buffer.smax - buffer.srate * handles.psd.sec2samp + 1):winspacing:max(1, buffer.smax + 1 - winlen);
else
    % calculate missing windows
    lastind = max(handles.psd.inds);
    ind = lastind + winspacing:winspacing:buffer.smax - winlen + 1;
end

if ~isempty(ind)
    % load ICA matrices
    W = evalin('base','pipeline.state.icaweights');
    sphere = evalin('base','pipeline.state.icasphere');
    
    % calculate ffts
    window = hanning(winlen);
    fftest = arrayfun(@(val) (W(handles.psd.curIC,:) * sphere * buffer.data{end}(:, mod((val:val + winlen - 1) - 1, buffer.pnts) + 1)), ...
        ind, 'uniformoutput', 0);
    fftest = fft(bsxfun(@times, cat(1, fftest{:})', window), buffer.srate);
    fftest = abs(fftest(2:min(ceil((buffer.srate + 1) / 2), handles.psd.fmax + 1), :));

    handles.psd.ffts = [handles.psd.ffts fftest];
    handles.psd.inds = [handles.psd.inds ind];
    
    % update plot
    if isempty(handles.psd.line) || ~isgraphics(handles.psd.line)
        handles.psd.line = plot(handles.axisInfo, ...
            1:size(handles.psd.ffts, 1), db(mean(handles.psd.ffts, 2)));
        grid(handles.axisInfo,'on');
        xlabel(handles.axisInfo,'Frequency (Hz)')
        ylabel(handles.axisInfo,'Power/Frequency (dB/Hz)')
        axis(handles.axisInfo,'tight')
        set(handles.axisInfo,'XTick',[0 10:10:size(handles.psd.ffts, 1)])
    else
        set(handles.psd.line, 'YData', db(mean(handles.psd.ffts, 2)))
        axis(handles.axisInfo,'tight')
        set(handles.axisInfo,'XTick',[0 10:10:size(handles.psd.ffts, 1)])
    end
    
    % save handles
    guidata(varargin{3}, handles)
end
end
end


function infoConverge(varargin)

% plot convergence statistics

% parse inputs
handle_learning_rate = varargin{3};
handle_panelInfo = varargin{4};

% load convergence statistics
learning_rate = evalin('base','learning_rate');

ylim_lr = minmax(learning_rate);
set(get(handle_learning_rate,'parent'),'YLim',ylim_lr,'YTick',linspace(ylim_lr(1),ylim_lr(2),5))

% update plots
set(handle_learning_rate,'YData',learning_rate)
set(handle_panelInfo,'Title','Adaptive Learning Rate')
end


function vis_topo(varargin)
% get the updated stream buffer
W = evalin('base','pipeline.state.icaweights');
sphere = evalin('base','pipeline.state.icasphere');

% get handles
handles = guidata(varargin{3});

% update topo plot
it = mod(get(varargin{1},'TasksExecuted')-1,handles.ntopo)+1;
% if it==1
%     updateICs(varargin{3}); end
hstr = ['axesIC' int2str(it)];
hand = get(handles.(hstr),'children');

try
    Winv = inv(W*sphere);
    
%     [map, cmin, cmax] = topoplotUpdate(Winv(:,handles.ics(it)), handles.chanlocs,'electrodes','off','gridscale',32);
    map = zeros(handles.topoNPixel,1);
    map(~handles.topoNaNMask) = handles.topoMat*Winv(:,handles.ics(it));
    maxabs = max(abs(map));
    cmin = -maxabs;
    cmax =  maxabs;
    map(handles.topoNaNMask) = NaN;
    map = reshape(map,sqrt(handles.topoNPixel),[]);
    
    surfInd = strcmp(get(hand,'type'),'surface');
    set(hand(surfInd),'CData',map);
    set(handles.(hstr),'CLim',[cmin cmax]);
end

% update name and buttons
lock = any(handles.lock==handles.ics(it));
reject = any(handles.reject==handles.ics(it));
set(handles.(['panelIC' int2str(it)]),'Title',['IC' int2str(handles.ics(it))])
set(handles.(['togglebuttonLock' int2str(it)]),'Value',lock,...
    'BackgroundColor',handles.color_lock*lock + handles.color_bg*(1-lock))
set(handles.(['togglebuttonReject' int2str(it)]),'Value',reject,...
    'BackgroundColor',handles.color_reject*reject + handles.color_bg*(1-reject))
end


% --- Outputs from this function are returned to the command line.
function varargout = REST_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

%% LSL In
% largely drawn from run_readlsl of BCILAB

% for refrence
function lslin(handles)
% load lsl
lib = lsl_loadlib(env_translatepath('dependencies:/liblsl-Matlab/bin'));

% find stream
result = lsl_resolve_bypred(lib, ['name=''' handles.streamName '''']);

% open inlet 
inlet = lsl_inlet(result{1});
% info = inlet.info();

% create online stream data structure in base workspace (using appropriate meta-data)
onl_newstream(handles.streamName, 'srate', 128, ...
    'chanlocs', {handles.chanlocs.labels}, 'buffer_len', 10);

% state variables for recursive least squares jitter correction
P = 1e10*eye(2);        % precision matrix (inverse covariance matrix of predictors)
w = [0 0]';             % linear regression coefficients [offset,slope]
lam = 2^(-1/(128 * 30)); % forget factor in RLS calculation
n = 0;                  % number of samples observed so far    
numeric_offset = [];    % time-stamp offset to keep numerics healthy; will be initialized with first measured time stamp

% start reading
onl_read_background(handles.streamName, @read_data, 20);

    % reads from inlet
    function results = read_data()
        [chunk, stamps] = inlet.pull_chunk();
        data_clock = inlet.time_correction([], 'median', 30);
        stamps = stamps + data_clock;
        stamps = update_regression(stamps);
        chunk = double(chunk);
        
        % this is the source of grief
        % taking only the last timestamp allows REST to run but the
        % timestamps are garbage
        % giving all the timestamps soft errors in onl_append. onl_append
        % seems to expect only one value so perhaps this idea of many
        % timestamps if not correct.
        try %#ok<TRYNC>
            stamps = stamps(end);
        end
        results = {chunk, stamps};
    end
    
    
    % perform RLS block update of regression coefficients
    % this is a regression from sample index onto timestamp of the sample
    function y = update_regression(y)
        if ~isempty(y)
            % sanitize numerics (all done relative to the first observed time stamp)
            if isempty(numeric_offset)
                numeric_offset = y(1); end
            y = y - numeric_offset;        
            % define predictor matrix (bias, sample index)
            X = [ones(1,length(y)); n + (1:length(y))];
            n = n + length(y);            
            % apply updates...
            for t=1:length(y)
                u = X(:,t);
                d = y(t);
                pi = u'*P;
                gam = lam + pi*u;
                k = pi'/gam;
                al = d - w'*u;
                w = w + k*al;
                Pp = k*pi;
                P = (1/lam)*(P-Pp);
            end            
            % predict y
            y = w'*X + numeric_offset;
        end
    end


end

%% LSL Out

function pushbuttonLSLOut_Callback(hObject, eventdata, handles)
% open a window with the following:
%   a pipeline filter selector
%   a stream name text box
%   a start/stop button
% the window should stop the output stream if closed.
% the window should be closed if REST is closed.

% create figure
% menu callback changes stream
% edit callback changes stream name
% button start locks menu/edit, changes button text, and starts timer
% button stop unlocks menu/edit, changes button text, and pauses timer

% create figure
fhandle = figure('toolbar','none','Menubar','none','Name','LSL Output', ...
    'position',[500 500 500 200], 'Resize','on', ...
    'DeleteFcn',{@closeFigLSLOut, get(hObject, 'parent')});

% save fhandle to handles structure as lslout cell array
if ~isfield(handles, 'lslout')
    handles.lslout{1} = fhandle;
    lslout_ind = 1;
else
    ind = cellfun(@isempty, handles.lslout);
    if ~isempty(ind)
        handles.lslout{ind} = fhandle;
        lslout_ind = ind;
    else
        handles.lslout{end + 1} = fhandle;
        lslout_ind = length(handles.lslout);
    end
end

% stream selector
hstream = uicontrol('style', 'popupmenu', 'string', get(handles.popupmenuEEG,'String'), ...
    'units', 'normalized', 'Position', [0.05 0.5 .4 .2]);
% stream selector label
uicontrol('style', 'text', 'string', {'Select data stream';'to broadcast'}, ...
   'units', 'normalized', 'Position', [0.05 0.7 .4 .25]);
% stream name
hname = uicontrol('style', 'edit', 'units', 'normalized', 'Position', [0.55 0.55 .4 .15]);
% stream name label
uicontrol('style', 'text', 'string', {'Enter name for'; 'broadcasted stream'},...
    'units', 'normalized', 'Position', [0.55 0.7 .4 .25])
% start/stop button
uicontrol('style', 'pushbutton', 'string', 'Start Broadcast', ...
    'units', 'normalized', 'Position', [0.25 .1 0.5 .3], ...
    'Callback', {@pushbuttonStartLslout_Callback, get(hObject, 'parent'), hstream, hname})


% THIS CAN BE USED TO GET FUNCTIONS IF WE CHANGE THE GUI
% % Gather =pipeline function names
% funs = get_pipeline_functions();
% funnames = cell(length(funs)-2,1);
% for it = 2:length(funs)
%     temp = arg_report('properties',funs{it});
%     funnames{it-1} = temp.name;
%     if iscell(funnames{it-1})
%         funnames{it-1} = funnames{it-1}{1}; end
% end

    function closeFigLSLOut(varargin)
        zhObject = varargin{3};
        % load handles
        zhandles = guidata(zhObject);
        zhandles.lslout{lslout_ind} = [];
        guidata(zhObject, zhandles);
    end

end

% largely from run_writestream
function pushbuttonStartLslout_Callback(button, evnt, hfig, hstream, hname)

persistent smax

if strcmp(get(button, 'string'), 'Start Broadcast')
    
    % disable uimenu and edit
    set(hstream, 'enable', 'off')
    set(hname, 'enable', 'off')
    
    % change button text
    set(button, 'string', 'Stop Broadcast')
    
    % set values
    handles = guidata(hfig);
    stream_name = get(hname, 'string');
    stream_ind = get(hstream, 'value');
    
    % load lsl library
    lib = lsl_loadlib(env_translatepath('dependencies:/liblsl-Matlab/bin'));

    % try to calculate a UID for the stream
    try
        if strcmp(opts.source_id,'model')
            uid = hlp_cryptohash({rmfield(model,'timestamp'),opts.predict_at,opts.in_stream,opts.out_stream});
        elseif strcmp(opts.source_id,'input_data')
            uid = hlp_cryptohash({model.source_data,opts.predict_at,opts.in_stream,opts.out_stream});
        else
            error('Unsupported SourceID option: %s',hlp_tostring(opts.source_id));
        end
    catch e
        disp('Could not generate a unique ID for the predictive model; the BCI stream will not be recovered automatically after the provider system had a crash.');
        hlp_handleerror(e);
        uid = '';
    end

    % describe the stream
    disp('Creating a new streaminfo...');
    info = lsl_streaminfo(lib, stream_name, 'EEG',length(handles.chanlocs), 128, 'cf_float32',uid);

    % create an outlet
    outlet = lsl_outlet(info);

    % sample after which data will be transmitted
    smax = evalin('base', [handles.bufferName '.smax']);

    % create timer
    lsloutTimer = timer('ExecutionMode','fixedRate', 'Name',[stream_name '_timer'], 'Period',1/20, ...
        'StartDelay', 0, 'TimerFcn', @send_samples);
    
    % save timer and set delButtonFcn
    set(button, 'UserData', lsloutTimer, 'DeleteFcn', @delButtonFcn);
    
    % start timer
    start(lsloutTimer)
    
else
    % stop and delete timer
    lsloutTimer = get(button, 'UserData');
    stop(lsloutTimer)
    delete(lsloutTimer)
    set(button, 'DeleteFcn', []);
    
    % delete outlet?
    
    % reenable uicontrols
    set(hstream, 'enable', 'on')
    set(hname, 'enable', 'on')
    
    % change button text
    set(button, 'string', 'Start Broadcast')
    
end


    function send_samples(varargin)
        
        % load handles
        zhandles = guidata(hfig);
        
        % load buffer
        buffer = evalin('base', zhandles.bufferName);
        
        if buffer.smax > smax

            % determine time
            % ???

            % if ica_cleaned, generate data from ica buffer
            if stream_ind > size(buffer.data, 1)
                p = evalin('base', 'pipeline');
                ind = setdiff(1:length(p.state.icaweights), zhandles.reject);
                chunk = (p.state.icasphere \ p.state.icaweights(ind, :)') ...
                    * buffer.data{end}(ind, mod((smax + 1:buffer.smax) - 1, buffer.pnts) + 1);
            % otherwise use data directly from buffer
            else
                chunk = buffer.data{stream_ind}(:, mod((smax + 1:buffer.smax) - 1, buffer.pnts) + 1);
            end

            % push samples
            outlet.push_chunk(chunk)

            % adjust smax
            smax = buffer.smax;
        end
    end
    
    % delete timer if figure closes
    function delButtonFcn(button, evnt)
        % stop and delete timer
        lsloutTimer = get(button, 'UserData');
        stop(lsloutTimer)
        delete(lsloutTimer)
    end


end


%% button functions

% --- Executes on selection change in popupmenuEEG.
function popupmenuEEG_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuEEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
eegTimer = timerfind('Name','eegTimer');
if strcmpi(get(handles.pushbuttonPause,'string'),'Pause'); % if running
    stop(eegTimer)
    temp = get(eegTimer,'UserData');
    set(eegTimer,'UserData',[temp(1) get(hObject,'value')]);
    start(eegTimer)
else % if paused
    temp = get(eegTimer,'UserData');
    set(eegTimer,'UserData',[temp(1) get(hObject,'value')]);
end
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuEEG contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuEEG
end


% --- Executes during object creation, after setting all properties.
function popupmenuEEG_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD,*DEFNU>
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

if flag_resume %#ok<*ALIGN>
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
colorbar('hide');

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
if length(handles.hmInd) == length(J)
    Jest = zeros(handles.nVertices,1);
    Jest(handles.hmInd) = J;
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
    scalp = geometricTools.localGaussianInterpolator(handles.headModel.channelSpace,sd.surfData(1).vertices,32/1000);
end

% create figure and dipole plot
fhandle = handles.headModel.plotDipoles([0 0 0],[0 0 0]);
set(fhandle.hFigure,'DeleteFcn',{@closeFigLoc,hObject},'name',['IC' num2str(handles.curIC)]);
hold(fhandle.hAxes,'on');

% dipole fit
[dipoles, ~, moments, rv, state] = dipole_fit(Winv(:,handles.curIC),handles.K,vertices);

% set scalp potentials
set(fhandle.hScalp,'FaceVertexCData',scalp*Winv(:,handles.curIC),'facecolor','interp')

% plot dipole
arrowlen = max(range(sd.surfData(3).vertices))/10;
if length(moments)/3==length(vertices)
    dipoles.moment = dipoles.moment/norm(dipoles.moment)*arrowlen;
    handles.figLoc.fixed_dip = false;
    arrows_dip = quiver3(dipoles.location(1),dipoles.location(2),dipoles.location(3), ...
        dipoles.moment(1),dipoles.moment(2),dipoles.moment(3),'ko','filled');
else
    handles.figLoc.fixed_dip = true;
    normals = geometricTools.getSurfaceNormals(sd.surfData(end).vertices,sd.surfData(end).faces,false);
    arrows_dip = quiver3(dipoles.location(1),dipoles.location(2),dipoles.location(3), ...
        normals(dipoles.L,1)*arrowlen,normals(dipoles.L,2)*arrowlen,normals(dipoles.L,3)*arrowlen,'ko','filled');
end
colormap(bipolar(512, 0.99))
maxabs = max(abs(vec(scalp*Winv(:,handles.curIC))));
caxis(fhandle.hAxes,[-maxabs maxabs]);

% create residual variance text
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
handles.figLoc.scalp = scalp;
handles.figLoc.arrowlen = arrowlen;
if handles.figLoc.fixed_dip
    handles.figLoc.normals = normals; end
guidata(hObject,handles);

end


function figLoc_update_dipolefit(varargin)
% get ica decomposition
W = evalin('base','pipeline.state.icaweights');
sphere = evalin('base','pipeline.state.icasphere');
Winv = inv(W*sphere);

% parse inputs
handles = guidata(varargin{3});
if handles.figLoc.fixed_dip
    normals = handles.figLoc.normals; end

% load surfaces
temp = load(handles.headModel.surfaces);
vertices = temp.surfData(3).vertices(handles.hmInd,:);

% update dipole
[dipoles, ~, ~, rv, handles.figLoc.state] = ...
    dipole_fit(Winv(:,handles.figLoc.IC),handles.K,vertices,handles.figLoc.state);

% update dipole plot
if ~handles.figLoc.fixed_dip
    dipoles.moment = dipoles.moment/norm(dipoles.moment)*handles.figLoc.arrowlen;
    set(handles.figLoc.arrows_dip, ...
        'XData',dipoles.location(1), ...
        'YData',dipoles.location(2), ...
        'ZData',dipoles.location(3), ...
        'UData',dipoles.moment(1), ...
        'VData',dipoles.moment(2), ...
        'WData',dipoles.moment(3));
else
    set(handles.figLoc.arrows_dip, ...
        'XData',dipoles.location(1), ...
        'YData',dipoles.location(2), ...
        'ZData',dipoles.location(3), ...
        'UData',normals(dipoles.L,1)*handles.figLoc.arrowlen, ...
        'VData',normals(dipoles.L,2)*handles.figLoc.arrowlen, ...
        'WData',normals(dipoles.L,3)*handles.figLoc.arrowlen);
end

% update scalp plot
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
stop(handles.pauseTimers(locTimerInd));
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
fhandle = figure('toolbar','none','Menubar','none','Name','IC Select','position',[1 1 temp(3:4)],'Resize','on','Colormap',colormap('jet'),'DeleteFcn',{@closeFigIC,hObject});
W = evalin('base','pipeline.state.icaweights');
sphere = evalin('base','pipeline.state.icasphere');
Winv = inv(W*sphere);

rowcols(2) = ceil(sqrt(handles.nic));
rowcols(1) = ceil(handles.nic/rowcols(2));
scaleMatTopo = [1 0 0 0;0 1 0 0;0 0 1 0;0 .2 0 .8];
buttonGap = .1;
scaleMatReject = [1 0 0 0;0 1 0 0;0 0 .5-buttonGap/2 0;.5+buttonGap/2 0 0 .2];
for it = 1:handles.nic
    h(it) = subaxis(rowcols(1),rowcols(2),it,'MR',.025,'ML',.025,'MT',.025,'MB',.025,'SH',0,'SV',0.02);
    tempPos = get(h(it),'Position');
    set(h(it),'position',get(h(it),'position')*scaleMatTopo)
    topoplotFast_LRBF(Winv(:,it),handles.chanlocs);
    title(['IC' int2str(it)])
    
    lock = any(handles.lock==it);
    reject = any(handles.reject==it);
    buttonLock(it) = uicontrol('Style', 'togglebutton', 'String', 'Lock',...
        'Units','normalize','Position', tempPos.*[1 1 .5-buttonGap/2 .2],...
        'Callback', {@lockIC,it,hObject},'Value',lock,...
        'BackgroundColor',handles.color_lock*lock + handles.color_bg*(1-lock));
    buttonReject(it) = uicontrol('Style', 'togglebutton', 'String', 'Reject',...
        'Units','normalize','Position', tempPos*scaleMatReject,...
        'Callback', {@rejectIC,it,hObject},'Value',reject,...
        'BackgroundColor',handles.color_reject*reject + handles.color_bg*(1-reject));
end
handles.figIC.buttonLock = buttonLock;
handles.figIC.buttonReject = buttonReject;
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
    set(button,'BackgroundColor',[0.5 1 0.5])
    if isfield(handles,'figIC')
        set(handles.figIC.buttonLock(ic),'Value',1,'BackgroundColor',handles.color_lock); end
else % turned lock off
    handles.lock(handles.lock==ic) = [];
    set(button,'BackgroundColor',handles.color_bg)
    if isfield(handles,'figIC')
        set(handles.figIC.buttonLock(ic),'value',0,'BackgroundColor',handles.color_bg); end
end
% save handles
guidata(hObject,handles);
% update ics to plot
updateICs(hObject,false)
end


function rejectIC(varargin)
ic = varargin{3};
button = varargin{1};
% load handles
if numel(varargin)>3
    hObject = varargin{4};
else
    hObject = get(button,'parent');
end
handles = guidata(hObject);
if get(button,'Value') % turned reject on
    handles.reject = sort([handles.reject ic]);
    set(button,'BackgroundColor',[1 0.5 0.5])
    if isfield(handles,'figIC')
        set(handles.figIC.buttonReject(ic),'Value',1,'BackgroundColor',handles.color_reject); end
else % turned reject off
    handles.reject(handles.reject==ic) = [];
    set(button,'BackgroundColor',handles.color_bg)
    if isfield(handles,'figIC')
        set(handles.figIC.buttonReject(ic),'value',0,'BackgroundColor',handles.color_bg); end
end
% save handles
guidata(hObject,handles);
end


function updateICs(hObject,flag_sort)
% load handles
handles = guidata(hObject);
if flag_sort
    % load info
    S = evalin('base','pipeline.state.icasphere');
    W = evalin('base','pipeline.state.icaweights');
    V = evalin('base','pipeline.state.Var');
    other = setdiff(1:handles.nic,[handles.lock]);
    % sort by residual variance
    meanvar = mean(pinv(W*S).^2).*V';
    [~,ind_lock] = sort(meanvar(handles.lock),'descend');
    [~,ind_other] = sort(meanvar(other),'descend');
    handles.ics = [handles.lock(ind_lock) other(ind_other)];
else
    handles.ics = [handles.lock setdiff(1:handles.nic,[handles.lock])];
end
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
%         set(handles.(['togglebuttonReject' int2str(ind)]),'Value',onFlag); end
%     % update fig
%     if isfield(handles,'figIC')
%         set(handles.figIC.buttonReject(ic),'Value',onFlag); end
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
        learning_rate = evalin('base','learning_rate');
        srate = evalin('base',[parseStreamName(handles.streamName) '.srate']);
        x = -(length(learning_rate)-1)/srate:1/srate:0;
        axes(handles.axisInfo)
        [handles.axisInfo] = plot(x,learning_rate);
        line1 = get(handles.axisInfo,'children');
        set(get(get(handles.axisInfo(1),'parent'),'XLabel'),'String','Time (seconds)')
        set(get(get(handles.axisInfo(1),'parent'),'YLabel'),'String','Learning Rate (dB)')
        set(infoTimer,'Period',1,'ExecutionMode','fixedRate', ...
            'TimerFcn',{@infoConverge,handles.axisInfo,handles.panelInfo},'StartDelay',0);
        axis(get(handles.axisInfo,'parent'),'tight')
        start(infoTimer)
    otherwise
        warning('REST: popupmenuInfo recieved a strange input')
end
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


% --- Executes on button press in togglebuttonReject8.
function togglebuttonReject8_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonReject8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rejectIC(hObject,[],handles.ics(8))
end


% --- Executes on button press in togglebuttonLock8.
function togglebuttonLock8_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lockIC(hObject,[],handles.ics(8))
end


% --- Executes on button press in togglebuttonReject7.
function togglebuttonReject7_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonReject7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rejectIC(hObject,[],handles.ics(7))
end


% --- Executes on button press in togglebuttonLock7.
function togglebuttonLock7_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lockIC(hObject,[],handles.ics(7))
end


% --- Executes on button press in togglebuttonReject6.
function togglebuttonReject6_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonReject6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rejectIC(hObject,[],handles.ics(6))
end


% --- Executes on button press in togglebuttonLock6.
function togglebuttonLock6_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lockIC(hObject,[],handles.ics(6))
end


% --- Executes on button press in togglebuttonReject5.
function togglebuttonReject5_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonReject5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rejectIC(hObject,[],handles.ics(5))
end


% --- Executes on button press in togglebuttonLock5.
function togglebuttonLock5_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lockIC(hObject,[],handles.ics(5))
end


% --- Executes on button press in togglebuttonReject4.
function togglebuttonReject4_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonReject4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rejectIC(hObject,[],handles.ics(4))
end


% --- Executes on button press in togglebuttonLock4.
function togglebuttonLock4_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lockIC(hObject,[],handles.ics(4))
end


% --- Executes on button press in togglebuttonReject3.
function togglebuttonReject3_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonReject3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rejectIC(hObject,[],handles.ics(3))
end


% --- Executes on button press in togglebuttonLock3.
function togglebuttonLock3_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lockIC(hObject,[],handles.ics(3))
end


% --- Executes on button press in togglebuttonReject2.
function togglebuttonReject2_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonReject2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rejectIC(hObject,[],handles.ics(2))
end


% --- Executes on button press in togglebuttonLock2.
function togglebuttonLock2_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonLock2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lockIC(hObject,[],handles.ics(2))
end


% --- Executes on button press in togglebuttonReject1.
function togglebuttonReject1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonReject1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rejectIC(hObject,[],handles.ics(1))
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
    % close figures
    if isfield(handles,'figIC')
        try
            close(handles.figIC.handle); end, end
    if isfield(handles,'figLoc')
        try
            close(handles.figIC.handle); end, end
    if isfield(handles, 'lslout')
        for it = 1:length(handles.lslout)
            try
                close(handles.lslout{it}); end, end, end

    % delete timers
    timerNames = {'eegTimer','oricaTimer','topoTimer','infoTimer','locTimer',[parseStreamName(handles.streamName) '_timer']};

    for it = 1:length(timerNames)
        temp = timerfind('name',timerNames{it});
        if isempty(temp)
            continue; end
        stop(temp)
        delete(temp)
    end
    if isfield(handles, 'playbackStream')
        stop(handles.playbackStream); end
end
end


% --- Executes on button press in pushbuttonSortICs.
function pushbuttonSortICs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSortICs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateICs(hObject,true)
end


% utility functions
% parse streamname
function streamnames = parseStreamName(streamnames)
if ~isvarname(streamnames)
    streamnames = streamnames(~ismember(streamnames,['-' ' ']));
end
end
