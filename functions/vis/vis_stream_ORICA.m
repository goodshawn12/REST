function [inlet, buffername] = vis_stream_ORICA(varargin)
% Display an LSL stream.
%
% Keyboard shortcuts:
%   [up arrow]   : increase the y scale of the time series
%   [down arrow] : decrease the y scale of the time series
%   [right arrow]: increase the displayed time range
%   [left arrow] : decrease the displayed time range
%   [page up]    : go up by one page of channels
%   [page down]  : go down by one page of channels
%
% In:
%   StreamName : Stream to display. The name of the stream that you would like to display.
%
%   TimeScale : Initial time scale in seconds. The time range of the display window;
%               can be changed with keyboard shortcuts (see help). Default=5
%
%   DataScale : Initial scale of the data. The scale of the data, in units between horizontal lines;
%               can be changed with keyboard shortcuts (see help). Default=150
%
%   ChannelRange : Channels to display. The channel range to display. Default=[1:32]
%
%   SamplingRate : Sampling rate for display. This is the sampling rate that is used for plotting, in Hz;
%                  for faster drawing. Default=100
%
%   RefreshRate : Refresh rate for display. This is the rate at which the graphics are updated, in Hz.
%                 Default=10
%
%   Rereference : Apply common-average re-referencing to the data. Useful for noisy EEG streams.
%                 Default=false
%
%   PageOffset : Channel page offset. Allows to flip forward or backward pagewise through the displayed channels.
%                Default=0
%
%   Position : Figure position. Allows to script the position at which the figures should appear.
%              This is a 4-element vector of the form [X-offset,Y-offset,Width,Height]
%              with all values in pixes.
%              Default=[]
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2012-07-10
%
%                                uses portions of vis_dataStreamViewer
%                                (c) 2012 by Tim Mullen


% make sure that everything is on the path and LSL is loaded
if ~exist('arg_define','file')
    addpath(genpath(fileparts(mfilename('fullpath')))); end
if ~exist('env_translatepath','file')
    % standalone case
    lib = lsl_loadlib();
else
    % if we're within BCILAB we want to make sure that the library is also found if the toolbox is compiled
    lib = lsl_loadlib(env_translatepath('dependencies:/liblsl-Matlab/bin'));
end

% handle input arguments
opts = arg_define(varargin, ...
    arg({'streamname','StreamName'},[],[],'LSL stream that should be displayed. The name of the stream that you would like to display.'), ...
    arg({'property','SelectionProperty'}, 'type',[],'Selection property. The selection criterion by which the desired device is identified on the net. This is a property that the desired device must have (e.g., name, type, desc/manufacturer, etc.'), ...
    arg({'value','SelectionValue'}, 'EEG',[],'Selection value. This is the value that the desired device must have for the selected property (e.g., EEG if searching by type, or Biosemi if searching by manufacturer).'), ...
    arg({'bufferrange','BufferRange'},360,[],'Maximum time range to buffer. Imposes an upper limit on what can be displayed.'), ...
    arg({'timerange','TimeRange'},5,[],'Initial time range in seconds. The time range of the display window; can be changed with keyboard shortcuts (see help).'), ...
    arg({'datascale','DataScale'},150,[],'Initial scale of the data. The scale of the data, in units between horizontal lines; can be changed with keyboard shortcuts (see help).'), ...
    arg({'channelrange','ChannelRange'},1:32,[],'Channels to display. The channel range to display.'), ...
    arg({'samplingrate','SamplingRate'},100,[],'Sampling rate for display. This is the sampling rate that is used for plotting; for faster drawing.'), ...
    arg({'refreshrate','RefreshRate'},10,[],'Refresh rate for display. This is the rate at which the graphics are updated.'), ...
    arg({'reref','Rereference'},false,[],'Common average reference. Enable this to view the data with a common average reference filter applied.'), ...
    arg_nogui({'pageoffset','PageOffset'},0,[],'Channel page offset. Allows to flip forward or backward pagewise through the displayed channels.'), ...
    arg_nogui({'position','Position'},[],[],'Figure position. Allows to script the position at which the figures should appear.'), ...
    arg_nogui({'figurehandles','FigureHandle'},[],[],'Handle to the desired figure.'), ...
    arg_nogui({'axishandles','AxisHandle'},[],[],'Handle to the desired axis.'));

% no arguments were passed? bring up GUI dialog
if isempty(varargin)
    opts = arg_guidialog;    
    if isempty(opts)
        return; end % -> user clicked cancel
end

% fix up some arguments
opts.bufferrange = max(opts.bufferrange,opts.timerange);

% init shared handles
[fig,ax,lines] = deal([]);

% streamnames = find_streams;
% if length(streamnames) == 1
%     assignin('base','streamname',streamnames); % save stream name in base workspace
% else
%     % if more than 2 (EEG) streams, pop up a GUI to select one.
%     selStream(streamnames); % result 'streamname' is saved in base workspace
% %     delete(h)
%     streamnames = evalin('base','streamname');
% end


% choose variable names to hold the stream's data (in the base workspace)
taken = evalin('base','whos(''lsl*'')');

% create a stream inlet
% inlet = create_inlet(opts);
inlet = [];

% create the stream data structure in the base workspace
% if ~isvarname(opts.streamname) opts.streamname = opts.streamname(~isspace(opts.streamname)); end
chunkname = genvarname(['lsl_' opts.streamname '_chunk'],{taken.name});
buffername = genvarname(['lsl_' opts.streamname '_stream'],{taken.name});

buffer = create_streambuffer(opts);
assignin('base', buffername, buffer);

% create the figure
create_figure(opts);

% set up a timer that reads from LSL
th = timer('Period', 1.0/opts.refreshrate,'ExecutionMode','fixedRate','TimerFcn',@on_timer,...
    'StartDelay',0.2,'Tag','lsl_visORICAst_timer',...
    'Name','eegTimer','UserData',1);

% th = timer('Period', 1.0/opts.refreshrate,'ExecutionMode','fixedRate','TimerFcn',@on_timer,...
%     'StartDelay',0.2,'Tag',['lsl_' genvarname(opts.streamname) '_timer'],...
%     'Name','eegTimer','UserData',0);
% start(th);


    % === nested functions (sharing some handles with each other) ===

    % create a new figure and axes
    function create_figure(opts)
        if ~isempty(opts.figurehandles)
            fig = opts.figurehandles;
            set(fig,'CloseRequestFcn','delete(gcbf)','KeyPressFcn',@(varargin)on_key(varargin{2}.Key))
        else
            fig = figure('Tag',['Fig' buffername],'Name',['LSL:Stream''' opts.streamname ''''],'CloseRequestFcn','delete(gcbf)', ...
                'KeyPressFcn',@(varargin)on_key(varargin{2}.Key));
        end
        if ~isempty(opts.position)
            set(fig,'Position',opts.position); end
        if ~isempty(opts.axishandles)
            ax = opts.axishandles;
            set(ax,'YDir','reverse')
        else
            ax = axes('Parent',fig, 'Tag','LSLViewer', 'YDir','reverse');
        end
    end

    function on_timer(varargin)
        try 
            % check if the buffer is still there
            if evalin('base',['exist(''' buffername ''',''var'')'])
                
                % check what to play
                timerdata = get(varargin{1},'userdata');
                handles = guidata(timerdata{1});
                plot_content = timerdata{2};
                
                % === update buffer contents (happens in the base workspace) ===
                
%                 % pull a new chunk from LSL
%                 assignin('base',chunkname,inlet.pull_chunk());
%                 
%                 % append it to the stream buffer
%                 evalin('base',['[' buffername '.smax,' buffername '.data(:,1+mod(' buffername '.smax:' buffername '.smax+size(' chunkname ',2)-1,' buffername '.pnts))] = deal(' buffername '.smax + size(' chunkname ',2),' chunkname ');']);
                
                % get the updated stream buffer
                stream = evalin('base',buffername);
                
                % reformat the stream buffer to contain only the current block that should be displayed
                samples_to_get = min(stream.pnts, round(stream.srate*stream.opts.timerange));
                if plot_content<=length(stream.data)
                    stream.data{plot_content} = stream.data{plot_content}(:, 1+mod(stream.smax-samples_to_get:stream.smax-1,stream.pnts));
                    [stream.nbchan,stream.pnts,stream.trials] = size(stream.data{plot_content});
                    stream.xmax = stream.smax/stream.srate;
                    stream.xmin = stream.xmax - (stream.pnts-1)/stream.srate;
                    plotdata = stream.data{plot_content}(:, round(1 : stream.srate/stream.opts.samplingrate : end));
                    if plot_content==length(stream.data)
                        W = evalin('base','pipeline.state.icaweights');
                        sphere = evalin('base','pipeline.state.icasphere');
                        stream.opts.datascale = stream.opts.datascale*mean(sqrt(mean((W*sphere).^2)));%component_scale/data_scale; %mean(abs(W*sphere)*ones(length(sphere),1));
                    end
                else % plot subtracted components
%                     stream.d= stream.data{plot_content}(:, 1+mod(stream.smax-samples_to_get:stream.smax-1,stream.pnts));
                    stream.data{end} = stream.data{end}(:, 1+mod(stream.smax-samples_to_get:stream.smax-1,stream.pnts));
                    plotica = stream.data{end}(:, round(1 : stream.srate/stream.opts.samplingrate : end));
%                     variance = evalin('base','pipeline.state.Var');
                    W = evalin('base','pipeline.state.icaweights');
                    sphere = evalin('base','pipeline.state.icasphere');
                    Winv = pinv(W*sphere);
                    [stream.nbchan,stream.pnts,stream.trials] = size(stream.data{end});
                    stream.xmax = stream.smax/stream.srate;
                    stream.xmin = stream.xmax - (stream.pnts-1)/stream.srate;
                    ind = setdiff(1:length(W),handles.reject);
                    plotdata = Winv(:,ind)*plotica(ind,:);
                end
                
                % find removed channels stage or pipeline they were removed
                flag_channel_removed = false;
                pipeline = evalin('base','pipeline');
                count = length(stream.data);
                index = [];
                while true
                    if ~isequal(pipeline.head,@flt_clean_channels) && ~isequal(pipeline.head,@flt_selchans)
                        try
                            pipeline = pipeline.parts{2}; % !!! make this more general
                            count = count-1;
                        catch
                            break
                        end
                    elseif isequal(pipeline.head,@flt_clean_channels)
                        index = pipeline.parts{2}.parts{end-1};
                        flag_channel_removed = true;
                        break
                    elseif isequal(pipeline.head,@flt_selchans)
                        plotica = arg_define(pipeline.parts, ...
                            arg_norep({'signal','Signal'}), ...
                            arg({'channels','Channels'}, [], [], 'Cell array of channel names to retain.','type','cellstr','shape','row'), ...
                            arg({'orderPreservation','OrderPreservation'}, 'query-order', {'query-order','dataset-order'}, 'Output channel order. The result will have its channels either in the order of the input set or in the order of the query list.'), ...
                            arg({'remove_selection','RemoveSelection'},false,[],'Remove selected channels.'), ...
                            arg({'find_closest','FindClosest'},false,[],'Find closest channels. This is for cases where the requested channels are not in the set.'), ...
                            arg({'relabel_to_query','RelabelToQuery'},false,[],'Relabel closest channels according to query. This only applies if FindClosest was true.'));
                        if plotica.remove_selection
                            index = find(ismember({stream.chanlocs.labels},plotica.channels));
                        else
                            [~,index] = setdiff({stream.chanlocs.labels},plotica.channels);
                            index = sort(index);
                        end
                        flag_channel_removed = true;
                        break
                    end
                end
                
                % determine channels and samples to display
                plotchans = stream.opts.channelrange + stream.opts.pageoffset*length(stream.opts.channelrange);
                if isempty(plotchans)
                    plotchans = 1:stream.nbchan;
                else
                    plotchans = intersect(1:stream.nbchan,plotchans);
                end
                plotdata = plotdata(plotchans,:);
                plottime = linspace(stream.xmin,stream.xmax,size(plotdata,2));
                if plot_content>=count
                    channels_remaining = setdiff(1:length(stream.chanlocs),index);
                    plotchans = channels_remaining(plotchans);
                end
                
                % re-reference
                if stream.opts.reref
                    plotdata = bsxfun(@minus,plotdata,mean(plotdata)); end
                                
                % zero-mean
                plotdata = bsxfun(@minus, plotdata, mean(plotdata,2));
                
                
                % arrange for plotting
                plotoffsets = (0:size(plotdata,1)-1)*stream.opts.datascale;
                plotdata = bsxfun(@plus, plotdata', plotoffsets);
                
                
                % === actual drawing ===
                
                % draw the block contents...
                if ~isempty(plotdata)
                    if ~exist('lines','var') || isempty(lines)                        
                        lines = plot(ax,plottime,plotdata);
                        title(ax,opts.streamname,'interpreter','none');
                        xlabel(ax,'Time (sec)','FontSize',12);
                        ylabel(ax,'Activations','FontSize',12);
                    else
                        for k=1:length(lines)
                            if k <= size(plotdata,2);
                                set(lines(k),'Ydata',plotdata(:,k),'visible','on','linestyle','-');
                                set(lines(k),'Xdata',plottime);
                            else
                                set(lines(k),'Ydata',zeros(size(plotdata(:,1))),'visible','off');
                                set(lines(k),'Xdata',plottime);
                            end
                        end
                    end
                
                    % update the axis limit and tickmarks
                    axis(ax,[stream.xmin stream.xmax -stream.opts.datascale size(plotdata,2)*stream.opts.datascale + stream.opts.datascale]);
                    if plot_content==length(stream.data)
                        set(ax, 'YTick',plotoffsets, 'YTickLabel',cellstr(int2str(plotchans')));
                    elseif plot_content<count
                        set(ax, 'YTick',plotoffsets, 'YTickLabel',{stream.chanlocs(plotchans).labels});
                    else
                        set(ax, 'YTick',plotoffsets, 'YTickLabel',{stream.chanlocs(plotchans).labels});
                    end
                end
                
                % show which channels are not in pipeline when applicable
                if flag_channel_removed && plot_content<count
                    set(lines(ismember(plotchans,index)),'linestyle','--')
                end
                
                drawnow;
            else
                try 
                    disp(['Deleting timer ' get(th,'Tag') '.']);
                catch e
                    disp('Deleting timer.');
                end
                % delete the timer
                warning off MATLAB:timer:deleterunning
                th = timerfindall;
                delete(th);
            end
        catch e
            if isempty(findobj('Tag',['Fig' buffername]))
                disp('Figure was closed.');
            else
                disp('An error occurred during the stream viewer update: ');
                hlp_handleerror(e);
            end
            warning off MATLAB:timer:deleterunning
            th = timerfindall;
            delete(th);
        end
    end

    function on_key(key)
        stream = evalin('base',buffername);
        switch lower(key)
            case 'uparrow'
                % decrease datascale
                stream.opts.datascale = stream.opts.datascale*0.9;
            case 'downarrow'
                % increase datascale
                stream.opts.datascale = stream.opts.datascale*1.1;
            case 'rightarrow'
                % increase timerange
                stream.opts.timerange = stream.opts.timerange*1.1;                
            case 'leftarrow'
                % decrease timerange
                stream.opts.timerange = stream.opts.timerange*0.9;                
            case 'pagedown'
                % shift display page offset down
                stream.opts.pageoffset = mod(stream.opts.pageoffset+1, ...
                    ceil(stream.nbchan/length(stream.opts.channelrange)));
            case 'pageup'
                % shift display page offset up
                stream.opts.pageoffset = mod(stream.opts.pageoffset-1, ...
                    ceil(stream.nbchan/length(stream.opts.channelrange)));
        end
        assignin('base',buffername,stream);
    end


    % === utility functions ===
    
    % find names of streams on the lab network...
    function names = find_streams
        streams = lsl_resolve_all(lib,3);
        names = cellfun(@(s)s.name(),streams ,'UniformOutput',false);
        if isempty(names)
            error('There is no stream visible on the network.'); end
    end


    % create an inlet to read from the stream with the given name
    function inlet = create_inlet(opts)
        % look for the desired device
        result = {};
        disp(['Looking for a stream with name = ' opts.streamname ' ...']);
        while isempty(result)
            result = lsl_resolve_byprop(lib,'name',opts.streamname); end
%             result = lsl_resolve_byprop(lib,opts.property,opts.value);end 
        % create a new inlet
        disp('Opening an inlet...');
        inlet = lsl_inlet(result{1},opts.bufferrange);
    end


    % create a new stream buffer in the base workspace
    function stream = create_streambuffer(opts)
        base_stream = evalin('base',opts.streamname);
        stream.srate = base_stream.srate;                                   % sampling rate in Hz
        stream.chanlocs = base_stream.chanlocs;                             % struct with per-channel meta-data
        stream.pnts = max(opts.bufferrange*stream.srate,100);               % number of data points in the buffer
        stream.nbchan = base_stream.nbchan;                                 % number of channels in the buffer
        stream.trials = 1;                                                  % number of segments in the buffer (always 1)
        pipe_len = find_pipeline_length();                                  % # of pipeline outputs
        stream.data = repmat({zeros(stream.nbchan,stream.pnts,stream.trials)},pipe_len,1);       % the circular buffer storage
        stream.smax = 0;                                                    % number of samples that have been written into the buffer so far (wrapping around)
        stream.opts = opts;                                                 % current display options for this stream
    end
    % create a new stream buffer in the base workspace
%     function stream = create_streambuffer(opts,info)
%         stream.srate = info.nominal_srate();                                % sampling rate in Hz
%         stream.chanlocs = struct('labels',derive_channel_labels(info));     % struct with per-channel meta-data
%         stream.pnts = max(opts.bufferrange*stream.srate,100);               % number of data points in the buffer
%         stream.nbchan = info.channel_count();                               % number of channels in the buffer
%         stream.trials = 1;                                                  % number of segments in the buffer (always 1)
%         stream.data = zeros(stream.nbchan,stream.pnts,stream.trials);       % the circular buffer storage
%         stream.smax = 0;                                                    % number of samples that have been written into the buffer so far (wrapping around)
%         stream.opts = opts;                                                 % current display options for this stream
%     end

    % derive a list of channel labels for the given stream info
    function channels = derive_channel_labels(info)
        channels = {};
        ch = info.desc().child('channels').child('channel');
        while ~ch.empty()
            name = ch.child_value_n('label');
            if name
                channels{end+1} = name; end
            ch = ch.next_sibling_n('channel');
        end
        if length(channels) ~= info.channel_count()
            disp('The number of channels in the steam does not match the number of labeled channel records. Using numbered labels.');
            channels = cellfun(@(k)['Ch' num2str(k)],num2cell(1:info.channel_count(),1),'UniformOutput',false);
        end
    end

    % find pipeline length
    function pipe_len = find_pipeline_length(p,pipe_len)
        if ~exist('p','var')
            p = evalin('base','pipeline'); end
        if ~exist('pipe_len','var')
            pipe_len = 0; end
        
        if p.subnodes
            for k = p.subnodes
                pipe_len = find_pipeline_length(p.parts{k},pipe_len);
            end
        end
        pipe_len = pipe_len + 1;
    end
    
end


