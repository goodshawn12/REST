function lslout(button, evnt, hfig, hstream, hname)

% if the 'Start Broadcast' button is pressed
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
        uid = '';
%         uid = hlp_cryptohash({button, evnt, hfig, hstream, hname});
    catch e
        disp('Could not generate a unique ID for the predictive model; the BCI stream will not be recovered automatically after the provider system had a crash.');
        hlp_handleerror(e);
        uid = '';
    end

    % describe the stream
    disp('Creating a new streaminfo...');
    info = lsl_streaminfo(lib, stream_name, 'EEG',length(handles.chanlocs), handles.srate, 'cf_float32', uid);

    % create an outlet
    outlet = lsl_outlet(info);

    % sample after which data will be transmitted
    smax = evalin('base', [handles.bufferName '.smax']);

    % create timer
    lsloutTimer = timer('ExecutionMode','fixedRate', 'Name',[stream_name '_timer'], 'Period',1/20, ...
        'StartDelay', 0, 'TimerFcn', @send_samples, 'UserData', smax);
    
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
        
        % load smax
        smax_local = get(varargin{1}, 'UserData');
        
%         persistent smax_local
%         if isempty(smax_local)
%             smax_local = smax;
%         end
        if buffer.smax > smax_local

            % determine time
            % ???

            % if ica_cleaned, generate data from ica buffer
            if stream_ind > size(buffer.data, 1)
                p = evalin('base', 'pipeline');
                ind = setdiff(1:length(p.state.icaweights), zhandles.reject);
                Winv = (p.state.icasphere \ p.state.icaweights');
                chunk =  Winv(:, ind) * buffer.data{end}(ind, mod((smax_local + 1:buffer.smax) - 1, buffer.pnts) + 1);
            % otherwise use data directly from buffer
            else
                chunk = buffer.data{stream_ind}(:, mod((smax_local + 1:buffer.smax) - 1, buffer.pnts) + 1);
            end

            % push samples
            outlet.push_chunk(chunk)

            % adjust smax
            set(varargin{1}, 'UserData', buffer.smax);
        end
    end
    
    % delete timer if figure closes
    function delButtonFcn(button, evnt)
        % stop and delete timer and outlet
        lsloutTimer = get(button, 'UserData');
        stop(lsloutTimer)
        delete(lsloutTimer)
        delete(outlet)
    end

end