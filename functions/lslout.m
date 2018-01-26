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
        uid = ['lsl_REST_' num2str(stream_ind) '_' stream_name];
%         uid = hlp_cryptohash({button, evnt, hfig, hstream, hname});
    catch e
        disp('Could not generate a unique ID for the predictive model; the BCI stream will not be recovered automatically after the provider system had a crash.');
        hlp_handleerror(e);
        uid = '';
    end
    
    % load buffer
    buffer = evalin('base', handles.bufferName);

    % describe the stream
    disp('Creating a new streaminfo...');
    if stream_ind <= length(buffer.data) + 1
        info = lsl_streaminfo(lib, stream_name, ...
            'EEG', size(buffer.data{min(stream_ind, length(buffer.data) - 1)}, 1), ...
            handles.srate, 'cf_float32', uid);
        if stream_ind ~= length(buffer.data)
            % channel labels
            chns = info.desc().append_child('channels');
            for label = {handles.chanlocs{min(stream_ind, length(buffer.data) - 1)}.labels}
                ch = chns.append_child('channel');
                ch.append_child_value('label',label{1});
                ch.append_child_value('unit','microvolts');
                ch.append_child_value('type','EEG');
            end
        end
    elseif stream_ind <= length(buffer.data) + 3
        info = lsl_streaminfo(lib, stream_name, ...
            'Parameters', size(buffer.data{end}, 1)^2, ...
            [], 'cf_double64', uid);
    else
        info = lsl_streaminfo(lib, stream_name, ...
            'Convergence', 1, handles.srate, 'cf_float32', uid);
    end

    % create an outlet
    outlet = lsl_outlet(info);

    % create timer
    lsloutTimer = timer('ExecutionMode','fixedRate', 'Name',[stream_name '_timer'], 'Period',1/20, ...
        'StartDelay', 0, 'TimerFcn', @send_samples, 'UserData', buffer.smax);
    
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
        zbuffer = evalin('base', zhandles.bufferName);
        
        % load smax
        smax_local = get(varargin{1}, 'UserData');
        
%         persistent smax_local
%         if isempty(smax_local)
%             smax_local = smax;
%         end
        if zbuffer.smax > smax_local

            % determine time labels for output stream
            % ???

            % if ica_cleaned, generate data from ica buffer
            if stream_ind == size(zbuffer.data, 1) + 1
                p = evalin('base', 'pipeline');
                ind = setdiff(1:length(p.state.icaweights), zhandles.reject);
                Winv = (p.state.icasphere \ p.state.icaweights');
                chunk =  Winv(:, ind) * zbuffer.data{end}(ind, mod((smax_local + 1:zbuffer.smax) - 1, zbuffer.pnts) + 1);
                
            % if icasphere
            elseif stream_ind == size(zbuffer.data, 1) + 2
%                 if zbuffer.ica.smax > smax_local
                    chunk = zbuffer.ica.icasphere(:);
%                 end
            % if icaweights
            elseif stream_ind == size(zbuffer.data, 1) + 3
%                 if zbuffer.ica.smax > smax_local
                    chunk = zbuffer.ica.icaweights(:);
%                 end
                
            % if normRn
            elseif stream_ind == size(zbuffer.data, 1) + 4
                chunk = zbuffer.ica.normRn(:, mod((smax_local + 1:zbuffer.smax) - 1, zbuffer.pnts) + 1);
                
            % otherwise use data directly from buffer
            else
                chunk = zbuffer.data{stream_ind}(:, mod((smax_local + 1:zbuffer.smax) - 1, zbuffer.pnts) + 1);
            end

            % push samples
            outlet.push_chunk(chunk)

            % adjust smax
            set(varargin{1}, 'UserData', zbuffer.smax);
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