% largely drawn from run_readlsl of BCILAB
function lslin(handles)
% load lsl
lib = lsl_loadlib(env_translatepath('dependencies:/liblsl-Matlab/bin'));

% find stream
result = lsl_resolve_bypred(lib, ['name=''' handles.streamName '''']);

% open inlet 
inlet = lsl_inlet(result{1});
% info = inlet.info();

% create online stream data structure in base workspace (using appropriate meta-data)
onl_newstream(parseStreamName(handles.streamName), 'srate', 128, ...
    'chanlocs', {handles.chanlocs.labels}, 'buffer_len', 10);

% state variables for recursive least squares jitter correction
P = 1e10*eye(2);        % precision matrix (inverse covariance matrix of predictors)
w = [0 0]';             % linear regression coefficients [offset,slope]
lam = 2^(-1/(128 * 30)); % forget factor in RLS calculation
n = 0;                  % number of samples observed so far    
numeric_offset = [];    % time-stamp offset to keep numerics healthy; will be initialized with first measured time stamp

% start reading
onl_read_background(parseStreamName(handles.streamName), @read_data, 20);

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


    % parse streamname
    function streamnames = parseStreamName(streamnames)
        if ~isvarname(streamnames)
            streamnames = streamnames(~ismember(streamnames,['-' ' ']));
        end
    end


end