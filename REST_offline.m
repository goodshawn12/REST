% Run the REST pipeline on pre-recorded data iteratively, but not online.
%
% Inputs:
%       dataset: Data to to be processed by REST. Can be either a path to
%           a *.set file or an EEGLAB EEG structure.
%       calibrationset: Data to use for calibration. Can be either a path 
%           to a *.set file, a cell array of inputs for io_loadset or the
%           output of io_loadset.
%       pipeline_desc: The processing pipeline BCILAB will run. Can be
%           either a cell array or structure accepted by flt_pipeline. If
%           omitted, the default REST pipeline is used.
%
% For ease of customization, the default pipeline_desc is copied here:
%     pipeline_desc = { ...
%         'selchans', {'channels', {'PO7'  'Packet Counter'  'TRIGGER'}, ...
%             'remove_selection', 1}, ...
%         'reref', {'ref_chn', {'A2'}}, ...
%         'fir', {'fspec', [0.5000 1 50 55], 'stopripple', -40, ...
%             'normalize_amplitude', 0}, ...
%         'repair_bursts', {'stddev_cutoff', 20, 'window_len', 0.5}, ...
%         'orica', {'onlineWhitening', {'blockSize', 16}, ... 
%             'options', {'blockSize', 16}, 'evalConvergence', 0.01} ...
%         'eyecatch', 'on', ...
%         'ica_reproject', {'ProjectionMatrix', '.icawinv', 'ComponentSubset', '.reject'}, ...
%         };
%     pipeline_desc = [{'FilterOrdering', strcat('flt_', pipeline_desc(1:2:end))}, pipeline_desc];
% 
% Example using default pipeline:
%   eegstr = 'data/Quick30_Shawn_EyeOpen.set';
%   calib = {eegstr, 'channels', 1:30, 'timerange', [0 30]};
%   REST_offline(eegstr, calib)
% 
% TODO:
%   * save more information about ic classifications across time
%   * save more information about ica decompositions across time
%   * figure out why saving to the buffer is so slow in onl_filtered_ORICA
%   * automatically setup paths



function buffer = REST_offline(dataset, calibrationset, pipeline_desc, rest_path)

% check inputs
assert(ischar(dataset) || isstruct(dataset), ...
    'calibrationset must be either a path to a *.set file or an EEGLAB EEG structure.')
assert(ischar(calibrationset) || iscell(calibrationset) || isstruct(calibrationset), ...
    'calibrationset must be either a path to a *.set file, a cell array of inputs for io_loadset, or the output of io_loadset.')
assert(~exist('pipeline_desc', 'var') || isempty(pipeline_desc) || iscell(pipeline_desc) || isstruct(pipeline_desc), ...
    'pipeline_desc must be either a cell array or structure of the format expected by flt_pipeline or left empty or excluded from the function call if the default REST pipeline is desired.')
assert(~exist('rest_path', 'var') || ischar(rest_path), ...
    'rest_path must be a path to the root folder of REST')

% setup environment
bcilab_path = which('bcilab.m');
if ~exist('rest_path', 'var')
    rest_path = '.'; end
if isempty(bcilab_path)
    current_path = pwd;
    addpath(fullfile(rest_path, 'dependencies', 'BCILAB'));
    bcilab
    close
    cd(current_path);
    addpath(genpath(rest_path));
end

% specify data
if ischar(calibrationset)
    calibration = io_loadset(calibrationset);
elseif iscell(calibrationset)
    calibration = io_loadset(calibrationset{:});
elseif isstruct(calibrationset)
    calibration = calibrationset;
end
if ischar(dataset)
    EEG = pop_loadset(dataset);
elseif isstruct(dataset);
    EEG = dataset;
end

% settings
block_size = EEG.srate / 10;
stream_name_base = 'REST_offline';

% create stream
alphabet = [48:57 65:90];
id = char(alphabet(randi(length(alphabet), 1, 8)));
stream_name =  [stream_name_base '_' id];
% stream_name = ['lsl_' stream_name_base '_stream'];
stream = onl_newstream(stream_name, 'srate', EEG.srate, ...
    'chanlocs', EEG.chanlocs);

% describe pipeline
if ~exist('pipeline_desc', 'var') || isempty(pipeline_desc)
    pipeline_desc = { ...
        'selchans', {'channels', {'PO7'  'Packet Counter'  'TRIGGER'}, ...
            'remove_selection', 1}, ...
        'reref', {'ref_chn', {'A2'}}, ...
        'fir', {'fspec', [0.5000 1 50 55], 'stopripple', -40, ...
            'normalize_amplitude', 0}, ...
        'repair_bursts', {'stddev_cutoff', 20, 'window_len', 0.5}, ... 
        'orica', {'onlineWhitening', {'blockSize', 16}, ... % does this need a simulated online switch to save ica matrices and times when they applied?
            'options', {'blockSize', 16}, 'evalConvergence', 0.01} ...
        'eyecatch', 'on', ...
        'ica_reproject', {'ProjectionMatrix', '.icawinv', 'ComponentSubset', '.reject'}, ...
        };
    pipeline_desc = [{'FilterOrdering', strcat('flt_', pipeline_desc(1:2:end))}, pipeline_desc];
end

% create buffer
create_streambuffer(pipeline_desc(3:end));

% run pipline on calibration data
cleaned_data = exp_eval(flt_pipeline(calibration, pipeline_desc{:}));

% initialize the pipeline for streaming data
pipeline = onl_newpipeline(cleaned_data, stream_name);
assignin('base', 'pipeline', pipeline)

% iterate over entire dataset
disp(sprintf('Processing dataset in chunks of %d points...', block_size))
start = 1;
while start < EEG.pnts
    
    % update indicies
    if start + 2 * block_size - 1 > EEG.pnts
        inds = start:EEG.pnts;
        start = EEG.pnts + 1;
    else
        inds = start:start + block_size - 1;
        start = start + block_size;
    end
    
    % update pipeline
    onl_append(stream_name, EEG.data(:, inds));
    onl_filtered_ORICA(nan, nan, stream_name);
    
end

% keep buffer
buffer = evalin('base', ['lsl_' stream_name '_stream']);

% clean base workspace
evalin('base', sprintf('clear temp pipeline %s %s_chunk_clr %s_range lsl_%s_stream', ...
    stream_name, stream_name, stream_name, stream_name))


    function stream = create_streambuffer(pipeline_desc)
        stream.stream_name = stream_name;
        stream.block_size = block_size;
        stream.srate = EEG.srate;                                    
        stream.pnts = EEG.pnts;                                            % number of data points in the buffer
        stream.nbchan = EEG.nbchan;                                 % number of channels in the buffer
        stream.trials = 1;                                                  % number of segments in the buffer (always 1)
        pipe_len = length(pipeline_desc)/2;                                  % # of pipeline outputs
        stream.smax = 0;                                                   % current display options for this stream
        stream.lsllib = lsl_loadlib(env_translatepath('dependencies:/liblsl-Matlab/bin'));
        
        stream.chanlocs{1} = EEG.chanlocs;
        stream.data{1} = zeros(stream.nbchan,stream.pnts,stream.trials);       % the circular buffer storage
        for it = 1:pipe_len
            % extend buffer
            stream.data(it + 1) = stream.data(it);
            % extend chanlocs
            stream.chanlocs(it + 1) = stream.chanlocs(it);
            % reduce size if necessary
            if strcmp('selchans', pipeline_desc{it * 2 - 1})
                ind_rm = find(strcmp('remove_selection', pipeline_desc{it * 2}));
                ind_chan = find(strcmp('channels', pipeline_desc{it * 2}));
                ind_chan_select = ismember({stream.chanlocs{it + 1}.labels}, pipeline_desc{it * 2}{ind_chan + 1});
                if isempty(ind_rm) || ~pipeline_desc{it * 2}{ind_rm + 1}
                    stream.data{it + 1} = stream.data{it + 1}(1:sum(ind_chan_select), :);
                    stream.chanlocs{it + 1} = stream.chanlocs{it + 1}(ind_chan_select);
                else
                    stream.data{it + 1}(1:sum(ind_chan_select), :) = [];
                    stream.chanlocs{it + 1}(ind_chan_select) = [];
                end
            elseif strcmp('reref', pipeline_desc{it * 2 - 1})
                ind_chan = find(strcmp('ref_chn', pipeline_desc{it * 2}));
                ind_chan_select = ismember({stream.chanlocs{it + 1}.labels}, pipeline_desc{it * 2}{ind_chan + 1});
                stream.data{it + 1}(1:sum(ind_chan_select), :) = [];
                stream.chanlocs{it + 1}(ind_chan_select) = [];
            end
            % create ica stream buffer
            if strcmp('orica', pipeline_desc{it * 2 - 1})
                c = size(stream.data{it + 1});
                len = ceil(c(2) / block_size);
                stream.ica.normRn = zeros(1, len);
                stream.ica.icasphere = zeros([c(1) c(1) len]);
                stream.ica.icaweights = zeros([c(1) c(1) len]);
            end
        end
        
        assignin('base', ['lsl_' stream_name '_stream'], stream);
        
    end

end
