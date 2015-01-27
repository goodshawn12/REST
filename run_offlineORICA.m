
% define parameters
opts.datapath = 'D:\Matlab Coding\VisEEG\data\';
opts.TrainingDataFile = '20150115_Calibration.set'; % Calibration
opts.PlaybackDataFile = '20150115_Experiment_raw_icainfo.set'; % Testing
opts.calibEpoch = [0 1]; % use 0-1 second calibration data
opts.winlen = 0; % If winlen=0, pull the most recent (variable length) chunk

opts.BCILAB_PipelineConfigFile = 'ORICA_pipeline_offine_cfg.mat'; % make sure this file doesn't have 'signal' entry
opts.lsl.StreamName          = 'EEGDATA';
opts.silence = true;

%% load calibration recording
calibData = exp_eval_optimized(io_loadset([opts.datapath opts.TrainingDataFile], ...
    'markerchannel',{'remove_eventchns',false}));

% select calibration data legnth
calibData = pop_select(calibData,'time',opts.calibEpoch);

% load playback data information
playbackData = exp_eval_optimized(io_loadset([opts.datapath opts.PlaybackDataFile]));
playbackLength = playbackData.pnts;
playbackSrate = playbackData.srate;
clear playbackData

%% define a pipeline
try    fltPipCfg = exp_eval(io_load([opts.datapath opts.BCILAB_PipelineConfigFile]));
catch, disp('-- no existing pipeline --'); fltPipCfg = {}; end

fltPipCfg = arg_guipanel('Function',@flt_pipeline,'Parameters',[{'signal',calibData} fltPipCfg],'PanelOnly',false);

if ~isempty(fltPipCfg)
    save(env_translatepath([opts.datapath opts.BCILAB_PipelineConfigFile]),'-struct','fltPipCfg'); end

%% apply the pipeline to calibration data
disp('-- Applying pipeline to calibration data (please wait) --');
cleaned_data = exp_eval(flt_pipeline('signal',calibData,fltPipCfg));


%% start streaming playback data 
turboRate = 1;
updateFrequency = 20 * turboRate;
run_readdataset_turbo('MatlabStream',opts.lsl.StreamName,                         ...
    'Dataset',io_loadset([opts.datapath opts.PlaybackDataFile],'markerchannel',{'remove_eventchns',false}), ...
    'UpdateFrequency', updateFrequency, ...
    'TurboRate', turboRate);

% initialize the pipeline for streaming data
pipeline = onl_newpipeline(cleaned_data,opts.lsl.StreamName);

% inspect streams
% gui_vis_filtered;

%% online processing of the pipeline
disp('-- Running pipeline --');
fprintf('Data processed (in secs), %f secs in total:\n',playbackLength/playbackSrate);
chunk_len = round(opts.winlen*cleaned_data.srate);

% run single pass over the playback data
storeInterval = 5; % sec
results = struct('time',[],'icaweights',{},'icasphere',{},'statIdx',[],'mir',[],'lambda',[]);
data_len = 0; dispMinute = 0; storeIdx = -1;
while data_len < playbackLength
    % grab data chunk and apply the pipeline 
    [eeg_chunk,pipeline] = onl_filtered(pipeline, chunk_len, opts.silence);

    data_len = data_len + eeg_chunk.pnts;
    
    if storeIdx < floor(data_len/playbackSrate)
        results(end+1).time = data_len; % sec
        results(end+1).icaweights = eeg_chunk.icaweights;
        results(end+1).icasphere = eeg_chunk.icasphere;
        results(end+1).statIdx = eeg_chunk.statIdx(1);
        results(end+1).mir = eeg_chunk.mir(1);
        results(end+1).lambda = eeg_chunk.lambda_k(1);
        fprintf([num2str(storeIdx+1),' ']);
        storeIdx = storeIdx + storeInterval;
    end
    
%     if dispMinute ~= floor(data_len/playbackSrate/60)
%         dispMinute = dispMinute + 1;
%         fprintf([num2str(dispMinute),'-']);
%     end
%     fprintf([num2str(eeg_chunk.pnts),' ']);
end
disp('Done');

eeg_chunk.icawinv = inv(eeg_chunk.icaweights * eeg_chunk.icasphere);
% pop_topoplot(eeg_chunk,0);
save([opts.datapath '20150115_result_decayW8B8_turbo2.mat'],'eeg_chunk')

%% visualize results





