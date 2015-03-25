% deal with some weird path 
rmpath('D:\Matlab Coding\VisEEG\BCILAB\dependencies\mexutils-2013-11-11\build-Jordan\');

% define parameters
opts.datapath = 'D:\Matlab Coding\VisEEG\data\';
opts.TrainingDataFile = 'SIM_NSTAT_3sess_16ch_3min.set'; %'20150115_Experiment_rmBadCH_AMICA.set'; % Calibration
opts.PlaybackDataFile = 'SIM_NSTAT_3sess_16ch_3min.set'; %'20150115_Experiment_rmBadCH_AMICA.set'; % Testing
opts.calibEpoch = [0 1]; % use 0-1 second calibration data
opts.winlen = 0; % If winlen=0, pull the most recent (variable length) chunk

opts.BCILAB_PipelineConfigFile = 'ORICA_pipeline_offine_sim_cfg_adapt.mat'; % make sure this file doesn't have 'signal' entry
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
turboRate = 3;
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
catData = []; 
sample = cleaned_data.pnts; lambda = cleaned_data.lambda_k; normRn = cleaned_data.normRn; normRnW = cleaned_data.normRnW; 
lambdaGrad = cleaned_data.lambdaGrad; normCovyy = cleaned_data.normCovyy; normCovfy = cleaned_data.normCovfy; 
covyy = cleaned_data.covyy; covfy = cleaned_data.covfy; RnW = cleaned_data.RnW; minNormRn = cleaned_data.minNormRn;
ratioOfNormRn = cleaned_data.ratioOfNormRn; columnWiseNormRn = cleaned_data.columnWiseNormRn;
while data_len < playbackLength
    % grab data chunk and apply the pipeline 
    [eeg_chunk,pipeline] = onl_filtered(pipeline, chunk_len, opts.silence);

    data_len = data_len + eeg_chunk.pnts;
    %     fprintf([num2str(eeg_chunk.pnts),' ']);
    %     catData = [catData, eeg_chunk.data];
    sample = [sample, data_len];
    lambda = [lambda, eeg_chunk.lambda_k];
    normRn = [normRn, eeg_chunk.normRn];
    minNormRn = [minNormRn, eeg_chunk.minNormRn];
    ratioOfNormRn = [ratioOfNormRn, eeg_chunk.ratioOfNormRn];
    normRnW = [normRnW, eeg_chunk.normRnW];
    lambdaGrad = [lambdaGrad, eeg_chunk.lambdaGrad];
    normCovyy = [normCovyy, eeg_chunk.normCovyy];
    normCovfy = [normCovfy, eeg_chunk.normCovfy];
    covyy = cat(3,covyy,eeg_chunk.covyy);
    covfy = cat(3,covfy, eeg_chunk.covfy);
    RnW = cat(3,RnW, eeg_chunk.RnW);
    columnWiseNormRn = [columnWiseNormRn, eeg_chunk.columnWiseNormRn];
    
    if storeIdx < floor(data_len/playbackSrate) && ~isempty(eeg_chunk.icaweights)
        results(end+1).time = data_len; % sec
        results(end).icaweights = eeg_chunk.icaweights;
        results(end).icasphere = eeg_chunk.icasphere;
        results(end).statIdx = eeg_chunk.statIdx(1);
        results(end).mir = eeg_chunk.mir(1);
        results(end).lambda = eeg_chunk.lambda_k(1);
        if isfield(eeg_chunk,'normRn'),     results(end).normRn = eeg_chunk.normRn(1); end
        if isfield(eeg_chunk,'normRnW'),    results(end).normRnW = eeg_chunk.normRnW(1); end
        if isfield(eeg_chunk,'lambdaGrad'), results(end).lambdaGrad = eeg_chunk.lambdaGrad(1); end
        if isfield(eeg_chunk,'normCovyy'),  results(end).normCovyy = eeg_chunk.normCovyy(1); end
        if isfield(eeg_chunk,'normCovfy'),  results(end).normCovfy = eeg_chunk.normCovfy(1); end
        
        fprintf([num2str(storeIdx+1),' ']);
        storeIdx = storeIdx + storeInterval;
    end
    
end
disp('Done');

% eeg_chunk.icawinv = inv(eeg_chunk.icaweights * eeg_chunk.icasphere);
% pop_topoplot(eeg_chunk,0);
% save([opts.datapath 'sim_3sess_16ch_3min_len1000.mat'],'results')

%% visualize results


