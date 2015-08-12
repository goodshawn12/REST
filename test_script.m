%% Run REST
clear
close all
delete(timerfind)

addpath(genpath(['.' filesep]));

playback = 1;

% load calibration data
data_location = ['data' filesep 'Demo_EmotivEPOC_EyeClose.set'];
opts.headModel = ['data' filesep 'head_models' filesep 'emotivHeadModel_file'];
opts.calibration_data = exp_eval_optimized(io_loadset(data_location, ...
        'markerchannel',{'remove_eventchns',false}));

opts.calibration_data = pop_select(opts.calibration_data,'time',[0 1]);

if playback
    if exist('playbackStream','var')
        try %#ok<TRYNC>
            stop(playbackStream)
            delete(playbackStream)
        end
    end
    EEG = pop_loadset(data_location);
    playbackStream = play_eegset_lsl(EEG,'REST_test_data','REST_test_markers',[],true);
end

% start REST
REST(opts)

% if playback
%     stop(playbackStream)
%     delete(playbackStream)
% end