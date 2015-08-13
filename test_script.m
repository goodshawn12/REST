%% Test REST

%% refresh workspace
% clear all will break bcilab and require it to restart as it uses global
% variables 
clear
close all
delete(timerfind)

%% set path and 
addpath(genpath(['.' filesep]));
% start bcilab or add it to the path as well

%% load data and create opts structure
% point to headModel
opts.headModel = ['data' filesep 'head_models' filesep 'emotivHeadModel_file'];

% load calibration data
data_location = ['data' filesep 'Demo_EmotivEPOC_EyeClose.set'];
opts.calibration_data = exp_eval_optimized(io_loadset(data_location, ...
        'markerchannel',{'remove_eventchns',false}));
opts.calibration_data = pop_select(opts.calibration_data,'time',[0 10]);

%% create LSL stream for online EEG data playback
playback = 1;
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

%% start REST
REST(opts)

%% delete playback stream
if playback
    stop(playbackStream)
    delete(playbackStream)
end