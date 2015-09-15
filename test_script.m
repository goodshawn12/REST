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

%% define opts structure
opts.customize_pipeline = true;
opts.config = 'Config_ORICA_quick20';

% load channel location file
load(['data' filesep 'chanlocs' filesep 'Quick20.mat']);
opts.chanlocs = chanlocs;

% point to headModel
% opts.headModel = ['data' filesep 'head_models' filesep 'emotivHeadModel_file'];
opts.headModel = ['data' filesep 'head_models' filesep 'quick20HeadModel'];

% (optional) path to calibration data and select time window
% opts.calibration_data = ['data' filesep 'Demo_EmotivEPOC_EyeClose.set'];
opts.calibration_data = ['data' filesep 'quick20_calib_EyeOpen.set'];
opts.calibration_window = [0,60]; % sec

% use playback data
opts.playback = 1;

%% start REST
REST(opts)

%% delete playback stream
% if playback
%     stop(playbackStream)
%     delete(playbackStream)
% end