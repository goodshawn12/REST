%% Test REST (go to 'REST' folder)
%% set path and initialize bcilab
bcilab_path = which('bcilab.m');
if isempty(bcilab_path)
    current_path = pwd;
    addpath(['dependencies' filesep 'BCILAB']); bcilab
    cd(current_path);
    addpath(genpath('./'));
end

%% refresh workspace
% clear all will break bcilab and require it to restart as it uses global variables 
close all hidden
clear
timer = timerfindall;
if ~isempty(timer)
    stop(timerfindall); delete(timerfindall); disp('Delete timers...'); end

%% define opts structure
% whether to customize pipeline 
opts.customize_pipeline = true;
opts.save_config = true;

% % % Emotiv - setting
% % (optional) define config file name
% opts.config = 'Config_ORICA_EmotivEPOC';
% 
% % point to headModel
% opts.headModel = ['data' filesep 'head_models' filesep 'emotivHeadModel_file'];
% 
% % (optional) path to calibration data and select time window
% opts.calibration_data = ['data' filesep 'Demo_EmotivEPOC_EyeOpen.set'];
% opts.calibration_window = [0,60]; % sec


% % Quick 20 - setting

% channel locations
load('data/chanlocs/Quick20.mat');
opts.chanlocs = chanlocs; 

% (optional) define config file name
opts.config = 'Config_ORICA_quick20';

% point to headModel
opts.headModel = ['data' filesep 'head_models' filesep 'quick20HeadModel'];

% (optional) path to calibration data and select time window
opts.calibration_data = ['data' filesep 'Quick20_Luca_calib_EyeOpen.set'];
opts.calibration_window = [0,60]; % sec


% (optional) load eyeCatch library
opts.libEyeCatch = load(['dependencies' filesep 'eyeCatch' filesep 'libEyeCatch.mat']);

% (optional) load IC_MARC model
opts.modIcMarc = load(['dependencies' filesep 'IC_MARC' filesep 'spatial2.mat']);

% use playback data
opts.playback = 1;

%% start REST
REST(opts)