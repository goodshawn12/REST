%% Test REST (go to 'REST' folder)
%% set path and initialize bcilab
rest_path = './';
bcilab_path = which('bcilab.m');
if isempty(bcilab_path)
    current_path = pwd;
    addpath(fullfile(rest_path, 'dependencies', 'BCILAB'));
    bcilab
    cd(current_path);
    addpath(genpath(rest_path));
end

%% refresh workspace
% clear all will break bcilab and require it to restart as it uses global variables 
close all hidden
clearvars -except rest_path
timer = timerfindall;
if ~isempty(timer)
    stop(timerfindall); delete(timerfindall);
    disp('Delete timers...');
end

%% define opts structure
% whether to customize pipeline 
opts.customize_pipeline = true;
opts.save_config = false;

% % % Emotiv - setting
% % (optional) define config file name
% opts.config = 'Config_ORICA_EmotivEPOC';
% 
% % point to headModel
% opts.headModel = fullfile(rest_path, 'data', 'head_models', 'emotivHeadModel_file');
% 
% % (optional) path to calibration data and select time window
% opts.calibration_data = fullfile(rest_path, 'data', 'Demo_EmotivEPOC_EyeOpen.set');
% opts.calibration_window = [0,60]; % sec


% % % Quick 20 - setting
% % channel locations
% load(fullfile(rest_path, 'data', 'chanlocs', 'Quick20.mat'));
% opts.chanlocs = chanlocs; 
% 
% % (optional) define config file name
% opts.config = 'Config_ORICA_quick20';
% 
% % point to headModel
% opts.headModel = fullfile(rest_path, 'data', 'head_models', 'quick20HeadModel');
% 
% % (optional) path to calibration data and select time window
% opts.calibration_data = fullfile(rest_path, 'data', 'Quick20_Luca_calib_EyeOpen.set');
% opts.calibration_window = [0,60]; % sec


% % Quick 30 - setting
% channel locations
load(fullfile(rest_path, 'data/chanlocs/Quick30.mat'));
opts.chanlocs = chanlocs; 

% (optional) define config file name
opts.config = 'Config_ORICA_quick30';

% (optional) path to calibration data and select time window
opts.calibration_data = fullfile(rest_path, 'data', 'Quick30_Shawn_EyeOpen.set');
opts.calibration_window = [0,133]; % sec

% (optional) load eyeCatch library
opts.libEyeCatch = load(fullfile(rest_path, 'dependencies', 'eyeCatch', 'libEyeCatch.mat'));

% (optional) load IC_MARC model
opts.modIcMarc = load(fullfile(rest_path, 'dependencies', 'IC_MARC', 'spatial2.mat'));

% use playback data
opts.playback = 1;

%% start REST
REST(opts)
