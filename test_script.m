%% Test REST (go to 'REST' folder)
%% set path and initialize bcilab
bcilab_path = which('bcilab.m');
if isempty(bcilab_path)
    current_path = pwd;
    addpath('dependencies\BCILAB'); bcilab
    cd(current_path);
    addpath(genpath('./'));
end

%% refresh workspace
% clear all will break bcilab and require it to restart as it uses global
% variables 
close all
delete(timerfind)
clear

%% define opts structure
% whether to customize pipeline 
opts.customize_pipeline = true;

% (optional) define config file name
opts.config = 'Config_ORICA_EmotivEPOC';

% point to headModel
opts.headModel = ['data' filesep 'head_models' filesep 'emotivHeadModel_file'];

% (optional) path to calibration data and select time window
opts.calibration_data = ['data' filesep 'Demo_EmotivEPOC_EyeClose.set'];
opts.calibration_window = [0,60]; % sec

% use playback data
opts.playback = 1;

%% start REST
REST(opts)
