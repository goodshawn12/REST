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
timer = timerfindall;
if ~isempty(timer)
    stop(timerfindall); delete(timerfindall); disp('Delete timers...'); end
close all
clear

%% define opts structure
% whether to customize pipeline 
opts.customize_pipeline = true;
opts.save_config = true;

% % Emotiv - setting
% (optional) define config file name
opts.config = 'Config_ORICA_EmotivEPOC';

% point to headModel
opts.headModel = ['data' filesep 'head_models' filesep 'emotivHeadModel_file'];

% (optional) path to calibration data and select time window
opts.calibration_data = ['data' filesep 'Demo_EmotivEPOC_EyeOpen.set'];
opts.calibration_window = [0,60]; % sec


% % Quick 20 - setting
%{
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
%}

% (optional) load eyeCatch library
opts.libEyeCatch = load(['dependencies' filesep 'eyeCatch' filesep 'libEyeCatch.mat']);

% (optional) load IC_MARC model
opts.modIcMarc = load(['dependencies' filesep 'IC_MARC' filesep 'spatial2.mat']);
if isfield(opts, 'modIcMarc')
    labels = {'Fp1', 'Fp2', 'F3', 'F4', 'C3', 'C4', 'P3', 'P4', 'O1', 'O2', 'F7', 'F8', 'T7', 'T8', 'P7', 'P8', 'Fz', 'Pz', 'FC1', 'FC2', 'CP1',...
    'CP2', 'FC5', 'FC6', 'CP5', 'CP6', 'F9', 'F10', 'TP9', 'TP10', 'Cz', 'Oz', 'F1', 'F2', 'C1', 'C2', 'P1', 'P2', 'AF3', 'AF4',...
    'FC3', 'FC4', 'CP3', 'CP4', 'PO3', 'PO4', 'F5', 'F6', 'C5', 'C6', 'P5', 'P6', 'AF7', 'AF8', 'FT7', 'FT8', 'TP7', 'TP8', 'PO7', ...
    'PO8', 'AFz', 'FCz', 'CPz', 'POz'};
    virtual_chanlocs_struct = struct('labels', lower(labels));
    virtual_chanlocs = pop_chanedit(virtual_chanlocs_struct, ...
        'lookup', 'standard-10-5-cap385.elp');
    for i=1:length(virtual_chanlocs)
        % set head radius to 9 cm
        virtual_chanlocs=pop_chanedit(virtual_chanlocs, 'changefield',{i 'sph_radius' '9'},'convert',{'sph2all'});
    end
    opts.virtual_chanlocs = virtual_chanlocs;
    
    opts.cdn_dipolefit = load('dipolfit_matrix'); % loads M100, clab
    
end


% use playback data
opts.playback = 1;

%% start REST
REST(opts)