%% Run visGUI

% addpath(genpath('dependencies'))
addpath(genpath('utility'))
addpath('vis')

Emotiv = 0;
Cognionics = 0;
test16 = 1;

% load channel locations
if Emotiv
    load data/chanlocs14
elseif Cognionics %#ok<*UNRCH>
    load data/chanlocs64
elseif test16
    load vis/artificial_data_16/chanlocs_test16.mat
end

% load calibration data
if Emotiv
    calibData = exp_eval_optimized(io_loadset('data/EmotivTrain_EyeClose_icainfo.set', ...
        'markerchannel',{'remove_eventchns',false}));
    hmObj = load('vis/emotivHeadModel');
    hmLoc = load('vis/emotivLFMDerivatives');
    calibData.headModel = hmObj.emotivHeadModel;
    calibData.localization = hmLoc;
elseif Cognionics
%     calibData = exp_eval_optimized(io_loadset('D:\Matlab Coding\VisEEG\data\20150115_Calibration.set', ...
%         'markerchannel',{'remove_eventchns',false}));
    calibData = exp_eval_optimized(io_loadset('/home/lpiontonachini/Dropbox/School/Research/VisEEG_local/20150115_Calibration.set', ...
        'markerchannel',{'remove_eventchns',false}));
    hmObj = load('vis/cognionicsHeadModel');
    hmLoc = load('vis/cognionicssLORETA');
    calibData.headModel = hmObj.cognionicsHeadModel;
    calibData.localization = hmLoc;
elseif test16
    calibData = exp_eval_optimized(io_loadset('vis/artificial_data_16/SIM_NSTAT_3sess_16ch_3min.set', ...
        'markerchannel',{'remove_eventchns',false}));
    hmObj = load('vis/artificial_data_16/hmTest16.mat');
    hmLoc = load('vis/artificial_data_16/test16LORETA.mat');
    calibData.headModel = hmObj.hmTest;
    calibData.localization = hmLoc;
end


% load 61-ch Flanker Task dataset
% calibData = exp_eval_optimized(io_loadset('D:\Matlab Coding\OnlineICA\BCILAB\userdata\Flanker_SH_cleaned_test_icainfo.set', ...
%     'markerchannel',{'remove_eventchns',false}));

% load 64-ch simulated dataset
% calibData = exp_eval_optimized(io_loadset('D:\Matlab Coding\OnlineICA\BCILAB\userdata\Colin27_Biosemi_1010_0_norot.set', ...
%     'markerchannel',{'remove_eventchns',false}));

% load 64-ch Cognionics Calibration Data
% calibData = exp_eval_optimized(io_loadset('D:\Matlab Coding\VisEEG\data\20141216_calib.set', ...
%     'markerchannel',{'remove_eventchns',false}));

% calibData.event = [];
% chanlocs = calibData.chanlocs;

% [filename, pathname] = uigetfile('*.set');

% select data
% calibData = set_selinterval(calibData, opts.calibEpoch);
% calibData = pop_loadset('data/EmotivTrain_EyeClose_icainfo.set');
calibData = pop_select(calibData,'time',[0 60]);

visORICA(calibData)