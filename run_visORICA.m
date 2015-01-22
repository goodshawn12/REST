%% Run visGUI

addpath(genpath('dependencies'))

Emotiv = 1;
Cognionics = 0;

% load channel locations
if Emotiv
    load data/chanlocs14
else
    load data/chanlocs64
end

% load calibration data
if Emotiv
    calibData = exp_eval_optimized(io_loadset('data/EmotivTrain_EyeClose_icainfo.set', ...
        'markerchannel',{'remove_eventchns',false}));
    hmObj = load('vis/emotivHeadModel');
    hmLoc = load('vis/emotivLFMDerivatives');
    calibData.headModel = hmObj.emotivHeadModel;
    calibData.localization = hmLoc;
else
    calibData = exp_eval_optimized(io_loadset('D:\Matlab Coding\VisEEG\data\20150115_Calibration.set', ...
        'markerchannel',{'remove_eventchns',false}));
%     calibData = exp_eval_optimized(io_loadset('data/EmotivTrain_EyeClose_icainfo.set', ...
%         'markerchannel',{'remove_eventchns',false}));
    hmObj = load('vis/cognionicsHeadModel');
    hmLoc = load('vis/cognionicssLORETA');
    calibData.headModel = hmObj.cognionicsHeadModel;
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