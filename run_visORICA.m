%% Run visGUI
addpath('vis')

% load channel locations
% load data/chanlocs14
load data/chanlocs64

% load emotive headModel data
% hmObj = load('vis/emotivHeadModel');
% hmLoc = load('vis/emotivLFMDerivatives');

% load cognionics headModel data
hmObj = load('vis/cognionicsHeadModel');
hmLoc = load('vis/cognionicssLORETA');


% load calibration data
calibData = exp_eval_optimized(io_loadset('data/EmotivTrain_EyeClose_icainfo.set', ...
    'markerchannel',{'remove_eventchns',false}));

% % attach headModel data to calibration data
% calibData.headModel = hmObj.emotivHeadModel;
% calibData.localization = hmLoc;


% load 61-ch Flanker Task dataset
% calibData = exp_eval_optimized(io_loadset('D:\Matlab Coding\OnlineICA\BCILAB\userdata\Flanker_SH_cleaned_test_icainfo.set', ...
%     'markerchannel',{'remove_eventchns',false}));

% load 64-ch simulated dataset
% calibData = exp_eval_optimized(io_loadset('D:\Matlab Coding\OnlineICA\BCILAB\userdata\Colin27_Biosemi_1010_0_norot.set', ...
%     'markerchannel',{'remove_eventchns',false}));

% load 64-ch Cognionics Calibration Data
% calibData = exp_eval_optimized(io_loadset('D:\Matlab Coding\VisEEG\data\20141216_calib.set', ...
%     'markerchannel',{'remove_eventchns',false}));
calibData = exp_eval_optimized(io_loadset('/home/lpiontonachini/Dropbox/School/Research/VisEEG_local/20150115_Experiment.set', ...
    'markerchannel',{'remove_eventchns',false}));
% attach headModel data to calibration data
calibData.headModel = hmObj.cognionicsHeadModel;
calibData.localization = hmLoc;

% calibData.event = [];
% chanlocs = calibData.chanlocs;

% [filename, pathname] = uigetfile('*.set');

% select data
% calibData = set_selinterval(calibData, opts.calibEpoch);
% calibData = pop_loadset('data/EmotivTrain_EyeClose_icainfo.set');
calibData = pop_select(calibData,'time',[0 60]);

visORICA(calibData)