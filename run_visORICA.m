%% Run visGUI

% load channel locations
load vis/chanlocs14

% load calibration data
calibData = exp_eval_optimized(io_loadset('data/EmotivTrain_EyeClose_icainfo.set', ...
    'markerchannel',{'remove_eventchns',false}));
% select data
% calibData = set_selinterval(calibData, opts.calibEpoch);
% calibData = pop_loadset('data/EmotivTrain_EyeClose_icainfo.set');
calibData = pop_select(calibData,'time',[0 1]);

visORICA(calibData)