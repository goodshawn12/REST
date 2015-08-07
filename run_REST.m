%% Run visGUI
clear
close all
delete(timerfind)

addpath(genpath('dependencies'))
addpath(genpath('utility'))
addpath('vis')

Emotiv = 1;
Cognionics = 0;
test16 = 0;
test64 = 0;

playback = 1;

% load channel locations
if Emotiv
    load data/chanlocs14
elseif Cognionics %#ok<*UNRCH>
    load data/chanlocs64
elseif test16
    load vis/artificial_data/16ch/chanlocs_test16.mat
elseif test64
    load vis/artificial_data/64ch/test64chanlocs.mat
end

% load calibration data
if Emotiv
    data_location = 'data/EmotivTrain_EyeClose_icainfo.set';
    opts.headModel = 'head_models/emotivHeadModel_file';
%     calibData.localization = hmLoc;
elseif Cognionics
%     data_location = '/home/lpiontonachini/Dropbox/School/Research/VisEEG_local/20150115_Experiment.set';
    data_location = 'data/20150115_Experiment.set';
    opts.headModel = 'head_models/cognionicsHeadModel_file';
%     calibData.localization = hmLoc;
elseif test16
    data_location = 'vis/artificial_data/16ch/SIM_NSTAT_3sess_16ch_3min.set';
    opts.headModel = 'head_models/sim16HeadModel_file';
%     calibData.localization = hmLoc;
elseif test64
    data_location = 'vis/artificial_data/64ch/SIM_STAT_64ch_10min.set';
    opts.headModel = 'head_models/sim64HeadModel_file';
%     calibData.localization = hmLoc;
end
opts.calibration_data = exp_eval_optimized(io_loadset(data_location, ...
        'markerchannel',{'remove_eventchns',false}));

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
opts.calibration_data = pop_select(opts.calibration_data,'time',[0 1]);
W = [];sphere = [];



if playback
    if exist('playbackStream','var')
        try
            stop(playbackStream)
            delete(playbackStream)
        end
    end
    EEG = pop_loadset(data_location);
    playbackStream = play_eegset_lsl(EEG,'REST_test_data','REST_test_markers',[],true);
end
REST(opts)

% if playback
%     stop(playbackStream)
%     delete(playbackStream)
% end