%% PIPELINE FOR RESONANCE_MEG SOURCE RECONSTRUCTION

clear all
% Ensure that all the toolboxes are added to matlab path
addpath('D:\Documents and settings\bgauthie\MNE-2.7.0-3106-Linux-i686\share\matlab'); % mne toolbox contains useful low level functions
PSpace  = 'E:\MEG\'; 
addpath ([PSpace 'SCRIPTS/Matlab_fieldtrip']);% add personal functions and copies of "private" subfolders
addpath ([PSpace 'SCRIPTS/Matlab_fieldtrip/priv']);
addpath ([PSpace 'SCRIPTS/Matlab_fieldtrip/priv2']);
addpath ([PSpace 'SCRIPTS/Matlab_fieldtrip/functions_pipeline']);
addpath ([PSpace 'SCRIPTS/Matlab_fieldtrip/functions_pipeline_2']);
addpath ([PSpace 'SCRIPTS/fieldtrip_corrections']);
% addpath ([PSpace 'SCRIPTS/Matlab_fieldtrip/my_toolbox']);
addpath 'E:\spm2\';
addpath 'D:\Documents and settings\bgauthie\FIELDTRIP\fieldtrip-20100420\fileio\'
addpath 'D:\Documents and settings\bgauthie\FIELDTRIP\fieldtrip-20100420\public'
%fieldtrip
ft_defaults

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subject                     = '20091110_sm';
% AnatMri                     = 'ptk_server_task_1/Examme070045-1530_001_Series000002_flip2.hdr'; 
Subject                     = '20100205_vvw';
AnatMri                     = 'ptk_server_task_4/Examvv100048-1620_001_Series000002_flip2.hdr';
ndataset                    = 6;
ncond                       = 4;
Condlist                    = {'cond50','cond100','cond200','cond400'};
EEGtag                      = 0;
MEGtag                      = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute PCA matrix correction for ECG and EOGv artifacts 
[Mgrads, Mmags, Meeg]       = PCA_matrix(PSpace,Subject,'_run1_raw_trans_sss.fif',MEGtag, EEGtag); % vvw
% [Mgrads, Mmags, Meeg]       = PCA_matrix(PSpace,Subject,'_run3_raw_trans_sss.fif',MEGtag, EEGtag); % seb

Neurospin_Epoching_fornutmeg(Subject,Mgrads,Mmags,Meeg,ndataset,Condlist,InputDataDir,OutputDataDir)



