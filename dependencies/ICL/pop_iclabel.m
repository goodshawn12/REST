function [EEG, LASTCOM] = pop_iclabel(EEG)

% get parameters
% inputdlg?

% run iclabel
EEG = iclabel(EEG);
LASTCOM = 'EEG = pop_iclabel(EEG);';

% visualize?
% pass to 