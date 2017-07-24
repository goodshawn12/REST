function [avgskew1s, avgskew15s, varskew1s, varskew15s] = avg_skew(varargin)
% avg_skew returns the average skewness in 1 second and 5 second intervals of each IC time series in an EEG
% data structure.
%
% Input:
% EEGLab data structure with additional fields icaact_filtered_resampled
% and virtual_srate. See the readme file for the toolbox for details on
% these fields.
%
% celleegs (optional): Cell array containing EEGLab data structures with
% the same ICA decomposition.
%
% Output:
% avgskew1s: Column vector containing, for each IC, the average skewness in one-second
% intervals of the IC's time series.
%
% avgskew15s: Column vector containing, for each IC, the average skewness in 15-second
% intervals of the IC's time series.
%
% varskew1s: Column vector containing, for each IC, the variance of skewness over one-second
% intervals of the IC's time series.
%
% varskew15s: Column vector containing, for each IC, the average skewness over 15-second
% intervals of the IC's time series.

% Copyright (C) 2013  Laura Froelich (laura.frolich@gmail.com)
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

arg_define(1,varargin, ...
    arg({'eeg', 'EEG'},[],[],'EEGLab struct containing EEG data.'),...
    arg('celleegs', [], [],...
    'EEGLab data structures in cell array.'));

if isempty(celleegs)
    celleegs{1} = eeg;
end

eeg_skew1s_all = [];
eeg_skew15s_all = [];

for ieeg = 1:length(celleegs)
    eeg = celleegs{ieeg};
    eeg.srate = eeg.virtual_srate;
    eeg.icaact = eeg.icaact_filtered_resampled;
    eeg.icaact = temporally_normalize_icdecomp(eeg.icaact);
    
    eeg_skew1s = eeg_skewness(eeg, 'interval', 1);
    eeg_skew1s_all=[eeg_skew1s_all reshape(eeg_skew1s, size(eeg_skew1s,1), size(eeg_skew1s,2)*size(eeg_skew1s,3))];
    
    eeg_skew15s = eeg_skewness(eeg, 'interval', 15);
    eeg_skew15s_all=[eeg_skew15s_all reshape(eeg_skew15s, size(eeg_skew15s,1), size(eeg_skew15s,2)*size(eeg_skew15s,3))];
end

avgskew1s = mean(eeg_skew1s_all,2);
varskew1s = var(eeg_skew1s_all, [],2);
avgskew15s = mean(eeg_skew15s_all,2);
varskew15s = var(eeg_skew15s_all,[],2);
end
