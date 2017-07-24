function  [kurt, varkurt]= kurtosis_temporal(varargin)
% kurtosis_temporal calculates the mean and variance of the kurtosis in
% one-second intervals of epochs in IC time series.
%
% Input:
% eeg: EEGLab data structure with additional fields
% virtual_srate and icaact_filtered_resampled. See
% the readme file for the toolbox for details on these fields.
%
% celleegs (optional): Cell array containing EEGLab data structures with
% the same ICA decomposition.
%
% Output:
% kurt: Vector of length equal to the number of ICs in the input eeg. The
% vector contains the means of the kurtosis calculated for one-second
% intervals for each IC time series.
%
% varkurt: Vector of length equal to the number of ICs in the input eeg.
% The vector contains the variance of the kurtosis calculated for
% one-second intervals for each IC time series.

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

K_all = [];
for ieeg = 1:length(celleegs)
    eeg = celleegs{ieeg};
    icaacts = eeg.icaact_filtered_resampled;
    icaacts = temporally_normalize_icdecomp(icaacts);
    interval = 1; % calculate over intervals of one second
    
    interval_points = ceil(interval*eeg.virtual_srate);
    interval_starts = 1:interval_points:size(eeg.icaact_filtered_resampled,2);
    n_intervals = length(interval_starts) -1; % subtract one to get the number of whole intervals
    
    num_epoch = size(icaacts, 3);
    n_ics = size(icaacts, 1);
    K = zeros(n_ics, n_intervals, num_epoch);
    
    for k=1:n_ics
        for j=1:num_epoch
            for i=1:(n_intervals)
                K(k, i, j) = kurtosis(...
                    icaacts(k, interval_starts(i):(interval_starts(i+1)-1), j));
            end
        end
    end
    
    K_all = [K_all reshape(K, size(K,1), size(K,2)*size(K,3))];
end
kurt = mean(K_all,2);
varkurt = var(K_all,[], 2);
