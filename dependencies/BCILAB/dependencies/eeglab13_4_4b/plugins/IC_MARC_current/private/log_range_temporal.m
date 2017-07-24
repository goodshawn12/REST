function [log_range, log_rangevar] = log_range_temporal(varargin)
% local_variance calculates the mean and variance of the variance in
% intervals of epochs in IC time series.
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
% log_range: Vector of length equal to the number of ICs in the input
% eeg. The vector contains the means of the logarithms of the ranges
% calculated in each one-second interval of each IC time series.
%
% log_rangevar: Vector of length equal to the number of ICs in the input
% eeg. The vector contains the variances of the logarithms of the ranges
% calculated in each one-second interval of each IC time series.

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

interval = 1; % the log(range()) is calculated over intervals of one second

interval_points = ceil(interval*eeg.virtual_srate);
lograngs = [];

for ieeg = 1:length(celleegs)
    eeg = celleegs{ieeg};
    icaacts = eeg.icaact_filtered_resampled;
    icaacts = temporally_normalize_icdecomp(icaacts);
    
    interval_starts = 1:interval_points:size(eeg.icaact_filtered_resampled,2);
    n_intervals = length(interval_starts)-1; % subtract one to get the number of whole intervals
    
    num_epoch = size(icaacts, 3);
    n_ics = size(icaacts, 1);
    lograngs_temp = zeros(n_ics, n_intervals, num_epoch);
    
    % split each epoch for each ic into intervals for which the variance is calculated.
    for k=1:n_ics
        for j=1:num_epoch
            for i=1:(n_intervals)
                lograngs_temp(k, i, j) = log(range(...
                    icaacts(k, interval_starts(i):(interval_starts(i+1)-1), j)));
            end
        end
    end
    
    lograngs = [lograngs reshape(lograngs_temp, size(lograngs_temp,1), size(lograngs_temp,2)*size(lograngs_temp,3))];
end

log_range = mean(lograngs,2);
log_rangevar = var(lograngs,[],2);
