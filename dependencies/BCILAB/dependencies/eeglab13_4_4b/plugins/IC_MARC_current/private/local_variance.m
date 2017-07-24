function [mean_variance, var_variance] = local_variance(varargin)
% local_variance calculates the mean and variance of the variance in
% intervals of epochs in IC time series.
%
% Input:
% eeg: EEGLab data structure with additional fields
% virtual_srate and icaact_filtered_resampled. See
% the readme file for the toolbox for details on these fields.
%
% interval (optional): The length in seconds of intervals over which the
% variance is to be calculated (default: 1).
%
% celleegs (optional): Cell array containing EEGLab data structures with
% the same ICA decomposition.
%
% Output:
% mean_variance: Vector of length equal to the number of ICs in the input
% eeg. The vector contains the means of the variances calculated in
% intervals of each IC time series.
%
% var_variance: Vector of length equal to the number of ICs in the input
% eeg. The vector contains the variances of the variances calculated in
% intervals of each IC time series.

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

arg_define([0 1], varargin, ...
    arg({'eeg', 'EEG'},[],[],'EEGLab struct containing EEG data.'),...
    arg('interval', 1, [],...
    'Interval in seconds in which variance is calculated.'),...
    arg('celleegs', [], [],...
    'EEGLab data structures in cell array.'));

if isempty(celleegs)
    celleegs{1} = eeg;
end

vars_all = [];
for ieeg = 1:length(celleegs)
    eeg = celleegs{ieeg};        
    eeg.icaact = temporally_normalize_icdecomp(eeg.icaact_filtered_resampled);
    
    num_epoch = size(eeg.icaact,3);
    n_ics = size(eeg.icaact, 1);
    interval_points = ceil(interval*eeg.virtual_srate);
    interval_starts = 1:interval_points:size(eeg.icaact,2);
    n_intervals = length(interval_starts)-1; % subtract one to get the number of whole intervals
    vars = zeros(n_ics, n_intervals, num_epoch);
   
    % split each epoch for each ic into intervals for which the variance is calculated.
    for k=1:n_ics
        for j=1:num_epoch
            for i=1:n_intervals
                vars(k, i, j) = var(...
                    eeg.icaact(k, interval_starts(i):(interval_starts(i+1)-1), j));
            end
        end
    end
    
    vars_all = [vars_all reshape(vars, size(vars,1), size(vars,2)*size(vars,3))];
end

mean_variance = mean(vars_all,2);
var_variance = var(vars_all,[],2);
end
