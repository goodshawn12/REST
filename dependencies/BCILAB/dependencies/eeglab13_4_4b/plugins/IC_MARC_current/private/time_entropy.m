function [meanent, varent] = time_entropy( varargin )
% time_entropy calculates, for each IC, the mean and variance of the
% entropy over time of the IC's time series.
%
% Input:
% eeg: EEGLab data structure with the additional fields
% virtual_srate and icaact_filtered_resampled. See
% the readme file for the toolbox for details on these fields.
%
% celleegs (optional): Cell array containing EEGLab data structures with
% the same ICA decomposition.
%
% Output:
% meanent: Vector which, for each IC, contains the mean of the entropy of
% the IC's time series. It is calculated by finding the entropy in
% one-second intervals of each epoch, and then taking the mean of those
% entropies.
%
% varent: Vector which, for each IC, contains the variance of the entropy
% of the IC's time series. It is calculated by finding the entropy in
% one-second intervals of each epoch, and then taking the variance of those
% entropies.

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
    arg('eeg', 'EEG'),...
    arg('celleegs', [], [],...
    'EEGLab data structures in cell array.'));

if isempty(celleegs)
    celleegs{1} = eeg;
end

ents_all = [];
for ieeg = 1:length(celleegs)
    eeg = celleegs{ieeg};
icaacts = eeg.icaact_filtered_resampled;
interval = 1; % the entropy is calculated over intervals of one second
interval_points = ceil(interval*eeg.virtual_srate);
interval_starts = 1:interval_points:size(icaacts,2);
n_intervals = length(interval_starts)-1;

num_epoch = size(icaacts,3);
n_ics = size(icaacts, 1);
ents = zeros(n_ics , n_intervals, num_epoch);
% split each epoch for each ic into intervals for which the variance is calculated.
for k=1:n_ics
    for j=1:num_epoch
        for m=1:n_intervals
            [y,xi]=ksdensity(icaacts(k, interval_starts(m):(interval_starts(m+1)-1), j),'width',0.1);
            d1=diff(xi); % binwidth
            y=y*d1(1)+eps;
            y=y/sum(y); % normalize so that the vector y sums to one, i.e. fulfills requirement of probability distribution
            ents(k, m, j)=-sum(y.*log(y));
        end
    end
end
ents_all = [ents_all reshape(ents, size(ents,1), size(ents,2)*size(ents,3))];
end

meanent = mean(ents_all, 2);
varent = var(ents_all, [], 2);
end
