function skew = eeg_skewness(eeg, varargin)
% eeg_skewness calculates the skewness in intervals of epochs for each IC.
%
% Input:
% eeg: EEGLab data structure.
%
% Options:
%   'interval': The length in seconds of intervals over which the skewness
%   is to be calculated (default: 15).
%
% Output:
% skew: 3D-matrix containing the calculated skewness. Each row contains the
% skewness calculated for one IC. The second dimension (columns in a 2d
% matrix) correspond to intervals. The third dimension contains the skewness
% calculated for each epoch. The matrix thus has dimensions
% (number of ICs)*(number of intervals per epoch)*(number of epochs)

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

opts = hlp_varargin2struct(varargin, 'interval', 15);

if(~isfield(opts, 'interval'))
    opts.interval = 15;
end

num_epoch = size(eeg.icaact,3);

interval_points = ceil(opts.interval*eeg.srate); % opts.interval is
% desired length of interval in seconds. The desired samples in each
% interval is opts.interval multiplied by the sampling rate.

interval_starts = 1:interval_points:size(eeg.icaact_filtered_resampled,2);

n_intervals = length(interval_starts)-1; % subtract one to get the number of whole intervals
n_ics = size(eeg.icaact, 1);
skew = zeros(n_ics, n_intervals, num_epoch);

% Calculate the skewness in the desired intervals in each epoch for each IC
for k=1:n_ics
    for j=1:num_epoch
        for i=1:n_intervals
            skew(k, i, j) = skewness(...
                eeg.icaact(k, interval_starts(i):(interval_starts(i+1)-1), j));
        end
    end
end
