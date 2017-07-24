function [meanmaxdiff, varmaxdiff]= max_first_derivative( varargin )
% max_amplitude calculates the mean and variance of the maximal amplitude
% in epochs of the IC time series.
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
% meanmaxdiff: Vector of length equal to the number of ICs in the input
% eeg. The vector contains the means over epochs of the maximal difference
% between two samples in an epoch.
%
% varmaxdiff: Vector of length equal to the number of ICs in the input
% eeg. The vector contains the variances over epochs of the maximal
% difference between two samples in an epoch.

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
    arg({'eeg', 'EEG'},[],[],'EEGLab data structure.'),...
    arg('celleegs', [], [],...
    'EEGLab data structures in cell array.'));

if isempty(celleegs)
    celleegs{1} = eeg;
end

maxdiffs_all = [];
for ieeg = 1:length(celleegs)
    eeg = celleegs{ieeg};
    icaacts = eeg.icaact_filtered_resampled;
    icaacts = temporally_normalize_icdecomp(icaacts);
    diffs = diff(icaacts,1,2);
    
    n_ics = size(diffs, 1);
    dim = size(diffs);
    if(length(dim)<3) % make sure that data is split into epochs such that the
        % max_epoch_var feature makes sense and is not just one
        
        num_epoch = max(floor(dim(2)/eeg.virtual_srate), 5); % make sure we get some epochs, if possible of one seconds duration each
        epochlength = floor(dim(2)/num_epoch);
        diffs = reshape(diffs(:, 1:num_epoch*epochlength), n_ics, epochlength,num_epoch);
    end
    dim = size(diffs);
    num_epoch = dim(3);
    
        maxdiffs = NaN(n_ics, num_epoch);
    for ic=1:n_ics
        for iepoch = 1:num_epoch
            maxdiffs(ic,iepoch) = max(abs(diffs(ic, :, iepoch)));
        end
    end
    maxdiffs_all = [maxdiffs_all maxdiffs];
end
        meanmaxdiff=mean(maxdiffs_all, 2);
        varmaxdiff=var(maxdiffs_all, [], 2);

end

