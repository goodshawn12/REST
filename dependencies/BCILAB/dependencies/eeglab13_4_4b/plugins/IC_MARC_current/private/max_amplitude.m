function [meanmaxampl, varmaxampl] = max_amplitude( varargin )
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
% meanmaxampl: Vector of length equal to the number of ICs in the input
% eeg. The vector contains the means of the maximal amplitudes in epochs.
%
% varmaxampl: Vector of length equal to the number of ICs in the input
% eeg. The vector contains the variances of the maximal amplitudes in
% epochs.

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

maxamps_all = [];
for ieeg = 1:length(celleegs)
    eeg = celleegs{ieeg};
    icaacts = temporally_normalize_icdecomp(eeg.icaact_filtered_resampled);
    n_ics = size(icaacts, 1);
    dim = size(icaacts);
    if(length(dim)<3) % make sure that data is split into epochs such that the
        % max_epoch_var feature makes sense and is not just one
        
        num_epoch = max(floor(dim(2)/eeg.virtual_srate), 5); % make sure we get some epochs, if possible of one seconds duration each
        epochlength = floor(dim(2)/num_epoch);
        icaacts = reshape(icaacts(:, 1:num_epoch*epochlength), n_ics, epochlength,num_epoch);
    end
    dim = size(icaacts);
    num_epoch = dim(3);
    
    maxamps = NaN(n_ics, num_epoch);
    for ic=1:n_ics
        for iepoch = 1:num_epoch
            maxamps(ic,iepoch) = max(abs(icaacts(ic, :, iepoch)));
        end
    end
    maxamps_all = [maxamps_all maxamps];
end

meanmaxampl=mean(maxamps_all, 2);
varmaxampl=var(maxamps_all, [], 2);
end
