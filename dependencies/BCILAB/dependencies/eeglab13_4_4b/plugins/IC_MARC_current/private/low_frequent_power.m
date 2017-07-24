function [delta_relative_pows, vardelta_relative_pows]= low_frequent_power(varargin)
% low_frequent_power calculates the mean and variance of the proportion
% power contributed by frequencies in the delta band (1Hz - 3Hz).
%
% Input:
% eeg: EEGLab data structure with additional fields
% virtual_srate and icaact_unfiltered. See
% the readme file for the toolbox for details on these fields.
%
% celleegs (optional): Cell array containing EEGLab data structures with
% the same ICA decomposition.
%
% Output:
% delta_relative_pows: Vector containing, for each IC, the mean proportion
% of power contributed by the delta band (1Hz-3Hz). The mean is
% taken over one-second intervals.
%
% vardelta_relative_pows: Vector containing, for each IC, the variance of the
% proportion of power contributed by the delta band (1Hz-3Hz).
% The variance is taken over one-second intervals.

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
    arg('eeg',[],[],'EEGLab data structure'),...
    arg('celleegs', [], [],...
    'EEGLab data structures in cell array.'));

if isempty(celleegs)
    celleegs = {eeg};
end

windowlength_seconds = 1; % window length in seconds. Must be a reciprocal of a whole number.
windowlength = ceil(windowlength_seconds*eeg.virtual_srate); % window lengths in samples

pxxalleegs = [];
for ieeg = 1:length(celleegs)
    eeg = celleegs{ieeg};
    icaacts = eeg.icaact_unfiltered;
    nepochs = size(icaacts,3);
    tslength = size(icaacts,2);
    n_one_second_intervals = fix((tslength-eeg.virtual_srate/2)/(eeg.virtual_srate-eeg.virtual_srate/2)) ;
    n_ics = size(icaacts,1);
    pxx = NaN(n_ics, n_one_second_intervals*nepochs);
    icaacts = temporally_normalize_icdecomp(icaacts);
    
    for ic=1:n_ics
        for epoch=1:nepochs
            pxxcols=((epoch-1)*n_one_second_intervals+1):epoch*n_one_second_intervals;
            [spectoutput, f] = spectrogram(icaacts(ic, :, epoch), windowlength,windowlength/2,[],eeg.virtual_srate);
            spectoutput = abs(spectoutput).^2;
            
            deltainds = find(f>=1, 1, 'first'):find(f<=3, 1, 'last'); % find the rows that correspond to frequencies in the delta band (1-3Hz)
            % in the matrix spectoutput.
            delta_relative_pow = mean(spectoutput(deltainds,:))./mean(spectoutput); % In each interval, get the proportion of mean power in frequencies
            % in the delta band (mean(spectoutput(deltainds,:))) to the mean power over all frequencies (mean(spectoutput)).
            pxx(ic,pxxcols) = mean(reshape(...
                delta_relative_pow(1:n_one_second_intervals/windowlength_seconds), 1/windowlength_seconds,n_one_second_intervals),1);
            
        end
    end
    pxxalleegs = [pxxalleegs pxx];
end

delta_relative_pows = mean(pxxalleegs,2);
vardelta_relative_pows = var(pxxalleegs,[],2);
end


