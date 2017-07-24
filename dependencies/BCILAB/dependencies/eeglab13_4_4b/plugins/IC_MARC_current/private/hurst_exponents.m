function [hurst1, hurst2, hurst3, varhurst1, varhurst2, varhurst3] = hurst_exponents(varargin)
% hurst_exponents calculates three estimates of the Hurst exponent for each
% one-second interval of each epoch for each IC. The means and variances of
% these Hurst exponent estimates are returned.
%
% Input:
% eeg: EEGLab data structure with additional fields
% icaact_filtered_resampled and virtual_srate. See the readme file for the
% toolbox for details on these fields.
%
% celleegs (optional): Cell array containing EEGLab data structures with
% the same ICA decomposition.
%
% Output:
% hurst1, hurst2, hurst3: The function wfbmesti in the Matlab wavelet
% toolbox returns three estimates. The outputs hurst1, hurst2, and hurst3
% contain, for each IC, the averages of the first, second, and third 
% outputs from wfbmesti in one-second intervals, respectively. Each of
% hurst1, hurst2, and hurst3 is a vector with length equal to the number of
% ICs.
%
% varhurst1, varhurst2, varhurst3: Vectors of length equal to the number of
% ICs. varhurst1, varhurst2, and varhurst3 contain the variances of Hurst 
% exponent estimates in one-second intervals.

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

interval = 1; % the hurst exponents are calculated over intervals of one second
hursts1_all = [];
hursts2_all = [];
hursts3_all = [];
for ieeg = 1:length(celleegs)
    eeg = celleegs{ieeg};
icaacts = eeg.icaact_filtered_resampled;
icaacts = temporally_normalize_icdecomp(icaacts);


interval_points = ceil(interval*eeg.virtual_srate);
interval_starts = 1:interval_points:size(eeg.icaact_filtered_resampled,2); 
n_intervals = length(interval_starts)-1; % subtract one to get the number of whole intervals

num_epoch = size(icaacts, 3);

n_ics = size(icaacts, 1);
hursts1 = zeros(n_ics, n_intervals, num_epoch);
hursts2 = zeros(n_ics, n_intervals, num_epoch);
hursts3 = zeros(n_ics, n_intervals, num_epoch);
% split each epoch in each ic time series into intervals over which the variance and mean are calculated.
for k=1:n_ics
    for j=1:num_epoch
        for m=1:n_intervals
            hurstests = wfbmesti(icaacts(k, interval_starts(m):(interval_starts(m+1)-1), j));
            
            hursts1(k, m, j) = hurstests(1);
            hursts2(k, m, j) = hurstests(2);
            hursts3(k, m, j) = hurstests(3);
        end
    end
end
hursts1_all = [hursts1_all reshape(hursts1, size(hursts1,1), size(hursts1,2)*size(hursts1,3))];
hursts2_all = [hursts2_all reshape(hursts2, size(hursts2,1), size(hursts2,2)*size(hursts2,3))];
hursts3_all = [hursts3_all reshape(hursts3, size(hursts3,1), size(hursts3,2)*size(hursts3,3))];
end

hurst1 = mean(hursts1_all, 2);
hurst2 = mean(hursts2_all, 2);
hurst3 = mean(hursts3_all, 2);

varhurst1 = var(hursts1_all, [], 2);
varhurst2 = var(hursts2_all, [], 2);
varhurst3 = var(hursts3_all, [], 2);

end
