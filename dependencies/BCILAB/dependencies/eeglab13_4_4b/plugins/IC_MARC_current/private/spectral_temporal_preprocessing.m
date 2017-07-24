function eeg= spectral_temporal_preprocessing(eeg)
% spectral_temporal_preprocessing calculates the mean of the logarithm of the band
% power in the following frequency bands, for each IC: theta (4-7Hz), alpha
% (8-13Hz), beta (14-20Hz), low gamma (21-30Hz), medium gamma (30-45Hz),
% power grid voltage (46-65Hz), and high gamma (66-80Hz)
%
% Input:
% eeg: EEGLab data structure. If the argument sp_icas is not supplied, the 
% data structure must have the additional fields
% virtual_srate and icaact_filtered_resampled. See
% the readme file for the toolbox for details on these fields.
%
% Output:
% eeg: The EEGLab data structure with additional fields 
% icaact_filtered_resampled, virtual_pnts, and icaact_unfiltered.  See
% the readme file for the toolbox for details on these fields.%
%
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

if isempty(eeg.icaact)
    eeg.icaact = eeg_getica(eeg);
end
icaact_filtered_resampled = eeg.icaact;
n_ics = size(icaact_filtered_resampled,1);
icaact_filtered_resampled = double(icaact_filtered_resampled);

P = 200;
Q = eeg.srate;
nepochs = size(icaact_filtered_resampled,3);
tslength = size(icaact_filtered_resampled(1,:, 1),2);
newlength = ceil(P/Q*tslength);
newicaact = NaN(n_ics, newlength, nepochs);
for ic = 1:n_ics
    for epoch=1:nepochs
        newicaact(ic,:, epoch) = resample(icaact_filtered_resampled(ic,:, epoch), P, Q);
    end
end
eeg.icaact_unfiltered = newicaact;
eeg.virtual_pnts = size(newicaact, 2);

rp=0.025; rs=30; locutoff=3; hicutoff=90;
eeg.icaact_filtered_resampled = eeg.icaact_unfiltered;
eeg = bandpass_icaactivation(eeg, rp, rs, locutoff, hicutoff);