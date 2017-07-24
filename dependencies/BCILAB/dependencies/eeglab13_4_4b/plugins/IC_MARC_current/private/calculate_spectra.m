function sp_icas = calculate_spectra(eeg)
% calculate_spectra calculates the power in the theta (4-7Hz), alpha
% (8-13Hz), beta (14-20Hz), low gamma (21-30Hz), medium gamma (30-45Hz),
% power grid voltage (46-65Hz), and high gamma (66-80Hz) bands for each IC
% in the EEGLab data structure given as input.
%
% Input:
% eeg: EEGLab data structure with additional fields
% icaact_filtered_resampled and virtual_srate. See the readme file for the
% toolbox for details on these fields.
%
% Output:
% sp_icas: Cell array containing a structure with the following two fields
% for each IC:
%   pxx: Matrix containing the power in each frequency band mentioned above
%   over time. The entry in the i'th row in the j'th column gives the power
%   calculated by spectrogram for the i'th frequency band for the j'th
%   one-second interval.
%
%   f: Vector containing the midpoints of the frequency ranges for which
%   power is calculated.

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

theta = 4:7; alpha=8:13; beta=14:20; gamma = 21:30; gammamed = 30:45; gamma_elec = 46:65; gammah = 66:80;

icaacts = eeg.icaact;
dim = size(icaacts);
if length(dim)==3
    nepochs=dim(3);
else
    nepochs = 1;
end

n_ics = size(eeg.icawinv,2);
sp_icas = cell(n_ics,1);

% To estimate low frequencies, longer time intervals are necessary. To
% capture characteristics of high-frequent activity, smaller time intervals
% are necessary. For this reason, we use different window lengths for
% calculating the spectrogram for different frequency ranges.
windowlength_seconds = 1; % window length in seconds
windowlength = ceil(windowlength_seconds*eeg.virtual_srate); % window lengths in samples

npnts = fix((size(icaacts,2)-windowlength/2)/(windowlength-windowlength/2));

for ic = 1:n_ics
    sp_ica.thetapow = NaN(1, npnts*nepochs);
    sp_ica.alphapow = NaN(1, npnts*nepochs);
    sp_ica.betapow = NaN(1, npnts*nepochs);
    sp_ica.gammapow = NaN(4, npnts*nepochs);
    for iepoch=1:nepochs        
        % calculate the spectrogram for one-second long windows
        [spectoutput f] = spectrogram(icaacts(ic, :, iepoch), windowlength,windowlength/2,[],eeg.virtual_srate);
        spectoutput = abs(spectoutput).^2;
        thetainds = find(f>=theta(1), 1, 'first'):find(f<=theta(end), 1, 'last');
        thetapow = mean(spectoutput(thetainds,:));
        sp_ica.thetapow(((iepoch-1)*npnts+1):(iepoch*npnts)) = thetapow;
        
        alphainds = find(f>=alpha(1), 1, 'first'):find(f<=alpha(end), 1, 'last');
        alphapow = mean(spectoutput(alphainds,:));
        sp_ica.alphapow(((iepoch-1)*npnts+1):(iepoch*npnts)) = alphapow;
        
        betainds = find(f>=beta(1), 1, 'first'):find(f<=beta(end), 1, 'last');
        betapow = mean(spectoutput(betainds,:));
        sp_ica.betapow(((iepoch-1)*npnts+1):(iepoch*npnts)) = betapow;
        
        gammainds = find(f>=gamma(1), 1, 'first'):find(f<=gamma(end), 1, 'last');
        gammamedinds = find(f>=gammamed(1), 1, 'first'):find(f<=gammamed(end), 1, 'last');
        gamma_elecinds = find(f>=gamma_elec(1), 1, 'first'):find(f<=gamma_elec(end), 1, 'last');
        gammahinds = find(f>=gammah(1), 1, 'first'):find(f<=gammah(end), 1, 'last');
        gammapow = mean(spectoutput(gammainds,:));
        gammamedpow = mean(spectoutput(gammamedinds,:));
        gamma_elecpow = mean(spectoutput(gamma_elecinds,:));
        gammahpow = mean(spectoutput(gammahinds,:));
        sp_ica.gammapow(:,((iepoch-1)*npnts+1):(iepoch*npnts)) = [gammapow; gammamedpow; gamma_elecpow; gammahpow];
    end
    sp_ica.f = [mean(theta) mean(alpha) mean(beta) mean(gamma) mean(gammamed) mean(gamma_elec) mean(gammah)];
    sp_icas{ic} = sp_ica;
end

