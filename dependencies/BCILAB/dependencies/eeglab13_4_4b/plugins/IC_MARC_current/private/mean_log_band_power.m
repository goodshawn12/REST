function [theta, alpha, beta, gamma, gammamed, gamma_elec, gammah,...
    vartheta, varalpha, varbeta, vargamma, vargammamed, vargamma_elec,...
    vargammah]= mean_log_band_power(varargin)
% mean_log_band_power calculates the mean of the logarithm of the band
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
% sp_icas (optional): Cell array as constructed by calculate_spectra
% (default: []). If not supplied, it is calculated.
%
% celleegs (optional): Cell array containing EEGLab data structures with
% the same ICA decomposition.
%
% Output:
% theta alpha beta gamma gammamed gamma_elec gammah: Vectors of lengths
% equal to the number of ICs in the input eeg. The vectors contain the
% means over time of the logarithms of the power in the frequency bands
% mentioned above.
%
% vartheta varalpha varbeta vargamma vargammamed vargamma_elec vargammah:
% Vectors of lengths equal to the number of ICs in the input eeg. The
% vectors contain the variances over time of the logarithms of the power in
% the frequency bands mentioned above.

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
    arg('sp_icas', [], [],...
    'Cell array of spectra for ic activations.'),...
    arg('celleegs', [], [],...
    'EEGLab data structures in cell array.'));

if(isempty(sp_icas))
    sp_icas = calculate_spectra_wrapper(eeg, 'celleegs', celleegs);
end

n_ics = size(eeg.icawinv,2);
theta = NaN(n_ics,1);
alpha = NaN(n_ics,1);
beta = NaN(n_ics,1);
gamma = NaN(n_ics,1);
gammamed = NaN(n_ics,1);
gamma_elec = NaN(n_ics,1);
gammah = NaN(n_ics,1);
vartheta = NaN(n_ics,1);
varalpha = NaN(n_ics,1);
varbeta = NaN(n_ics,1);
vargamma = NaN(n_ics,1);
vargammamed = NaN(n_ics,1);
vargamma_elec = NaN(n_ics,1);
vargammah = NaN(n_ics,1);

for ic=1:n_ics
    sp_ica = sp_icas{ic};
    theta(ic) = mean(log(sp_ica.thetapow));
    alpha(ic) = mean(log(sp_ica.alphapow));
    beta(ic) = mean(log(sp_ica.betapow));
    gamma(ic) = mean(log(sp_ica.gammapow(1,:)));
    gammamed(ic) = mean(log(sp_ica.gammapow(2,:)));
    gamma_elec(ic) = mean(log(sp_ica.gammapow(3,:)));
    gammah(ic) = mean(log(sp_ica.gammapow(4,:)));
    
    vartheta(ic) = var(log(sp_ica.thetapow));
    varalpha(ic) = var(log(sp_ica.alphapow));
    varbeta(ic) = var(log(sp_ica.betapow));
    vargamma(ic) = var(log(sp_ica.gammapow(1,:)));
    vargammamed(ic) = var(log(sp_ica.gammapow(2,:)));
    vargamma_elec(ic) = var(log(sp_ica.gammapow(3,:)));
    vargammah(ic) = var(log(sp_ica.gammapow(4,:)));    
end


end


