function [meanentropy, varentropy]= spectral_entropy(varargin)
% spectral_entropy calculates, for each IC, the entropy of the mean power 
% over the following bands: theta (4-7Hz), alpha
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
% meanentropy: Vector of length equal to the number of ICs in the input 
% eeg. The vector contains the means over time of the entropy of the power 
% over the frequency bands mentioned above. It is calculated by finding the
% entropy of the power distribution for each time window (column of
% sp_ica.pxx), and then taking the mean over these entropies. This is done
% for each IC.
%
% varentropy: Vector of length equal to the number of ICs in the input 
% eeg. The vector contains the variances over time of the entropy of the 
% power over the frequency bands mentioned above. It is calculated by 
% finding the entropy of the power distribution for each time window 
% (column of sp_ica.pxx), and then taking the variance over these 
% entropies. This is done for each IC.

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
meanentropy = NaN(n_ics,1);
varentropy = NaN(n_ics,1);
for ic=1:n_ics
    sp_ica = sp_icas{ic};
    allfreqs = [sp_ica.thetapow; sp_ica.alphapow; sp_ica.betapow; sp_ica.gammapow];
    nwindows = size(allfreqs,2);
    ents = NaN(1, nwindows);
    for window=1:nwindows
        y=allfreqs(:,window)/sum(allfreqs(:,window));
        ents(window)=-sum(y.*log(y));
    end
    meanentropy(ic) = mean(ents);
    varentropy(ic) = var(ents);
end
