function eegout = check_ica_presence(eeg, varargin)
% checks that an ica decomposition is present. If not, it runs one. Also, inverse weights are calculated
% if not already present, as is the activation of the ic's over time
%
% Input:
% eeg: EEGLab data structure.
%
% force_preprocessing: Force preprocessing steps to be performed even if
% unnecessary. If no temporal or spectral features are requested, it is not
% necessary to downsample and filter time series. However, if the outputs
% produced during preprocessing are desired, force_preprocessing can be set
% to true. (default: false)
%
% Output:
% eegout: The EEGlab data structure given as input, with IC fields filled
% in.
%
% ica_options: Options to pass to the IC method if no IC decomposition is
% attached to the eeg data structure given as input. (default: [])
%
% ica_method: IC decomposition method. If no IC decomposition is
% attached to the eeg data structure given as input, one will be
% calculated. (default: 'fast_ica')

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


arg_define(0,varargin,...
    arg('do_icaact', true, [],...
    'Whether or not to calculate IC activations.'));

eegout = eeg;

if isempty(eeg.icaact) || isempty(eeg.icawinv) || isempty(eeg.icaweights)

if(isempty(eegout.icawinv) && isempty(eegout.icaact)  && isempty(eegout.icaweights))
    error('check_ica_presence.m: No ica decomposition associated with dataset, run ICA before calling IC_MARC or ic_feature_extraction.m')
end

if(isempty(eegout.icawinv) && ~isempty(eegout.icaweights))
    if(~isempty(eegout.icasphere))
        eegout.icaweights = eegout.icaweights*eegout.icasphere;
        eegout.icasphere = eye(size(eegout.weights,2));
    end
        eegout.icawinv = pinv(eegout.icaweights);
end

if(isempty(eegout.icaweights) && ~isempty(eegout.icawinv))
        eegout.icaweights = pinv(eegout.icawinv);
end

% the case of icaact being non-empty and weights being empty should not arise, and has not been
% accounted for here
if(isempty(eegout.icaact) && ~(isempty(eegout.icaweights) || ~isempty(eegout.icawinv)))
    warning('check_ica_presence.m: EEG dataset does not contain IC weight or inverse weight matrices (but does have IC activations)')
end

if isempty(eegout.icaact) 
    if ~isempty(eegout.data)
    eegout.icaact = eeg_getica(eegout);
    else
       warning('check_ica_presence.m: EEG dataset does not contain data. IC activations not recalculated.') 
    end
end

end
end
