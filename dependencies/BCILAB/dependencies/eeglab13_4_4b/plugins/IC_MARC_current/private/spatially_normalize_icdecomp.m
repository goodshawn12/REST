function eeg = spatially_normalize_icdecomp(eeg)
% Standardizes IC scalp maps to have mean zero and variance one.
%
% Input:
% eeg: EEGLab data structure.

% Output:
% eegout: The EEGlab data structure given as input, with standardized IC
% scalp maps (columns of eeg.icawinv have mean zero and variance one).

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

if isempty(eeg.icaact) && isempty(eeg.icawinv) && isempty(eeg.icaweights)
    error('no ica decomposition present in eeg struct, cannot normalize non-existent spatial map')
end

eeg.icawinv = zscore(eeg.icawinv);
