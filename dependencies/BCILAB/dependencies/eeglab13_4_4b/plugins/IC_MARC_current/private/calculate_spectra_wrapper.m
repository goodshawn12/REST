function sp_icas = calculate_spectra_wrapper(eeg, varargin)
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


arg_define(0,varargin,...
    arg('celleegs', [], [],...
    'EEGLab data structures in cell array.'));

if isempty(celleegs)
    celleegs = {eeg};
end

nics = size(celleegs{1}.icawinv,2);
sp_icas = cell(nics, 1);

for ic=1:nics
    sp_icas{ic}.thetapow = [];
    sp_icas{ic}.alphapow = [];
    sp_icas{ic}.betapow = [];
    sp_icas{ic}.gammapow = [];
end
    
for ieeg = 1:length(celleegs)
    eeg = celleegs{ieeg};
    eeg.icaact = temporally_normalize_icdecomp(eeg.icaact_filtered_resampled);
    sp_icascur = calculate_spectra(eeg);
    for ic=1:nics
       sp_icacur = sp_icascur{ic}; 
       sp_icas{ic}.thetapow = [sp_icas{ic}.thetapow sp_icacur.thetapow];
       sp_icas{ic}.alphapow = [sp_icas{ic}.alphapow sp_icacur.alphapow];
       sp_icas{ic}.betapow = [sp_icas{ic}.betapow sp_icacur.betapow];
       sp_icas{ic}.gammapow = [sp_icas{ic}.gammapow sp_icacur.gammapow];
    end
    sp_icas{ic}.f = sp_icacur.f;
end