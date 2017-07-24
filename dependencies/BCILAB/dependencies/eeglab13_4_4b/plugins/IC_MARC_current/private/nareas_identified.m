function ideds = nareas_identified(eeg)
% ndipoles_identified counts the number of labels of anatomical areas
% associated with the dipole fit of each IC in the EEGLab data structure
% given as input.
%
% Input:
% eeg: EEGLab data structure with a dipole fit and labels for the dipoles
% added by the function add_dipfit_labels.
%
% Output:
% ideds: Vector of length equal to the number of ICs in the input
% eeg. The vector contains the number of anatomical areas labeled for each
% IC via the dipole fit.

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

ncomps = size(eeg.icawinv, 2);
ideds = zeros(ncomps, 1);

for i = 1:ncomps
    labs = regexp(eeg.dipfit.label{i}, ',', 'split');%strsplit(eeg.dipfit.label{i}, ',');
    ided = 0; % id'ed, number of identified IC locations
    for j = 1:length(labs)
        if(~strcmp('*', labs(j)))
            ided = ided+1;
        end
    end
    ideds(i) = ided;
end
