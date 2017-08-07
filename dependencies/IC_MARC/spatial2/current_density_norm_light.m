function cdn = current_density_norm_light(virtual_topography, virtual_chanlocs)
% current_density_norm calculates the current density norm as described in 
% Automatic Classification of Artifactual ICA-Components for Artifact 
% Removal in EEG Signals by Irene Winkler, Stefan Haufe and Michael 
% Tangermann http://www.behavioralandbrainfunctions.com/content/7/1/30
%
% Input:
% eeg: EEGLab data structure with additional fields
% virtual_topography and virtual_chanlocs. See the readme file for the
% toolbox for details on these fields.
%
% Output:
% cdn: Vector of current density norms for the ICs in the EEGLab data
% structure.

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

icawinv = virtual_topography';
icaweights = pinv(icawinv); % Why pinv? Make sure icawinv exists?

A = icawinv;
load('dipolfit_matrix'); % loads M100, clab
channel_labels={virtual_chanlocs.labels};
[dummy, idx_M, idx_IC] = intersect(lower(clab), channel_labels);

cdn = NaN(size(icawinv,2),1);
for ic = 1:size(icawinv,2)
    cdn(ic) = log(norm(M100(:,idx_M)*A(idx_IC,ic)));
end
end
