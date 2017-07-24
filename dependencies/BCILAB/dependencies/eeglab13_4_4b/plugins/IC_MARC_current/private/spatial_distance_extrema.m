function spatial_dist_extrema = spatial_distance_extrema(varargin)
% spatial_dist_extrema calculates the log of the 2-norm of the distance
% between the two electrodes with the highest and lowest activations
% as described in 
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
% spatial_dist_extrema: A vector containing the log of the 2-norm of the
% distance between the electrodes with highest and lowest activation in the
% scalp map of each IC in the EEGLab data structure given as input.

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
    arg('eeg', 'EEG'));

topography = eeg.virtual_topography;
eeg.chanlocs = eeg.virtual_chanlocs;

nics = size(eeg.icawinv,2);
spatial_dist_extrema = NaN(nics,1);
for ic =1:nics
    channel_activation = topography(ic,:);
    [dummy maxindex] = max(channel_activation);
    [dummy minindex] = min(channel_activation);
    maxlocation = [eeg.chanlocs(maxindex).X eeg.chanlocs(maxindex).Y eeg.chanlocs(maxindex).Z];
    minlocation = [eeg.chanlocs(minindex).X eeg.chanlocs(minindex).Y eeg.chanlocs(minindex).Z];
    
    spatial_dist_extrema(ic) = log(norm(maxlocation-minlocation,2));
end

end

