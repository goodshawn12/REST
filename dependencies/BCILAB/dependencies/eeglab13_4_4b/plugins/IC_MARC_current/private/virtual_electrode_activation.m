function activations = virtual_electrode_activation(phis, thetas, sigmas, chanlocs, scalpmap, sph_radius)
% virtual_electrode_activation estimates the activation of electrodes in
% arbirtrary locations given a scalp map for an IC and the locations of the
% electrodes in that scalp map.
%
% Input:
% phis: Longitudes of points for which a virtual activation is desired.
%
% thetas: Latitudes of points for which a virtual activation is desired.
%
% sigmas: Standard devations of Gaussian kernels used to weight the
% original electrode activations in the calculation of the virtual
% electrode activations.
%
% chanlocs: Locations of original electrodes. Must be a struct containing
% the fields sph_phi and sph_theta. These fields should contain the
% longitudes (sph_phi) and latitudes (sph_theta) of the original 
% electrodes.
%
% scalpmap: Scalp map (topography, i.e. a column of the field icawinv in
% the EEGLab data structure) for a single IC.
%
% sph_radius: Radius of the head.
%
% Output:
% activations: Topography of IC on the virtual electrode cap, i.e.
% activations of the specified virtual electrodes on this IC's scalp map.

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


all_electrodes_phi = cell2mat({chanlocs.sph_phi});
all_electrodes_theta = cell2mat({chanlocs.sph_theta});
nelectrodes = length(all_electrodes_theta);
nvirtual_electrodes = length(phis);

activations = NaN(nvirtual_electrodes, 1);
for ivirtual_electrode = 1:nvirtual_electrodes
    distances = spherical_distance(repmat(phis(ivirtual_electrode), 1, nelectrodes),...
        repmat(thetas(ivirtual_electrode), 1, nelectrodes),...
        all_electrodes_phi, all_electrodes_theta, sph_radius);
    distances = distances.^2;
    exp_distances = exp(-distances/(2*sigmas(ivirtual_electrode)^2));
    exp_distances=exp_distances/sum(exp_distances);
    activations(ivirtual_electrode) = sum(scalpmap.*exp_distances);
    
end
end

