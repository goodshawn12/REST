function distance = spherical_distance(sph_phi1, sph_theta1, sph_phi2, sph_theta2, sph_radius)
% spherical_distance calculates the distance between two points along a
% sphere.
% Input:
% sph_phi1, sph_phi2: The longitudes of two points. Can also be vectors
% giving longitudes of several points between which the spherical distance
% is to be calculated.
%
% sph_theta1, sph_theta2: The azimuths (also referred to as latitudes) of 
% the two points. Can also be vectors
% giving latitudes of several points between which the spherical distance
% is to be calculated.
%
% sph_radius: The radius of the sphere on which the two points lie.
%
% Output:
% distance: The distance(s) between the two (or several) points along the
% sphere.

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

lambdaf = sph_phi1/180*pi; % longitude
lambdas = sph_phi2/180*pi; % longitude
phif = sph_theta1/180*pi; %latitude
phis = sph_theta2/180*pi; %latitude
delta_lambda = lambdaf-lambdas;
cosphif = cos(phif);
cosphis = cos(phis);
sinphif = sin(phif);
sinphis = sin(phis);

% http://en.wikipedia.org/wiki/Great-circle_distance
central_angle = atan2(sqrt((cosphif.*sin(delta_lambda)).^2+(cosphis.*sinphif -sinphis.*cosphif.*cos(delta_lambda)).^2),...
    (sinphis.*sinphif+cosphis.*cosphif.*cos(delta_lambda)));

distance = (sph_radius*central_angle);
