function [front, post, leftarea, rightarea] = scalpmap_features_light(topog, chanlocs, virtual_chanlocs)
% scalpmap_features calculates the degree of activation in various areas of
% the scalp for each IC in the EEGLab data structure given as input. It
% also calculates a discontinuity measure of spatial activity and measures
% of the difference in activation between the two eye areas and the frontal
% and posterior areas of the scalp.
%
% Input:
% eeg: EEGLab data structure.
%
% Output:
% lateral_eyes: The absolute value of the difference between lefteye and
% righteye, for each IC.
%
% vertical_polarity: The absolute value of the difference between front and
% post, for each IC.
%
% lefteye: The sum of activations of electrodes weighted by a Gaussian
% kernel with center in Fp1 and standard deviation one. Calculated for each
% IC in the EEGLab data structure given as input.
%
% righteye: The sum of activations of electrodes weighted by a Gaussian
% kernel with center in Fp2 and standard deviation one. Calculated for each
% IC in the EEGLab data structure given as input.
%
% front: The sum of activations of electrodes weighted by a Gaussian
% kernel with center in AFz and standard deviation two. Calculated for each
% IC in the EEGLab data structure given as input.
%
% cent: The sum of activations of electrodes weighted by a Gaussian
% kernel with center in Cz and standard deviation two. Calculated for each
% IC in the EEGLab data structure given as input.
%
% post: The sum of activations of electrodes weighted by a Gaussian
% kernel with center in POz and standard deviation two. Calculated for each
% IC in the EEGLab data structure given as input.
%
% leftarea: The sum of activations of electrodes weighted by a Gaussian
% kernel with center in C5 and standard deviation two. Calculated for each
% IC in the EEGLab data structure given as input.
%
% rightarea: The sum of activations of electrodes weighted by a Gaussian
% kernel with center in C4 and standard deviation two. Calculated for each
% IC in the EEGLab data structure given as input.
%
% abs_med: Absolute value of the median of the scalp map for each IC.

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

% 356 long struct containing locations of electrodes that we use as centers of scalp areas for 
% which features are calculated

% for each of the areas left_eye, right_eye, frontal, central, posterior,
% left, and right, get the theta and phi coordinates of that area by
% finding the position in chanlocs that has the correct label. Since
% chanlocs is 356 long, make a struct containing the desired label repeated
% 356 times and find the match between the constructed struct and chanlocs. 

frontal = cellfun(@strcmp, {virtual_chanlocs.labels}, repmat({'afz'}, 1,64));
frontal_theta = virtual_chanlocs(frontal).sph_theta; frontal_phi = virtual_chanlocs(frontal).sph_phi;
posterior = cellfun(@strcmp, {virtual_chanlocs.labels}, repmat({'poz'}, 1,64));
posterior_theta = virtual_chanlocs(posterior).sph_theta; posterior_phi = virtual_chanlocs(posterior).sph_phi;
left = cellfun(@strcmp, {virtual_chanlocs.labels}, repmat({'c5'}, 1,64));
left_theta = virtual_chanlocs(left).sph_theta; left_phi = virtual_chanlocs(left).sph_phi;
right = cellfun(@strcmp, {virtual_chanlocs.labels}, repmat({'c4'}, 1,64));
right_theta = virtual_chanlocs(right).sph_theta; right_phi = virtual_chanlocs(right).sph_phi;

thetas = [frontal_theta, posterior_theta, left_theta, right_theta];
phis = [frontal_phi, posterior_phi, left_phi, right_phi];
sigmas = repmat(2, 1,4); % the standard deviation of the Gaussian kernel for each scalp area

n_ics = size(topog, 1);
front=zeros(n_ics,1);
post=zeros(n_ics,1);
leftarea=zeros(n_ics,1);
rightarea=zeros(n_ics,1);
for ic=1:n_ics
    activations = virtual_electrode_activation(phis, thetas, sigmas, chanlocs, topog(ic, :), virtual_chanlocs(1).sph_radius);
    front(ic)=abs(activations(1));
    post(ic)=abs(activations(2));
    leftarea(ic)=abs(activations(3));
    rightarea(ic)=abs(activations(4));
end

end
