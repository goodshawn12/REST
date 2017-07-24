function [lateral_eyes, vertical_polarity,...
    lefteye, righteye, front, cent, post, leftarea, rightarea, abs_med] = scalpmap_features(varargin)
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

arg_define(1,varargin, ...
    arg('eeg', 'EEG'));

eeg = spatially_normalize_icdecomp(eeg);
topog = eeg.icawinv';

filename = mfilename('fullpath');
folder_path = filename(1:(strfind(filename, 'scalpmap_features')-1));

% 356 long struct containing locations of electrodes that we use as centers of scalp areas for 
% which features are calculated
chanlocs = readlocs([ folder_path 'standard-10-5-cap385.elp']);

% for each of the areas left_eye, right_eye, frontal, central, posterior,
% left, and right, get the theta and phi coordinates of that area by
% finding the position in chanlocs that has the correct label. Since
% chanlocs is 356 long, make a struct containing the desired label repeated
% 356 times and find the match between the constructed struct and chanlocs. 
left_eye = cellfun(@strcmp, {chanlocs.labels}, repmat({'Fp1'}, 1,356));
left_eye_theta = chanlocs(left_eye).sph_theta; left_eye_phi = chanlocs(left_eye).sph_phi;
right_eye = cellfun(@strcmp, {chanlocs.labels}, repmat({'Fp2'}, 1,356));
right_eye_theta = chanlocs(right_eye).sph_theta; right_eye_phi = chanlocs(right_eye).sph_phi;
frontal = cellfun(@strcmp, {chanlocs.labels}, repmat({'AFz'}, 1,356));
frontal_theta = chanlocs(frontal).sph_theta; frontal_phi = chanlocs(frontal).sph_phi;
central = cellfun(@strcmp, {chanlocs.labels}, repmat({'Cz'}, 1,356));
central_theta = chanlocs(central).sph_theta; central_phi = chanlocs(central).sph_phi;
posterior = cellfun(@strcmp, {chanlocs.labels}, repmat({'POz'}, 1,356));
posterior_theta = chanlocs(posterior).sph_theta; posterior_phi = chanlocs(posterior).sph_phi;
left = cellfun(@strcmp, {chanlocs.labels}, repmat({'C5'}, 1,356));
left_theta = chanlocs(left).sph_theta; left_phi = chanlocs(left).sph_phi;
right = cellfun(@strcmp, {chanlocs.labels}, repmat({'C4'}, 1,356));
right_theta = chanlocs(right).sph_theta; right_phi = chanlocs(right).sph_phi;

thetas = [left_eye_theta, right_eye_theta, frontal_theta, central_theta, posterior_theta, left_theta, right_theta];
phis = [left_eye_phi, right_eye_phi, frontal_phi, central_phi, posterior_phi, left_phi, right_phi];
sigmas = [1 1 repmat(2, 1,5)]; % the standard deviation of the Gaussian kernel for each scalp area

n_ics = size(eeg.icawinv, 2);
lateral_eyes=zeros(n_ics,1);
vertical_polarity=zeros(n_ics,1);
lefteye=zeros(n_ics,1);
righteye=zeros(n_ics,1);
front=zeros(n_ics,1);
cent=zeros(n_ics,1);
post=zeros(n_ics,1);
leftarea=zeros(n_ics,1);
rightarea=zeros(n_ics,1);
for ic=1:n_ics
    activations = virtual_electrode_activation(phis, thetas, sigmas, eeg.chanlocs, topog(ic, :), eeg.chanlocs(1).sph_radius);
    lefteye(ic)=abs(activations(1));
    righteye(ic)=abs(activations(2));
    front(ic)=abs(activations(3));
    cent(ic)=abs(activations(4));
    post(ic)=abs(activations(5));
    leftarea(ic)=abs(activations(6));
    rightarea(ic)=abs(activations(7));
    lateral_eyes(ic) =abs(activations(1)-activations(2));
    vertical_polarity(ic) = abs(activations(3)-activations(5));
end

abs_med = abs(median(topog,2));

end
