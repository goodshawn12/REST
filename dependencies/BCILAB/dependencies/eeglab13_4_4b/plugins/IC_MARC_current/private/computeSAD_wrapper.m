function [SAD, var_front, var_back, mean_front, mean_back] = computeSAD_wrapper(eeg, varargin)
% Wrapper for computeSAD_variable_ics. computeSAD_wrapper takes an EEGLab
% data structure (must have additional fields virtual_topography and 
% virtual_chanlocs) as input and calls computeSAD_variable_ics with the
% appropriate variables from the EEGLab data structure.
%
% Input:
% eeg: EEGLab data structure with additional fields
% virtual_topography and virtual_chanlocs. See the readme file for the
% toolbox for details on these fields.
%
% Output:
% SAD: Spatial Average Difference, the difference between activation in
% frontal and posterior areas of the scalp. Calculated as abs(mean(frontal
% area)) - abs(mean(posterior area)).
%
% var_front: The variance in activation in the frontal part of the scalp.
%
% var_back: The variance in activation in the posterior part of the scalp.
%
% mean_front: The mean activation in the frontal part of the scalp.
%
% mean_back: The mean activation in the posterior part of the scalp.

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


[SAD,var_front,var_back,mean_front,mean_back] = ...
    computeSAD_variable_ics(eeg.virtual_topography,eeg.virtual_chanlocs,size(eeg.virtual_topography, 2), size(eeg.virtual_topography, 1));

mean_front = abs(mean_front);
mean_back = abs(mean_back);

