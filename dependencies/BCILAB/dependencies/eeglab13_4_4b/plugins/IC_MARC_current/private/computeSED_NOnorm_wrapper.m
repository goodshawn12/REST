function [SED medie_left medie_right] = computeSED_NOnorm_wrapper(eeg, varargin)
% Wrapper for computeSED_NOnorm_variable_ics. computeSED_NOnorm_wrapper 
% takes an EEGLab data structure (must have additional fields 
% virtual_topography and virtual_chanlocs) as input and calls 
% computeSED_NOnorm_variable_ics with the appropriate variables from the 
% EEGLab data structure.
%
% Input:
% eeg: EEGLab data structure with additional fields
% virtual_topography and virtual_chanlocs. See the readme file for the
% toolbox for details on these fields.
%
% Output:
% SED: Spatial Eye Difference. The difference between activation in the two
% eye areas. Calculated as abs(mean(left)-mean(right)).
%
% medie_left: Mean activation in the left eye area.
%
% medie_right: Mean activation in the right eye area.

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


[SED,medie_left,medie_right] = ...
    computeSED_NOnorm_variable_ics(eeg.virtual_topography,eeg.virtual_chanlocs,size(eeg.virtual_topography, 2), size(eeg.virtual_topography, 1));

medie_left = abs(medie_left);
medie_right = abs(medie_right);
end
