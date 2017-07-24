function log_range = log_range_spatial(eeg, varargin)
% log_range_spatial calculates the logarithm of the range of values in the
% topographies associated with ICs.
%
% Input:
% eeg: EEGLab data structure with additional field virtual_topography. See 
% the readme file for the toolbox for details on this field.
%
% Output:
% log_range: Vector of length equal to the number of ICs in the input 
% eeg. The vector contains the logarithms of the ranges of values in the
% virtual topography for each IC.

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

log_range = log(range(eeg.virtual_topography, 2));
