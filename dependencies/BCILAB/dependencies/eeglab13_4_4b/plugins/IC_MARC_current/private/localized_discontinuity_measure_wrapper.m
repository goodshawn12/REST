function general_discontinuity = localized_discontinuity_measure_wrapper(eeg, varargin)
% localized_discontinuity_measure_wrapper is a wrapper for the function
% localized_discontinuity_measure. localized_discontinuity_measure_wrapper
% takes as input an EEGLab data structure and calls 
% localized_discontinuity_measure with the inputs needed by that function.
%
% Input:
% eeg: EEGLab data structure with additional fields virtual_chanlocs and
% virtual_topography. See the readme file for the toolbox for details on
% these fields.
%
% Output:
% general_discontinuity: Column vector containing, for each IC, a measure 
% of discontinuity based on the scalp map.

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

general_discontinuity = localized_discontinuity_measure(eeg.virtual_topography,eeg.virtual_chanlocs,size(eeg.icawinv, 2));
