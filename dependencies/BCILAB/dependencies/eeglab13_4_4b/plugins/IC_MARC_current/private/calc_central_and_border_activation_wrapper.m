function [f_cent,f_border] = calc_central_and_border_activation_wrapper(varargin)
% calc_central_and_border_activation_wrapper calls the function
% calc_central_and_border_activation for each IC in the EEGLab data
% structure given as input.
%
% Input:
% eeg: EEGLab data structure with additional fields virtual_chanlocs and
% virtual_topography. See the readme file for the toolbox for details on
% these fields.
%
% spatially_normalized_eeg (optional): EEGLab data structure in which the 
% inverse ICA weights are standardized (the columns of icawinv have mean 
% zero and variance one). If not supplied, it is calculated. 
%
% Output:
% f_cent: Column vector containing central activation on the scalp for each
% IC.
%
% f_border: Column vector containing border activation on the scalp for
% each IC.

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
    arg({'eeg', 'EEG'},[],[],'EEGLab struct containing EEG data.'), ...
    arg('spatially_normalized_eeg', [], [],...
    'EEG struct with ICA activation normalized.')...
    );

if isempty(spatially_normalized_eeg)
   eeg = spatially_normalize_icdecomp(eeg); 
else
    eeg = spatially_normalized_eeg;
end

topography = eeg.virtual_topography;
eeg.chanlocs = eeg.virtual_chanlocs;

f_cent = NaN(size(eeg.icawinv,2),1); 
f_border = NaN(size(eeg.icawinv,2),1);
for ic = 1:size(eeg.icawinv,2)
    [f_cent(ic),f_border(ic)] = calc_central_and_border_activation(topography(ic,:) ,eeg.chanlocs);
end

