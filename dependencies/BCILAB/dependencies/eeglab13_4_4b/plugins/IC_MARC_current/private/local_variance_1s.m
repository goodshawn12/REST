function [mean_vars1s, var_vars1s] = local_variance_1s(varargin)
% local_variance_1s calls local_variance with the argument interval set to
% 1.
%
% Input:
% eeg: EEGLab data structure with additional fields
% virtual_srate and icaact_filtered_resampled. See
% the readme file for the toolbox for details on these fields.
%
% celleegs (optional): Cell array containing EEGLab data structures with
% the same ICA decomposition.
%
% Output:
% mean_vars1s: Vector of length equal to the number of ICs in the input 
% eeg. The vector contains the means of the variances calculated for 
% one-second intervals for each IC time series.
%
% var_vars1s: Vector of length equal to the number of ICs in the input 
% eeg. The vector contains the variances of the variances calculated for 
% one-second intervals for each IC time series.

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
    arg({'eeg', 'EEG'},[],[],'EEGLab data structure.'),...
    arg('celleegs', [], [],...
    'EEGLab data structures in cell array.'));

if isempty(celleegs)
    celleegs{1} = eeg;
end

[mean_vars1s, var_vars1s] = local_variance(eeg, 'interval', 1, 'celleegs', celleegs);
end

