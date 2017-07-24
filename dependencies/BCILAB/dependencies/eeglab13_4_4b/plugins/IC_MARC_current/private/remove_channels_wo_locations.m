function [eeg, missing_channel_locations] = remove_channels_wo_locations(varargin)
% remove_channels_wo_locations returns the EEGLab data structure given as input,
% but without channels that did not have locations.
%
% Input:
% EEGLab data structure.
%
% Output:
% eeg: The same data structure that was given as input, but without channels
% that did not have locations or not used in the ICA decomposition. These 
% channels have been removed by the function pop_select from EEGLab.
%
% missing_channel_locations: Vector of removed channels (original indices).

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
    arg({'eeg'},[],[],'EEGLab struct containing EEG data.'));

x_chanlocs = {eeg.chanlocs.X};
missing_channel_locations = find(cellfun(@isempty, x_chanlocs));
if(~isempty(missing_channel_locations))
    eeg = pop_select(eeg, 'nochannel', missing_channel_locations);
end
