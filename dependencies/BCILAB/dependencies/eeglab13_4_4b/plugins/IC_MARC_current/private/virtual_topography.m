function eeg = virtual_topography(varargin)
% virtual_topography adds the fields virtual_chanlocs, virtual_nbchan, and
% virtual_topography to the input EEGLab data structure. The field icawinv
% in the input is standardized.
%
% Input:
% eeg: EEGLab data structure. All channels must have associated locations,
% otherwise an error results.
%
% Output:
% eeg: The EEGLab data structure that was given as input, with the 
% additional fields virtual_chanlocs, virtual_nbchan, and 
% virtual_topography. The virtual channel locations field (virtual_chanlocs)
% consists of the locations in the standard 64 electrode EEG cap and have
% the standard labels. Correspondingly, the number of virtual channels
% (virtual_nbchan) is 64. The virtual topography (virtual_topography) is
% calculated for each IC by summing, for each virtual electrode location,
% all original electrode activations weighted by Gaussian kernel with
% standard deviation 0.5cm and center in the virtual electrode location. A
% head radius of 9cm is used. The field virtual_topography is then a matrix
% with rows containing the IC's topographies, and each column corresponding
% to a virtual electrode.
% The field icawinv is standardized to have columns with mean zero and 
% variance one. 

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
    arg('eeg', 'EEG', [], 'EEGLab struct.'));

eeg=remove_channels_wo_locations(eeg);


eeg = spatially_normalize_icdecomp(eeg);
topog = eeg.icawinv';

labels = {'Fp1', 'Fp2', 'F3', 'F4', 'C3', 'C4', 'P3', 'P4', 'O1', 'O2', 'F7', 'F8', 'T7', 'T8', 'P7', 'P8', 'Fz', 'Pz', 'FC1', 'FC2', 'CP1',...
'CP2', 'FC5', 'FC6', 'CP5', 'CP6', 'F9', 'F10', 'TP9', 'TP10', 'Cz', 'Oz', 'F1', 'F2', 'C1', 'C2', 'P1', 'P2', 'AF3', 'AF4',...
'FC3', 'FC4', 'CP3', 'CP4', 'PO3', 'PO4', 'F5', 'F6', 'C5', 'C6', 'P5', 'P6', 'AF7', 'AF8', 'FT7', 'FT8', 'TP7', 'TP8', 'PO7', ...
'PO8', 'AFz', 'FCz', 'CPz', 'POz'};

chanlocs = struct('labels', lower(labels));
virtual_chanlocs=pop_chanedit(chanlocs, ...
    'lookup', 'standard-10-5-cap385.elp');

for i=1:length(virtual_chanlocs)
    % set head radius to 9 cm
    virtual_chanlocs=pop_chanedit(virtual_chanlocs, 'changefield',{i 'sph_radius' '9'},'convert',{'sph2all'});
end

thetas = cell2mat({virtual_chanlocs.sph_theta});
phis =  cell2mat({virtual_chanlocs.sph_phi});
sigmas = repmat(0.5, 64,1);
head_radius=9;

n_ics = size(eeg.icawinv,2);
virtual_topog = NaN(n_ics,64);

for ic=1:n_ics
activations = virtual_electrode_activation(phis, thetas, sigmas, eeg.chanlocs, topog(ic, :), head_radius);
virtual_topog(ic,:) = activations;
end

eeg.virtual_topography = zscore(virtual_topog,[],2);
eeg.virtual_chanlocs = virtual_chanlocs;
eeg.virtual_nbchan = size(labels);
end
