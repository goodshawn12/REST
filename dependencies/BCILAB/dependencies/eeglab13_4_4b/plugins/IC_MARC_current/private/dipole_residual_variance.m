function [residual_variance, eeg] = dipole_residual_variance(varargin)
% residual_variance returns the residual variance of the dipole fit. That
% is, the percentage of variation in scalp activation that is not explained
% by the dipole fit.
%
% Input:
% eeg: EEGLab data structure with additional fields
% virtual_topography, virtual_chanlocs, and icaact_filtered_resampled. See
% the readme file for the toolbox for details on these fields.
%
% besapath: The path on the computer to BESA files in EEGLab plug-in
% dipfit.
%
% Output:
% residual_variance: The residual variance of the dipole fit.

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


filename = mfilename('fullpath');
folder_path = filename(1:(strfind(filename, 'dipole_residual_variance')-1));

arg_define(1,varargin, ...
    arg('eeg', 'EEG'), ...
    arg('besapath', folder_path, [],...
    'Path to BESA files.'));

eegtemp=eeg;
eegtemp.icawinv = eegtemp.virtual_topography';
nics = size(eegtemp.icawinv, 2);
if nics < size(eegtemp.icawinv, 1) % to avoid dipfit error: "inconsistent number of channels in trial 1"
    eegtemp.icawinv = [eegtemp.icawinv, repmat(eegtemp.icawinv(:,1), 1,...
        size(eegtemp.icawinv, 1)-nics)];
end
eegtemp.icaweights = pinv(eegtemp.icawinv);
eegtemp.chanlocs = eegtemp.virtual_chanlocs;
eegtemp.icasphere = eye(length(eegtemp.chanlocs));
eegtemp.icachansind = 1:length(eegtemp.chanlocs);
% the field data is not used for anything here, but must be present to pass
% various eeg_checkset calls.
eegtemp.data  = randn(size(eegtemp.icawinv,1), size(eegtemp.icaact_filtered_resampled,2), size(eegtemp.icaact_filtered_resampled,3));
eegtemp.nbchan =eegtemp.virtual_nbchan;

% if a dipole fit is not already attached to the eeg struct.
if(~isfield(eeg, 'dipfit') || isempty(eeg.dipfit))
    eegtemp = pop_select(eegtemp, 'channel', eegtemp.icachansind);
    eegtemp = pop_dipfit_settings(eegtemp ...
        ,'mrifile', [besapath 'avg152t1.mat']...
        ,'hdmfile',[besapath 'standard_BESA.mat'] ...
        ,'coordformat', 'Spherical', 'coord_transform', [0 0 -10 0 0 0 894.4470 1020 894.4470] ...
        ,'chanfile', [besapath 'standard-10-5-cap385.elp']...
        );
    eegtemp = pop_multifit(eegtemp, 1:size(eegtemp.icaact_filtered_resampled, 1), 'threshold', 100);

    eeg.dipfit = eegtemp.dipfit; % add the fitted dipole fit to the eeg struct given as input
    eeg.dipfit.model = eegtemp.dipfit.model(1:nics); 
end

residual_variance = cell2mat({eeg.dipfit.model.rv});
