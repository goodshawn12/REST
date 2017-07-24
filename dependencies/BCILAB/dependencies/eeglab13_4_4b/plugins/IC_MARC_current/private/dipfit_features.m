function [xcoord, ycoord, zcoord, ndipole_labels, eeg] = dipfit_features(varargin)
% dipfit_features calculates a number of features (described under output) based on 
% the dipole fit.
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
% xcoord: The X-coordinate of the dipole fit.
%
% ycoord: The Y-coordinate of the dipole fit.
%
% zcoord: The Z-coordinate of the dipole fit.
%
% ndipole_labels: The number of anatomical labels that match the dipole
% fit. The labels are added by the function add_dipfit_labels.
%
% eeg: The EEGlab data structure given as input, with the dipole fit found 
% in this function.

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
folder_path = filename(1:(strfind(filename, 'dipfit_features')-1));

arg_define(1,varargin, ...
    arg('eeg', 'EEG'), ...
    arg('besapath', folder_path, [],...
    'Path to BESA files.'));

eegtemp=eeg; % make a copy of the eeg struct that was given as input that
% can be modified freely without affecting the eeg struct that will be
% returned
eegtemp.icawinv = eegtemp.virtual_topography';
nics = size(eegtemp.icawinv, 2);
if nics < size(eegtemp.icawinv, 1) % to avoid dipfit error: "inconsistent number of channels in trial 1"
    eegtemp.icawinv = [eegtemp.icawinv, repmat(eegtemp.icawinv(:,1), 1,...
        size(eegtemp.icawinv, 1)-nics)]; % the time it takes to obtain the 
    % dipole fit is independent of the number of ICs being fitted according
    % the DIPFIT and these extra fits are thrown away at the end of this
    % function.
end

eegtemp.icaweights = pinv(eegtemp.icawinv); % this is not used for dipole 
% fitting, so the redundant ICs added above do not influence the fit.
% (EEG.icawinv is used in the dipole fit, see the field data.topo set in
% line 129 in eeglab2fieldtrip.m)
eegtemp.chanlocs = eegtemp.virtual_chanlocs;
eegtemp.icasphere = eye(length(eegtemp.chanlocs));
eegtemp.icachansind = 1:length(eegtemp.chanlocs);
% the fields data and icaact are not used for anything here, but must 
% be present to pass various eeg_checkset and fieldtrip calls.
eegtemp.data  = randn(eegtemp.virtual_nbchan(2), eegtemp.pnts, eegtemp.trials);
eegtemp.icaact = randn(64, eegtemp.pnts, eegtemp.trials);
eegtemp.nbchan =eegtemp.virtual_nbchan(2);
eegtemp.srate = eeg.virtual_srate;


% if a dipole fit is not already attached to the eeg struct.
if(~isfield(eeg, 'dipfit') || isempty(eeg.dipfit)) || ~isfield(eeg.dipfit, 'model')
    eeglab_options
    oldoption_computeica = option_computeica;
    pop_editoptions('option_computeica', 1);
    
    
    eegtemp = pop_select(eegtemp, 'channel', eegtemp.icachansind);
    eegtemp = pop_dipfit_settings(eegtemp ...
        ,'mrifile', [besapath 'avg152t1.mat']...
        ,'hdmfile',[besapath 'standard_BESA.mat'] ...
        ,'coordformat', 'Spherical', 'coord_transform', [0 0 -10 0 0 0 894.4470 1020 894.4470] ...
        ,'chanfile', [besapath 'standard-10-5-cap385.elp']...
        );
    
    eegtemp = pop_multifit(eegtemp, 1:size(eegtemp.icawinv, 2), 'threshold', 100);
    
    eeg.dipfit = eegtemp.dipfit; % add the fitted dipole fit to the eeg struct given as input
    eeg.dipfit.model = eegtemp.dipfit.model(1:nics);
    pop_editoptions('option_computeica', oldoption_computeica);
end

eegtemp = add_anatomical_labels(eegtemp);
xcoord = NaN(nics,1);
ycoord = NaN(nics,1);
zcoord = NaN(nics,1);

eegtemp.icawinv = eegtemp.icawinv(:, 1:nics); % throw away the redundant ICs in beginning of function
% the number of labels (anatomical brain areas) given to each IC
ndipole_labels = nareas_identified(eegtemp);
ndipole_labels = reshape(ndipole_labels, length(ndipole_labels),1);

dipole_locs = {eeg.dipfit.model.posxyz};
for ic=1:nics
    dip = dipole_locs{ic};
    xcoord(ic) =  dip(1, 1);
    ycoord(ic) =  dip(1, 2);
    zcoord(ic) =  dip(1, 3);
end
end
