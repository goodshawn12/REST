function [EEG, com] = ICMARC_interface(EEG, feature_set, force_preprocessing)
% ICMARC_interface() calls ICMARC_run() to obtain predicted classes of the
% ICs in EEG and adds these in the field reject.classtype to the the EEG 
% data structure given as input.
%
% Input:
%   EEG: current dataset structure or structure array.
%
%   feature_set: String containing the name of the the feature set to use 
%   to classify ICs. Possible values are 'established_spatial_features', 
%   'spatial2', and 'established_features'. Default:
%   'established_spatial_features'.
%   Possible values:
%   'established_feature_set': Calculate the features in the established
%       feature set consisting of 14 features as described in the
%       accompanying paper.
%   'established_spatial_features': Calculate the features in the
%       established spatial feature set (described in the
%       accompanying paper.)
%   'spatial2': Calculate a set of spatial features without the
%   computationally demanding dipole features. This feature set was not
%   optimized systematically, but by some fiddling by hand.
%
%   force_preprocessing: Boolean. If true, force preprocessing 
%   steps to be performed even if unnecessary. If no temporal or spectral 
%   features are requested, it is not necessary to downsample and filter 
%   time series. However, if the outputs
%   produced during preprocessing are desired, force_preprocessing can be set
%   to true. (default: false)
%
% Output:
%   EEG: The EEGLab struct given as input with the fields
%   EEG.reject.classtype is modified to contain the predicted classes for
%   each IC.

% Copyright (C) 2014 Laura Froelich
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

if nargin < 1
    help ICMARC_interface;
    return;
end


disp(' ')
disp (['Running ICMARC on dataset ' strrep(EEG.filename, '.set', '') '.set'])
[EEG, predclass, predprob] = ICMARC_run(EEG, feature_set, force_preprocessing);

EEG.reject.classtype = predclass;
EEG.reject.probabilities = predprob;
disp('Done classifying ICs, EEG.reject.classtype has been updated.')

% return the command used to invoke this function
% -------------------------

com = ['ICMARC_interface(' inputname(1) ', '''  feature_set ''', ' num2str(force_preprocessing) ');'];

return;