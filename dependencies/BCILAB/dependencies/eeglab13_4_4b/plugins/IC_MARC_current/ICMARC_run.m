function [EEG, predclass, predprob] = ICMARC_run(EEG, feature_set, force_preprocessing)
% ICMARC_run() calls ic_feature_extraction() to calculate features of the
% ICs in the input EEG dataset and applies a pre-trained classifier to the
% calculated features to predict the classes of the ICs.
%
% Input:
%   EEG: current dataset structure or structure array.
%
%   feature_set: String containing the name of the the feature set to use 
%   to classify ICs. Possible values are 'established_spatial_features', 
%   'spatial2', and 'established_features'.
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
%   force_preprocessing: Force preprocessing steps to be performed even if
%   unnecessary. If no temporal or spectral features are requested, it is not
%   necessary to downsample and filter time series. However, if the outputs
%   produced during preprocessing are desired, force_preprocessing can be set
%   to true. (default: false)
%
% Output:
%   predclass: Vector containing the classes of ICs
%   predprob: Matrix containing probabilities for each IC belonging to each
%   class. Rows contain class membership probabilities for each IC such
%   that each row adds up to one. Columns correspond to the following
%   classes: blink, neural, heart, lateral eye, muscle, and mixed (in that
%   order)
%

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
    help ICMARC_run;
    error('ICMARC_run: The first argument must be given.')
end

if nargin < 2
    feature_set = 'established_spatial_features';
end

switch feature_set
    case 'established_spatial_features'
        loadfile = 'spatial_established_features.mat';
    case 'spatial2'
        loadfile = 'spatial2.mat';
    case 'established_features'
        loadfile = 'established_features.mat';
    case 'established_feature_set'
        loadfile = 'established_features.mat';
    otherwise
   error(['ICMARC_run: the feature set requested is not a valid option. '...
       'Invoke ic_feature_extraction directly to calculate other feature '...
       'sets than those possible through ICMARC_run.'])
end

savedvars = load(loadfile);
mu = savedvars.mu;
mod = savedvars.mod;
sigma = savedvars.sigma;

[ic_feats, EEG] = ic_feature_extraction(EEG, {feature_set}, 'force_preprocessing', force_preprocessing);


featsnorm = ic_feats - repmat(mu, size(ic_feats,1),1);
featsstand = featsnorm./repmat(sigma, size(ic_feats,1),1);
predprob = mnrval(mod, featsstand);
[temp1, predclass] = max(predprob, [], 2); % using temp1 instead of ~ for 
% backward compatibility with older versions of Matlab
end

