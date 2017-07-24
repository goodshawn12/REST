function [EEG,com] = pop_ICMARC_interface(EEG, varargin)
% pop_ICMARC_interface() offers a graphical user interface to receive
% inputs for ICMARC_interface() if all the inputs required by
% ICMARC_interface() are not given.
%
% Input:
%   EEG: Current dataset structure or structure array.
%
% Optional input, if both are not given, a window will pop up to ask for
% their values.
%   feature_set: (Optional) string containing the name of the the feature set to use 
%   to classify ICs. Possible values are 'established_spatial_features', 
%   'spatial2', and 'established_features'. Default:
%   'established_spatial_features'.
%   Possible values:
%   'established_features': Calculate the features in the established
%       feature set consisting of 14 features as described in the
%       accompanying paper.
%   'established_spatial_features': Calculate the features in the
%       established spatial feature set (described in the
%       accompanying paper.)
%   'spatial2': Calculate a set of spatial features without the
%   computationally demanding dipole features. This feature set was not
%   optimized systematically, but by some fiddling by hand.
%
%   force_preprocessing: (Optional) boolean. If true, force preprocessing 
%   steps to be performed even if unnecessary. If no temporal or spectral 
%   features are requested, it is not necessary to downsample and filter 
%   time series. However, if the outputs
%   produced during preprocessing are desired, force_preprocessing can be set
%   to true. (default: false)
%
% Output:
%   EEG: The input EEGLab struct, with the field EEG.reject.classtype
%   added.
%
%   com: Command used to invoke this function
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



com = ''; % this initialization ensures that the function will return something
% if the user presses the cancel button

% display help if not enough arguments
% ------------------------------------
if nargin < 1
    help pop_ICMARC_interface;
    return;
end;

if nargin >= 2
       feature_set = varargin{1};
      if nargin >= 3
          force_preprocessing = varargin{2};
      else
          force_preprocessing = [];
      end
else
    feature_set = [];
    force_preprocessing = [];
end

feature_sets ={'established_spatial_features','spatial2','established_features'};
if nargin == 1 || (isempty(feature_set) && isempty(force_preprocessing))
    	txtfeatureset1= 'Feature set';
        
        uilist = { { 'style'   'text'     'string'    txtfeatureset1 } ...
               { 'style'   'listbox'     'string'    'established_spatial_features|spatial2|established_features' } ...
               { } ...
               { 'style'   'text'     'string'    'Force temporal and spectral preprocessing' } ...
               { 'style'   'checkbox' 'string'    '' } { } ...
               { } };
           
           uigeom = { [1.5 1] [1] [1.55 0.2 0.8] [1] };
           
    guititle = 'pop_ICMARC_interface';
    
    result = inputgui( uigeom, uilist, 'pophelp(''pop_ICMARC_interface'')', guititle, [], 'normal');
	if isempty(result), return; end;
	feature_set = feature_sets{result{1}};
    force_preprocessingstr = result{2};
    force_preprocessing = fastif(force_preprocessingstr==1, true, false);
end

EEG = ICMARC_interface(EEG, feature_set, force_preprocessing);
% return the string command
% -------------------------
com = ['pop_ICMARC_interface(' inputname(1) ', '''  feature_set ''', ' num2str(force_preprocessing) ');'];

return;
