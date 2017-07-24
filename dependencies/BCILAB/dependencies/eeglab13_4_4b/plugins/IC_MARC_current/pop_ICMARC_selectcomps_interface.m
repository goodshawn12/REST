function [EEG,com] = pop_ICMARC_selectcomps_interface(EEG, varargin)
% pop_ICMARC_selectcomps_interface() offers a graphical user interface to receive
% inputs for pop_ICMARC_selectcomps() if all the inputs required by
% pop_ICMARC_selectcomps() are not given.
%
% Input:
%   EEG: Current dataset structure or structure array.
%
% Optional input, if both are not given, a window will pop up to ask for
% their values.
%   componentstoplot: (Optional) string containing the component numbers to
%   plot.
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
       componentstoplot = varargin{1};
end
       
if nargin == 1 || isempty(componentstoplot)
    	txtcomponents= 'Components to plot';
        

    uilist = { { 'style'   'text'     'string'    txtcomponents } ...
               { 'style'   'edit'     'string'    ['1:' num2str(size(EEG.icawinv,2))]} ...
               { } ...
               };
           
           uigeom = { [1.5 1] [1]};
           
    guititle = 'pop_ICMARC_selectcomps_interface';
    
    result = inputgui( uigeom, uilist, 'pophelp(''pop_ICMARC_selectcomps_interface'')', guititle, [], 'normal');
	if isempty(result), return; end;
	componentstoplot = str2num(result{1});
end



EEG = pop_ICMARC_selectcomps(EEG, componentstoplot);
% return the string command
% -------------------------
com = ['pop_ICMARC_selectcomps_interface(' inputname(1) ', [' num2str(componentstoplot) ']);'];

return;
