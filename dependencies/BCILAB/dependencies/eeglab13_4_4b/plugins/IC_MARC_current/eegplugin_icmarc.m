% eegplugin_icmarc()
%
% Usage:
%   >> eegplugin_icmarc( fig, try_strings, catch_strings);
%
% Inputs:
%
%   fig            - [integer]  EEGLAB figure
%   try_strings    - [struct] "try" strings for menu callbacks.
%   catch_strings  - [struct] "catch" strings for menu callbacks.
%
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


function eegplugin_icmarc(fig, try_strings, catch_strings)
toolsmenu = findobj(fig, 'tag', 'tools');
plotmenu = findobj(fig, 'tag', 'plot');

% create command to run IC_MARC
comrunicmarc = '[EEG,LASTCOM]= pop_ICMARC_interface(EEG); eeglab redraw;';
comrunicmarc_final = [try_strings.no_check comrunicmarc];
%comrunicmarc_final = [comrunicmarc_final 'LASTCOM = ''' comrunicmarc_plot ''';'];
comrunicmarc_final = [comrunicmarc_final catch_strings.store_and_hist ]; 

% create command to plot scalp maps with classes
complot = '[EEG, LASTCOM] = pop_ICMARC_selectcomps_interface(EEG);';
complotfinal = [try_strings.no_check complot];
%complotfinal = [complotfinal 'LASTCOM = ''' complotfinal ''';'];
complotfinal = [complotfinal catch_strings.store_and_hist ];

% create menus
% ------------
submenutools = uimenu( toolsmenu, 'label', 'ICMARC', 'separator', 'on');
submenuplot = uimenu( plotmenu, 'label', 'ICMARC', 'separator', 'on');
uimenu( submenutools, 'Label', 'Run ICMARC', 'CallBack', comrunicmarc_final, ...
    'userdata', 'epoch:on;continuous:on;startup:off;study:off;chanloc:on');

uimenu( submenuplot, 'Label', 'Plot scalpmaps with classes'  , 'CallBack', complotfinal, ...
    'userdata', 'epoch:on;continuous:on;startup:off;study:off;chanloc:on', 'interruptible', 'on');
end