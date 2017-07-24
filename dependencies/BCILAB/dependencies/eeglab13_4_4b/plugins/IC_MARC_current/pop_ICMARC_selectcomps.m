% pop_ICMARC_selectcomps() - Display component scalp maps, their classifications,
%                  rejection status. This function is a modified version of
%                  pop_selectcomps_ADJ in the ADJUST plugin. The function pop_selectcomps_ADJ is
%                  based on pop_selectcomps.
%
% Usage:
%   >> [EEG, com] = pop_ICMARC_selectcomps( EEG, compnum, fig );
%
% Inputs:
%   EEG: Current dataset structure or structure array
%   compnum: Vector of component numbers to plot
%   fig: Handle to figure in which scalp map(s) should be plotted
%
% Outputs:
%   EEG: Output dataset with updated rejected components
%   com: Command used to invoke function
%
%
% Modified by Laura Froelich.

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


function [EEG,com] = pop_ICMARC_selectcomps(EEG, compnum, fig )
classtype = {'blink', 'neural', 'heart', 'lat. eye', 'muscle', 'mixed'};
COLREJ = '[1 0.6 0.6]';
COLACC = '[0.75 1 0.75]';
COLART = '[1 0 0]';
PLOTPERFIG = 35;

com = [ 'pop_selectcomps(' inputname(1) ', ' vararg2str(compnum) ');' ];
if nargin < 1
    help pop_ICMARC_selectcomps;
    return;
end;


if isempty(EEG.icaact)
    disp('Warning: EEG.icaact missing! Recomputed using eeg_getica for plotting purposes');
    % Next instruction: see eeg_checkset
    EEG.icaact = eeg_getica(EEG);
end;

if nargin < 2
    promptstr = { 'Components to plot:' };
    initstr   = { [ '1:' int2str(size(EEG.icaweights,1)) ] };
    
    result = inputdlg2(promptstr, 'Reject comp. by map -- pop_selectcomps',1, initstr);
    if isempty(result), return; end;
    compnum = eval( [ '[' result{1} ']' ]);
    
    if length(compnum) > PLOTPERFIG
        ButtonName=questdlg2(strvcat(['More than ' int2str(PLOTPERFIG) ' components so'],'this function will pop-up several windows'), ...
            'Confirmation', 'Cancel', 'OK','OK');
        if ~isempty( strmatch(lower(ButtonName), 'cancel')), return; end;
    end;
    
end;

if ~isfield(EEG.reject,'classtype') || isempty(EEG.reject.classtype)
    supergui( 'geomhoriz', { 1 1 1 }, 'uilist', { ...
        {'style', 'text', 'string', 'Warning: Classtype (EEG.reject.classtype) not found,'},...
        {'style', 'text', 'string', ' initializing with class ''mixed'' for all components'},...
        {'style', 'pushbutton', 'string', 'OK' 'callback' 'close(gcbf);'}...
        });
    warning('Classtype (EEG.reject.classtype) not found, initializing with class ''mixed'' for all components')
    EEG.reject.classtype = ones(1,size(EEG.icawinv,2))*6;
end
origclasstypefield = EEG.reject.classtype;
origrejfield = EEG.reject.gcompreject;


predclasses = EEG.reject.classtype;

fprintf('Drawing figure...\n');
currentfigtag = ['selcomp' num2str(rand)]; % generate a random figure tag

% Call this function several times if there are too many components to plot
% in one figure
if length(compnum) > PLOTPERFIG
    for index = 1:PLOTPERFIG:length(compnum)
        EEG = pop_ICMARC_selectcomps(EEG, compnum([index:min(length(compnum),index+PLOTPERFIG-1)]));
    end;
    
    com = [ 'pop_ICMARC_selectcomps(' inputname(1) ', ' vararg2str(compnum)  ');' ];
    return;
end;

if ~isfield(EEG.reject, 'gcompreject') || isempty(EEG.reject.gcompreject)
    EEG.reject.gcompreject = ones(1, size(EEG.icawinv,2));
    EEG.reject.gcompreject(EEG.reject.classtype==2) = 0;
    EEG.reject.gcompreject(EEG.reject.classtype==6) = 0;
end;

if exist('icadefs', 'file');
    icadefs;
end
if ~exist('BACKCOLOR', 'var')
    BACKCOLOR = [0.8 0.8 0.8];
end;
if ~exist('GUIBUTTONCOLOR', 'var')
    GUIBUTTONCOLOR   = [0.8 0.8 0.8];
end

% set up the figure
% -----------------
column =ceil(sqrt( length(compnum) ))+1;
rows = ceil(length(compnum)/column);
if ~exist('fig', 'var') % 12/12/2014: the second argument, 'var', was added to the function call
    fighandle = figure('name', [ 'Reject components by map - pop_ICMARC_selectcomps() (dataset: ' EEG.setname ')'],...
        'tag', currentfigtag, 'numbertitle', 'off', 'color', BACKCOLOR);
    set(fighandle,'MenuBar', 'none');
    pos = get(fighandle,'Position');
    set(fighandle,'Position', [pos(1) 20 800/7*column 600/5*rows]);
    incx = 120;
    incy = 110;
    sizewx = 100/column;
    if rows > 2
        sizewy = 90/rows;
    else
        sizewy = 80/rows;
    end;
    pos = get(gca,'position'); % plot relative to current axes
    hh = gca;
    q = [pos(1) pos(2) 0 0];
    s = [pos(3) pos(4) pos(3) pos(4)]./100;
    axis off;
end;

% figure rows and columns
% -----------------------
if EEG.nbchan > 64
    disp('More than 64 electrodes: electrode locations not shown');
    plotelec = 0;
else
    plotelec = 1;
end;

count = 1;
for ri = compnum
    if exist('fig', 'var') % 12/12/2014: the second argument, 'var', was added to the function call
        button = findobj('parent', fig, 'tag', ['comp' num2str(ri)]);
        if isempty(button)
            error( 'pop_ICMARC_selectcomps(): figure does not contain the component button');
        end;
    else
        button = [];
    end;
    
    if isempty( button )
        % compute coordinates
        % -------------------
        X = mod(count-1, column)/column * incx-10;
        Y = (rows-floor((count-1)/column))/rows * incy - sizewy*1.3;
        
        % plot the head
        % -------------
        if ~strcmp(get(fighandle, 'tag'), currentfigtag);
            disp('Aborting plot');
            return;
        end;
        ha = axes('Units','Normalized', 'Position',[X Y sizewx sizewy].*s+q);
        if plotelec
            topoplot( EEG.icawinv(:,ri), EEG.chanlocs, 'verbose', ...
                'off', 'style' , 'fill', 'chaninfo', EEG.chaninfo);
        else
            topoplot( EEG.icawinv(:,ri), EEG.chanlocs, 'verbose', ...
                'off', 'style' , 'fill','electrodes','off', 'chaninfo', EEG.chaninfo);
        end;
        axis square;
        
        % plot the button
        % ---------------
        tagnameprop = ['comp' num2str(ri)];
        tagname = tagnameprop;
        while exist(tagname, 'var')
            tagname = [tagnameprop '_' num2str(randn(1))];
        end
        button = uicontrol(fighandle, 'Style', 'pushbutton', 'Units','Normalized', 'Position',...
            [X Y+sizewy sizewx sizewy*0.25].*s+q, 'tag', tagname);
        %command = sprintf('pop_ICMARC_prop( %s, 0, %d, %3.15f, %d);', ...
        %    inputname(1), ri, button, EEG.reject.classtype(ri));
        %set( button, 'callback', command );
        set(button, 'callback', {@callpop_ICMARC_prop, EEG,  0, ri,...
            button, EEG.reject.classtype(ri), tagname})
    end;
    
    % MODIFY BUTTON COLOR and class string: ARTIFACT IC?
    %set( button, 'backgroundcolor', eval(fastif(EEG.reject.gcompreject(ri), COLREJ,COLACC)),...
    %    'string', [int2str(ri) ': ' classtype{predclasses(ri)}]);
    set( button, 'backgroundcolor', eval(fastif(EEG.reject.gcompreject(ri), COLREJ,COLACC)),...
        'string', [int2str(ri) ': ' classtype{predclasses(ri)}]);
       
    drawnow;
    count = count +1;
end;

setappdata(fighandle, 'origclasstypefield', origclasstypefield);
setappdata(fighandle, 'origrejfield', origrejfield);
setappdata(fighandle, 'compnum', compnum); % the component number plotted in the figure

% draw the bottom button
% ----------------------
if ~exist('fig', 'var') % 13/12/2014: second argument added to "exist" function call
    uicontrol(fighandle, 'Style', 'pushbutton', 'string', 'Cancel', 'Units','Normalized',...
        'backgroundcolor', GUIBUTTONCOLOR, 'Position',[20 -10  15 sizewy*0.25].*s+q,...
        'callback', {@startcancelcommand, fighandle} );
    
    uicontrol(fighandle, 'Style', 'pushbutton', 'string', 'Help', 'Units','Normalized',...
        'backgroundcolor', GUIBUTTONCOLOR, 'Position',[40 -10  15 sizewy*0.25].*s+q,...
        'callback', 'pophelp(''pop_ICMARC_selectcomps'');');
    
    command = '[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET); eegh(''[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);''); close(gcf)';
    
    uicontrol(fighandle, 'Style', 'pushbutton', 'string', 'OK', 'Units','Normalized', 'backgroundcolor', GUIBUTTONCOLOR, ...
        'Position',[60 -10  15 sizewy*0.25].*s+q, 'callback',  command);
end;

end

function EEG=startcancelcommand(hObj,evnt, fighandle)
global EEG;
origclasstypefield = getappdata(fighandle, 'origclasstypefield');
origrejfield = getappdata(fighandle, 'origrejfield');
compnum = getappdata(fighandle, 'compnum');
EEG.reject.classtype(compnum) = origclasstypefield(compnum);
EEG.reject.gcompreject(compnum) = origrejfield(compnum);
close(gcf); fprintf('Operation cancelled\n')
end

function callpop_ICMARC_prop(hobj, evnt, EEG, typecomp, numcompo, winhandle, predclass, tagname)
pop_ICMARC_prop(EEG, typecomp, numcompo, winhandle, predclass, tagname);
end
