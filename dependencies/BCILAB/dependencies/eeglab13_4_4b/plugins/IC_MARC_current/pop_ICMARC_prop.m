% pop_ICMARC_prop()
% Modified version of pop_prop_ADJ in the ADJUST plugin.
%
%
% Usage:
%   >> com = pop_ICMARC_prop(EEG, typecomp, numcompo, winhandle, predclass)
%
% Inputs
%   EEG        - current dataset structure or structure array
%   typecomp   - [0|1] compute electrode property (1) or component
%                property (0). Default is 1.
%   numcompo   - channel or component number
%   winhandle  - if this parameter is present or non-NaN, buttons for the
%                rejection of the component are drawn. If
%                non-zero, this parameter is used to backpropagate
%                the color of the rejection button.
%
% Output:
%   com: command used to invoke function
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


function com = pop_ICMARC_prop(EEG, typecomp, numcompo, winhandle, predclass, tagname)

global EEG;
classtype = {'blink', 'neural', 'pulse', 'lat. eye', 'muscle', 'mixed', 'electrode'};
com = '';
if nargin < 1
    help pop_ICMARC_prop;
    return;
end;

if nargin == 1
    typecomp = 1;
end;
if typecomp == 0 & isempty(EEG.icaweights)
    error('No ICA weights recorded for this set, first run ICA');
end;
if nargin == 2
    promptstr    = { fastif(typecomp,'Channel number to plot:','Component number to plot:') ...
        'Spectral options (see spectopo help):' };
    inistr       = { '1' '''freqrange'', [2 50]' };
    result       = inputdlg2( promptstr, 'Component properties - pop_ICMARC_prop()', 1,  inistr, 'pop_ICMARC_prop');
    if size( result, 1 ) == 0 return; end;
    
    numcompo   = eval( [ '[' result{1} ']' ] );
    spec_opt   = eval( [ '{' result{2} '}' ] );
end;

% plotting several component properties 
% -------------------------------------
if length(numcompo) > 1
    for index = numcompo
        pop_ICMARC_prop(EEG, typecomp, numcompo, winhandle, predclass(index));
    end;
    com = sprintf('pop_prop( %s, %d, [%s]);', inputname(1), typecomp, int2str(numcompo));
    return;
end;

if numcompo < 1 | numcompo > EEG.nbchan
    error('Component index out of range');
end;

% assumed input is numcompo
% -------------------------
try, icadefs;
catch,
    BACKCOLOR = [0.8 0.8 0.8];
    GUIBUTTONCOLOR   = [0.8 0.8 0.8];
end;
basename = [fastif(typecomp,'Channel ', 'Component ') int2str(numcompo) ];

fh = figure('name', ['pop_ICMARC_prop() - ' basename ' properties'], 'color', BACKCOLOR, 'numbertitle', 'off', 'visible', 'off');
pos = get(gcf,'Position');
set(gcf,'Position', [pos(1) pos(2)-500+pos(4) 500 500], 'visible', 'on');
pos = get(gca,'position'); % plot relative to current axes
hh = gca;
q = [pos(1) pos(2) 0 0];
s = [pos(3) pos(4) pos(3) pos(4)]./100;
axis off;

% plotting topoplot
% -----------------
h = axes('Units','Normalized', 'Position',[-10 65 40 35].*s+q);

%topoplot( EEG.icawinv(:,numcompo), EEG.chanlocs); axis square;

if typecomp == 1 % plot single channel locations
    topoplot( numcompo, EEG.chanlocs, 'chaninfo', EEG.chaninfo, ...
        'electrodes','off', 'style', 'blank', 'emarkersize1chan', 12); axis square;
else             % plot component map
    topoplot( EEG.icawinv(:,numcompo), EEG.chanlocs, 'chaninfo', EEG.chaninfo, ...
        'shading', 'interp', 'numcontour', 3); axis square;
end;
%basename = [fastif(typecomp,'Channel ', 'IC') int2str(numcompo), ': ' classtype(predclass) ];
basename = [fastif(typecomp,'Channel ', 'IC') int2str(numcompo)];
% title([ basename fastif(typecomp, ' location', ' map')], 'fontsize', 14);
title(basename, 'fontsize', 12);

% plotting erpimage
% -----------------
hhh = axes('Units','Normalized', 'Position',[45 67 48 33].*s+q); %era height 38
eeglab_options;
if EEG.trials > 1
    % put title at top of erpimage
    axis off
    hh = axes('Units','Normalized', 'Position',[45 67 48 33].*s+q);
    EEG.times = linspace(EEG.xmin, EEG.xmax, EEG.pnts);
    if EEG.trials < 6
        ei_smooth = 1;
    else
        ei_smooth = 3;
    end
    if typecomp == 1 % plot component
        offset = nan_mean(EEG.data(numcompo,:));
        erpimage( EEG.data(numcompo,:)-offset, ones(1,EEG.trials)*10000, EEG.times , ...
            '', ei_smooth, 1, 'caxis', 2/3, 'cbar','erp');
    else % plot channel
        if option_computeica
            offset = nan_mean(EEG.icaact(numcompo,:));
            erpimage( EEG.icaact(numcompo,:)-offset, ones(1,EEG.trials)*10000, EEG.times , ...
                '', ei_smooth, 1, 'caxis', 2/3, 'cbar','erp', 'yerplabel', '');
        else
            
            icaacttmp = (EEG.icaweights(numcompo,:) * EEG.icasphere) ...
                * reshape(EEG.data, EEG.nbchan, EEG.trials*EEG.pnts);
            %                     icaacttmp = (EEG.icaweights(numcompo,:) * EEG.icasphere) ...
            %                                    * EEG.data(EEG.nbchan, EEG.trials*EEG.pnts);
            offset = nan_mean(icaacttmp);
            erpimage( icaacttmp-offset, ones(1,EEG.trials)*10000, EEG.times, ...
                '', ei_smooth, 1, 'caxis', 2/3, 'cbar','erp', 'yerplabel', '');
        end;
    end;
    axes(hhh);
    title(sprintf('%s activity \\fontsize{10}(global offset %3.3f)', basename, offset), 'fontsize', 12);
else
    
    % put title at top of erpimage
    EI_TITLE = 'Continous data';
    axis off
    hh = axes('Units','Normalized', 'Position',[45 62 48 38].*s+q);
    ERPIMAGELINES = 200; % show 200-line erpimage
    %while size(EEG.data,2) < ERPIMAGELINES*EEG.srate
    %    ERPIMAGELINES = round(0.9 * ERPIMAGELINES);
    %end
    ERPIMAGELINES = ceil(size(EEG.data,2)/EEG.srate);
    if ERPIMAGELINES > 2   % give up if data too small
        if ERPIMAGELINES < 10
            ei_smooth = 1;
        else
            ei_smooth = 3;
        end
        erpimageframes = floor(size(EEG.data,2)/ERPIMAGELINES);
        erpimageframestot = erpimageframes*ERPIMAGELINES;
        eegtimes = linspace(0, erpimageframes-1, EEG.srate/1000);
        if typecomp == 1 % plot component
            offset = nan_mean(EEG.data(numcompo,:));
            erpimage( reshape(EEG.data(numcompo,1:erpimageframestot),erpimageframes,ERPIMAGELINES)-offset, ones(1,ERPIMAGELINES)*10000, eegtimes , ...
                EI_TITLE, ei_smooth, 1, 'caxis', 2/3, 'cbar');
        else % plot channel
            if option_computeica
                offset = nan_mean(EEG.icaact(numcompo,:));
                erpimage( ...
                    reshape(EEG.icaact(numcompo,1:erpimageframestot),erpimageframes,ERPIMAGELINES)-offset, ...
                    ones(1,ERPIMAGELINES)*10000, eegtimes , ...
                    EI_TITLE, ei_smooth, 1, 'caxis', 2/3, 'cbar','yerplabel', '');
            else
                %                     icaacttmp = reshape(EEG.icaweights(numcompo,:) * EEG.icasphere) ...
                %                                      * reshape(EEG.data, erpimageframes, ERPIMAGELINES);
                
                icaacttmp = EEG.icaweights(numcompo,:) * EEG.icasphere ...
                    *EEG.data(:,1:erpimageframes*ERPIMAGELINES);
                
                offset = nan_mean(icaacttmp);
                erpimage( icaacttmp-offset, ones(1,ERPIMAGELINES)*10000, eegtimes, ...
                    EI_TITLE, ei_smooth, 1, 'caxis', 2/3, 'cbar', 'yerplabel', '');
            end;
        end
    else
        axis off;
        text(0.1, 0.3, [ 'No erpimage plotted' 10 'due to too little continuous data']);
    end;
    axes(hhh);
end;

% plotting spectrum
% -----------------

if exist('winhandle', 'var')
    h = axes('units','normalized', 'position',[10 25 85 25].*s+q);
    %h = axes('units','normalized', 'position',[5 10 95 35].*s+q); %%%
    %CHANGE!
else
    h = axes('units','normalized', 'position',[10 15 85 30].*s+q);
    %h = axes('units','normalized', 'position',[5 0 95 40].*s+q); %%%
    %CHANGE!
end;

%h = axes('units','normalized', 'position',[45 5 60 40].*s+q);
try
    eeglab_options;
    %next instr added for correct function! Andrea
    option_computeica=1;
    if typecomp == 1
        %[spectra freqs] = spectopo( EEG.data(numcompo,:), EEG.pnts, EEG.srate, spec_opt{:},'freqrange', [0 45] );
        [spectra freqs] = spectopo( EEG.data(numcompo,:), EEG.pnts, EEG.srate, 'freqrange', [0 45] );
    else
        if option_computeica
            
            %[spectra freqs] = spectopo( EEG.icaact(numcompo,:), EEG.pnts, EEG.srate, 'mapnorm', EEG.icawinv(:,numcompo), spec_opt{:}, 'freqrange', [0 45]);
            % CONTROL ADDED FOR CONTINUOUS DATA
            if size(EEG.data,3)==1
                EEG.icaact = EEG.icaweights*EEG.icasphere*EEG.data;
            end
            [spectra freqs] = spectopo( EEG.icaact(numcompo,:), EEG.pnts, EEG.srate, 'mapnorm', EEG.icawinv(:,numcompo),  'freqrange', [0 45]);
        else
            if exist('icaacttmp')~=1,
                
                icaacttmp = (EEG.icaweights(numcompo,:)*EEG.icasphere)*reshape(EEG.data, EEG.nbchan, EEG.trials*EEG.pnts);
            end;
            
            %[spectra freqs] = spectopo( icaacttmp, EEG.pnts, EEG.srate, 'mapnorm', EEG.icawinv(:,numcompo), spec_opt{:} ,'freqrange', [0 45]);
            [spectra freqs] = spectopo( icaacttmp, EEG.pnts, EEG.srate, 'mapnorm', EEG.icawinv(:,numcompo), 'freqrange', [0 45]);
        end;
    end;
    set(gca,'fontsize',8);
    % set up new limits
    % -----------------
    %freqslim = 50;
    %set(gca, 'xlim', [0 min(freqslim, EEG.srate/2)]);
    %spectra = spectra(find(freqs <= freqslim));
    %set(gca, 'ylim', [min(spectra) max(spectra)]);
    
    %tmpy = get(gca, 'ylim');
    %set(gca, 'ylim', [max(tmpy(1),-1) tmpy(2)]);
    set( get(gca, 'ylabel'), 'string', 'Power 10*log_{10}(\muV^{2}/Hz)', 'fontsize', 8);
    set( get(gca, 'xlabel'), 'string', 'Frequency (Hz)', 'fontsize', 8);
    title('Activity power spectrum', 'fontsize', 12);
catch
    axis off;
    text(0.1, 0.3, [ 'Error: no spectrum plotted' 10 ' make sure you have the ' 10 'signal processing toolbox']);
end;

% ----------------------------------------------------------------
% plotting IC properties
% -----------------

if ~exist('winhandle', 'var')
    h = axes('units','normalized', 'position',[3 2 95 10].*s+q);
else
    h = axes('units','normalized', 'position',[3 0 95 10].*s+q);
end;

axis off



% -----------------------------------------------------------------


% display buttons
% ---------------

if exist('winhandle', 'var')
    COLREJ = '[1 0.6 0.6]';
    COLACC = '[0.75 1 0.75]';
    % CANCEL button
    % -------------
    h  = uicontrol(gcf, 'Style', 'pushbutton', 'backgroundcolor', GUIBUTTONCOLOR, 'string', 'Cancel', 'Units','Normalized','Position',[-15 -10 15 6].*s+q, 'callback', 'close(gcf);');
    
    % VALUE button
    % -------------
%    hval  = uicontrol(gcf, 'Style', 'pushbutton', 'backgroundcolor', GUIBUTTONCOLOR, 'string', 'Values', 'Units','Normalized', 'Position', [5 -10 15 6].*s+q);
    
    % REJECT button
    % -------------
    status = EEG.reject.gcompreject(numcompo);
    hr = uicontrol(gcf, 'Style', 'pushbutton', 'backgroundcolor',...
        fastif(EEG.reject.gcompreject(numcompo)==1, num2str(COLREJ),num2str(COLACC)), ...
        'string', fastif(EEG.reject.gcompreject(numcompo)==1,...
        'REJECT', 'ACCEPT'), 'Units','Normalized', 'Position', [27 -10 15 6].*s+q, 'userdata', status,...
        'tag', 'rejstatus');
    command = [ 'set(gcbo, ''userdata'', get(gcbo, ''userdata'')~=1);' ...
        'if get(gcbo, ''userdata'')==1,' ...
        ...%'     set( gcbo, ''string'', ''REJECT'', ''userdata'',  true ,''backgroundcolor'',' COLREJ ' );' ...
        '     set( gcbo, ''string'', ''REJECT'',''backgroundcolor'',' COLREJ ' );' ...
        'else ' ...
        '     set( gcbo, ''string'', ''ACCEPT'', ''backgroundcolor'',' COLACC ');' ...
        'end;' ];
    set( hr, 'callback', command);
    
    
    % Class type button
    % -------------
    hclass = uicontrol(gcf, 'Style', 'popupmenu', 'backgroundcolor',GUIBUTTONCOLOR, ...
        'string', classtype', 'Units','Normalized', 'Position', [47 -10 20 6].*s+q, 'Value', EEG.reject.classtype(numcompo),...
        'tag', 'classstatus', 'userdata', EEG.reject.classtype(numcompo));
    command =[ 'set(gcbo, ''backgroundcolor'',',...
        '[' num2str(GUIBUTTONCOLOR(1)) ',' num2str(GUIBUTTONCOLOR(2)) ',' num2str(GUIBUTTONCOLOR(3)) ']' ');'...
        'set(gcbo, ''userdata'', get(gcbo, ''Value''));'];
    set( hclass, 'callback', command);
    
    % HELP button
    % -------------
    h  = uicontrol(gcf, 'Style', 'pushbutton', 'backgroundcolor', GUIBUTTONCOLOR, 'string', 'HELP', 'Units','Normalized', 'Position', [70 -10 15 6].*s+q, 'callback', 'pophelp(''pop_ICMARC_prop'');');
    
    % OK button
    % ---------
    
    command = [ 'global EEG;' ...
        'tmpstatus = get( findobj(''parent'', gcbf, ''tag'', ''rejstatus''), ''userdata''); ' ...
        'EEG.reject.gcompreject(' num2str(numcompo) ') = tmpstatus; '... 
        'classtmpstatus = get( findobj(''parent'', gcbf, ''tag'', ''classstatus''), ''userdata'');' ...
        'EEG.reject.classtype(' num2str(numcompo) ') = classtmpstatus;'];    
    if winhandle ~= 0
        %command = [ command ...
        %    sprintf('if tmpstatus set(findobj(''tag'', tagname), ''backgroundcolor'', %s); else set(findobj(''tag'', tagname), ''backgroundcolor'', %s); end;', ...
        %    COLREJ, COLACC)];
        command = [command ...
            'set(findobj(''tag'',  ''' tagname...
            ''' ), ''backgroundcolor'', fastif(tmpstatus == 1, '...
            num2str(COLREJ) ', ' num2str(COLACC) ' ));'];
            %'if tmpstatus==1 set(findobj(''tag'', ''' tagname...
            %' '' ), ''backgroundcolor'','  num2str(COLREJ)...
            %'); else set(findobj(''tag'',''' tagname ...
            %'''), ''backgroundcolor'','  num2str(COLACC) '); end;'];
        
        command = [command ...
            'set(findobj(''tag'',  ''' tagname '''), ''string'', ' '[num2str(' num2str(numcompo)...
            ') '': '' ' ' fastif(classtmpstatus<7, fastif(classtmpstatus<4, fastif(classtmpstatus<3, fastif(classtmpstatus==1, '' '...
             classtype{1} ' '' , '' ' classtype{2} ' ''), '' ' classtype{3} ' '') , fastif(classtmpstatus<6, fastif(classtmpstatus==4, '' '...
             classtype{4} ' '' , '' ' classtype{5} ' ''), '' ' classtype{6} ' '')), '' ' classtype{7} ' '') ' ']' ' );'];
            
    end;
    command = [ command 'close(gcf); clear tmpstatus classtmpstatus' ];
    h  = uicontrol(gcf, 'Style', 'pushbutton', 'string', 'OK', 'backgroundcolor', GUIBUTTONCOLOR, 'Units','Normalized', 'Position',[90 -10 15 6].*s+q, 'callback', command);
    
    
    % MODIFICA
    %com = sprintf('pop_prop( %s, %d, %d, 0, %s);', inputname(1), typecomp, numcompo, vararg2str( { spec_opt } ) );
    com =  sprintf('pop_ICMARC_prop( %s, %d, %d, %d);',...
        inputname(1), typecomp, numcompo, predclass);
    
else
    %com = sprintf('pop_prop( %s, %d, %d, NaN, %s);', inputname(1), typecomp, numcompo, vararg2str( { spec_opt } ) );
    com = sprintf('pop_ICMARC_prop( %s, %d, %d, %d);',...
        inputname(1), typecomp, numcompo, predclass);
    
    
end;
return;
end

function out = nan_mean(in)

nans = find(isnan(in));
in(nans) = 0;
sums = sum(in);
nonnans = ones(size(in));
nonnans(nans) = 0;
nonnans = sum(nonnans);
nononnans = find(nonnans==0);
nonnans(nononnans) = 1;
out = sum(in)./nonnans;
out(nononnans) = NaN;
end
