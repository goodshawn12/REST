function fcm_showcorr
% Called by nut_results_viewer to display correlation plots.

SCALE=1;

global rivets beam corrs

if isempty(rivets), error('NUT_RESULTS_VIEWER is not open.'), end

%corrs.fig=findobj('Tag','fcmcorr');
if isempty(corrs)
    corrs.fig=figure('Name','Correlation','NumberTitle','off','Units','normalized', ...
        'position',[0.55    0.10    0.15    0.20],'color','w','Tag','fcmcorr', ...
        'CloseRequestFcn','clear global corrs; closereq','InvertHardcopy','off');
    corrs.showpatnr=false;
    corrs.linecol='b';
    corrs.bback=false;
else
    figure(corrs.fig);
end

load(beam.corr.imagingdata);
if exist('pop','var')   
    imdata = pop.s(:,rivets.currMEGvoxelindex,rivets.timeselect,rivets.freqselect);
% elseif exist('allsubj','var')   % legacy compatibilty
%     imdata = squeeze(allsubj(rivets.currMEGvoxelindex,:,rivets.timeselect));
else
    error('Invalid imagingdata file.')
end
corrs.axis=gca;
cla

corrs.hp = plot(imdata,SCALE*beam.corr.behavdata,'ks');

good = find(isfinite(beam.corr.behavdata))';
corrs.ht = zeros(length(good),1);
for k=1:length(good)
    corrs.ht(k)=text(imdata(good(k)),SCALE*beam.corr.behavdata(good(k)),[' ' int2str(beam.corr.nr(good(k)))]);
end
if ~corrs.showpatnr, set(corrs.ht,'visible','off'); end

minx=min(imdata(good));
maxx=max(imdata(good));
if minx==maxx, minx=minx-.5; maxx=maxx+.5; end
miny=min(beam.corr.behavdata);
maxy=max(beam.corr.behavdata);
axis([minx-0.2*(maxx-minx) maxx+0.2*(maxx-minx) miny-0.2*(maxy-miny) maxy+0.2*(maxy-miny)]);
set(corrs.axis,'linewidth',2,'fontweight','bold')

warning('off','MATLAB:polyfit:RepeatedPointsOrRescale')
corrs.hl=lsline;
warning('on','MATLAB:polyfit:RepeatedPointsOrRescale')
set(corrs.hl,'color',corrs.linecol,'linewidth',2);

if corrs.bback, blackback; end

xlabel('Neural');
ylabel('Behavioral');
if isfield(pop,'R')
    sr = find(pop.R.roi2voxel_tfm(rivets.currMEGvoxelindex,:),1);
    if ~isempty(sr)
        title(pop.R.roilabel{pop.R.goodroi(sr)},'Interpreter','none')
    end
end
    
% Define the context menues
axmenu = uicontextmenu;
set(corrs.axis,'UIContextMenu',axmenu);
uimenu(axmenu, 'Label','Show Patient Numbers', 'Callback', 'global corrs; corrs.showpatnr=true; set(corrs.ht,''visible'',''on'')');
uimenu(axmenu, 'Label','Hide Patient Numbers', 'Callback', 'global corrs; corrs.showpatnr=false; set(corrs.ht,''visible'',''off'')');
uimenu(axmenu, 'Label','Invert Colors', 'Callback', @blackback);
uimenu(axmenu, 'Label','Uninvert Colors', 'Callback', @whiteback);

linemenu = uicontextmenu;
set(corrs.hl,'UIContextMenu',linemenu);
uimenu(linemenu, 'Label','Blue', 'Callback', 'global corrs; corrs.linecol=''b''; set(corrs.hl,''color'',''b'')');
uimenu(linemenu, 'Label','Black', 'Callback', 'global corrs; corrs.linecol=''k''; set(corrs.hl,''color'',''k'')');


function blackback(dum1,dum2)
global corrs;
corrs.bback=true; 
set(gcf,'color','k')
set(corrs.axis,'color','k','xcolor','w','ycolor','w');
set([corrs.hp;corrs.ht],'color','w');
switch corrs.linecol
    case 'b'
        set(corrs.hl,'color','y')
    case 'k'
        set(corrs.hl,'color','w')
end

function whiteback(dum1,dum2)
global corrs;
corrs.bback=false; 
set(gcf,'color','w')
set(corrs.axis,'color','w','xcolor','k','ycolor','k');
set([corrs.hp;corrs.ht],'color','k');
set(corrs.hl,'color',corrs.linecol)


