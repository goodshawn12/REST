function vis_topo(varargin)
% handles:
% handles.ics(1-8) has which ICs to plot

% get the updated stream buffer
Winv = evalin('base','Winv');
Winv = Winv + rand(size(Winv));
handles = varargin{3};

it = mod(get(varargin{1},'TasksExecuted')-1,8)+1;
hstr = ['axesIC' int2str(it)];
[map, cmin, cmax] = topoplotUpdate(Winv(:,handles.ics(it)), handles.chanlocs,'electrodes','off','gridscale',32);
hand = get(handles.(hstr),'children');
set(hand(end),'CData',map);
set(handles.(hstr),'CLim',[cmin cmax]);
