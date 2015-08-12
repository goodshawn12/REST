function vis_topo(varargin)
% handles:
% handles.ics(1-8) has which ICs to plot

% get the updated stream buffer
W = evalin('base','W');
sphere = evalin('base','sphere');
Winv = inv(W*sphere);
handles = varargin{3};

it = mod(get(varargin{1},'TasksExecuted')-1,8)+1;
hstr = ['axesIC' int2str(it)];
hand = get(handles.(hstr),'children');
[map, cmin, cmax] = topoplotUpdate(Winv(:,handles.ics(it)), handles.chanlocs,'electrodes','off','gridscale',32);
set(hand(end),'CData',map);
set(handles.(hstr),'CLim',[cmin cmax]);
