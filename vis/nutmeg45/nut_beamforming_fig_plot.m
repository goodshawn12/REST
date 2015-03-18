function [lh,lh2,newlh,fnum1,fnum2] = nut_beamforming_fig_plot

lh = findall(0,'tag','nut_beamforming_fig');
lh2 = findall(lh,'type','axes');
c = get(lh,'colormap');

fnum1 = figure;
set(fnum1,'colormap',c);

fnum2 = figure;
set(fnum2,'colormap',c);

pos = get(fnum1,'position');
pos(2:4) = [max([20 pos(2)-100]) 470 320];
set(fnum1,'position',pos); % change size

pos = get(fnum2,'position');
pos(1:4) = [max([20 pos(2)-100]) pos(2) 350 320];
set(fnum2,'position',pos); % change size

newlh = copyobj(lh2(2),fnum1);
newlh(2) = copyobj(lh2(1),fnum2);

pos = get(newlh(1),'position');
pos(2) = 40;
set(newlh(1),'position',pos);

pos = get(newlh(2),'position');
pos(1:2) = [40 40];
set(newlh(2),'position',pos);

return
