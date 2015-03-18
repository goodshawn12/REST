% rebalances colorscale on plots such that center is at 0
clim=caxis;
caxis([-1 1]*max(abs(clim)));