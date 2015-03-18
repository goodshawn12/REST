function tsout=nut_smooth(ts,width)
% function tsout=nut_smooth(ts,width)
% 
% Intention to replicate smooth.m from curve-fitting-toolbox
% ts: time x locations
% width: number of samples to one side to average in

[nt,nk]=size(ts);
if nt<6
    error('not enough data points for smoothing')
end

tsout=ts;
for ii=(width+1):(nt-width)
    tsout(ii,:)=mean(ts((ii-width):(ii+width),:),1);
end

