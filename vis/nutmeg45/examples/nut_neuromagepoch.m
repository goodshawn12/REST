function nuts = nut_neuromagepoch(nuts, codenumber,t)
jj = find(nuts.meg.markersnmag.codes==11);

timeind = t(1)*nuts.meg.srate:t(2)*nuts.meg.srate;
timewin = zeros(length(nuts.meg.markersnmag.latencies{jj}),length(timeind));
nuts.meg.data2 = zeros(length(timewin),size(nuts.meg.data,2),length(nuts.meg.markersnmag.latencies{jj}),'single');

for ii=1:length(nuts.meg.markersnmag.latencies{jj})
    timewin(ii,:) = double(nuts.meg.markersnmag.latencies{jj}(ii))+timeind - nuts.meg.srate*nuts.meg.latency(1);
    nuts.meg.data2(:,:,ii) = nuts.meg.data(timewin(ii,:),:);
end

nuts.meg.latency = 1000*(t(1):(1/nuts.meg.srate):t(2))';