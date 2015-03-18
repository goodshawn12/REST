function thresh = nut_find_thresh(vol)
% THRESH = NUT_FIND_THRESH(VOL)
%
% This function findd the threshold for the given 3d image.
%

h = waitbar(0,'Please wait, finding optimal threshold...');
for k=1:size(vol,3)
	% using sobel filter for edge detection
	[BW,thresh(k)] = edge(vol(:,:,k),'sobel');
	waitbar(k/(size(vol,3)),h);
end
%finding the maximum threshold  
thresh=max(thresh);
close(h);
return
