function out = plot_bfield(coil_coord,A,phi)

% draw image showing strength of magnetic field
% code taken from nut_beamforming_gui.m on March 23, 2006
%
% A = (N x 1) mixing matrix
% coil_coord = (M x 3) matrix of coil coordinates
% goodchannels = (1 x N) vector of indices of coil_coord (1 <= integer <= M)
%
% outputs 2-D mapping of sensor coordinates, for, e.g., later manipulation
% of sensor locations

   goodchannels = 1:size(coil_coord,1);

% normalize coil_coord using all channels
max_value1 = max(abs(coil_coord(:,1)));
max_value2 = max(abs(coil_coord(:,2)));
coil_coord(:,2) = coil_coord(:,2)*max_value1/max_value2; % make extent of first and second columns equal
max_value2 = max(abs(coil_coord(:,2)));
coil_coord(:,3) = coil_coord(:,3) - min(coil_coord(:,3)); % make min value 0 and max value 1
coil_coord(:,3) = coil_coord(:,3)/max(coil_coord(:,3));

% convert 3-d coil_coord into 2-d variable, out
out(:,1) = coil_coord(:,1).*exp(-coil_coord(:,3));
out(:,2) = coil_coord(:,2).*exp(-coil_coord(:,3));

% normalize to uniformly fill a circle
out(:,1) = out(:,1)*max_value1/max(abs(out(:,1)));
out(:,2) = out(:,2)*max_value2/max(abs(out(:,2)));

% imagesc(out(:,1),out(:,2),A)


% adjust for size of border
border = 4;
resolution = 100;
out = out*resolution/(resolution + border);

% initialize
L = resolution + 2*border;
radius = 0.95*(L+1)/2; % radius of circle that represents head
image_matrix = zeros(L); % create empty image

% colormap
% first row should be [0 0 0] (black) and second row should be [1 1 1], white
c = zeros(4096+2,3); c(2,:) = ones(1,3);
c((1:2048)+2,1) = (1:-1/2047:0)'; c((2049:4096)+2,3) = (0:1/2047:1)';
cmap_length = size(c,1);

% row indices where sensors lie (in 2-d image space)
row = round((resolution-1)*(out(:,1)+max_value1)/(2*max_value1)) + 1 + border;
col = round((resolution-1)*(out(:,2)+max_value2)/(2*max_value2)) + 1 + border;

% ensure sensors lie in circle (that represents the head) - fix needed for BTI sensors
sensor_radii = sqrt(((L+1)/2-row).^2 + ((L+1)/2-col).^2);
if max(sensor_radii) > radius
   out = out*0.95*radius/max(sensor_radii);
   row = round((resolution-1)*(out(:,1)+max_value1)/(2*max_value1)) + 1 + border;
   col = round((resolution-1)*(out(:,2)+max_value2)/(2*max_value2)) + 1 + border;
end

% define indices
sensor_ndx_good = (col(goodchannels)-1).*(size(image_matrix,1)) + row(goodchannels);
sensor_ndx = (col-1).*(size(image_matrix,1)) + row;

if size(A,1) >= 2
   Arange = max(max(A)) - min(min(A));
   maxvalue = max(max(A)) - Arange*0.1;
   minvalue = min(min(A)) + Arange*0.1;
   
   % create image (1 is black, 2 is white, the rest are defined in the colormap above)
   meandata = A(:,1); % ignore other columns
   f = find(meandata < minvalue); meandata(f) = minvalue; % introduce user-specified clipping
   f = find(meandata > maxvalue); meandata(f) = maxvalue;
   f = find(meandata < 0); f2 = find(meandata >= 0);
   meandata(f) = meandata(f)/(2*abs(minvalue)); % -0.5 to 0 are for negative A values
   meandata(f2) = meandata(f2)/(2*abs(maxvalue)); % 0 to 0.5 are for positive A values
   image_matrix(sensor_ndx_good) = A(:,1); % output takes values between -0.5 and 0.5
   image_matrix = fliplr(flipud(image_matrix));

   if(1)
   % smooth image
   maxnum = max(max(image_matrix));
   minnum = min(min(image_matrix));
   b = exp(-[-2:0.1:2].^2)'*exp(-[-2:0.1:2].^2);
   image_matrix = filter2(b,image_matrix);
   f = find(image_matrix < 0); f2 = find(image_matrix > 0);
   minnum_new = min(min(image_matrix));
   if minnum_new == 0, minnum_new = eps; end
   image_matrix(f) = image_matrix(f)*abs(minnum/minnum_new); % restore previous minimum
   image_matrix(f2) = image_matrix(f2)*maxnum/max(max(image_matrix)); % restore previous maximum
   end
end

% compute x and y indices of image
xstep = 2*max_value1/(size(image_matrix,1) - 2*border - 1);
ystep = 2*max_value2/(size(image_matrix,2) - 2*border - 1);
xndx = (-max_value1 - 2*xstep):xstep:(max_value1 + 2*xstep);
yndx = (-max_value2 - 2*ystep):ystep:(max_value2 + 2*ystep);
   
% adjust values of amplitude map
%image_matrix = image_matrix + 0.5; % output takes values between 0 and 1

% adjust values for size of colormap
%image_matrix = round((cmap_length-3)*image_matrix) + 3; % avoid first 2 rows!

% find indices of region lying outside of circle
f2 = [];
for i = 1:L
   f = find(((1:L)-(L+1)/2).^2 + (i-(L+1)/2)^2 > radius^2);
   f2 = [f2 f+L*(i-1)];
end

% remove portion outside of circle
outsidecircle = f2;
% image_matrix(f2) = NaN; % makes area outside of circle white


%lh = findobj('name','Plot B Field');
%if isempty(lh)
%    lh = figure;
%    set(lh,'name','Plot B Field')
%    pos = get(lh,'position');
%    set(lh,'position',[pos(1:2) 460 420])
%end
% set(lh,'Colormap',c);
% figure(lh); clf

if(min(image_matrix(:)) < 0)
    clim = [-max(abs(image_matrix(:))) max(abs(image_matrix(:)))];
else
    clim = [0 max(abs(image_matrix(:)))];
end

%bf_image = imagesc(xndx,yndx,image_matrix,[-clim clim]); colorbar
%bf_image = imagesc(xndx,yndx,image_matrix,[0 clim]); colorbar
%bf_image = imagesc(xndx,yndx,image_matrix); colorbar
bf_image = imagesc(xndx,yndx,image_matrix,clim); colorbar




set(gca,'XTick',[],'XTickLabel','','YTick',[],'YTickLabel','');
if(0)
a = axis; text(a(2)+0.04*(a(2)-a(1)),a(3)+0.44*(a(4)-a(3)),' Right','Rotation',270);
title('Anterior'); ylabel('  Left'); xlabel('Posterior'); % title and axis labels must occur after text command above
end

if(1)
hold on
if(exist('phi','var'))
    quiver(-out(:,2),-out(:,1),phi(:,1),phi(:,2),0.4,'w')
else
    plot(-out(goodchannels,2),-out(goodchannels,1),'.k');
    badchannels = 1:size(coil_coord,1); badchannels(goodchannels) = [];
    plot(-out(badchannels,2),-out(badchannels,1),'.','Color',[0.2 0.2 0.2]);
end
end


