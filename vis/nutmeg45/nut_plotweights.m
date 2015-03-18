function nut_plotweights(weights)
% draw image showing strength of magnetic field

global nuts bolts beam rivets

figh = rivets.W_fig;
c = zeros(4096+2,3); c(2,:) = ones(1,3);
c((1:2048)+2,1) = (1:-1/2047:0)'; c((2049:4096)+2,3) = (0:1/2047:1)';
set(figh,'Colormap',c);

coil_coord = nuts.meg.sensorCoord(nuts.meg.goodchannels,:,1);
bolts.plot_bfield.max_value1 = max(abs(coil_coord(:,1)));
bolts.plot_bfield.max_value2 = max(abs(coil_coord(:,2)));
coil_coord(:,2) = coil_coord(:,2)*bolts.plot_bfield.max_value1/bolts.plot_bfield.max_value2; % make extent of first and second columns equal
bolts.plot_bfield.max_value2 = max(abs(coil_coord(:,2)));
coil_coord(:,3) = coil_coord(:,3) - min(coil_coord(:,3)); % make min value 0 and max value 1
coil_coord(:,3) = coil_coord(:,3)/max(coil_coord(:,3));

% convert 3-d coil_coord into 2-d variable, out
out(:,1) = coil_coord(:,1).*exp(-coil_coord(:,3));
out(:,2) = coil_coord(:,2).*exp(-coil_coord(:,3));

% normalize to uniformly fill a circle
out(:,1) = out(:,1)*bolts.plot_bfield.max_value1/max(abs(out(:,1)));
out(:,2) = out(:,2)*bolts.plot_bfield.max_value2/max(abs(out(:,2)));

% adjust for size of border
bolts.plot_bfield.border = 4;
bolts.plot_bfield.resolution = 100;
bolts.plot_bfield.out = out*bolts.plot_bfield.resolution/(bolts.plot_bfield.resolution + bolts.plot_bfield.border);


% initialize
L = bolts.plot_bfield.resolution + 2*bolts.plot_bfield.border;
radius = 0.95*(L+1)/2; % radius of circle that represents head
image_matrix = zeros(L); % create empty image
cmap_length = size(get(figh,'colormap'),1);

% row indices where sensors lie (in 2-d image space)
row = round((bolts.plot_bfield.resolution-1)*(bolts.plot_bfield.out(:,1)+bolts.plot_bfield.max_value1)/(2*bolts.plot_bfield.max_value1)) + 1 + bolts.plot_bfield.border;
col = round((bolts.plot_bfield.resolution-1)*(bolts.plot_bfield.out(:,2)+bolts.plot_bfield.max_value2)/(2*bolts.plot_bfield.max_value2)) + 1 + bolts.plot_bfield.border;

% ensure sensors lie in circle (that represents the head) - fix needed for BTI sensors
sensor_radii = sqrt(((L+1)/2-row).^2 + ((L+1)/2-col).^2);
if max(sensor_radii) > radius
   bolts.plot_bfield.out = bolts.plot_bfield.out*0.95*radius/max(sensor_radii);
   row = round((bolts.plot_bfield.resolution-1)*(bolts.plot_bfield.out(:,1)+bolts.plot_bfield.max_value1)/(2*bolts.plot_bfield.max_value1)) + 1 + bolts.plot_bfield.border;
   col = round((bolts.plot_bfield.resolution-1)*(bolts.plot_bfield.out(:,2)+bolts.plot_bfield.max_value2)/(2*bolts.plot_bfield.max_value2)) + 1 + bolts.plot_bfield.border;
end

% define indices
% sensor_ndx_good = (col(nuts.meg.goodchannels)-1).*(size(image_matrix,1)) + row(nuts.meg.goodchannels);
sensor_ndx = (col-1).*(size(image_matrix,1)) + row;

% if size(meg,2) >= 2
%    minvalue = min(weights);
   maxvalue = max(abs(weights));

   if(0)
      f = find(weights < minvalue); normdata(f) = minvalue; % introduce user-specified clipping
      f = find(weights > maxvalue); normdata(f) = maxvalue;
   end
   normdata = weights;
      f = find(normdata < 0); f2 = find(weights >= 0);
      
      normdata(f) = normdata(f)/(2*abs(maxvalue)); % -0.5 to 0 are for negative meg values
      normdata(f2) = normdata(f2)/(2*abs(maxvalue)); % 0 to 0.5 are for positive meg values
      image_matrix(sensor_ndx) = normdata; % output takes values between -0.5 and 0.5

      image_matrix = fliplr(flipud(image_matrix));

   % smooth image
   maxnum = max(max(image_matrix));
   minnum = min(min(image_matrix));
   b = exp(-[-.5:0.05:.5].^2)'*exp(-[-.5:0.05:.5].^2);
   image_matrix = filter2(b,image_matrix);
   f = find(image_matrix < 0); f2 = find(image_matrix > 0);
%    if get(handles.nut_amplitude_figure,'Value') == 1 % amplitude map
      minnum_new = min(min(image_matrix));
      if minnum_new == 0, minnum_new = eps; end
      image_matrix(f) = image_matrix(f)*abs(minnum/minnum_new); % restore previous minimum
%    end
   image_matrix(f2) = image_matrix(f2)*maxnum/max(max(image_matrix)); % restore previous maximum
% end

% compute x and y indices of image
xstep = 2*bolts.plot_bfield.max_value1/(size(image_matrix,1) - 2*bolts.plot_bfield.border - 1);
ystep = 2*bolts.plot_bfield.max_value2/(size(image_matrix,2) - 2*bolts.plot_bfield.border - 1);
xndx = (-bolts.plot_bfield.max_value1 - 2*xstep):xstep:(bolts.plot_bfield.max_value1 + 2*xstep);
yndx = (-bolts.plot_bfield.max_value2 - 2*ystep):ystep:(bolts.plot_bfield.max_value2 + 2*ystep);
   
% adjust values of amplitude map
% if get(handles.nut_amplitude_figure,'Value') == 1 % amplitude map
   image_matrix = image_matrix + 0.5; % output takes values between 0 and 1
% end
   
% adjust values for size of colormap
image_matrix = round((cmap_length-3)*image_matrix) + 3; % avoid first 2 rows!

% find indices of region lying outside of circle
f2 = [];
for i = 1:L
   f = find(((1:L)-(L+1)/2).^2 + (i-(L+1)/2)^2 > radius^2);
   f2 = [f2 f+L*(i-1)];
end

% remove portion outside of circle
image_matrix(f2) = 2; % makes area outside of circle white

%c = colormap;set(handles.nut_beamforming_fig,'Colormap',c);
%lh = get(handles.nut_select_axes,'Children');
%lh = handles.sensor_locations;
%if isempty(lh) || ~strcmp(lower(get(lh(end),'Type')),'image') % if the last line handle is not of type image
handles.bf_image = image(xndx,yndx,image_matrix,'Parent',rivets.W_ax);
% guidata(hObject,handles);
if(0)
    set(handles.nut_select_axes,'XTick',[],'XTickLabel','','YTick',[],'YTickLabel','');
    a = axis; text(a(2)+0.04*(a(2)-a(1)),a(3)+0.44*(a(4)-a(3)),' Right','Rotation',270);
    title('Anterior'); ylabel('  Left'); xlabel('Posterior'); % title and axis labels must occur after text command above
end
%lh = get(handles.nut_select_axes,'Children');

   % plot circles for the sensor locations, if ismember then the circles must be redrawn
% if ismember(hObject,[handles.nut_select_button handles.nut_deselect_button handles.nut_redeem handles.nut_renounce])
%    %delete(lh(1:end-5)); % delete previous circles (the last 2 children are the text "Right" and the image)
   if isfield(handles,'sensor_locations')
      delete(handles.sensor_locations);
   end
   hold(rivets.W_ax,'on');
   
   % fast way, but can't assign unique context menus
   %       h=plot(-bolts.plot_bfield.out	(:,2),-bolts.plot_bfield.out	(:,1),'ow');
   
   % slow way -- gives you context menus
   lh_sensor_locations = zeros(size(nuts.meg.goodchannels,2),1);
   for i=1:size(nuts.meg.goodchannels,2)%(bolts.plot_bfield.out	,1)
      sensorui=uicontextmenu('Parent',rivets.W_fig);
      item1=uimenu(sensorui,'Label',nuts.meg.sensor_labels{nuts.meg.goodchannels(i)},'Enable','off');
      togglecmd = ['nut_beamforming_gui(''nut_toggle_channel'',gcbo,[],guidata(gcbf),' num2str(i) ')'];
%       if(ismember(i,nuts.meg.goodchannels))
         item2=uimenu(sensorui,'Label',num2str(weights(i,:)));
         % this is too slow to make the ButtonDownFcn practical...
         %  plot(-bolts.plot_bfield.out(i,2),-bolts.plot_bfield.out(i,1),'ow','UIContextMenu',sensorui,'ButtonDownFcn',togglecmd);
         lh_sensor_locations(i) = plot(rivets.W_ax,-bolts.plot_bfield.out(i,2),-bolts.plot_bfield.out(i,1),'ow','UIContextMenu',sensorui);
%       else
% %          item2=uimenu(sensorui,'Label','Select','Callback',['nut_beamforming_gui(''nut_toggle_channel'',gcbo,[],guidata(gcbf),' num2str(i) ')'],'Tag','nut_sensor_locations');
%          % this is too slow to make the ButtonDownFcn practical...
%          % plot(-bolts.plot_bfield.out(i,2),-bolts.plot_bfield.out(i,1),'o','Color',[0.2 0.2 0.2],'UIContextMenu',sensorui,'ButtonDownFcn',togglecmd);
%          lh_sensor_locations(i) = plot(-bolts.plot_bfield.out(i,2),-bolts.plot_bfield.out(i,1),'o','Color',[0.2 0.2 0.2],'UIContextMenu',sensorui);
%       end
   end
   
   hold(rivets.W_ax,'off');
   handles.sensor_locations = lh_sensor_locations;
%    guidata(hObject,handles);
% end

return
