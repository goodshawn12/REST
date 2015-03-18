function nut_plotbeamtf(MEGvoxelindex)
global beam rivets ndefaults
handles = guidata(rivets.fig);
rivets.tf_ax = handles.nut_ts_axes;

axes(rivets.tf_ax);
tmp=unique(beam.bands(:));
set(rivets.tf_ax,'YTick',tmp);
if ~isinteger(tmp)
    tmp=round(tmp*10)/10;
end
set(rivets.tf_ax,'YTickLabel',num2str(tmp));
if get(handles.nut_check_logscale,'Value')
    set(rivets.tf_ax,'YScale','log')           %,'YTick',[4 12 30 60 120 180 220]);
else
    set(rivets.tf_ax,'YScale','linear')
end
% set(rivets.tf_ax,'YDir','normal');
tfh = zeros(rivets.timeselect,rivets.freqselect);
% if isfield(beam,'sa')
%     if(size(beam.sa{1},1)==1)
%         beam.sa{1} = beam.sa{1}';
%     end
% else
if(size(beam.s{1},1)==1)
    beam.s{1} = beam.s{1}';
end
% end
cla(rivets.tf_ax);
    
for timebin=1:size(beam.timewindow,1)
    for freqbin=1:size(beam.bands,1)
        tfh(timebin,freqbin) = patch([beam.timewindow(timebin,1) beam.timewindow(timebin,1) beam.timewindow(timebin,2) beam.timewindow(timebin,2)],[beam.bands(freqbin,[1 2]) beam.bands(freqbin,[2 1])],rivets.s(MEGvoxelindex,timebin,freqbin),'Parent',rivets.tf_ax,'HitTest','off');
        if(~ndefaults.tfbf.border)
            set(tfh(timebin,freqbin),'LineStyle','none');
        end
    end
end

time = str2double(get(handles.nut_time_text,'String'));
rivets.timeselect = dsearchn(beam.timepts,time);

%% hack to bring selected patch to foreground with complete red box
delete(tfh(rivets.timeselect,rivets.freqselect));
for timebin=rivets.timeselect
    for freqbin=rivets.freqselect
        tfh(timebin,freqbin) = patch([beam.timewindow(timebin,1) beam.timewindow(timebin,1) beam.timewindow(timebin,2) beam.timewindow(timebin,2)],[beam.bands(freqbin,[1 2]) beam.bands(freqbin,[2 1])],rivets.s(MEGvoxelindex,timebin,freqbin),'Parent',rivets.tf_ax,'HitTest','off');
    end
end
set(tfh(rivets.timeselect,rivets.freqselect),'EdgeColor',[1 0 0],'LineStyle','-','LineWidth',2);
%%

axis(rivets.tf_ax,[min(beam.timewindow(:)) max(beam.timewindow(:)) 0 max(beam.bands(:))]);
caxis(rivets.tf_ax,[rivets.scalemin rivets.scalemax]);
rivets.cbar = colorbar('peer',rivets.tf_ax);
if ( isfield(beam,'labels') && isfield(beam.labels,'colorbar') )
    cblab = beam.labels.colorbar;
else
    cblab = 'dB';
end
% dbH = title(rivets.cbar,cblab,'FontSize',14,'FontWeight','bold');
set(rivets.tf_ax,'YDir','normal');
set(rivets.tf_ax,'ylim',[beam.bands(1,1) beam.bands(end,end)])
drawnow;        % This avoids stupid matlab bug in version 2003a
axes(rivets.tf_ax); 
if isfield(beam,'labels')
    xlabel(beam.labels.xaxis)
    ylabel(beam.labels.yaxis)
else
    xlabel('Time (ms)'); 
    ylabel('Frequency (Hz)');
end
