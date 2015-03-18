function nut_view_beamforming_activation(time, whichblob, threshold, addblob)
% NUT_VIEW_BEAMFORMING_ACTIVATION(TIME, WHICHBLOB, THRESHOLD)
%
% Display pre-computed beamformer activations in the SPM MR display window

global rivets beam st;
if nargin<4
    addblob=false;
end

% retrieve slice -- if time is a vector, average over given time points
% FIXME: need workaround for kludgy findobj method
handles = guihandles(rivets.fig);

if ~isfield(rivets,'scalemax')
    rivets.scalemax=rivets.maxblob*rivets.scalefactor;
    rivets.scalemin=rivets.minblob*rivets.scalefactor;
end

if isfield(handles,'nut_orientations_menu')         % if called by nut_results_viewer
    switch get(handles.nut_orientations_menu,'Value')
        case 1
            timeslice = rivets.scalefactor*squeeze(mean(rivets.s_beam(:,time),2));
        case 2
            timeslice = rivets.scalefactor*squeeze(mean(rivets.s_II(:,time),2));
        case 3
            timeslice = rivets.scalefactor*squeeze(mean(beam.s{rivets.tsflag}(:,time),2));
        case {4,5,6,7}
            timeslice = rivets.s_beam(:,time);
    end
    % if(get(handles.nut_statsenabled_box,'Value'))
    %     timesliceII = rivets.scalefactor*squeeze(mean(rivets.s_II(:,time),2));
    %     select = abs(timesliceII)>rivets.thresh;
    % else
    select = (timeslice > threshold(:,1)) | (timeslice < threshold(:,2));
    % end
elseif isfield(handles,'nut_orientated_menu')  % if called by nut_results_viewer
    if rivets.timefmode
        timeslice = rivets.scalefactor*squeeze(mean(rivets.s(:,time,rivets.freqselect),2));
        % I introduced new format of rivets.threshold in nut_timef_viewer on 6/12/07 to
        % avoid problems with statistical thresholds and peak finding functions. (Adrian)
        % Now, threshold is a logical array (voxels*time*[pos neg]), true=show.
        select = find(any(squeeze(threshold(:,rivets.timeselect,:)),2));
    else
        orientcontents=get(handles.nut_orientated_menu,'String');
        orienttype=orientcontents{get(handles.nut_orientated_menu,'Value')};
        switch orienttype
            case {'Power All Orientations','Multi Time Series'}
                timeslice = rivets.scalefactor*squeeze(mean(rivets.s(:,time,1),2));
            case {'Amplitude: Constant Orientation','Amplitude: Instant Orientation','Amplitude: Mob Orientation','Amplitude: Original Orientation'}
                timeslice = rivets.scalefactor*squeeze(mean(rivets.s(:,time,1,:),2));
%             case 'Multi Time Series'
%                 timeslice = rivets.scalefactor*squeeze(mean(beam.sa{rivets.tsflag}(:,time),2));
        end
        % if(get(handles.nut_statsenabled_box,'Value'))
        %     timesliceII = rivets.scalefactor*squeeze(mean(rivets.s_II(:,time),2));
        %     select = abs(timesliceII)>rivets.thresh;
        % else
%         select = (timeslice > threshold(:,1)) | (timeslice < threshold(:,2));
        select = find(any(squeeze(threshold(:,rivets.timeselect,:)),2));
        % end
    end
else % if called by nut_timef_viewer
    timeslice = rivets.scalefactor*squeeze(mean(rivets.s(:,time,rivets.freqselect),2));
    % I introduced new format of rivets.threshold in nut_timef_viewer on 6/12/07 to
    % avoid problems with statistical thresholds and peak finding functions. (Adrian)
    % Now, threshold is a logical array (voxels*time*[pos neg]), true=show.
    select = find(any(squeeze(threshold(:,rivets.timeselect,:)),2));
end

if addblob
%     spm_orthviews('rmblobs',1);
    spm_orthviews('addcolouredblobs',1,rivets.voxelsblob(select,:)',timeslice(select),rivets.blob2mri_tfm,[0 0 1]);
%     st.vols{1}.blobs{2}.max = rivets.maxblob*rivets.scalefactor;
%     st.vols{1}.blobs{2}.min = rivets.minblob*rivets.scalefactor;
    st.vols{1}.blobs{2}.max = 1000;
    st.vols{1}.blobs{2}.min = 0;
else
spm_orthviews('rmblobs',1);
if(find(select))  % if voxels are above threshold...
    switch(whichblob)
        case 1
            % "hot" opaque blob, with colorbar
            spm_orthviews('rmblobs',1);
            spm_orthviews('addblobs',1,rivets.voxelsblob(select,:)',timeslice(select),rivets.blob2mri_tfm);
%             if rivets.timefmode
%             st.vols{1}.blobs{1}.max = rivets.maxblob*rivets.scalefactor;
%             st.vols{1}.blobs{1}.min = rivets.minblob*rivets.scalefactor;
%             else
            st.vols{1}.blobs{1}.max = rivets.scalemax;
            st.vols{1}.blobs{1}.min = rivets.scalemin;
%             end

            % we shouldn't have to do this -- even after above hack, we need to force
            % the colorbar to reflect the full range as well
            image([0 1],[st.vols{1}.blobs{1}.min st.vols{1}.blobs{1}.max],[1:64]' + 64,'Parent',st.vols{1}.blobs{1}.cbar);
            set(st.vols{1}.blobs{1}.cbar,'YDir','normal','XTickLabel',[]);
        case 2
            % with last field set to [1 0 0], yields a translucent red blob
            spm_orthviews('rmblobs',1)
            spm_orthviews('addcolouredblobs',1,rivets.voxelsblob(select,:)',timeslice(select),rivets.blob2mri_tfm,[1 0 0]);
%             if rivets.timefmode
%             st.vols{1}.blobs{1}.max = rivets.maxblob*rivets.scalefactor;
%             st.vols{1}.blobs{1}.min = rivets.minblob*rivets.scalefactor;
%             else
            st.vols{1}.blobs{1}.max = rivets.scalemax;
            st.vols{1}.blobs{1}.min = rivets.scalemin;
%             end
        case 3  % second translucent blob
            % with last field set to [1 1 0], yields a translucent yellow blob
            spm_orthviews('rmblobs',1)
            spm_orthviews('addcolouredblobs',1,rivets.voxelsblob(select,:)',timeslice(select),rivets.blob2mri_tfm,[1 0 0]);
            timeslice = rivets.scalefactor*squeeze(mean(beam.s_perp(:,time),2));
            select = abs(timeslice)>threshold;
            spm_orthviews('addcolouredblobs',1,rivets.voxelsblob(select,:)',timeslice(select),rivets.blob2mri_tfm,[0 1 0]);
%             if rivets.timefmode
%             st.vols{1}.blobs{1}.max = rivets.maxblob*rivets.scalefactor;
%             st.vols{1}.blobs{1}.min = rivets.minblob*rivets.scalefactor;
%             st.vols{1}.blobs{2}.max = rivets.maxblob*rivets.scalefactor;
%             st.vols{1}.blobs{2}.min = rivets.minblob*rivets.scalefactor;
%             else
            st.vols{1}.blobs{1}.max = rivets.scalemax;
            st.vols{1}.blobs{1}.min = rivets.scalemin;
            st.vols{1}.blobs{2}.max = rivets.scalemax;
            st.vols{1}.blobs{2}.min = rivets.scalemin;
%             end
        case 4 % difference s_II and s_perp -- this is an elaborate and hopefully temporary hack to difference two s_beams
            % "hot" opaque blob, with colorbar
            spm_orthviews('rmblobs',1);
            timeslice = rivets.scalefactor*(squeeze(mean(rivets.s_II(:,time),2)) - squeeze(mean(beam.s_perp(:,time),2)));
            spm_orthviews('addblobs',1,beam.voxelsblob(select,:)',timeslice(select),rivets.blob2mri_tfm);
            st.vols{1}.blobs{1}.max = max(abs(timeslice));
            st.vols{1}.blobs{1}.min = -max(abs(timeslice));

            % we shouldn't have to do this -- even after above hack, we need to force
            % the colorbar to reflect the full range as well
            image([0 1],[st.vols{1}.blobs{1}.min st.vols{1}.blobs{1}.max],[1:64]' + 64,'Parent',st.vols{1}.blobs{1}.cbar);
            set(st.vols{1}.blobs{1}.cbar,'YDir','normal','XTickLabel',[]);
    end
end
end

warning('off','MATLAB:FINITE:obsoleteFunction'); % disable complaints of ne Matlab versions about SPM2
spm_orthviews('redraw');
