function nut_import_neuromag(megdata_filename);
% imports *.fif Neuromag data
% and optionally events *eve.fif file

global nuts ndefaults

mnecode=which('fiff_setup_read_raw');
if isempty(mnecode)
    error('Please place the MNE Matlab software in your Matlab path. If you need to download it, please see http://nmr.mgh.harvard.edu/martinos/userInfo/data/MNE_register/index.php');
end

% testing if raw or averaged fiff data
[ fid, tree ] = fiff_open(megdata_filename);
[ info, meas ] = fiff_read_meas_info(fid,tree);
FIFF = fiff_define_constants();
raw = fiff_dir_tree_find(meas,FIFF.FIFFB_RAW_DATA);

if isempty(raw)  % for averaged data *.fif file
    avdata=fiff_read_evoked_all(megdata_filename);
    if length(avdata.evoked)>1
        whichepoch = input(['Which epoch to use of the ' num2str(length(avdata.evoked)) ' available?  ']);
    end
    megtmp=avdata.evoked(whichepoch).epochs';
    latency=avdata.evoked(whichepoch).times;
    [info]=fiff_read_meas_info(megdata_filename);
    %             sensor_labels=info.ch_names;
    jj=1;
    for ii=1:length(info.chs)
        if info.chs(ii).kind==1 % MEG channel type
            coil_coord(jj,:)=info.chs(ii).loc;
            meg(:,jj)=megtmp(:,ii);
            nuts.meg.sensor_labels{jj}=info.chs(ii).ch_name;
            jj=jj+1;
        elseif info.chs(ii).kind==2 % EEG channel type
            coil_coord(jj,:)=info.chs(ii).loc;
            ind(jj)=ii;
            %                     meg(:,jj)=megtmp(:,ii);
            nuts.meg.sensor_labels{jj}=info.chs(ii).ch_name;
            jj=jj+1;
            nuts.meg.type{jj}='eeg';
        end
    end

else % for raw *.fif data
    blah=fiff_setup_read_raw(megdata_filename);
    [meg,latency] = fiff_read_raw_segment(blah);
    nuts.meg.latency = latency';   % transpose for ease of use later

    jj=1;
    for ii=1:length(info.chs)
        if info.chs(ii).kind==1 % MEG channel type
            coil_coord(jj,:)=info.chs(ii).loc;
            ind(jj)=ii;
            %                     meg(:,jj)=megtmp(:,ii);
            nuts.meg.sensor_labels{jj}=info.chs(ii).ch_name;
            jj=jj+1;
            nuts.meg.type{jj}='meg';
        elseif info.chs(ii).kind==2 % EEG channel type
            coil_coord(jj,:)=info.chs(ii).loc;
            ind(jj)=ii;
            %                     meg(:,jj)=megtmp(:,ii);
            nuts.meg.sensor_labels{jj}=info.chs(ii).ch_name;
            jj=jj+1;
            nuts.meg.type{jj}='eeg';            
        end
    end
%     meg=meg(ind,:)';
    kill = setdiff(1:length(info.chs),ind);
    meg(kill,:)=[];
    meg = meg';

end

if size(nuts.meg.sensor_labels,1) > size(nuts.meg.sensor_labels,2)
    nuts.meg.sensor_labels = nuts.meg.sensor_labels'; % must be a row vector
end

nuts.meg.sensorOrient = nan(size(coil_coord,1),1);  % for now this will be meaningless
% coil_coord2
% coil_coord2=zeros(size(coil_coord));
nuts.meg.srate = info.sfreq;

if ndefaults.meg.single
    nuts.meg.data = meg;
else
    %     nuts.meg.data=double(nuts.meg.data);
    nuts.meg.data=double(meg);
end
nuts.meg.sensorCoord(:,:,1) = coil_coord;
% if(exist('coil_coord2','var'))
%     nuts.meg.sensorCoord(:,:,2) = coil_coord2;
% end

nuts.meg.goodchannels = 1:size(nuts.meg.data,2);

%% events files loading
[pathname,filebase,fileext]=fileparts(megdata_filename);
eventfile = fullfile(pathname,[filebase '-eve' fileext]);

if exist(eventfile,'file')
    markerdelay = 0.038 * nuts.meg.srate   %%%% NOTE: this delay is only for NeuroSpin's projection system as of May 28, 2009!!!!
    disp('markerdelay correct for NeuroSpin lab only');
    eventlist = mne_read_events(eventfile);
    markerlist = unique(eventlist(:,3));
    for ii=1:length(markerlist)
        markerlatencies{ii}=eventlist(find(eventlist(:,3)==markerlist(ii)),1) + markerdelay;
    end
    nuts.meg.markersnmag.latencies = markerlatencies;
    nuts.meg.markersnmag.codes = markerlist;
end
% FIXME this should be called here
if(0)
    nuts=nut_neuromagepoch(nuts,codenumber,t);
end

%% Do we need to deal with lsc here?  is that the first three
% coords of info.chs(*).loc?
% nuts.meg.lsc = lsc;



%%% lsc now specified in coregistration for CTF, we still do this line up
%%% in BTi and KIT code though.
% nuts.meg.lsc = lsc;
% if ~exist('sensor_labels','var')
%     leftchans = find(coil_coord(:,2) >= 0);
%     rightchans = find(coil_coord(:,2) < 0);
%     for i = 1:length(leftchans)
%         channel_string = int2str(i);
%         sensor_labels{leftchans(i)} = ['L' '0'*ones(1,4-length(channel_string)) channel_string];
%     end
%     for i = 1:length(rightchans)
%         channel_string = int2str(i);
%         sensor_labels{rightchans(i)} = ['R' '0'*ones(1,4-length(channel_string)) channel_string];
%     end
% end


% needed for CTF data only
%     nuts.meg.Gcoef = 0;
%     nuts.meg.grad_order = 0;



