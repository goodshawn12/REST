function nut_import_BTi_UCSF(megdata_filename)

% this is based on UCSF's way of converting BTi data,
% unfort not general to all BTi data

global nuts ndefaults

[sensor_filename, sensor_path]=uigetfile('*.sensor','Select MEG sensor locations...');
if isequal(sensor_filename,0)|isequal(sensor_path,0)
    return;
else
    sensor_filename = [sensor_path sensor_filename];
end

[lsc_filename, lsc_path]=uigetfile('*.A;*.B','Select MEG local sphere center...');
if isequal(lsc_filename,0)|isequal(lsc_path,0)
    return;
else
    lsc_filename = [lsc_path lsc_filename];
end


[events_filename, events_path]=uigetfile('*.*','Events file (press cancel if none)');
if isequal(events_filename,0)|isequal(events_path,0)
    [rawdata srate no_points start_time]=nut_read_BTi(megdata_filename);
else
    events=load([events_path events_filename]);
    prompt   = 'Which events?';
    title    = 'Pick your poison.';
    lines = 1;
    eventcodes = unique(events(:,2)'); % figure out which codes are in file
    eventcodes(eventcodes==99)=[];  % 99s are bad trials, get rid of them
    def{1} = int2str(eventcodes);
    answer   = inputdlg(prompt,title,lines,def);
    if (isempty(answer)) msgbox('MEG loading cancelled. Freak.');return; end;
    eventselect = str2num(answer{1});

    select = find(ismember(events(:,2),eventselect));   % select given codes
    [rawdata srate no_points start_time]=nut_read_BTi(megdata_filename,select);
end

switch lsc_filename(end)
    case 'A'
        kept_channels = 1:37;
        %badchannels = 28;   % channel 28 is busted, yo
        readsensors = 37;
    case 'B'
        kept_channels = (size(rawdata,1)-37):(size(rawdata,1)-1);
        %badchannels = [27 34];
        readsensors = 74;
    otherwise
        error('you''re nothing but a dirty junkie.');
end

[a,b,c,d,e,f,g,h,i,j,k,l,m]=...
    textread(sensor_filename,'%f%f%f%n%n%n%n%n%n%n%n%s%s',readsensors,...
    'headerlines',1);
coil_coord = 1000*[a b c];   % .sensors file has coords in meter, we want mm
normal_direction = [d e f];
if ~ndefaults.meg.btitype
    sensor_labels=l;  %nice to see sensor label (i.e. 'A24') if only 37 channels,
    %versus left/right labels
end

switch lsc_filename(end)
    case 'A'
        coil_coord = coil_coord(1:37,:);
        normal_direction = normal_direction(1:37,:);
    case 'B'
        coil_coord = coil_coord(38:74,:);
        normal_direction = normal_direction(38:74,:);
    otherwise
        error('you''re nothing but a dirty junkie.');
end

% location of outer magnetometer in gradiometer pair is 50.5 mm away
coil_coord2=coil_coord+normal_direction*50.5;

lsc = load(lsc_filename);  % given in meters
lsc = 1000*lsc; % we want millimeters
latency = 1000*((0:(no_points-1))/srate + start_time);  % in ms

% we have enough memory (and reasons) to leave data unaveraged at this point
% meg = squeeze(mean(rawdata,3))';   % average into ERPs for each channel

% meg = rawdata(:,kept_channels,:);  
meg = permute(rawdata(kept_channels,:,:),[2 1 3]);  % discard unkept channels
clear rawdata


%%%%%% EXPERIMENTAL CODE: uses individual trials to compute covariance matrix
covmethod = false;
if(covmethod)
    disp('Using trial-by-trial covariance computation (EXPERIMENTAL!!!)');
    rawdata = rawdata(kept_channels,:,:);
    %rawdata(badchannels,:,:) = [];
    dc = mean(rawdata(:,1:zeropoint,:),2);
    rawdata = rawdata - repmat(dc,[1 size(rawdata,2) 1]);

    save('megcovar.mat','rawdata');
end
%%%%% END covariance experiment


if ndefaults.meg.single
    nuts.meg.data = meg;
else
    %     nuts.meg.data=double(nuts.meg.data);
    nuts.meg.data=double(meg);
end

%%% lsc now specified in coregistration for CTF, we still do this line 
%%% in BTi and KIT code though.
nuts.meg.lsc = lsc;
nuts.meg.sensorCoord(:,:,1) = coil_coord;
if(exist('coil_coord2','var'))
    nuts.meg.sensorCoord(:,:,2) = coil_coord2;
end
if ~exist('sensor_labels','var')
    leftchans = find(coil_coord(:,2) >= 0);
    rightchans = find(coil_coord(:,2) < 0);
    for i = 1:length(leftchans)
        channel_string = int2str(i);
        sensor_labels{leftchans(i)} = ['L' '0'*ones(1,4-length(channel_string)) channel_string];
    end
    for i = 1:length(rightchans)
        channel_string = int2str(i);
        sensor_labels{rightchans(i)} = ['R' '0'*ones(1,4-length(channel_string)) channel_string];
    end
end
if size(sensor_labels,1) > size(sensor_labels,2)
    sensor_labels = sensor_labels'; % must be a row vector
end
nuts.meg.sensor_labels = sensor_labels;
nuts.meg.type=repmat({'meg'},size(nuts.meg.data,2),1);

nuts.meg.sensorOrient = normal_direction;

if(exist('Gcoef','var'))  % store the good bits for synthetic gradient correction
    nuts.meg.Gcoef = Gcoef;
    nuts.meg.grad_order = grad_order;
    nuts.meg.refSensorCoord(:,:,1) = ref_coord;
    nuts.meg.refSensorCoord(:,:,2) = ref_coord2;

    % refGradiometerChans = find(all(ref_coord2,2)); % find reference coils that are gradiometers, i.e., ref_coord2(chan,:) ~= [0 0 0])
    refMagnetometerChans = find(~any(ref_coord2,2)); % find reference coils that are magnetometers, i.e., ref_coord2(chan,:) == [0 0 0])
    nuts.meg.refSensorCoord(refMagnetometerChans,:,2) = NaN;  % magnetometers don't have second coil, so set to NaN

    nuts.meg.refSensorOrient = ref_direction;
    nuts.meg.ref_sensor_labels = ref_sensor_labels;
else
    nuts.meg.Gcoef = 0;
    nuts.meg.grad_order = 0;
end

nuts.meg.goodchannels = 1:size(nuts.meg.data,2);
nuts.meg.srate = srate;
nuts.meg.latency = latency';   % transpose for ease of use later
