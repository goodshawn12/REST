function nut_import_4D_BTi(megdata_filename)

% requires MSI>>Matlab functions from Eugene Kronberg
% available from http://biomag.wikidot.com/msi-matlab
if(isempty(which('pdf4D')))
    error('To read 4D/BTi data, please install the MSI>>Matlab functions, available from http://biomag.wikidot.com/msi-matlab');
end

global nuts 
pdf = pdf4D(megdata_filename);
hdr = pdf.header;

% convert these header values to double, since the math/indexing screws up
% on ints, etc.
hdr.header_data.total_epochs = double(hdr.header_data.total_epochs);
hdr.header_data.sample_period = double(hdr.header_data.sample_period);
hdr.epoch_data{1}.pts_in_epoch = double(hdr.epoch_data{1}.pts_in_epoch);


numchans = hdr.header_data.total_chans;

cht=channel_type(pdf,1:numchans);
megchans = strmatch('meg',[cht{:}]);
nuts.meg.type = cht(megchans);
nuts.meg.sensor_labels = channel_label(pdf,megchans);

chp = channel_position(pdf,megchans);

nuts.meg.lsc = 10*s_fit(pdf)';
switch(size([chp(1).position],2))
    case 1
        nuts.meg.sensorCoord(:,:,1) = 1000*[chp(:).position]';
        nuts.meg.sensorOrient(:,:,1) = [chp(:).direction]';
    case 2 % coordinate pairs only for gradiometer systems
        sensorCoord = 1000*[chp(:).position]';
        sensorOrient = 1000*[chp(:).direction]';
        nuts.meg.sensorCoord(:,:,1) = sensorCoord(:,1:2:end);
        nuts.meg.sensorCoord(:,:,2) = sensorCoord(:,2:2:end);
        nuts.meg.sensorOrient(:,:,1) = sensorOrient(:,1:2:end);
        nuts.meg.sensorOrient(:,:,2) = sensorOrient(:,2:2:end);
end


if(0) % event trigger info in m4d file
% MSI.TotalEvents: 1
% MSI.Events: 1
% MSI.EventCodes: 0
% MSI.TrigEventCount: 312
% MSI.TrigEvents: 6089,7932,9379,10826,12443,13777,15608,17078,18321,20096,21894,23567,25409,27252,28688,30067,31887,33266,35120,36499,37879,39179,40524,42378,43700,44955,46674,47985,49421,50981,52710,54078,55842,57300,59052,60477,62172,63834,65202,66841,68684,70515,72132,73941,75625,77230,79028,80871,82137,83685,85347,86794,88287,89575,91226,92978,94730,96290,97534,98789,100677,102214,103808,105605,106849,108624,110093,111552,113270,114717,116515,117815,119250,121025,122314,123557,125389,127096,128792,130103,131833,133438,135292,136682,137926,139294,140888,142402,143691,145240,146495,148179,149886,151288,152905,154431,156093,157754,159507,160942,162762,164560,165871,167669,176136,177402,179030,180805,182376,184027,185734,187102,188888,190504,192189,193964,195863,197231,198632,200170,201719,203460,205291,207055,208716,210254,211814,213623,215284,216618,217998,219592,221344,222723,224509,225877,227200,229042,230342,231755,233044,234819,236504,238357,240053,241466,243117,244485,246384,247955,249459,251313,252986,254331,256083,257836,259350,260843,262685,264291,265964,267388,269084,270983,272215,273515,275222,276737,278376,279982,281677,283554,284888,286414,287692,289127,290970,292304,293627,295096,296826,298442,300262,301551,302908,304524,306118,307837,309668,310979,312325,313715,315298,317129,318565,320374,322194,323652,325212,326693,328027,329587,330989,332843,334561,336426,337953,339196,347449,349223,350964,352581,354073,355825,357080,358878,360291,362156,364010,365627,367040,368928,370522,372376,373755,375111,376841,378424,379927,381691,383432,384867,386574,388282,389525,390780,392182,393889,395720,397258,398659,399937,401644,403238,404719,406618,408517,410088,411581,413129,414972,416849,418454,419958,421653,423179,424626,426458,428244,430098,431443,433195,434462,435863,437412,438814,440238,441911,443178,444715,446026,447666,448954,450367,451667,453273,454776,456189,458032,459886,461491,462814,464555,466002,467732,469179,471078,472502,474119,475566,476888,478426,479907,481772,483683,485288,486656,488114,489550,491370,492806,494501,496050,497904,499181,500945,502844,504596,506021,507807,509525,510995
% MSI.TrigEventCodes: 4206,4216,4196,4216,4216,4206,4216,4206,4216,4206,4216,4196,4296,4196,4206,4216,4216,4196,4206,4196,4216,4196,4216,4296,4216,4196,4196,4196,4196,4196,4196,4206,4206,4196,4196,4196,4216,4206,4196,4196,4206,4206,4206,4206,4216,4216,4206,4296,4206,4206,4196,4196,4206,4196,4206,4216,4196,4196,4206,4206,4216,4196,4206,4206,4206,4216,4196,4196,4216,4216,4196,4216,4196,4216,4216,4206,4196,4206,4206,4196,4196,4206,4206,4196,4196,4216,4206,4196,4216,4206,4206,4216,4206,4216,4296,4216,4206,4206,4196,4216,4216,4216,4206,4206,4216,4196,4296,4196,4206,4206,4216,4216,4206,4216,4206,4216,4196,4206,4196,4216,4196,4196,4216,4206,4196,4196,4196,4206,4196,4206,4216,4206,4216,4216,4216,4206,4196,4296,4196,4206,4206,4216,4216,4216,4196,4206,4206,4216,4196,4206,4196,4216,4196,4206,4216,4206,4206,4206,4296,4196,4196,4216,4196,4216,4206,4206,4216,4206,4216,4196,4206,4206,4206,4206,4216,4296,4196,4216,4206,4206,4216,4216,4196,4206,4216,4196,4196,4206,4216,4196,4196,4216,4216,4216,4206,4216,4216,4206,4216,4196,4206,4206,4216,4196,4206,4196,4206,4216,4206,4196,4196,4196,4206,4196,4206,4196,4196,4196,4206,4296,4196,4216,4196,4206,4216,4206,4206,4196,4206,4216,4216,4206,4196,4216,4206,4216,4196,4216,4206,4206,4216,4216,4196,4196,4216,4216,4216,4216,4216,4216,4206,4206,4216,4216,4216,4196,4206,4216,4206,4296,4196,4196,4206,4196,4196,4216,4196,4196,4196,4196,4196,4216,4216,4196,4196,4206,4196,4196,4216,4216,4196,4206,4196,4296,4196,4196,4216,4206,4206,4196,4216,4196,4216,4196,4206,4206,4206,4216,4216,4216,4296,4216,4216,4206,4206,4216,4216,4196,4196,4206,4216,4206

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


nuts.meg.srate = 1/hdr.header_data.sample_period;
nuts.meg.latency = 1000*((0:hdr.epoch_data{1}.pts_in_epoch-1)*hdr.header_data.sample_period + ind2lat(pdf,1))';   % transpose for ease of use later
nuts.meg.data = zeros(size(nuts.meg.latency,1),size(megchans,1),hdr.header_data.total_epochs,'single');


for ii=1:hdr.header_data.total_epochs;
    % trigger around chosen event code???
    %lat = lat2ind(pdf,ii,[nuts.meg.latency(1)-1 nuts.meg.latency(end)+1]/1000);
    lat = [1 hdr.epoch_data{1}.pts_in_epoch]+(ii-1)*hdr.epoch_data{1}.pts_in_epoch;
    nuts.meg.data(:,:,ii) = read_data_block(pdf,lat,megchans)'; %pdf4D code
    nut_progress(ii,hdr.header_data.total_epochs,2);
end



if(exist('Gcoef','var'))  % store the good bits for synthetic gradient/noise correction
    error('not ready yet -- awaiting details of 4D/BTi noise correction')
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

% for magnetometer systems
disp('NOTE: assuming all sensors are magnetometers! (Do you have gradiometers? Please contact Sarang Dalal with details about your system.');
nuts.meg.chanmixMtx{1}=eye(size(megchans,1));
