function nut_importmegcov(megdata_fullpath,lsc_fullpath);
% NUT_IMPORTMEG
% currently supports CTF, 37/74-channel BTi Magnes, and KIT
% 
% work on input /outpur relations
% do we want to output R? save it in nuts? i dunno, i'll decide tomorrow

global nuts ndefaults;

if exist('megdata_fullpath','var')
    [megdata_path,megdata_filename,megdata_ext]=fileparts(megdata_fullpath);
    megdata_path = [megdata_path '/'];  % to maintain consistency with uigetfile output
    megdata_filename = [megdata_filename megdata_ext];
else
    [megdata_filename, megdata_path]=uigetfile('*.bin;*.meg4;*.txt','Select MEG data file (BTi,CTF,KIT)...');
    if isequal(megdata_filename,0)|isequal(megdata_path,0)
        return;
    end
end

% compare path (sans trailing '/') to .ds -- CTF data if matched
if(strcmp(megdata_path((end-3):(end-1)),'.ds'))
    megsys = 'CTF'
    megfolder = megdata_path(1:(end-1));
elseif(strcmp(megdata_filename((end-2):end),'.ds'))
    megsys = 'CTF'
    megfolder = megdata_fullpath;
else  % otherwise it's something else
    [crap,morecrap,ext] = fileparts(megdata_filename);
    switch(ext)
        case '.bin'
            megsys = 'BTi'
        case '.txt'
            megsys = 'KIT'
        case '.mat'
            megsys = 'matfile'
        otherwise
            errordlg('Sorry, unknown data type.');
    end
    megdata_filename = fullfile(megdata_path,megdata_filename);
end

[crap,megfile,morecrap] = fileparts(megdata_filename);
clear crap morecrap;

%%%% FIXME FIXME FIXME: make following 'cases' into function calls
%%%% [meg, latency, sensors, lsc] = bti_read(path,file) or ctf_read or kit_read;


switch megsys
    case 'BTi',
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

        meg = squeeze(mean(rawdata,3))';   % average into ERPs for each channel
        meg = meg(:,kept_channels);   % discard unkept channels


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
        nuts.meg.lsc = lsc;


    case 'KIT'
        meg=load('-ASCII',megdata_filename);

        % drop last 10 characters from selected filename, so we can construct filenames of parameter files
        kit_basename = megdata_filename(1:(end-10));

        sensor_filename = [kit_basename 'sns.txt'];

        sensor2head = kit_matrix([kit_basename 'matrix.txt']);   % read sensor2head transform matrix
        [lsc_kit r] = kit_sphere([kit_basename 'sphere.txt']);

        sensor2head(4,:) = [0 0 0 1]; % augment transformation matrix to make 4x4


        sensor_filename

        fid=fopen(sensor_filename);

        lineL = fgetl(fid);
        lineL = fgetl(fid);
        lineL = fgetl(fid);

        %  amdmax=192;
        % [channelnum,channeltype,coil_th,coil_ph,coilsize,baseline]=deal(zeros(amdmax,1));
        % coil_coord = zeros(amdmax,3);

        ichannel = 0;
        while 1
            lineL = fgetl(fid);
            if(~ischar(lineL)),break,end;
            [crap,channeltype] = strread(lineL,'%f%s',1);
            if(strcmp(channeltype,'AxialGradioMeter'))
                ichannel = ichannel+1;
                [channelnum(ichannel,1),channeltype(ichannel,1),coil_coord(ichannel,1),coil_coord(ichannel,2),coil_coord(ichannel,3),coil_th(ichannel,1),coil_ph(ichannel,1),coilsize(ichannel,1),baseline(ichannel,1)] = strread(lineL,'%f%s%f%f%f%f%f%f%f',1);
            end
        end

        fclose(fid);

        meg = meg(:,channelnum+1);  % discard bogus channels

        baseline=baseline(1);  % all baseline values should be the same; just keep first one

        coil_th=coil_th*pi/180;  % convert to radians
        coil_ph=coil_ph*pi/180;  % convert to radians
        normal_direction = [sin(coil_th).*cos(coil_ph) sin(coil_th).*sin(coil_ph) cos(coil_th)];
        coil_coord2 = coil_coord + normal_direction.*baseline;

        %%%% convert from KIT sensor coords to head coords -- should already be in millimeters
        coil_coord = nut_coordtfm(coil_coord,sensor2head);
        coil_coord2 = nut_coordtfm(coil_coord2,sensor2head);

        % convert KIT coord convention to BTi/CTF
        % (x_ctf = -y_KIT; y_ctf = x_KIT; z_ctf = z_KIT)
        lsc = [-lsc_kit(2) lsc_kit(1) lsc_kit(3)];
        nuts.meg.lsc=lsc;
        coil_coord = [-coil_coord(:,2) coil_coord(:,1) coil_coord(:,3)];
        coil_coord2 = [-coil_coord2(:,2) coil_coord2(:,1) coil_coord2(:,3)];

        normal_direction = coil_coord2-coil_coord;
        normal_direction = normal_direction./repmat(nut_rownorm(normal_direction),[1 3]);
        % normal_direction3 = nut_coordtfm(normal_direction,[sensor2head(:,1:3) [0;0;0;1]]);


        % WE NEED TO PROMPT USER FOR LATENCY INFO
        prompt   = 'Enter starting latency and sampling rate (e.g., [-100 256] for 256 Hz starting at -100 ms)';
        title    = 'Input latency info';
        lines = 1;
        def{1}   = num2str('-100 1000');

        answer   = inputdlg(prompt,title,lines,def);
        if (isempty(answer)) errordlg('MEG loading cancelled. Freak.');return; end;
        timeinfo = str2num(answer{1});

        start_time = timeinfo(1);
        srate = timeinfo(2);

        no_points = size(meg,1);
        latency = (1000*(0:(no_points-1))/srate) + start_time;  % in ms
    case 'CTF'
        currentdir = pwd;  % remember current directory so we can get back to it later

        [crap, filename, ext] = fileparts(megfolder);
        clear crap;

        ctf = ctf_read_res4(megfolder);
        if ctf.setup.number_trials > 3
            numtrials = ctf.setup.number_trials;
        else
            % if numtrials is 1, 2, or 3, it's 99% likely to be an averaged
            % dataset, in which case only the first "trial" is of importance
            numtrials = 1;
        end

        %%%%% pick channels to read
        kept_channels = ctf.sensor.index.meg_sens;
        %kept_channels=1:length(ctf.sensor.info); % useful for grabbing all channels (not just MEG channels)

        % setting up "meg" might run out of memory, so let's do it before we read in our big-ass file
        % meg = zeros(length(kept_channels), ctf.setup.number_samples, length(trials));  % channels x time x trials

        [tempcell{1:length(kept_channels)}]=deal(ctf.sensor.info(kept_channels).grad_order_no);
        grad_order = max(cell2mat(tempcell));
        clear tempcell;
        
        R = zeros(length(kept_channels));  % set up covariance matrix

        if(grad_order > 0)
            disp('Synthetic gradient correction detected.');
            sprintf('Order = %d',grad_order)

            [tempcell{1:length(ctf.sensor.index.meg_ref)}]=deal(ctf.sensor.info(ctf.sensor.index.meg_ref).location);  % in cm
            refLocations = [tempcell{:}]*10;  % convert to millimeters for compatibility
            ref_coord = refLocations(:,1:2:end)';  % odd columns (closer sensor in gradiometer pair)
            ref_coord2 = refLocations(:,2:2:end)';  % even columns (farther sensors in gradiometer pair)

            [tempcell{1:length(ctf.sensor.index.meg_ref)}]=deal(ctf.sensor.info(ctf.sensor.index.meg_ref).orientation);
            refOrientations = [tempcell{:}];
            ref_direction = refOrientations';  % reject even columns (duplicates)
            
            [ref_sensor_labels{1:length(ctf.sensor.index.meg_ref)}] = deal(ctf.sensor.info(ctf.sensor.index.meg_ref).label);

            %strip off weird '-2202' on the labels
            for i=1:length(ref_sensor_labels)
                ref_sensor_labels{i} = strrep(ref_sensor_labels{i},'-2202','');
            end


            %%% get synthetic gradient coefs
            COEFS = true;
            for trialidx=1:numtrials
                ctf = ctf_read(megfolder,kept_channels,'all',trialidx,COEFS);
                R = R + ctf.data*ctf.data';
            end
            Gcoef = cell2mat({ctf.sensor.info(kept_channels).Gcoef}');

            clear tempcell
        else
            disp('No synthetic gradient correction detected.');
            for trialidx=1:numtrials
                ctf = ctf_read(megfolder,kept_channels,'all',trialidx);
                R = R + ctf.data*ctf.data';
            end
        end
        
        R = R/numtrials;


        [tempcell{1:length(kept_channels)}]=deal(ctf.sensor.info(kept_channels).location);  % in cm
        sensorLocations = [tempcell{:}]*10;  % convert to millimeters for compatibility
        coil_coord = sensorLocations(:,1:2:end)';  % odd columns (closer sensor in gradiometer pair)
        coil_coord2 = sensorLocations(:,2:2:end)';  % even columns (farther sensors in gradiometer pair)

        [tempcell{1:length(kept_channels)}]=deal(ctf.sensor.info(kept_channels).orientation);
        sensorOrientations = [tempcell{:}]';

        clear tempcell;

        normal_direction = sensorOrientations;
        start_time = ctf.setup.start_sec;

        latency = ctf.setup.time_msec';
        srate=ctf.setup.sample_rate;

        % get local sphere origin + head radius from default.hdm -- we will want to
        % compute our own from the mesh later...
        % should be on lines 33-36; dlmread indexes from 0 -> rows 32-35
        %         if(~exist('lsc_fullpath','var'))
        %             if(exist(fullfile(megfolder,'default.hdm'),'file'))
        %                 lsc_fullpath = fullfile(megfolder,'default.hdm');
        %             elseif(exist(fullfile(pwd,'default.hdm'),'file'))
        %                 lsc_fullpath = fullfile(pwd,'default.hdm');
        %             else
        %                 [lsc_filename, lsc_path]=uigetfile('*.hdm','Select CTF file with lsc...');
        %                 if isequal(lsc_filename,0)|isequal(lsc_path,0)
        %                     lsc_fullpath = [];
        %                 else
        %                     lsc_fullpath = fullfile(lsc_path,lsc_filename);
        %                 end
        %             end
        %         end
        if(exist(fullfile(megfolder,'localSpheres.hdm'),'file'))
            nut_load_fiducials(fullfile(megfolder,'localSpheres.hdm'))
            disp('you are using multiple spheres for lead field since there is a localSpheres.hdm in your .ds directory')
            disp('be sure to load the appropriate MRI in the Coregistration Tool')

        elseif((isfield(nuts,'meg') && ~isfield(nuts.meg,'lsc')) | ~isfield(nuts,'meg'))
            nut_load_fiducials
        end



        % if no file specified after all this, come up with a best guess for sphere origin
%         if isempty(lsc_fullpath)
%             global st;
%             if(isfield(nuts,'mesh'))  %%% guess lsc from mean of mesh coords (not very accurate)
%                 headsurf_vx = reshape(nuts.coreg.mesh,size(nuts.coreg.mesh,1)*size(nuts.coreg.mesh,2),3);
%                 headsurf = nut_voxels2mm(headsurf_vx);
%                 lsc_mri = mean(headsurf);
%                 lsc = nut_mri2meg(lsc_mri); % convert to meg coords since user is more likely to be familiar with lsc values in MEG coord system
%                 r_head = mean(nut_rownorm(headsurf - repmat(lsc_mri,[length(headsurf) 1])));
%             else  % no mesh, so use default guess of (0,0,50) mm
%                 lsc = [0 0 50];
%             end
% 
%             prompt   = 'Gimme a local sphere origin (MEG coords, in mm)';
%             title    = 'Input lsc';
%             lines = 1;
%             def{1}   = num2str([lsc]);
% 
%             answer   = inputdlg(prompt,title,lines,def);
%             if (isempty(answer)) msgbox('What''s wrong with you?! You really need to either specify a head model file or enter a local sphere origin.');return; end;
%             lsc = str2num(cell2mat(answer));
% 
%         end

        sensor_labels = ctf.sensor.label;


    case 'matfile'
        nut_importmat(megdata_filename);
    otherwise  % user cancelled
        return
end


% if(isfield(nuts,'meg'))
%     nuts.meg = []; % remove any old MEG info from nuts
% end

% remove any lead field since we've obviously made it invalid by loading a new MEG dataset
if(isfield(nuts,'Lp'))
    nuts = rmfield(nuts,'Lp');
end

nuts.meg.filename = megfile;
%nuts.meg.pathname = megdata_fullpath;
if ~isempty(gcbf)  %will be empty during batching
    handles = guihandles(gcbf);
    if isfield(handles,'nut_megfile')
        set(handles.nut_megfile,'String',megfile);
    end
end

% nuts.meg.data = meg;
nuts.meg.Ract = R;

%%% lsc now specified in coregistration for CTF, we still do this line up
%%% in BTi and KIT code though.
% nuts.meg.lsc = lsc;
nuts.meg.sensorCoord(:,:,1) = coil_coord;
nuts.meg.sensorCoord(:,:,2) = coil_coord2;
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


if isfield(handles,'nut_beamforming_button')  %sometimes this will be called from activation viewer, and so don't want to call nut_enabler
    nut_enabler;
end




function meg2norm=kit_matrix(matrix_name)
%read KIT sensor coord to head coord transformation matrix
disp('point conversion info file')
%      [filename, pathname, filterindex] = uigetfile('*matrix.txt', 'conversion info file');
%       matrix_name=[pathname,filename];
fid=fopen(matrix_name);

str = fgetl(fid);

while  length(str)~=length('[Transform Matrix]')|str~='[Transform Matrix]'
    str = fgetl(fid);
end
str1 = fgetl(fid);
str2 = fgetl(fid);
str3 = fgetl(fid);

ku1=findstr(str1,'=')+1;
ku246=findstr(str1,'*')-1;, ku2=ku246(1);, ku4=ku246(2);, ku6=ku246(3);
ku357=findstr(str1,'+')+1;, ku3=ku357(1);, ku5=ku357(2);, ku7=ku357(3);
ku8=findstr(str1,'[')-1;

mat11=str2num( str1(ku1:ku2) );, mat12=str2num( str1(ku3:ku4) );,
mat13=str2num( str1(ku5:ku6) );, mat14=str2num( str1(ku7:ku8) );

ku1=findstr(str2,'=')+1;
ku246=findstr(str2,'*')-1;, ku2=ku246(1);, ku4=ku246(2);, ku6=ku246(3);
ku357=findstr(str2,'+')+1;, ku3=ku357(1);, ku5=ku357(2);, ku7=ku357(3);
ku8=findstr(str2,'[')-1;

mat21=str2num( str2(ku1:ku2) );, mat22=str2num( str2(ku3:ku4) );,
mat23=str2num( str2(ku5:ku6) );, mat24=str2num( str2(ku7:ku8) );


ku1=findstr(str3,'=')+1;
ku246=findstr(str3,'*')-1;, ku2=ku246(1);, ku4=ku246(2);, ku6=ku246(3);
ku357=findstr(str3,'+')+1;, ku3=ku357(1);, ku5=ku357(2);, ku7=ku357(3);
ku8=findstr(str3,'[')-1;

mat31=str2num( str3(ku1:ku2) );, mat32=str2num( str3(ku3:ku4) );,
mat33=str2num( str3(ku5:ku6) );, mat34=str2num( str3(ku7:ku8) );

meg2norm=[mat11 mat12 mat13 mat14;mat21 mat22 mat23 mat24;mat31 mat32 mat33 mat34];

fclose(fid);


function [lsc,r]= kit_sphere(spferef_name)
% read KIT sphere origin
disp('point sphere info file');
%     [filename, pathname, filterindex] = uigetfile('*sphere.txt', 'sphere info file');
%     spferef_name=[pathname,filename];
fid=fopen(spferef_name);
%
for i=1:4
    str = fgetl(fid);
end
str = fgetl(fid);
idx1 = findstr(str,'center');
idx2 = findstr(str,',');
idx3 = findstr(str,')');
idx4 = findstr(str,'radius');
idx5 = findstr(str,'[');

center_x = str2num(str(idx1(1,1)+8:idx2(1,2)-1));
center_y = str2num(str(idx2(1,2)+1:idx2(1,3)-1));
center_z = str2num(str(idx2(1,3)+1:idx3(1,1)-1));
r = str2num(str(idx4(1,1)+7:idx5(1,1)-1));

% fprintf('lsc=[%f,%f,%f] \n head_radius=%f \n',center_x,center_y,center_z,radius);
fclose(fid);
lsc=[center_x center_y center_z];   % should be in mm

