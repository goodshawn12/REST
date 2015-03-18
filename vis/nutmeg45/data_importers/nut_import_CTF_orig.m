function nut_import_ctf(megfolder)

global nuts ndefaults

currentdir = pwd;  % remember current directory so we can get back to it later

[crap, filename, ext] = fileparts(megfolder);
clear crap;

ctf = ctf_read_res4(megfolder);
if ctf.setup.number_trials > 3
    trials = 1:ctf.setup.number_trials;

    %            warndlg('Please load averaged data. Do you really think you can load all those trials without running out of memory??? I mock thee. Begone!');
else
    % if number of trials is 1, 2, or 3, dataset is most likely
    % averaged, and only first "trial" contains relevant info
    trials = 1:1;
end


% setting up "meg" might run out of memory, so let's do it before we read in our big-ass file
%         meg = zeros(length(kept_channels), ctf.setup.number_samples, length(trials));  % channels x time x trials

%%%%% pick channels to read
kept_channels = ctf.sensor.index.meg_sens;
% kept_channels = union(ctf.sensor.index.meg_sens,ctf.sensor.index.eeg_sens);
%kept_channels=1:length(ctf.sensor.info); % useful for grabbing all channels (not just MEG channels)

[tempcell{1:length(kept_channels)}]=deal(ctf.sensor.info(kept_channels).grad_order_no);
grad_order = max(cell2mat(tempcell));
clear tempcell;

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
    ctf = ctf_read(megfolder,kept_channels,'all',trials,COEFS);
    Gcoef = cell2mat({ctf.sensor.info(kept_channels).Gcoef}');

    clear tempcell
else
    disp('Hmmm... doesn''t look like the file was saved with synthetic gradient correction.');
    ctf = ctf_read(megfolder,kept_channels,'all',trials);
end
%         if iscell(ctf.data)
%             for ii=1:length(ctf.data)
%                 meg(:,:,ii)=cell2mat(ctf.data(ii,:));
%             end
% %         meg = cell2mat(ctf.data);   % take first (and hopefully only) "trial" and keep in meg matrix
%         else
nuts.meg.data = ctf.data;   % take first (and hopefully only) "trial" and keep in meg matrix
%         end

[tempcell{1:length(kept_channels)}]=deal(ctf.sensor.info(kept_channels).location);  % in cm
sensorLocations = [tempcell{:}]*10;  % convert to millimeters for compatibility
clear tempcell;

nuts.meg.sensorCoord(:,:,1) = sensorLocations(:,1:2:end)';  % odd columns (closer sensor in gradiometer pair)
nuts.meg.sensorCoord(:,:,2) = sensorLocations(:,2:2:end)';  % even columns (farther sensors in gradiometer pair)

[tempcell{1:length(kept_channels)}]=deal(ctf.sensor.info(kept_channels).orientation);
sensorOrientations = [tempcell{:}]';
clear tempcell;
nuts.meg.sensorOrient = sensorOrientations;

start_time = ctf.setup.start_sec;

nuts.meg.srate = ctf.setup.sample_rate;
nuts.meg.latency = ctf.setup.time_msec;


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


%         if(exist(fullfile(megfolder,'localSpheres.hdm'),'file'))
%             nut_load_fiducials(fullfile(megfolder,'localSpheres.hdm'))
%             disp('you are using multiple spheres for lead field since there is a localSpheres.hdm in your .ds directory')
%             disp('be sure to load the appropriate MRI in the Coregistration Tool')
%
if(exist(fullfile(megfolder,'default.hdm'),'file'))
    nut_load_fiducials(fullfile(megfolder,'default.hdm'))
    disp('Using the default.hdm in your .ds directory...')
    disp('Be sure to load the appropriate MRI in the Coregistration Tool and verify fiducials')

    % i admit i'm a little embarrassed about this. but it works, so just go with it, okay?
    % no it doesn't always!
    global coreg
    if (~isempty(nuts.coreg) & isempty(coreg))
        %then leave nuts.coreg alone, since no fiducial info was found in this default.hdm
    else
        nuts.coreg = coreg;
    end
    clear global coreg

elseif((isfield(nuts,'meg') && ~isfield(nuts.meg,'lsc')) | ~isfield(nuts,'meg'))
    nut_load_fiducials

    global coreg
    nuts.coreg = coreg;
    clear global coreg
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


if ndefaults.meg.single
    nuts.meg.data=single(nuts.meg.data);
else
    nuts.meg.data=double(nuts.meg.data);
end

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
nuts.meg.sensor_labels = ctf.sensor.label;
if size(nuts.meg.sensor_labels,1) > size(nuts.meg.sensor_labels,2)
    nuts.meg.sensor_labels = nuts.meg.sensor_labels'; % must be a row vector
end


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
nuts.meg.keepref_ts = 0; % new CTF code allows keeping reference channels; nutmeg now expects this flag

