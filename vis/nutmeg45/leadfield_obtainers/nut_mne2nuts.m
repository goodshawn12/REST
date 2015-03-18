function nuts = nut_mne2nuts(nuts,forwfile)
% %% OLD: function nut_mne2nuts(megfile,forwfile)
% converts MGH MNE's data format into nutmeg session structure
% nuts = nut_fcdc2nuts(megfile, forwfile,covfile)
% megfile = file containing MEG data (FIFF format)
% forwfile = file containing forwards solution, i.e., the lead field
% covfile = file containing covariance (but more importantly, also
% containing source locations for forwfile)

% global nuts
mnecode=which('fiff_setup_read_raw');
if isempty(mnecode)
    error('Please place the MNE Matlab software in your Matlab path. If you need to download it, please see http://nmr.mgh.harvard.edu/martinos/userInfo/data/MNE_register/index.php');
end

load_data=0;  % now call nut_import_Neuromag(megfile);
if(load_data)
    [ fid, tree ] = fiff_open(megfile);
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
                sensor_labels{jj}=info.chs(ii).ch_name;
                jj=jj+1;
            end
        end

    else % for raw *.fif data
        blah=fiff_setup_read_raw(megdata_filename);
        [meg,latency] = fiff_read_raw_segment(blah);

        jj=1;
        for ii=1:length(info.chs)
            if info.chs(ii).kind==1 % MEG channel type
                coil_coord(jj,:)=info.chs(ii).loc;
                ind(jj)=ii;
                %                     meg(:,jj)=megtmp(:,ii);
                sensor_labels{jj}=info.chs(ii).ch_name;
                jj=jj+1;
            end
        end
        meg=meg(ind,:)';


    end
    normal_direction=zeros(size(coil_coord,1),1);  % for now this will be meaningless
    coil_coord2=zeros(size(coil_coord));
    srate=info.sfreq;
    nuts.meg.srate = srate;
    % data.hdr.nChans
    % data.hdr.nSamples
    % data.hdr.nSamplesPre
    % data.hdr.nTrials
end


fwd = mne_read_forward_solution(forwfile); % ,force_fixed,surf_ori,include,exclude);
fwd.source_ori; % =1 if one-component (normal to cortex); =2 for 3-component
fwd.nsource;
fwd.nchan;
fwd.coord_frame;
switch(fwd.source_ori)
    % should be size (fwd.nchan x fwd.nsource) for source_ori=1, [fwd.nsource x (3*fwd.nsource)] for source_ori=3;
    case 1
        nuts.meg.Lp(:,1,:) = fwd.sol.data;
    case 212324324320978318905789017
        warning('need to verify order of sources here...');
        nuts.meg.Lp(:,1,:) = fwd.sol.data(:,1:fwd.nsource);
        nuts.meg.Lp(:,2,:) = fwd.sol.data(:,fwd.nsource + (1:fwd.nsource));
        nuts.meg.Lp(:,3,:) = fwd.sol.data(:,2*fwd.nsource + (1:fwd.nsource));
    case 2
        % other possibility for case 2
        nuts.Lp(:,1,:) = fwd.sol.data(:,1:3:(3*fwd.nsource));
        nuts.Lp(:,2,:) = fwd.sol.data(:,2:3:(3*fwd.nsource));
        nuts.Lp(:,3,:) = fwd.sol.data(:,3:3:(3*fwd.nsource));
end

% fwd.mri_head_t.trans   Neuromag-MRI transform
% nuts.voxels = nut_nmag2meg(1000*fwd.source_rr);
nuts.voxels =  nut_nmag2meg(fwd.source_rr*1000);

% nuts.voxels = cov.source_rr;
% warning('SSD: need to check units of cov.source_rr');
% nuts.voxelsize = [5 5 5];
% warning('SSD: need to properly calculate nuts.voxelsize');

% fwd structure:
%methods int32 Has the solution been computed using MEG data (1), EEG data (2), or both (3).
%source_ori int32 Has the solution been computed for the current component normal to the cortex only (1) or all three source orientations (2)
%nsource int32 Total number of source space points.
%nchan int32 Number of channels.
%coord_frame int32 Coordinate frame in which the locations and orientations are expressed.
%source_nn double(*,3) The source orientations. Number of rows is either nsource (?xed source orientations) or 3*nsource (all source orientations).


if(load_data)
    nuts.meg.sensorCoord = data.grad.pnt*1000;
    nuts.meg.Gcoef = 0;
    nuts.meg.grad_order = 0;
    % Gcoef and grad_order should always be zero for Neuromag since there is no reference
    % sensor correction (?), or if there is, it's already been accounted for
    % within the lead field that MNE has already created
    % nuts.meg.sensorOrient = data.grad.ori;
    nuts.meg.latency = data.time{1}'*1000;
    nuts.meg.goodchannels = [1:size(nuts.meg.data,2)];

    nuts.coreg.mripath = which('blank.img');
    nuts.coreg.meg2mri_tfm = eye(4);
end


if(0)  % assign original names
    nuts.meg.sensor_labels = data.grad.label;
    % data.grad.label has chosen sensors as ordered in raw file
    % data.label has chosen sensors, user-specified order
    % data.hdr.label contains all channels from raw file
else   % assign labels based on left/right position

    leftchans = find(nuts.meg.sensorCoord(:,2) >= 0);
    rightchans = find(nuts.meg.sensorCoord(:,2) < 0);
    for i = 1:length(leftchans)
        channel_string = int2str(i);
        nuts.meg.sensor_labels{leftchans(i)} = ['L' '0'*ones(1,4-length(channel_string)) channel_string ' (' data.grad.label{leftchans(i)} ')'];
    end
    for i = 1:length(rightchans)
        channel_string = int2str(i);
        nuts.meg.sensor_labels{rightchans(i)} = ['R' '0'*ones(1,4-length(channel_string)) channel_string ' (' data.grad.label{rightchans(i)} ')'];
    end
end

% selected channels



% lead field
% nuts.meg.lsc
% nuts.meg.lsc_sensor_labels
