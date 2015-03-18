function nuts=nut_importmeg(megdata_fullpath)
% NUT_IMPORTMEG
% currently supports CTF, 37/74-channel BTi Magnes, and KIT

global nuts ndefaults;


if exist('megdata_fullpath','var')
    [megdata_path,megdata_filename,megdata_ext]=fileparts(megdata_fullpath);
    megdata_path = [megdata_path '/'];  % to maintain consistency with uigetfile output
    megdata_filename = [megdata_filename megdata_ext];

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
            case '.fif'
                megsys = 'Neuromag'
            case '.mat'
                megsys = 'matfile'
            case '.avg'                 % NUTEEG mod
                megsys = 'Neuroscan'
            case '.eeg'
                megsys = 'Neuroscan'
            otherwise
                errordlg('Sorry, unknown data type.');
        end
        megdata_filename = fullfile(megdata_path,megdata_filename);
    end

else
    nutmegpath = fileparts(which('nutmeg'));
    loadpath = fullfile(nutmegpath,'data_importers');
    pluginfiles = dir(fullfile(loadpath,'nut_*.m'));
    for k = 1:length(pluginfiles)
        [p pluginname e] = fileparts(pluginfiles(k).name);
        % remove 'nut_' and replace other underscores with a space
        importers{k} = strrep(strrep(pluginname, 'nut_',''),'_',' ');
    end
    megsys_ind=menu('Import which data type?',importers);
    megsys=pluginfiles(megsys_ind).name(12:end-2);

    if(0) % for now..
    [megdata_filename, megdata_path]=uigetfile('*.bin;*.meg4;*.txt;*.fif','Select MEG data file (BTi,CTF,KIT,Neuromag)...');
    if isequal(megdata_filename,0)|isequal(megdata_path,0)
        return;
    end
    end

end

% get rid of any old information that might confuse things later...
% remove any lead field since we've obviously made it invalid by loading a new MEG dataset
% if(isfield(nuts,'Lp'))
%     nuts = rmfield(nuts,'Lp');
% end

%get rid of previously loaded MEG data since channel dimensions might be different
if isfield(nuts,'meg')
    %    if isfield(nuts.meg,'lsc')
    %       tempmeg.lsc=nuts.meg.lsc; % hang on to sphere center, if present -- wait, WHY?
    %       if isfield(nuts.meg,'lsc_sensor_labels')
    %           tempmeg.lsc_sensor_labels=nuts.meg.lsc_sensor_labels;
    %       end
    nuts=rmfield(nuts,'meg')
    %       nuts.meg.lsc=tempmeg.lsc;
    %       if isfield(tempmeg,'lsc_sensor_labels')
    %           nuts.meg.lsc_sensor_labels=tempmeg.lsc_sensor_labels;
    %       end
    %       clear tempmeg
    %    end
end

%%%% FIXME FIXME FIXME: make following 'cases' into function calls
%%%% [meg, latency, sensors, lsc] = bti_read(path,file) or ctf_read or kit_read;

% nuts.meg.filename=megfile;

nuts.meg.system=megsys;
switch megsys
    case 'BTi_UCSF', % special UCSF BTI-export format
        %             [megdata_filename, megdata_path]=uigetfile('*.bin;*.meg4;*.txt;*.fif','Select MEG data file (BTi,CTF,KIT,Neuromag)...');
        [megdata_filename, megdata_path]=uigetfile('*.bin','Select MEG BTi data file');
        if isequal(megdata_filename,0)|isequal(megdata_path,0)
            return;
        end
        [crap,megfile,morecrap] = fileparts(megdata_filename);
        clear crap morecrap;
        nuts.meg.filename=megfile;
        nut_import_BTi_UCSF(fullfile(megdata_path,megdata_filename));
    case '4D_BTi', % uses MSI>>Matlab functions
        %             [megdata_filename, megdata_path]=uigetfile('*.bin;*.meg4;*.txt;*.fif','Select MEG data file (BTi,CTF,KIT,Neuromag)...');
        [megdata_filename, megdata_path]=uigetfile('*','Select 4D/BTi MEG data file');
        if isequal(megdata_filename,0)|isequal(megdata_path,0)
            return;
        end
        [crap,megfile,morecrap] = fileparts(megdata_filename);
        % obtain parent directory, which is more descriptive for 4D data
        [crap,parentdir] = fileparts(megdata_path(1:(end-1)));
        clear crap morecrap
        nuts.meg.filename=[parentdir '-' megdata_filename];
        nut_import_4D_BTi(fullfile(megdata_path,megdata_filename));
    case 'KIT'
        [megdata_filename, megdata_path]=uigetfile('*.txt','Select KIT MEG data file');
        if isequal(megdata_filename,0)|isequal(megdata_path,0)
            return;
        end
        [crap,megfile,morecrap] = fileparts(megdata_filename);
        clear crap morecrap;
        nuts.meg.filename=megfile;
        nut_import_KIT(fullfile(megdata_path,megdata_filename));
    case 'CTF'
        if ~exist('megdata_fullpath','var')
            [megdata_filename, megdata_path]=uigetfile('*.meg4','Select CTF MEG data file');
            if isequal(megdata_filename,0)|isequal(megdata_path,0)
                return;
            end
        end
        [crap,megfile,morecrap] = fileparts(megdata_filename);
        clear crap morecrap;
        if(strcmp(megdata_path((end-3):(end-1)),'.ds'))
            megfolder = megdata_path(1:(end-1));
        else(strcmp(megdata_filename((end-2):end),'.ds'))
            megfolder = megdata_fullpath;
        end
        nuts.meg.filename=megfile;
        nut_import_CTF(megfolder,ndefaults.ctf.gradcorr);
    case 'CTF_orig'
        if ~exist('megdata_fullpath','var')
            [megdata_filename, megdata_path]=uigetfile('*.meg4','Select CTF MEG data file');
            if isequal(megdata_filename,0)|isequal(megdata_path,0)
                return;
            end
        end
        [crap,megfile,morecrap] = fileparts(megdata_filename);
        clear crap morecrap;
        if(strcmp(megdata_path((end-3):(end-1)),'.ds'))
            megfolder = megdata_path(1:(end-1));
        else(strcmp(megdata_filename((end-2):end),'.ds'))
            megfolder = megdata_fullpath;
        end
        nuts.meg.filename=megfile;
        nut_import_CTF_orig(megfolder);       
    case 'via_FieldTrip_fileio'
        if ~exist('megdata_fullpath','var')
            [megdata_filename, megdata_path]=uigetfile('*','Select MEG/EEG data file');
            if isequal(megdata_filename,0)|isequal(megdata_path,0)
                return;
            end
        end
        nut_import_via_FieldTrip_fileio(fullfile(megdata_path,megdata_filename)); 
        nuts.meg.filename=megdata_filename;
    case 'Cartool_EPH'
        datapath = uigetdir(pwd,'Select Folder containing all EPH files of current subject and condition.');
        if isequal(datapath,0), return, end
        nut_import_Cartool_EPH(datapath);
    case 'Neuromag'
        [megdata_filename, megdata_path]=uigetfile('*.fif','Select Neuromag MEG data file');
        if isequal(megdata_filename,0)|isequal(megdata_path,0)
            return;
        end
        [crap,megfile,morecrap] = fileparts(megdata_filename);
        clear crap morecrap;
        nuts.meg.filename=megfile;
        nut_import_Neuromag(fullfile(megdata_path,megdata_filename));
    case 'BrainProducts'
        [megdata_filename, megdata_path]=uigetfile('*.eeg;*.dat','Select BrainProducts EEG data file');
        if isequal(megdata_filename,0)|isequal(megdata_path,0)
            return;
        end
        [crap,megfile,morecrap] = fileparts(megdata_filename);
        clear crap morecrap;
        nuts.meg.filename=megfile;
        nut_import_BrainProducts(fullfile(megdata_path,megdata_filename));
    case 'matfile'
        nut_importmat(megdata_filename);
    case 'Neuroscan'
        [megdata_filename, megdata_path]=uigetfile('*.eeg;*.avg','Select Neuroscan EEG data file');
        [crap,megfile,morecrap] = fileparts(megdata_filename);
        nuts.meg.filename = megfile;
        clear crap morecrap;
        nut_import_Neuroscan('filename',strcat(megdata_path,megdata_filename));
        nut_import_eeg_coords;
    case 'FieldTrip'
        ft_defaults; % this sets proper paths for FieldTrip
        nut_import_FieldTrip;

    otherwise  % user cancelled
        return
end

%nuts.meg.pathname = megdata_fullpath;
if ~isempty(gcbf)  %will be empty during batching
    handles = guihandles(gcbf);
    if isfield(handles,'nut_megfile') && isfield(nuts.meg,'filename')
        set(handles.nut_megfile,'String',nuts.meg.filename);
    end
end

if exist('handles') %will not exist during batching
    if isfield(handles,'nut_beamforming_button')  %sometimes this will be called from activation viewer, and so don't want to call nut_enabler
        nut_enabler;
    end
end


if(0)  %this is all inside relevant importers now
    
    % if iscell(meg) %current issue of new import meg, should remove this eventually
    %     nuts.meg.data=cell2mat(meg);
    %     nuts.meg.data=reshape(nuts.meg.data,length(meg),size(meg{1},1),size(meg{1},2));
    %     nuts.meg.data=permute(nuts.meg.data,[2 3 1]);
    % else
    %     nuts.meg.data=meg;
    % end
    if ndefaults.meg.single
        nuts.meg.data = meg;
    else
        %     nuts.meg.data=double(nuts.meg.data);
        nuts.meg.data=double(meg);
    end

    %%% lsc now specified in coregistration for CTF, we still do this line up
    %%% in BTi and KIT code though.
    % nuts.meg.lsc = lsc;
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
end


