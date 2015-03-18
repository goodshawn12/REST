function nut_import_BrainProducts(eegdata_filename,stimonset,start_end_samples,srate,sphererad);
% calls Fieldtrip code to load brainvision EEG
% or loads in .dat file if you've already converted using BrainVision

global nuts 

fieldtripcode=which('ft_read_data');
if isempty(fieldtripcode)
    error('Please download, and place in your Matlab path, Fieldtrip code');
end

[dum,eegname,eegext] = fileparts(eegdata_filename);
clear crap;

nuts.meg.filename=eegname;

switch eegext
    case '.eeg'
        data = ft_read_data(eegdata_filename);
    case '.dat'
        data=load(eegdata_filename);
end
data=data';
tot_chan=size(data,2);

if ~exist([eegdata_filename(1:end-4) '.vhdr'])
    [hdr_filename, hdr_path]=uigetfile('*.vhdr','Is there a header file .vhdr for this .dat file?');
end

if exist([eegdata_filename(1:end-4) '.vhdr']) | hdr_filename~=0
    if exist([eegdata_filename(1:end-4) '.vhdr'])
        bp_hdr = read_brainvision_vhdr([eegdata_filename(1:end-4) '.vhdr']);
    elseif hdr_filename~=0
        bp_hdr = read_brainvision_vhdr([hdr_path hdr_filename]);
    end
    nuts.meg.srate=bp_hdr.Fs;
    ecg=find(strcmp(bp_hdr.label,'ECG'));
    eog=find(strcmp(bp_hdr.label,'EOG'));
    noteeg_channels=union(ecg,eog);
    eeg_channels=setdiff([1:tot_chan],noteeg_channels);
    nuts.meg.sensor_labels=bp_hdr.label(eeg_channels);
else
    if ~exist('noteeg_channels','var')
        nch=inputdlg('Which are not EEG channels? (i.e. EOG, ECG etc)');
        noteeg_channels=str2num(nch{1});
    end
    eeg_channels=setdiff([1:tot_chan],noteeg_channels);
    if ~isfield(nuts.meg,'srate')
        if ~exist('srate','var')
            sr=inputdlg('What is sample rate? (Hz)');
            srate=str2num(sr{1});
        end
        nuts.meg.srate=srate;
    end
    if ~isfield(nuts.meg,'sensor_labels')
        for ii=1:length(eeg_channels)
            nuts.meg.sensor_labels{ii}=['E' num2str(ii)];
        end
    end
end




% if ~exist('numtrials','var')
%     ntrials=inputdlg('How many trials?');
%     numtrials=str2num(ntrials{1});
% end

if ~exist('stimonset','var')
    [events_filename, events_path]=uigetfile('*.mat','Load events file converted from *.vmrk to *.mat?');
    if isequal(events_filename,0)|isequal(events_path,0)
        disp('exiting out...');
        return;
    end
    tmp=load([events_path events_filename]);
    stimonset=struct2array(tmp);
end

if ~exist('start_end_samples','var')
    tl=inputdlg('How many samples before and after stimonset per trial?');
    triallength=str2num(tl{1});
end


ttrial=-triallength(1)/nuts.meg.srate:1/nuts.meg.srate:triallength(2)/nuts.meg.srate;
nuts.meg.latency=1000*ttrial';

numtrials=length(stimonset);
for ii=1:numtrials
    eegt(:,:,ii)=data(stimonset(ii)-triallength(1):(stimonset(ii)+triallength(2)),eeg_channels);
    eegtrm(:,:,ii)=eegt(:,:,ii)-repmat(mean(squeeze(eegt(:,:,ii)),1),size(eegt,1),1);
end

nuts.meg.data=eegtrm;




nuts.meg.grad_order=0;
nuts.meg.type=repmat({'eeg'},size(nuts.meg.data,2),1);

if ~exist('sphererad','var')
        sphr=inputdlg('What is subject specific head radius (mm) (roughly 88mm for girls, 92mm for boys)');
        nuts.meg.sphererad=str2num(sphr{1});
end

disp('note: remember to set Local Sphere Centre in Coreg Tool');

[sensor_filename, sensor_path]=uigetfile('*.mat','Load Sensor Coordinates, such as from Polhemus?');
if isequal(sensor_filename,0)|isequal(sensor_path,0)
    return;
end

tmp=load([sensor_path sensor_filename]);
sensors=struct2array(tmp);
nuts.meg.sensorCoord=1000*sensors; % true for Polhemus




