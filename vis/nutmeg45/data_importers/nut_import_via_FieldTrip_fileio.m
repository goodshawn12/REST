function nut_import_via_FieldTrip_fileio(eegdata_filename)
% calls Fieldtrip code to load EEG data.

global nuts 

fieldtripcode=which('ft_preprocessing');
if isempty(fieldtripcode)
    error('FieldTrip not found -- Please download and place the top-level FieldTrip directory in your Matlab path');
end
ft_defaults

if ~exist(eegdata_filename,'file'), error('Cannot locate file %s',eegdata_filename), end

% Options GUI
MT=nut_maketrial_gui;
if isempty(MT) % if user clicked cancel
    return
end
domt=(isfield(MT,'trialdef') || isfield(MT,'trialfun')); 
dora=isfield(MT,'rema');
    
[eegpath,eegfile,eegext] = fileparts(eegdata_filename);
syst = ft_filetype(eegdata_filename);
if strcmpi(syst,'unknown'), error('Cannot read this file with FieldTrip.'), end
hdr  = ft_read_header(eegdata_filename,'headerformat',syst);

% prepare fieldtrip config
cfg.dataset     = eegdata_filename;
cfg.channel = {'EEG' 'MEG'};
% load in small portion of data so we have complete channel info, etc.
% (ft_read_header does not seem to be sufficient for finding EEG channels...)
cfg.trl = [1 1 0]; % defines a single sample to read in
datatmp = ft_preprocessing(cfg);
if isempty(datatmp.label), datatmp.label=hdr.label; end
cfg.channel = ft_channelselection('gui',datatmp.label); % select desired channels via GUI
cfg = rmfield(cfg,'trl');
%cfg.channel     = {'all'};
% % Try to find EEG channels. FT_CHANNELSELECT with setting 'eeg' seems to
% % bug as on Aug 01 2011.
% switch lower(strtok(syst,'_'))   
%     case 'biosemi'
%         cfg.channel{end+1} = '-Status';
%         R=find(strncmp('EXG',hdr.label,3));
%         for k=1:length(R)
%             cfg.channel{end+1} =  ['-' hdr.label{R(k)}];
%         end
%     % Add other systems here if needed
%     otherwise
%         R=[find(strncmpi('EOG',hdr.label,3));find(strncmpi('ECG',hdr.label,3));find(strncmpi('EKG',hdr.label,3))];
%         for k=1:length(R)
%             cfg.channel{end+1} =  ['-' hdr.label{R(k)}];
%         end
% end

if domt
    if isfield(MT,'trialdef')
        cfg.trialdef=MT.trialdef;
        cfg.trialdef.eventtype ='gui';
        cfg.trialdef.eventvalue='gui';
    else
        cfg.trialfun=MT.trialfun;
    end
    cfg=ft_definetrial(cfg);
elseif dora
    d=dir(fullfile(eegpath,[eegfile '.mrk']));
    if length(d)~=1
        [fMRK, pMRK] = uigetfile('*.mrk;*.txt','Select Cartool marker file');
        if isequal(fMRK,0), return, end
    else
        fMRK=d.name;
        pMRK=eegpath;
    end 
    cfg.trl = getCartoolEpochs(fullfile(pMRK,fMRK));
end
cfg

% select desired channels
%goodlab=ft_channelselection('gui',data.label); 

ssdtest = 0;
if(ssdtest)
%    cfg.bpfilter = 'yes';
%    cfg.bpfreq = [6 90];
    cfg.demean = 'yes';
    cfg.dftfilter = 'yes';
end

% Import via FieldTrip
data=ft_preprocessing(cfg);   % Not really preprocessing, just convenient way of importing data.
% hdr  = ft_read_header(eegdata_filename,'headerformat',hdrformat);
% data = ft_read_data(eegdata_filename,'header',hdr,'dataformat',hdrformat);

% goodchannels = find(ismember(data.label,goodlab))';
numcells = length(data.trial);
 for k=1:numcells
%     data.trial{k}=data.trial{k}(goodchannels,:).';
     data.trial{k}=data.trial{k}';
 end

% Rearrange if requested
if dora
    fprintf('Rearranging to %d sec artifact free trials...\n',MT.rema.secpertrial)
    % find number of available trials
    lens = zeros(1,numcells);
    for tt=1:numcells
        lens(tt)=data.time{tt}(end);
    end
    TOTALTRIALNUMBER = sum(floor(lens./MT.rema.secpertrial))
    clear lens

    % rearrange to trials of MT.rema.secpertrial duration
    numchan=length(data.label);
    lentrial=MT.rema.secpertrial*data.fsample;
    try, D=zeros(lentrial,numchan,TOTALTRIALNUMBER); end   % try to preallocate for speed, but may not be possible because of memory issues.
    countA=1;
    for ee=1:numcells;
        currtrlen= length(data.time{ee});
        index=1;
        while( (index+lentrial-1<=currtrlen) && (ee<=numcells) )      % (countA<=TOTALTRIALNUMBER)
            D(:,:,countA) = data.trial{ee}( index:index+lentrial-1 , : );
            index=index+lentrial;
            countA=countA+1;
        end
        data.trial{ee}=[]; % liberate memory
    end
    data=rmfield(data,'trial');
    if countA<TOTALTRIALNUMBER, error('There seems to be a problem in your code.'), end
    latency  = 1000 .* [0:1/data.fsample:MT.rema.secpertrial-1/data.fsample]';
        
else
    disp('Rearranging to event-locked trials...')
    D = cat(3,data.trial{:,:});
    data=rmfield(data,'trial');
    latency = data.time{1}.*1000;
end

% put data into nuts.meg structure
% nuts.meg=[];

nuts.meg.data           = D;
nuts.meg.srate          = data.fsample;
nuts.meg.filename       = eegfile;
nuts.meg.latency        = latency(:);
nuts.meg.system         = strtok(syst,'_');
nuts.meg.grad_order     = 0;

if(isfield(data,'grad') & any( strcmp(ft_chantype(data.label),'meg')) || any(strcmp(ft_chantype(data.label),'meggrad')))
    % FIXME: not sure what to do if we load in MEG and EEG simultaneously...

    nuts.meg.eegflag        = false;

    % need to find gradiometer channels (including in reference sensors)
    % which need two entries in rawsensor_labels, since each coil needs to
    % have its own entry
    % find(strcmp('refgrad',ft_chantype(data.grad.label)));
    [label1,label2]=find(data.grad.tra == 1);
    nuts.meg.rawsensor_labels{1}=data.grad.label{label1(1)};
    for ii=2:length(label2)
        if(label1(ii) == label1(ii-1))
            nuts.meg.rawsensor_labels{ii-1}=[data.grad.label{label1(ii)} 'a'];
            nuts.meg.rawsensor_labels{ii}=[data.grad.label{label1(ii)} 'b'];
        else
            nuts.meg.rawsensor_labels{ii}=data.grad.label{label1(ii)};
        end
    end
        
    
    nuts.meg.sensor_labels  = data.label;
    if(str2num(ft_version) > 7000) % not sure when this switch was made...
        data.grad = ft_convert_units(data.grad,'mm');
        nuts.meg.sensorCoord = data.grad.coilpos;
        nuts.meg.sensorOrient = data.grad.coilori;
    else % old FT convention
        nuts.meg.sensorCoord = data.grad.pnt*1000;
        nuts.meg.sensorOrient = data.grad.ori;
    end
        
    % these should be "real" EEG/MEG channels (not ref or individual coils)
    for ii=1:size(data.label,1)
        labelmap(ii)=strmatch(data.label{ii},data.grad.label,'exact');
    end
    
    nuts.meg.chanmixMtx{1} = data.grad.tra(labelmap,:);
    nuts.meg.goodchannels = 1:length(nuts.meg.sensor_labels);

%    nuts.meg.data = D(:,labelmap,:);
%
else % not MEG -> EEG
    nuts.meg.eegflag        = true;
	nuts.meg.sensor_labels  = datatmp.label;
    nuts.meg.goodchannels   = find(ismember(datatmp.label,data.label));
    if isfield(data.hdr,'elec') && isfield(data.hdr.elec,'pnt')
        nuts.meg.sensorCoord = data.hdr.elec.pnt*1000;
    end
    
    aw=questdlg('Do you want to average reference your data?','NUTMEG question','Yes','No','Yes');
    if strcmp(aw,'Yes')
        nuts.meg=nut_eegref(nuts.meg,'AVG');
    end
    
end

clear data


% disp('note: remember to set Local Sphere Centre in Coreg Tool');
if(~isfield(nuts.meg,'sensorCoord'))
    aw=menu('Do you want to load electrode coordinates?','Individual digitized coordinates','Template coordinates aligned to MNI template brain','Coordinates previously aligned to individual MRI','None');
    if aw<4
        nut_import_eeg_coords(aw);
    end
end


%-----------------------------
function timerange = getCartoolEpochs(fn)

[idx,idxst,mark]=textread(fn,'%d%d%q','headerlines',1);

badidx=strmatch('bad',mark);
startidx=strmatch('start',mark);
stopidx=strmatch('stop',mark);

repidx=find(diff(startidx)<2);
if ~isempty(repidx)
    error('Marker "start" at sample %d has no "stop".',idx(startidx(repidx(1))))
end
repidx=find(diff(stopidx)<2);
if ~isempty(repidx)
    error('Marker "stop" at sample %d is repeated 2 times.',idx(stopidx(repidx(1))))
end
repidx=find(ismember(badidx+1,startidx));
if ~isempty(repidx)
    error('"Start" marker after "Bad" at sample %d',idx(badidx(repidx(1))))
end
repidx=find(ismember(badidx-1,stopidx));
if ~isempty(repidx)
    error('"Stop" marker before "Bad" at sample %d',idx(badidx(repidx(1))))
end

numseg = length(badidx) + length(startidx);
start = zeros(numseg,1);
stop = zeros(numseg,1);

nummark=length(idx);
counter = 0;
for k=1:nummark
    switch lower(mark{k})
        case 'bad'
            stop(counter)=idx(k);
            counter=counter+1;
            start(counter)=idxst(k);
        case 'start'
            counter=counter+1;
            start(counter)=idx(k);
        case 'stop'
            stop(counter)=idx(k);
    end
end
timerange = [start+1 stop+1 zeros(counter,1)];





