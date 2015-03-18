function nut_import_CTF(megfolder,gradcorr,trials,keepref_ts)
% megfolder:  required: data folder path name
% gradcorr:   optional: desired order gradient correction (default 3rd order)
% trials:     optional: vector of trials to load (default 'all')
% readref_ts: optional: keep time series of reference channels (default no)

global nuts ndefaults

if ~exist('keepref_ts','var')
    nuts.meg.keepref_ts = 0;
else
    nuts.meg.keepref_ts = keepref_ts;
end

currentdir = pwd;  % remember current directory so we can get back to it later

[crap, filename, ext] = fileparts(megfolder);
clear crap;

ctf = nut_readCTFds(megfolder);
% below is obsolete way of detecting pre-averaged data
% if ctf.res4.no_trials > 3
%     trials = 1:ctf.res4.no_trials;
%     % warndlg('Please load averaged data. Do you really think you can load all those trials without running out of memory??? I mock thee. Begone!');
% else
%     % if number of trials is 1, 2, or 3, dataset is most likely
%     % averaged, and only first "trial" contains relevant info
%     warning('assuming averaged data');
%     trials = 1;
% end

%%%%% pick channels to read
channeltypes = [ctf.res4.senres.sensorTypeIndex];
% channeltype codes:
%  0 = reference magnetometers
%  1 = reference gradiometers
%  5 = MEG channels
%  9 = EEG channels
%megchans = find((channeltypes == 0) | (channeltypes == 1) | (channeltypes == 5) | (channeltypes == 9));

% for the moment, only read MEG channels (including reference channels)
coilselect = find((channeltypes == 0) | (channeltypes == 1) | (channeltypes == 5));
refselect = find((channeltypes == 0) | (channeltypes == 1));
MEGindex = find(channeltypes == 5);

% coilselect(1) should be first of MEG and MEGREF sensors used (igonring EOG, STIM etc)
coilselectn=coilselect-coilselect(1)+1;
MEGindexn=MEGindex-coilselect(1)+1;
refselectn=refselect-coilselect(1)+1;

if exist('trials','var')
    dataList=trials;
else
    dataList=[];
end
if ndefaults.meg.single
    nuts.meg.data=getCTFdata(ctf,dataList,coilselect,'fT','single');
else
    nuts.meg.data=getCTFdata(ctf,dataList,coilselect,'fT');
end



%% synthetic coefficient / gradient handling

%warning('we should make it possible to set desired synthetic gradient order here; currently reads whatever is stored')
% desiredbalancing = 'NONE';
% desiredbalancing = [3];
if exist('gradcorr','var')
    switch gradcorr
        case 0
            desiredbalancing = ['NONE';'NONE']; % standard 0th gradient setup
        case 1
            desiredbalancing = ['G1BR';'G1BR']; % standard 1st gradient setup
        case 2
            desiredbalancing = ['G2BR';'G1BR']; % standard 2nd gradient setup
        case 3
            desiredbalancing = ['G3BR';'G1BR']; % standard 3rd gradient setup
        otherwise
            error('gradcorr, if set, must be 0, 1, 2 or 3');
    end
else
    desiredbalancing = ['G3BR';'G1BR']; % standard 3rd gradient setup
    warning(['SARANG has forced default GRADIENT CORRECTION : ' desiredbalancing(1,:) ' !!!']);
end

[nuts.meg.data,ctf]=setCTFDataBalance(nuts.meg.data,ctf,desiredbalancing,'fT',coilselect);


% 100th channel arbitrarily chosen to read grad_order_no (MEG channels should all be same)
nuts.meg.grad_order = ctf.res4.senres(100).grad_order_no;



% we don't need to save Gcoef anymore, since we have nuts.meg.chanmixMtx

switch(nuts.meg.grad_order)
    case 0
        balancing = 'NONE';
    case 1
        balancing = 'G1BR';
    case 2
        balancing = 'G2BR';
    case 3
        balancing = 'G3BR';
end

switch(nuts.meg.grad_order)
    case {1,2,3}
        [alphaMEG,MEGindex_chk,MEGbalanceindex] = getCTFBalanceCoefs(ctf,balancing,'fT');
        Gcoef(:,MEGbalanceindex-coilselect(1)+1) = alphaMEG';
        warning('this still needs to be compared with old Gcoef');
        if nuts.meg.keepref_ts
            Gcoefall=zeros(length(MEGindexn),length(refselectn));
            Gcoefall(:,MEGbalanceindex-coilselect(1)+1)=alphaMEG';
            toHighergrad=eye(length(coilselect));
            toHighergrad(MEGindexn,refselectn)=-Gcoefall;
        end
    case 0
        if nuts.meg.keepref_ts
            [alphaMEG,MEGindex_chk,MEGbalanceindex] = getCTFBalanceCoefs(ctf,'G3BR','fT');
            Gcoef(:,MEGbalanceindex-coilselect(1)+1) = alphaMEG';
            Gcoefall=zeros(length(MEGindexn),length(refselectn));
            Gcoefall(:,MEGbalanceindex-coilselect(1)+1)=alphaMEG';
            toHighergrad=eye(length(coilselect));
            toHighergrad(MEGindexn,refselectn)=-Gcoefall;
            warning('this still needs to be compared with old Gcoef');
        else
            MEGbalanceindex = [];
        end
end


nuts.meg.chanmixMtx=[];
% nuts.meg.chanmixMtx{1}=zeros(length(nuts.meg.sensor_labels),length(nuts.meg.rawsensor_labels));
jj=1;
for ii=1:length(coilselect)
    % these are sensor labels corresponding to data channels
    nuts.meg.sensor_labels{ii}=parse_sensor_label(ctf.res4.chanNames(coilselect(ii),:));
    
    nuts.meg.sensorCoord(jj,:) = ctf.res4.senres(coilselect(ii)).pos(:,1)'*10; % convert from cm to mm;
    nuts.meg.sensorOrient(jj,:) = ctf.res4.senres(coilselect(ii)).ori(:,1)';
    
    [isrefchan,alphaMEGidx]=ismember(coilselect(ii),MEGbalanceindex);
    if(isrefchan)
        Gcoef_vec(ii,jj) = -1;
    end
    
    [isMEGchan,alphaMEGidx]=ismember(coilselect(ii),MEGindex);
    if(isMEGchan & (nuts.meg.grad_order > 0) & ~nuts.meg.keepref_ts)
        Gcoef_coils = Gcoef*Gcoef_vec;
        nuts.meg.chanmixMtx{1}(ii,1:size(Gcoef_coils,2))=Gcoef_coils(alphaMEGidx,:);
    elseif (isMEGchan & (nuts.meg.grad_order > 0) & nuts.meg.keepref_ts)
        % nothing, but just note that chanmixMtx{2} will equal toHighergrad
    end
    
    if(size(ctf.res4.senres(coilselect(ii)).pos,2) == 2) % if it's a gradiometer pair, create 2nd entry for outer coil
        % create new entry and unique label for each coil, with record in chanmixMtx
        %        nuts.meg.rawsensor_labels{jj} = [deblank(ctf.res4.chanNames(coilselect(ii),:)) 'A'];
        nuts.meg.rawsensor_labels{jj} = [nuts.meg.sensor_labels{ii} 'A'];
        nuts.meg.chanmixMtx{1}(ii,jj)=1;
        %         hichanmixMtx(ii,jj)=1;
        
        jj = jj+1;
        nuts.meg.rawsensor_labels{jj} = [nuts.meg.sensor_labels{ii} 'B'];
        nuts.meg.chanmixMtx{1}(ii,jj)=1;
        
        nuts.meg.sensorCoord(jj,:) = ctf.res4.senres(coilselect(ii)).pos(:,2)'*10; % convert from cm to mm;
        nuts.meg.sensorOrient(jj,:) = ctf.res4.senres(coilselect(ii)).ori(:,2)';
        
        if(isrefchan)
            Gcoef_vec(ii,jj) = -1;
        end
        
        if(isMEGchan & (nuts.meg.grad_order > 0) & ~nuts.meg.keepref_ts)
            Gcoef_coils = Gcoef*Gcoef_vec;
            nuts.meg.chanmixMtx{1}(ii,1:size(Gcoef_coils,2))= Gcoef_coils(alphaMEGidx,:);
        elseif (isMEGchan & (nuts.meg.grad_order > 0) & nuts.meg.keepref_ts)
            % nothing, but just note that chanmixMtx{2} will equal toHighergrad
        end
    else
        % create label for simple magnetometer
        nuts.meg.rawsensor_labels{jj} = nuts.meg.sensor_labels{ii};
        nuts.meg.chanmixMtx{1}(ii,jj)=1;
    end
    jj = jj+1;
end



if nuts.meg.keepref_ts
    % save below, even if not applied if nuts.meg.grad_order=0.
    % if nuts.meg.grad_order=0, it won't get multiplied in for lead field computation
    nuts.meg.chanmixMtx{2}=toHighergrad;
    % example how to use this.  if nuts.meg.grad_order=0, then
    %     if ndefaults.meg.single=1
    %         for ii=1:size(nuts.meg.data,3)
    %             data_hiorder(:,:,ii)=[single(nuts.meg.chanmixMtx{2}*double(nuts.meg.data(:,:,ii)'))]';
    %         end
    %     else
    %         for ii=1:size(nuts.meg.data,3)
    %             data_hiorder(:,:,ii)=[nuts.meg.chanmixMtx{2}*nuts.meg.data(:,:,ii)']';
    %         end
    %     end
else
    % remove reference channels from final derived channel list, if not
    % treating reference as normal channel
    channeltypes = [ctf.res4.senres.sensorTypeIndex];
    nuts.meg.data(:,refselect-refselect(1)+1,:)=[];
    nuts.meg.chanmixMtx{1}(1:length(refselect),:)=[];
    nuts.meg.sensor_labels(1:length(refselect))=[];
    %     nuts.meg.chanmixMtx{1}(refselectn,:)=[];
    %     nuts.meg.sensor_labels(refselectn)=[];
end


%start_time = ctf.setup.start_sec;
nuts.meg.srate = ctf.res4.sample_rate;

nuts.meg.latency = 1000*((1:ctf.res4.no_samples)' - ctf.res4.preTrigPts - 1)/ctf.res4.sample_rate;


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




nuts.meg.goodchannels = 1:size(nuts.meg.data,2);



function sensorName = parse_sensor_label(temp)

% sensorName = parse_sensor_label(temp)
% parse sensor label names

temp = deblank(temp');
temp(temp>127) = 0;
temp(temp<0) = 0;
temp = strtok(temp,char(0));
temp = strtok(temp,'-');
sensorName = char(temp)';

return


