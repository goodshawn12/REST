function meg = nut_eegref(meg,REF,domagic,external)
% NUT_EEGREF  rereferences EEG data
%   nuts.meg = nut_eegref(nuts.meg,refname,¦domagic¦,¦numnoneeg¦)
%     (parameters in ¦¦ are optional)
%
% nuts.meg     NUTMEG data structure with EEG (sic!) data.
% refname      Name of reference channel or 'AVG' (for average reference).
% domagic      Average referenced data produces very high condition numbers
%              for the beamformer. This can be resolved when designating
%              one of the average referenced channels as bad. To do the 
%              magic, set domagic to true. Default is false. 
% numnoneeg    If the data has non-EEG channels that should not be used for
%              rereferencing, indicate the number of these extra channels
%              which must be located after all EEG channels in nuts.meg.data.
%              Default is 0.


if nargin<2, help nut_eegref, return, end
if nargin<3, domagic=false; end

if isfield(meg,'reference') && strcmpi(meg.reference,'deref')
    error('Your data has been preprocessed to be reference free.')
end

hadext = (nargin>3 && external>0);
if hadext
    E = meg.data(:,end-external+1:end,:);
    meg.data = meg.data(:,1:end-external,:);
end

% We delete bad channels in the meg.data array.
if length(meg.goodchannels) < size(meg.data,2)
    meg.data=meg.data(:,meg.goodchannels,:);
end

% Rereference
if strcmpi(REF,'avg')
    reference = mean(meg.data, 2);
    meg.reference = 'AVG';
    meg.referenceidx = meg.goodchannels;
    if domagic
        newbadchan=meg.goodchannels(1);
        meg.data(:,newbadchan,:)=[];
        meg.goodchannels(newbadchan)=[];
    end
elseif ( strcmpi(REF,'cz') && strcmpi(meg.system,'Biosemi') )       % Biosemi has Cz as the first channel (A1)
    meg.referenceidx =  find(ismember(1:3,meg.goodchannels),1);
    if isempty(meg.referenceidx), error('All electrodes close to Cz were designated as bad channels.'), end
    reference = meg.data(:,1,:);
    meg.reference = meg.sensor_labels{meg.referenceidx};
elseif strcmpi(REF,'none')
    meg.reference = 'none';
    return
else        % user defined REF electrode label.
    meg.referenceidx = strmatch(REF,meg.sensor_labels,'exact');
    if isempty(meg.referenceidx), error('Chosen reference does not exist.'), end
    if ~ismember(meg.referenceidx,meg.goodchannels), error('Chosen reference was designated as bad channel.'), end
    meg.reference = REF;
    reference = meg.data(:,meg.referenceidx,:);
end

for k=1:size(meg.data,2)
    meg.data(:,k,:) = meg.data(:,k,:) - reference;
end

if hadext
    meg.data = cat(2,meg.data,E);
end


