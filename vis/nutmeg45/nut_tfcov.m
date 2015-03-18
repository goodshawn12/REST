function [Ra,Rc]=nut_tfcov(meg,active,control,trialselect)
% NUT_TFCOV  calculates sensor covariance matrix for spatial beamformer.
%    [Ra,Rc]=nut_tfcov(meg,active,control,[trial])
% optionally specify specific trial(s) for, e.g., Wilcoxon tfbf

numtrial=size(meg.data,3);

% if nargin<3     % old syntax of this function had active and control concatenated.
%     control=active(end,:);
%     active(end,:)=[];
% end

% if (exist('trialselect','var'))
%     error('trialselect input is deprecated. please create meg.markers yourself')
% end
if(exist('trialselect','var'))
    % if we are selecting specific trials, mark designated trial as "good"
    [meg.markers.active,meg.markers.control] = deal(nan(numtrial,1));
    if iscell(trialselect)
        meg.markers.active(trialselect{1})=1;
        meg.markers.control(trialselect{2})=1;
    else
        meg.markers.active(trialselect)=1;
        meg.markers.control(trialselect)=1;
    end
elseif (~isfield(meg,'markers') || isempty(meg.markers))
% if (~isfield(meg,'markers') || isempty(meg.markers))
    disp('you have no meg.markers set.  All trials used')
    [meg.markers.active,meg.markers.control] = deal(zeros(numtrial,1)); % '0' is used because isfinite
end

isdual=~isempty(control);
if ~isdual, meg.markers.control=[]; end


%for ii = 1:size(timewins,1)
%    timepts(ii,:) = dsearchn(meg.latency,timewins(ii,1)):dsearchn(meg.latency,timewins(ii,2));
%end

% bpdata = nut_filter2(meg.data,'firls','bp',4,freqband(1),freqband(2),meg.srate,1);
% figure;subplot(1,2,1);plot(bpdata(:,200,1));

%lowmem = true;
%if(lowmem)  % this method might make more sense even with lots of RAM
Ra = zeros(size(meg.data,2),size(meg.data,2),size(active,1));
if isdual
    Rc = zeros(size(meg.data,2),size(meg.data,2),size(control,1));
else
    Rc=[];
end

goodactivetrials=find(isfinite(meg.markers.active))';    % take care of contrasts with different conditions in different trials
goodcontroltrials=find(isfinite(meg.markers.control))';
numactivetrials=length(goodactivetrials)
numcontroltrials=length(goodcontroltrials)

for trial=goodactivetrials
    curractive  = active  + meg.markers.active(trial);   % take care of timewins that are relative to markers
    for ii = 1:size(active,1)
        timepts = dsearchn(meg.latency,curractive(ii,1)):dsearchn(meg.latency,curractive(ii,2));
        Ra(:,:,ii) = Ra(:,:,ii) + cov(double(meg.data(timepts,:,trial)));
        % need double here to preserve numerical accuracy
        % meg.data is on the order of 1e-15 -> cov(meg.data) gets
        % close to realmin('single')
        % alternatively we can convert to femtoTesla, but this requires massive
        % nutmeg-wide changes (including lead field code)
    end
end
    
for trial=goodcontroltrials
    currcontrol = control + meg.markers.control(trial);
    for jj = 1:size(control,1)
        timepts = dsearchn(meg.latency,currcontrol(jj,1)):dsearchn(meg.latency,currcontrol(jj,2));
        Rc(:,:,jj) = Rc(:,:,jj) + cov(meg.data(timepts,:,trial));
%             R(:,:,ii) = R(:,:,ii) + meg.data(timepts(ii,:),:,trial)'*meg.data(timepts(ii,:),:,trial);
    end
end

Ra = Ra./numactivetrials;
if isdual
    Rc = Rc./numcontroltrials;
end

% if nargout==1   % old syntax
%     Ra=cat(3,Ra,Rc);
%     clear Rc
% end
    
%else
%    R = zeros(size(meg.data,2),size(meg.data,2),size(meg.data,3));
%    
%    for trial=1:size(meg.data,3)
%        %             covactive(:,:,trial) = bpdata(timepts,:,trial)'*bpdata(timepts,:,trial);
%        %             covcontrol(:,:,trial) = bpdata(controlwin,:,trial)'*bpdata(controlwin,:,trial);
%        R(:,:,trial) = cov(meg.data(timepts,:,trial));
%    end
%    for ii = 1:size(timewins,1)
%        R(:,:,ii) = mean(R(:,:,ii),3);
%    end
%end


% subplot(1,2,2);imagesc(R(:,:,1));
% cond(R(:,:,1))