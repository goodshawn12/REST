function nut_preprocessing_defaults
global nuts

nuts.preprocessing.baseline = 1;
nuts.preprocessing.notch = 1;
nuts.preprocessing.bpf = 0;
nuts.preprocessing.avecov = 0;
nuts.preprocessing.bpf_low_cutoff = 1;
nuts.preprocessing.bpf_high_cutoff = 100;
nuts.preprocessing.stepsize = 25;
% nuts.preprocessing.beamformertype = 1;
nuts.preprocessing.denoisertype = 'No denoising';
% nuts.preprocessing.beamformertype = 'Eigenspace Vector Beamformer';
nuts.preprocessing.beamformertype = 'LCMV Scalar Beamformer';
nuts.preprocessing.invtype = 'tikhonov';
nuts.preprocessing.signalspace = 3;

ptimeinterval = [nuts.meg.latency(1) -1];
nuts.preprocessing.bf_ptimeinterval = ptimeinterval; % time interval for control period
nuts.preprocessing.cs_ptimeinterval = ptimeinterval; % time interval for control period
ptimepts(1)=dsearchn(nuts.meg.latency,ptimeinterval(1));
ptimepts(2)=dsearchn(nuts.meg.latency,ptimeinterval(2));
% set initial cs_intervaltime (based on maximum energy time point)
% data = nut_filter(nuts.meg.data,nuts.meg.latency,nuts.meg.srate,nuts.preprocessing.baseline,nuts.preprocessing.notch,nuts.preprocessing.bpf,nuts.preprocessing.bpf_low_cutoff,nuts.preprocessing.bpf_high_cutoff);
% data = nut_filter2(nuts.meg.data,'firls','notch',50,nuts.preprocessing.bpf_low_cutoff,nuts.preprocessing.bpf_high_cutoff,nuts.meg.srate,1);
data = nut_filter2(mean(nuts.meg.data,3),'firls','baseline_only',4,nuts.preprocessing.bpf_low_cutoff,nuts.preprocessing.bpf_high_cutoff,nuts.meg.srate,ptimepts);
[jnk,pos] = max(sum(data.^2,2)); % find maximum energy index
timeinterval = nuts.meg.latency([pos pos]); % convert from indices to time points
timeinterval = round(timeinterval*10)/10;
nuts.preprocessing.cs_timeinterval = timeinterval; % time interval for channel selection

% set initial bf_intervaltime (based on maximum energy time point)

timeinterval = [0 nuts.meg.latency(end)];
%timeinterval = round(pos + [-0.20 0.20]*length(nuts.meg.latency)); % left/right vertical lines are +-20 percent of data length
%if timeinterval(1) < 1, timeinterval(1)=1; end
%if timeinterval(2) > length(nuts.meg.latency), timeinterval(2) = length(nuts.meg.latency); end
%timeinterval = nuts.meg.latency(timeinterval); % convert from indices to time points
%timeinterval = round(timeinterval*10)/10;
%if timeinterval(1) < nuts.meg.latency(1), timeinterval(1) = nuts.meg.latency(1); end
%if timeinterval(2) > nuts.meg.latency(end), timeinterval(2) = nuts.meg.latency(end); end
nuts.preprocessing.bf_timeinterval = timeinterval; % time interval for beamforming tool

