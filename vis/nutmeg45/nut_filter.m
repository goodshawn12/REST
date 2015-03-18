function data = nut_filter(data,timewindow,fs,baseline,notch,bpf,bpf_low_cutoff,bpf_high_cutoff)
% Used to filter time series (input 'data' can be 2 or 3-dimensional, output 'data' is always 2-dimensional)
% Program is called by nut_results_viewer.m and nut_beamforming_gui.m

global nuts % NUTEEG mod
if isfield(nuts.meg,'eegflag') & nuts.meg.eegflag==1
    disp('Subtracting mean from EEG data...');
    data = data - repmat(mean(data,2),[1 size(data,2) 1]);
end

butter_filt = 1; % if 1 use butterworth.m, else use firpm.m
filtorder = 2;

K = size(data,3);

if bpf_high_cutoff < bpf_low_cutoff, error('High cutoff frequency is less than low cutoff frequency'); end

% apply activation_notch filter
if notch == 1
elseif notch == 2
   if butter_filt == 1
      [b,a] = butter(filtorder,2*[47 53]/fs,'stop'); % 50 Hz notch filter
   else
      f = 0:0.05:1; if rem(length(f),2)~=0, f(end)=[]; end
      z = ones(1,length(f)); [val,pos] = min(abs(fs*f/2 - 50)); z(pos) = 0;
      b = firpm(10,f,z); % linear phase 50 Hz notch filter (10th order)
      a = 1;
   end
   if K == 1
      data = filtfilt(b,a,data);
   else
      for k = 1:K
         data(:,:,k) = filtfilt(b,a,squeeze(data(:,:,k)));
      end
   end
elseif notch == 3
   if butter_filt == 1
      [b,a] = butter(filtorder,2*[57 63]/fs,'stop'); % 60 Hz notch filter
   else
      f = 0:0.05:1; if rem(length(f),2)~=0, f(end)=[]; end
      z = ones(1,length(f)); [val,pos] = min(abs(fs*f/2 - 60)); z(pos) = 0;
      b = firpm(10,f,z); % linear phase 60 Hz notch filter (10th order)
      a = 1;
   end
   if K == 1
      data = filtfilt(b,a,data);
   else
      for k = 1:K
         data(:,:,k) = filtfilt(b,a,squeeze(data(:,:,k)));
      end
   end
else
   error('Invalid value for notch filter')
end

% apply band pass filter
if bpf == 0
elseif bpf == 1
   if butter_filt == 1
      if bpf_low_cutoff == 0
         [b,a] = butter(filtorder,2*bpf_high_cutoff/fs);
      elseif bpf_high_cutoff == 0
         [b,a] = butter(filtorder,2*bpf_low_cutoff/fs,'high');
      else
         [b,a] = butter(filtorder,2*[bpf_low_cutoff bpf_high_cutoff]/fs);
         %[b,a] = cheby1(filtorder,0.5,2*[bpf_low_cutoff bpf_high_cutoff]/fs);
         %[b,a] = cheby2(filtorder,20,2*[bpf_low_cutoff bpf_high_cutoff]/fs);
      end
   else
      f = 0:0.001:1; if rem(length(f),2)~=0, f(end)=[]; end
      z = zeros(1,length(f));
      [val,pos1] = min(abs(fs*f/2 - bpf_low_cutoff));
      [val,pos2] = min(abs(fs*f/2 - bpf_high_cutoff));
      z(pos1:pos2) = hanning(pos2-pos1+1);
      %z(pos1:pos2) = 1;
      b = firpm(10,f,z); % 10th order linear phase bpf
      a = 1;
   end
   if K == 1
      meandata = repmat(mean(data),size(data,1),1);
      data = data - meandata; % remove mean prior to applying filter (to reduce transients)
      for k = 1:size(data,2)
         data(:,k) = filtfilt(b,a,data(:,k));
      end
      data = data + meandata; % add the mean back
   else
      for k = 1:K
         tmp = squeeze(data(:,:,k));
         meandata = repmat(mean(tmp),size(data,1),1);
         tmp = tmp - meandata; % remove mean prior to applying filter (to reduce transients)
         tmp = filtfilt(b,a,tmp) + meandata;
         data(:,:,k) = tmp;
      end
   end
else
   error('Invalid value for band pass filter')
end

% apply linear detrending
if baseline == 3
   for i = 1:size(data,2)
      if K == 1
         if data(end,i) - data(1,i) ~= 0
            data(:,i) = data(:,i) - (data(1,i):((data(end,i)-data(1,i))/(size(data,1)-1)):data(end,i))';
         end
      else
         for k = 1:K
            if data(end,i,k) - data(1,i,k) ~= 0
               data(:,i,k) = data(:,i,k) - (data(1,i,k):((data(end,i,k)-data(1,i,k))/(size(data,1)-1)):data(end,i,k))';
            end
         end
      end
   end
end

% apply activation_baseline correction
if baseline == 2 | baseline == 3
   if (timewindow(1) < 0)  % assuming there is a prestimulus period (some t < 0)
      f = find(timewindow < 0);
   else  % if there is no prestimulus period (all t >= 0)
      f = 1:length(timewindow);
   end
   if K == 1
      for i = 1:size(data,2)
         data(:,i) = data(:,i) - mean(data(f,i)); % remove dc component of pre-stim data
         %data(:,i) = data(:,i) - mean(data(1:10,i)); % subtract first 10 values of each component
      end
   else
      for k = 1:K
         for i = 1:size(data,2)
            data(:,i,k) = data(:,i,k) - mean(data(f,i,k)); % remove dc component of pre-stim data
            %data(:,i,k) = data(:,i,k) - mean(data(1:10,i,k)); % subtract first 10 values of each component
         end
      end
   end
end

if K == 1
else
   data_trial_mean = zeros(size(data,1),size(data,2));
   for k = 1:K
      data_trial_mean = data_trial_mean + squeeze(data(:,:,k));
   end
%    data = data_trial_mean;
   data = data_trial_mean/K;
   clear data_trial_mean
end

return

