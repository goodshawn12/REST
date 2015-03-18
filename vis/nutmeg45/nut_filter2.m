function data = nut_filter2(data,filtclass,filttype,filtorder,lofreq,hifreq,fs,baseline,eegflag)
% Used to filter time series -- should replace nut_filter eventually
% TODO: need to adjust parameters in nut_results_viewer
% 
% nut_filter2(data,filtclass,filttype,filtorder,lofreq,hifreq,fs,baseline)
% data = [samples x channels] or [samples x channels x trials]
% filtclass = 'butter', 'cheby1', 'cheby2', 'ellip', 'firpm', 'firls'
% filttype = 'notch', 'bp' (bandpass), 'low', 'high', or 'baseline_only'
% filtorder = filter order, set to 'auto' for automatic.
% lofreq = low cutoff (Hz)  (0 for lowpass)
% hifreq = high cutoff (Hz) (inf for highpass)
% fs = sampling rate (Hz)
% baseline = 0 (off) or 1 (subtract mean of whole interval)
%            or, e.g., [1 100] (subtract mean of first 100 points)
% eegflag = nuts.meg.eegflag
%
% Program is called by nut_results_viewer.m and nut_beamforming_gui.m

if (strcmp(filttype,'baseline_only') & ( (length(baseline)==1 & baseline==0) | isempty(baseline)))
    error('Crackhead! you have baseline_only selected as filttype, but baseline is not set right')
end
if isempty(filttype), filttype=''; end

global ndefaults
if isempty(ndefaults), nut_defaults; end
% global nuts ndefaults   % NUTEEG mod
% if isfield(nuts.meg,'eegflag') & nuts.meg.eegflag==1
if nargin>8 && eegflag
   disp('Subtracting mean from EEG data...');
   data = data - repmat(mean(data,2),[1 size(data,2) 1]);
end

% apply activation_baseline correction
if ~isempty(baseline)
    if length(baseline) == 2
        f = baseline(1):baseline(2);
    else
        f = 1:size(data,1);
    end
    try
        meandata = repmat(mean(data(f,:,:),1),[size(data,1) 1 1]);
        data = data - meandata; % remove mean prior to applying filter (to reduce transients)
    catch ME
        if ( strcmpi(ME.identifier,'MATLAB:nomem') && all(baseline~=0) )
            for k2=1:size(data,2)
                for k3=1:size(data,3)
                    meandata = mean(data(f,k2,k3),1);
                    data(:,k2,k3) = data(:,k2,k3)-meandata;
                end
            end
        else
            rethrow(ME)
        end
    end
    if all(baseline~=0)
        clear meandata
    end
end

if ~isempty(hifreq) && ((hifreq>0 && hifreq<lofreq))
    error('Crackhead! High cutoff frequency is less than low cutoff frequency!');
end

switch(filttype)
    case 'notch'
        [b,a] = filtdef(filtclass,filtorder,lofreq,hifreq,fs,'stop',size(data,1));
    case 'baseline_only'
        return
    otherwise
        [b,a] = filtdef(filtclass,filtorder,lofreq,hifreq,fs,[],size(data,1));
end


% data forced to be double before going into filtfilt;
% some filter parameters results in NaN or other screwy results with single
for jj=1:size(data,3)
    for ii=1:size(data,2)
        data(:,ii,jj) = filtfilt(b,a,double(data(:,ii,jj)));
    end
end

if(baseline == 0 & ~ndefaults.postbc) % no final baseline correction, so add the mean back
    data = data + meandata;
elseif (baseline ~= 0 & ndefaults.postbc)
    % Baseline correct after filtering
    if length(baseline) == 2
        f = baseline(1):baseline(2);
    else
        f = 1:size(data,1);
    end
    
    meandata = repmat(mean(data(f,:,:),1),size(data,1),1);
    data = data - meandata;
end


%------------------------
function [b,a]=filtdef(filtclass,filtorder,lofreq,hifreq,fs,filttype,ndat)

if(~exist('filttype','var') || isempty(filttype))
    if(isempty(lofreq) || lofreq==0 )
        lofreq = [];
        filttype='low';
    elseif(isempty(hifreq) || hifreq==Inf || hifreq==0)
        hifreq = [];
        filttype='high';
    else
        filttype='bandpass';
    end
end

wn = 2*[lofreq hifreq]/fs;
if ischar(filtorder) || isempty(filtorder) || filtorder==0
    [filtorder,wn] = filtord(filtclass,filttype,lofreq,hifreq,fs,ndat);
end

switch(filtclass)
    case 'butter'
        [b,a] = butter(filtorder,wn,filttype);
    case 'cheby1'
        [b,a] = cheby1(filtorder,0.5,wn,filttype);
    case 'cheby2'
        [b,a] = cheby2(filtorder,20,wn,filttype);
    case 'ellip'
        [b,a] = ellip(filtorder,0.5,20,wn,filttype);
    case {'firpm','firls'}
        if ~iscell(wn)
            if strcmp(filttype,'stop')
                stopband=1;
            else
                stopband=0;
            end
            f = 0:0.001:1;
            if rem(length(f),2)~=0
                f(end)=[];
            end
            if stopband
                z = ones(1,length(f));
            else
                z = zeros(1,length(f));
            end
            if(isfinite(lofreq))
                [val,pos1] = min(abs(fs*f/2 - lofreq));
            else   % not sure when this should kick in....
                [val,pos2] = min(abs(fs*f/2 - hifreq));
                pos1=pos2;
            end
            if(isfinite(hifreq))
                [val,pos2] = min(abs(fs*f/2 - hifreq));
            else
                pos2 = length(f);
            end
            %       z(pos1:pos2) = hanning(pos2-pos1+1);
            if stopband
                z(pos1:pos2) = 0;
            else
                z(pos1:pos2) = 1;
            end
            a = 1;
            switch(filtclass)
                case 'firpm'
                    b = firpm(filtorder,f,z);
                case 'firls'
                    b = firls(filtorder,f,z);
            end
        else
            b = firpm(filtorder,wn{:});
        end
    otherwise
        error('Unknown filter class.')
end

%----------------------------------
function [filtorder,wn] = filtord(filtclass,filttype,lofreq,hifreq,fs,ndat)
        
% MINFREQ=5;
if isempty(hifreq), hifreq=0; end
if isempty(lofreq), lofreq=0; end
nyq = fs/2;

% Transition width of filter
if lofreq<fs/40 && lofreq>0           % make softer edges when filtering low frequencies
    trans_bw = lofreq/2;
elseif hifreq<fs/40 && hifreq > 0
    trans_bw = hifreq/2;
else
    trans_bw = 1;                           
end

% % Ripple in the passband
rp=0.5;
% if ((lofreq < MINFREQ) && (lofreq > 0)) || ((hifreq < MINFREQ) && (hifreq > 0))
%     rp=0.01; 
% else
%     rp=0.0025;
% end
% 
% % Ripple in the stopband
rs=20;
% if ((lofreq < MINFREQ) && (lofreq > 0)) || ((hifreq < MINFREQ) && (hifreq > 0))
%     rs=30; 
% else
%     rs=40;
% end

switch filttype
    case 'bandpass'
        ws=[(lofreq-trans_bw)/nyq (hifreq+trans_bw)/nyq];
        wp=[lofreq/nyq hifreq/nyq];
    case 'low'      
        ws=(hifreq+trans_bw)/nyq;
        wp=hifreq/nyq;
    case 'high'
        ws=(lofreq-trans_bw)/nyq;
        wp=(lofreq)/nyq;
    case 'stop'
        ws=[lofreq/nyq hifreq/nyq];
        wp=[(lofreq-trans_bw)/nyq (hifreq+trans_bw)/nyq];
end

switch filtclass
    case 'butter'
        [filtorder,wn] = buttord(wp,ws,rp,rs);
    case 'cheby1'
        [filtorder,wn] = cheb1ord(wp,ws,rp,rs);
    case 'cheby2'
        [filtorder,wn] = cheb2ord(wp,ws,rp,rs);
    case 'ellip'
        [filtorder,wn] = ellipord(wp,ws,rp,rs);
    case 'firls'    % taken from brainstorm filter
        if lofreq>0,
            filtorder = 3*fix(fs/lofreq);
        elseif hifreq>0,
            filtorder = 3*fix(fs/hifreq);
        end
        if filtorder < 15
            filtorder = 15;
        end
        wn = 2*[lofreq hifreq]/fs;
        if filtorder > floor((ndat-1)/3)
          filtorder = floor(ndat/3) - 1;
        end
%         if filtorder > floor((ndat)/3)
%           filtorder = floor(ndat/3);
%         end
    case 'firpm'
        switch filttype 
            case 'high' 
                fcuts = [ws wp];
                mags = [0 1];
                devs = [0.01 0.05];
            case 'bandpass'              
                fcuts = [ws(1) wp(1:2) ws(2)];
                mags = [0 1 0];
                devs = [0.01 0.05 0.01];
            case 'low'
                fcuts = [wp ws];
                mags = [1 0];
                devs = [0.05 0.01];    
            case 'stop'
                fcuts = [wp(1) ws(1:2) wp(2)];
                mags = [1 0 1];
                devs = [0.05 0.01 0.05];
        end
        [filtorder,fo,ao,w] = firpmord(fcuts,mags,devs,fs);
        wn={fo ao w};
end
filtorder
       
