function [sources,t,sim]=nut_create_sim_sources(sim);
% function [sources,t]=nut_create_sim_sources(sim);
%
% sim.nf = number of orientations of lead field per voxel location (normally 1, 2, or 3)
% sim.numdip = number of desired source locations
% sim.numint = number of desired interference locations (interfering brain noise)
% sim.numtime = number of total time points for time course (prestim and poststim)
% sim.srate = sample rate in Hz
% sim.prestimpercent = percent of sim.numtime that is prestim time
% sim.hifreq = max high frequency of source time course (may be vector)
% sim.lofreq = min low frequency of source time course (may be vector)
% sim.lofreqinter
% sim.hifreqinter
% sim.sourcetype = 'gauss' (drawn from gaussian distribution), 'sin' (sinusoidal), 'laplacian' (drawn from laplacian), 'mog' (mixture of gaussians)
% sim.sourcedamping = 0 or 1, if dampen ends of time course with hanning window
% sim.inter = 'gauss' (drawn from gaussian distribution), 'sin' (sinusoidal), 'laplacian' (drawn from laplacian), 'none' or 'real' (no additional interference despite sim.numint)
% sim.repeat = pick between 0 through 9 for 10 diff variations of 'rand'
% sim.numtrial = number of trials
% sim.trialweight = matrix numtrialx(size(sources,1)), amplitude scaling 
%                           for each trial and for each source
% sim.trialphaseconst    % =1 same phase for brain source for each trial for sinusoid. =0 rand phase per trial
% 
% extracted from kenny's 'create_lead_field' which of course does more than
% just create a lead field.  

% user-defined constants for brain sources
% bs_distribution = 'Sinusoidal'; % use 'Sinusoidal' (which has bimodal, sub-Gaussian distribution), 'Laplacian', or 'MOG'
bs_Hanning = 1; % use 1 to apply Hanning window to brain sources, use 0 to remove pre-stimulus signal from brain sources
bs_correlated_envelope = 0; % use 0 to use random peak locations of Hanning window of sources of interest, 1 otherwise (not used when bs_Hanning is false)
% bs_corr = 0; % -1 < bs_corr < 1, use 0 for uncorrelated brain sources

% user-defined constants for brain noises
% bn_distribution = 'Sinusoidal'; % use 'Sinusoidal' (which has bimodal, sub-Gaussian distribution), 'Laplacian', or 'Gaussian'
%bn_distribution = 'Gaussian'; %%%%%%%%%%%%
% bn_corr = 0; % -1 < bn_corr < 1, use 0 for uncorrelated brain noises
numsource=sim.nf*sim.numdip;
numinter=sim.numint;
numtotal = numsource + numinter;
N=sim.numtime;

sim.seed.sourceloc=13;sim.seed.eta=17;sim.seed.sourcegauss=19;sim.seed.noise=23;sim.seed.interts=29;sim.seed.interloc=31;sim.seed.interorient=37;sim.seed.sourcefreqphase=41;sim.seed.interfreqphase=43;sim.seed.sourcewindow=47;
seednames=fieldnames(sim.seed);
for ii=1:length(seednames)
    seedval=getfield(sim.seed,char(seednames(ii)));
    sim.seed=setfield(sim.seed,char(seednames(ii)),seedval+sim.repeat);
end

fs=sim.srate/1000; 
T = N/fs; % total length in ms
t = ((-sim.prestimpercent*T):T/(N-1):(T*(1-sim.prestimpercent)))-.5*T/(N-1); % time in ms, 0 corresponds to stimulus onset
% t = ((-sim.prestimpercent*T):T/(N-1):(T*(1-sim.prestimpercent))); % time in ms, 0 corresponds to stimulus onset
pre = find(t < 0);
post = find(t >= 0);
Npre=length(pre);
Npost=length(post);
freqrange=sim.hifreq-sim.lofreq;

% creates brain sources
switch sim.sourcetype
    case 'gauss'
        if sim.seed.sourcegauss
            randn('state',sim.seed.sourcegauss);
        else
            randn('state',sum(100*clock));
        end
        tsources =[zeros(nf*sim.numdip,Npre) randn(nf*sim.numdip,Npost)];
    case 'sin'
        tsources = zeros(numsource,N,sim.numtrial);
        rand('state',sim.seed.sourcefreqphase);
        for ii = 1:numsource
            wtmp(ii) = 2*pi*(sim.lofreq(ii)+freqrange(ii)*rand); % frequencies between 5-30 Hz
            if sim.trialphaseconst
                phaseuse=2*pi*rand;
                tsources(ii,:,:) = repmat(sin(t*1e-3*wtmp(ii) + phaseuse)',1,sim.numtrial);
            else
                phaseuse=2*pi*rand(sim.numtrial,1);
                for jj=1:sim.numtrial
                    tsources(ii,:,jj) = sin(t*1e-3*wtmp(ii) + phaseuse(jj));
                end
            end
        end
    case 'laplacian'
        tsources = laplace(numsource,N);
    case 'mog'
        z = -3:0.1:3;
        nu = [50 100 1000];
        w = [0.4 0.2 0.4];
        mu = [0 0 0];
        for ii = 1:length(w)
            fz = w(ii)/(sqrt(2*pi)*1/nu(ii))*exp(-0.5*nu(ii)*(z-mu(ii)).^2);
        end
        fz = fz/sum(fz)/0.1;

        tsources = reshape(arb_pdf(fz,z,N*numsource),numsource,N);
    otherwise
        error('Invalid value for params.sourcetype')
end

% creates/applies Hanning window for brain sources
if sim.sourcedamping
   rand('state',sim.seed.sourcewindow);
   hanning_window = ones(numsource,N);
   if bs_correlated_envelope
%       source_peak = 2*T/8*rand*ones(1,numsource) + T/6; % single envelope for all brain sources (3T/8 + T/8)
      source_peak = 1000/sim.srate*(2*Npost/5*rand*ones(1,numsource) + 4*Npost/15); % single envelope for all brain sources (3T/8 + T/8)
   else
%       source_peak = 2*T/8*rand(1,numsource) + T/6; % must stay within post stimulus window (3T/8 + T/8)
      source_peak = 1000/sim.srate*(2*Npost/5*rand(1,numsource) + 4*Npost/15); % must stay within post stimulus window (3T/8 + T/8)      
   end
%    window_size = 0.5*ones(1,numsource);  % 0.5 of total N, when 5/8 post was fixed
   window_size = 0.8*ones(1,numsource);  % 0.8 of Npost
   for ii = 1:numsource
       [val,pos] = min(abs(t - source_peak(ii)));
       %       N2 = round(window_size(i)*N); % length of non-zero portion of window (relative to overall length)
       N2 = round(window_size(ii)*Npost); % length of non-zero portion of window (relative to poststim length)
       N2 = N2 + 1 - rem(N2,2); % make N2 an odd number
       ndx = pos - 0.5*(N2-1);
       tmp = hanning(N2)';

       if ndx < 1
           hanning_window(ii,:) = [tmp(2-ndx:end) zeros(1,N-N2-ndx+1)];
       elseif ndx > N-N2+1
           hanning_window(ii,:) = [zeros(1,ndx-1) tmp(1:N-ndx+1)];
       else
           hanning_window(ii,:) = [zeros(1,ndx-1) tmp zeros(1,N-N2-ndx+1)];
       end
   end
   tsources = tsources.*hanning_window;
else
    tsources(:,pre,:) = 0; % brain sources are active only in post-stimulus
end

for ii=1:sim.numtrial
    if size(tsources,3)==1
        sources(:,:,ii)=tsources.*repmat(sim.trialweight(ii,1:size(tsources,1))',[1 size(tsources,2)]);
    elseif size(tsources,3)==sim.numtrial
        sources(:,:,ii)=tsources(:,:,ii).*repmat(sim.trialweight(ii,1:size(tsources,1))',[1 size(tsources,2)]);
    end
end

% creates additional neural sources (brain noise)
switch sim.inter
    case 'sin'
        freqrangeinter=sim.hifreqinter-sim.lofreqinter;
        sources = [sources; zeros(numinter,N,sim.numtrial)];
        rand('state',sim.seed.interfreqphase);
        for ii = (numsource + 1):numtotal
            for jj=1:sim.numtrial
                sources(ii,:,jj) = sim.trialweight(jj,ii)*sin(t*1e-3*2*pi*(sim.lofreqinter(ii-numsource)+freqrangeinter(ii-numsource)*rand) + 2*pi*rand);
            end
        end
    case 'laplacian'
        sources = [sources; laplace(numinter,N,sim.numtrial)];
    case 'gauss'
        randn('state',sim.seed.interts);
        sources = [sources; randn(numinter,N,sim.numtrial)];
    case {'none','real'}
        numtotal=numsource;
    otherwise
        error('Invalid value for sim.inter')
end

for ii = 1:numtotal
   sources(ii,:,:) = squeeze(sources(ii,:,:)) - repmat(squeeze(mean(sources(ii,pre,:),2)),1,N)';
end

