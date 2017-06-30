function [signal, state] = flt_orica(varargin)
% This function implements Online Recursive Indepnedent Component Analysis (ORICA) as a plug-in function compatible with BCILAB for online processing of live data stream.
% The ORICA processing pipeline (including online RLS whitening) was propsoed by [1], with the ORICA algorithm adapted from [2].
% It also implements three strategies for setting up the forgetting factors, i.e. the weight applied to new data in the recursive update, as described in [3].
%
% Advanced learning options:
%       - Adaptive Forgetting Factor : 
%           1. 'cooling': lambda = (lambda_0 / t ^ gamma). Larger gamma corresponds to faster convergence speed but unstable learning.
%           2. 'constant': Memory e-folding time at steady state (in sec). Set to this constant value when lambda decreases below the threshold.
%           3. 'adaptive': lambda = lambda - DecayRate*lambda + UpperBound*Gain*lambda^2, Gain(z) ~ tanh((NSI/NSI_{min} - TransBandCenter) / TransBandWidth).
%
% Author:   Sheng-Hsiou (Shawn) Hsu, 2014-15 SCCN/INC/UCSD
%           Tim Mullen, 2014, SCCN/INC/UCSD
%
% References: 
% [1] S-H. Hsu, T. Mullen, T-P Jung, and G. Cauwenberghs, "Online recursive independent component analysis for real-time source separation of high-density EEG," in IEEE EMBS, 2014.
% [2] Akhtar, Mu. T., Jung, T.-P., Makeig, S., & Cauwenberghs, G. (2012). Recursive independent component analysis for online blind source separation. 2012 IEEE International Symposium on Circuits and Systems, (6), 2813?2816. doi:10.1109/ISCAS.2012.6271896
% [3] S-H. Hsu, L. Pion-Tonachini, T-P. Jung, and G. Cauwenberghs, "Tracking non-stationary EEG sources using adaptive online recursive independent component analysis," in IEEE EMBS, 2015.

if ~exp_beginfun('filter') return; end

%% define input arguments
declare_properties('name','ORICA','experimental',true,'precedes',{'flt_fir','flt_iir'}, 'follows',{'flt_reref'}, 'independent_trials',false, 'independent_channels',false);

arg_define([0 1],varargin,...
    arg_norep({'signal','Signal'}), ...
    arg_subtoggle({'onlineWhitening','OnlineWhitening'},[],...
    {...
        arg({'blockSize','BlockSize'},8,[1 Inf],'L_{white}: Block size for online whitening block update. Suggested value: 4 to 8, depending on the noise of the data. Larger the noise, larger the block size.'), ...
        arg({'centerdata','CenterData'},true,[],'Remove data mean.'), ...
    }, 'Run online RLS whitening prior to ORICA. Suggested for online processing.'), ...
    arg_sub({'options','OnlineICA'},{}, ...
    { ...
        arg({'blockSize','BlockLen','BlockSize'},1,[1 Inf],'L_{ICA}: Block size for online ICA block update. Guideline: if signal is relatively stationary increasing blockSize will speed up runtime without sacrificing too much performance.'), ...
        arg({'nsub','NumSubgaussian'},0,[0 Inf],'Number of subgaussian sources in EEG signal. EEG brain sources are usually supergaussian. Subgaussian sources are motstly artifact or noise.'),...
        arg({'timeperm','TimePerm'},true,[],'Shuffle data order to reduce temporal correlation.'), ...
    },'Parameters and options for online recursive ICA.'), ...
    arg_subswitch({'adaptiveFF','AdaptiveFF'},'cooling',{ ...
        'cooling', { ...
            arg({'gamma','FFDecayRate','Gamma'},0.6,[0 Inf],sprintf('Forgetting factor decay rate.')), ...
            arg({'lambda_0','FFInitialValue'},0.995,[0 Inf],'Forgetting factor initial value.'), ...
            arg({'tau_const','ConstTauAtSteadyState'},3,[0 Inf],sprintf('Memory e-folding time at steady state (in sec). Set to this constant value when lambda decreases below the threshold.'))}, ...
        'constant', { ...
            arg({'tau_const','ConstTauSteadyState'},3,[0 Inf],sprintf('Memory e-folding time at steady state (in sec). Set to this constant value when lambda decreases below the threshold.'))}, ...
        'adaptive', { ...
            arg({'decayRateAlpha','DecayRate'},0.02,[0 Inf],sprintf('Decay rate')), ...
            arg({'upperBoundBeta','UpperBound'},0.001,[0 Inf],sprintf('Upper bound for lambda changes')), ...
            arg({'transBandWidthGamma','TransBandWidth'},1,[0 Inf],sprintf('Transition band width')), ...
            arg({'transBandCenter','TransBandCenter'},5,[1 Inf],sprintf('Transition band center')), ...
            arg({'lambdaInitial','InitialFF'},0.1,[0 1],sprintf('Forgetting factor initial value.'))}, ...
    },'Choices for forgetting factors.'), ...
    arg_subtoggle({'evalConvergence','EvalConvergence'},[],...
        {...
            arg({'leakyAvgDelta','LeakyAvgNSI'},0.01,[0 1],sprintf('Leaky average NSI')), ...
            arg({'leakyAvgDeltaVar','LeakyAvgVar'},1e-3,[0 1],sprintf('Leaky average variance of source activity')), ...
        }, 'Evaluate convergence such as Non-Stationarity Index (NSI).'), ...
    arg({'nlfunc','Nonlinearity'},[],[],'Nonlinear ICA function. Optional. For instance @(y)= -2*tanh(y) corresponds to the hyperbolic secant supergaussian prior density assumed in Infomax ICA (default)','type','expression'), ...
    arg_norep({'state','State'},unassigned));
        
%% initialize the filter
[nChs,nPts] = size(signal.data);
    
% handle special case of no data
if nPts == 0
    fprintf('No input data.\n');
    return;
end

% initialize state structure (only the first time)
if ~exist('state','var') || isempty(state)

    state.icaweights    = eye(nChs);                    % should use infinity for convergence
    state.lambda_k    	= zeros(1,options.blockSize);   % readout lambda

    % buffer for handling online data
    state.buffer             = zeros(nChs,options.blockSize-1);
    state.bufferIdx          = 0;

    if onlineWhitening.arg_selection
        state.icasphere = eye(nChs);
    end

    switch adaptiveFF.arg_selection
        case 'cooling'
            state.gamma         = adaptiveFF.gamma;         % store forgetting rate factor: adapt.gamma
            state.lambda_0      = adaptiveFF.lambda_0;      % store forgetting rate factor: adapt.lambda_0
            state.lambda_const  = 1-exp(-1/(adaptiveFF.tau_const*signal.srate)); % steady state constant lambda
            state.counter       = 0; % time index counter, used to keep track of time for computing lambda
        case 'constant'
            state.lambda_const  = 1-exp(-1/(adaptiveFF.tau_const*signal.srate));  % steady state lambda
        case 'adaptive'
            state.decayRateAlpha      = adaptiveFF.decayRateAlpha;
            state.upperBoundBeta      = adaptiveFF.upperBoundBeta;
            state.transBandWidthGamma = adaptiveFF.transBandWidthGamma;
            state.transBandCenter     = adaptiveFF.transBandCenter;
            state.lambda_k            = adaptiveFF.lambdaInitial;
            state.minNormRn     = []; 
        otherwise
            disp('Method not specified.')
    end
    
    if evalConvergence.arg_selection
        state.Var           = [];
        state.Rn            = [];
        state.normRn        = [];
    end
    
    % sign of kurtosis for each component: true(supergaussian), false(subgaussian)
    state.kurtsign      = ones(nChs,1) > 0;      % store kurtosis sign for each channels
    if options.nsub ~= 0
        state.kurtsign(1:options.nsub) = false;
    end 
    
end


% handle case of blockSize larger than number of samples 
if nPts < options.blockSize || state.bufferIdx > 0
    % concatenate buffer data to signal.data
    if state.bufferIdx + nPts < options.blockSize
        state.buffer(:,(state.bufferIdx+1):(state.bufferIdx+nPts)) = signal.data;
        state.bufferIdx = state.bufferIdx + nPts;
        fprintf('Concatenate data.\n');
        signal.icaact = state.icaweights * state.icasphere * signal.data;
        return;
    else
        signal.data = [state.buffer(:,1:state.bufferIdx) signal.data];
        state.bufferIdx = 0;
        [nChs,nPts] = size(signal.data);
        fprintf('Concatenate and process data.\n');
    end
end

%% online RLS whitening
if onlineWhitening.arg_selection
    state = dynamicWhitening(signal.data, state, options, onlineWhitening, adaptiveFF);
else
    if isempty(signal.icasphere)
        state.icasphere = eye(nChs);
    else
        state.icasphere = signal.icasphere;
    end
end

% whiten / sphere the data
Mixtures = state.icasphere * signal.data;
    
%% online recursive ICA 
if options.timeperm % shuffling data to help eliminate temporal correlation
    permIdx = randperm(nPts);
else
    permIdx = 1:nPts;
end   

% divide chunk data into blocks for batch update
% note: chunk is the size of online buffer while block is the user-defined size for block-update ICA. Usually chunk > block.
numsplits = floor(nPts/options.blockSize);

% apply online recursive ICA algorithm, dynamicORICA(), to each block with size options.blockSize
for bi = 0 : numsplits-1
    
    dataRange = 1 + floor(bi*nPts/numsplits) : min(nPts, floor((bi+1)*nPts/numsplits));

    % update weight matrix using online recursive ICA block update rule
    state = dynamicOrica(Mixtures(:, permIdx(dataRange)), state, options, dataRange, ...
                        adaptiveFF, evalConvergence);

    % store lambda information
    signal.lambda_k(dataRange) = state.lambda_k;

    % store convergence matrix for sanity check
    if evalConvergence.arg_selection
        signal.normRn(dataRange) = repmat(state.normRn,1,length(dataRange));
    end

end % for each block

% compute ica activation and store icaweights
signal.icaact = state.icaweights * Mixtures;
signal.icaweights = state.icaweights;
signal.icasphere = state.icasphere;

exp_endfun;
end


function state = dynamicWhitening(blockdata, state, options, sphere, adaptiveFF)

% Usage: 
% Compute whitening (sphering) matrix in blockwise RLS fashion. Store
% whitening matrix in signal.icasphere and return whitened signal.data. Use
% adpative forgetting rate, i.e. exponential decay. 

% Reference:
% Zhu X.L., Zhang X.D, Ye J.I. "Natural gradient-based recursive
% least-sqrares algorithm for adaptive blind source separation", Science in
% China, 2004. 

% handle input
if nargin < 5
    adaptiveFF.arg_selection = 'cooling'; end

% initialize
[nChs,nPts] = size(blockdata);

% remove data mean
if sphere.centerdata
    rowmeans = mean(blockdata,2);
    blockdata = bsxfun(@minus,blockdata,rowmeans);
end

% divide chunk data into blocks for batch update
numsplits    = floor(nPts/sphere.blockSize);

% online RLS whitening - block update for sphere matrix
for bi = 0 : numsplits-1

    dataRange = 1+floor(bi*nPts/numsplits) : min(nPts,floor((bi+1)*nPts/numsplits));

    % define adaptive forgetting rate: lambda
    switch adaptiveFF.arg_selection
        case 'cooling' % cooling FF: lambda_0 / (t+t_0)^gamma + lambda_ss
            lambda = genCoolingFF(state.counter+dataRange, state.gamma, state.lambda_0);
            if lambda(1) < state.lambda_const
                lambda = repmat(state.lambda_const,1,nPts); 
            end
        case 'constant' % constant FF
            lambda = repmat(state.lambda_const,1,nPts);
        case 'adaptive' % Murata adaptive FF:
            lambda = repmat(state.lambda_k(end),1,nPts); % using previous adaptive lambda_k from adaptiveOrica
        otherwise
            disp('Method not specified.')
    end
        
    % update sphere matrix using online RLS whitening block update rule
    xPreWhite = state.icasphere * blockdata(:, dataRange); % pre-whitened data 
    lambda_avg = 1 - lambda(ceil(end/2));    % median lambda
    QWhite = lambda_avg/(1-lambda_avg) + trace(xPreWhite' * xPreWhite) / length(dataRange);
    state.icasphere = 1/lambda_avg * (state.icasphere - xPreWhite * xPreWhite' / length(dataRange) / QWhite * state.icasphere);

end

end


function state = dynamicOrica(blockdata, state, options, dataRange, adaptiveFF, evalConvergence,nlfunc)
% This is an auxiliary function for flt_orica, implementing Online Recursive Independent Component Analysis. 
% It contains online ICA weight update with options (constant, cooling, and adaptive) for the forgetting factors.

% handle input
if nargin < 5
    adaptiveFF.arg_selection = 'cooling'; end
if nargin < 6
    evalConvergence = true; end
if nargin < 7
    nlfunc = []; end

% initialize
[nChs, nPts] = size(blockdata);
f            = zeros(nChs, nPts);

% compute source activation using previous weight matrix
y = state.icaweights * blockdata;

% choose nonlinear functions for super- vs. sub-gaussian 
if isempty(nlfunc)
    f(state.kurtsign,:)  = -2 * tanh(y(state.kurtsign,:));                        % Supergaussian
    f(~state.kurtsign,:) = tanh(y(~state.kurtsign,:)) - y(~state.kurtsign,:);     % Subgaussian
else
    f = nlfunc(y);
end

% compute Non-Stationarity Index (normRn) and variance of source dynamics (Var)
if evalConvergence.arg_selection
    modelFitness = eye(nChs)+y*f'/nPts;
    variance = blockdata.*blockdata;
    if isempty(state.Rn)
        state.Rn = modelFitness;
        state.Var = sum(variance,2)/(nPts-1);
    else
        state.Rn = (1-evalConvergence.leakyAvgDelta)*state.Rn + evalConvergence.leakyAvgDelta*modelFitness; % !!! this does not account for block update!
        state.Var = (1-evalConvergence.leakyAvgDeltaVar)^nPts*state.Var ...
            + sum(evalConvergence.leakyAvgDeltaVar*bsxfun(@times,variance, ...
            (1-evalConvergence.leakyAvgDeltaVar).^(nPts-1:-1:0)),2);
    end
    state.normRn = norm(state.Rn,'fro');
end

% compute the forgetting rate
switch adaptiveFF.arg_selection
    case 'cooling'
        state.lambda_k = genCoolingFF(state.counter+dataRange, state.gamma, state.lambda_0);
        if state.lambda_k(1) < state.lambda_const
            state.lambda_k = repmat(state.lambda_const,1,nPts); 
        end
        state.counter = state.counter + nPts;
    case 'constant'
        state.lambda_k = repmat(state.lambda_const,1,nPts);
    case 'adaptive'
        if isempty(state.minNormRn)
            state.minNormRn = state.normRn; 
        end
        state.minNormRn = max(min(state.minNormRn, state.normRn),1);
        ratioOfNormRn = state.normRn/state.minNormRn;
        state.lambda_k = genAdaptiveFF(dataRange,state.lambda_k,state.decayRateAlpha,state.upperBoundBeta,state.transBandWidthGamma,state.transBandCenter,ratioOfNormRn);
    otherwise
        disp('Method not specified.')
end

% update weight matrix using online recursive ICA block update rule
lambda_prod = prod(1./(1-state.lambda_k));
Q = 1 + state.lambda_k .* (dot(f,y,1)-1);
state.icaweights = lambda_prod * (state.icaweights - y * diag(state.lambda_k./Q) * f' * state.icaweights);

% orthogonalize weight matrix 
[V,D] = eig(state.icaweights * state.icaweights');
state.icaweights = V/sqrt(D)/V * state.icaweights; 

end


% cooling forgetting factor: 
% lambda = lambda_0 / sample^gamma
function lambda = genCoolingFF(t,gamma,lambda_0)
    lambda = lambda_0 ./ (t .^ gamma);
end


% adaptive forgetting factor:
% lambda = lambda - DecayRate*lambda + UpperBound*Gain*lambda^2
% Gain(z) ~ tanh((z/z_{min} - TransBandCenter) / TransBandWidth)
function lambda = genAdaptiveFF(dataRange,lambda,decayRateAlpha,upperBoundBeta,transBandWidthGamma,transBandCenter,ratioOfNormRn)
    gainForErrors = upperBoundBeta*0.5*(1+tanh((ratioOfNormRn-transBandCenter)/transBandWidthGamma));
    f = @(n) (1+gainForErrors).^n * lambda(end) - decayRateAlpha*((1+gainForErrors).^(2*n-1)-(1+gainForErrors).^(n-1))/gainForErrors*lambda(end)^2;
    lambda = f(1:length(dataRange));    
end



