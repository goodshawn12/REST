function [signal, state] = flt_arica(varargin)
% This function implemented an online recursive rule for ICA, ORICA[1].
% The recursive rule is a fixed point solution derived from infomax 
% learning rule. In another point of view, this function is an adaptive
% non-linear recursive least-square (RLS) filter.

% Author:   Sheng-Hsiou Hsu, 2013-2014, SCCN/INC/UCSD
%           Tim Mullen, 2013-2014, SCCN/INC/UCSD
% References: 
% [1] Akhtar, Mu. T., Jung, T.-P., Makeig, S., & Cauwenberghs, G. (2012). Recursive independent component analysis for online blind source separation. 2012 IEEE International Symposium on Circuits and Systems, (6), 2813?2816. doi:10.1109/ISCAS.2012.6271896

if ~exp_beginfun('filter') return; end

%% define input arguments
% has its own highpass filter, sometimes applied on re-referenced data
declare_properties('name','ARICA','experimental',true,'precedes',{'flt_fir','flt_iir'}, 'follows',{'flt_reref'}, 'independent_trials',false, 'independent_channels',false);

arg_define([0 1],varargin,...
    arg_norep({'signal','Signal'}), ...
    arg_sub({'options','Parameters'},{}, ...
    { ...
        arg({'blockSize','BlockLen','BlockSize'},1,[1 Inf],'Block size for weights batch update. Guideline: if signal is relatively stationary increasing blockSize will speed up runtime without sacrificing too much performance'), ...
        arg({'numPass','NumOfPasses'},1,[0 20],'Number of passes over the same chunk data. ORICA obtains a new icaweight at each pass.'), ...
        arg({'timeperm','TimePerm'},true,[],'Shuffle data order'), ...
        arg({'nsub','NumSubgaussian'},0,[0 Inf],'Number of subgaussian sources in EEG signal. EEG brain sources are usually supergaussian. Subgaussian sources are motstly artifact or noise.'),...
    },'Options for implementation features.'), ...
    arg_subtoggle({'internalSphere','InternalSphere'},[],...
    {...
        arg({'blockSize','BlockSize'},[],[1 Inf],'Block size for whitening block update'), ...
        arg({'centerdata','CenterData'},true,[],'Remove data mean'), ...
        arg({'numpcs','NumPCs'},[],[],'Number of principal components to retain. If empty, default to number of channels'), ...
    }, 'Run RLS whitening in ORICA.'), ...
    arg_subswitch({'adaptiveFF','AdaptiveFF'},'nonadapt_cooling',{ ...
        'nonadapt_cooling', { ...
            arg({'gamma','FFDecayRate','Gamma'},0.6,[0 Inf],sprintf('Forgetting factor decay rate. \nThis applies to the learning rule: \nlambda = (lambda_0 / (t + t_0) ^ gamma) + lambda_ss. \nlarger gamma -> faster convergence but may be unstable. Gamma=0 will disable decay of forgetting factor.')), ...
            arg({'lambda_0','FFInitialValue'},0.995,[0 Inf],'Forgetting factor initial value.'), ...
            arg({'t_0','TimeOffset'},0,[0 Inf],sprintf('Initial time offset (in sec). Suggested values 0~1. Default: 0.')), ...
            arg({'tau_ss','TauSteadyState'},Inf,[0 Inf],sprintf('Memory e-folding time at steady state (in sec). \nThis determines the steady-state value of the forgetting factor (lambda) as per the equation\ntau_ss = -1 / log(1-lambda_ss). \nIn this way the contribution to W_t of past weights W_k where (k = t - tau_ss) will be (1/e)*W_k. \nThis value should be chosen to be the approximate time window for local stationarity of the process/weights. Good choices for EEG are usually in the data_range of 1-3 seconds.')), ...
            arg({'tau_const','ConstTauSteadyState'},Inf,[0 Inf],sprintf('Memory e-folding time at steady state (in sec). Set to this constant value when lambda decreases below the threshold.'))}, ...
        'nonadapt_const', { ...
            arg({'tau_const','ConstTauSteadyState'},3,[0 Inf],sprintf('Memory e-folding time at steady state (in sec). Set to this constant value when lambda decreases below the threshold.'))}, ...
        'nonadapt_exp', { ...
            arg({'gamma','ExpDecayRate','Gamma'},0.0384,[0 Inf],sprintf('Forgetting factor decay rate. \nThis applies to the exponentially attenuated FF:  \nlambda_0 * exp(-gamma*(t-t_0)), t>t_0. \nlarger gamma -> faster convergence but may be unstable.')), ...
            arg({'lambda_0','FFInitialValue'},0.995,[0 Inf],sprintf('Forgetting factor initial value.')), ...
            arg({'t_0','TimeOffset'},5,[0 Inf],sprintf('Initial time offset (in sec). Default: 5.'))}, ...
        'adapt_Murata_block', { ...
            arg({'decayRateAlpha','DecayRate'},0.02,[0 Inf],sprintf('Larger alpha leads to faster but unstable rate changes.')), ...
            arg({'upperBoundBeta','UpperBound'},0.005,[0 Inf],sprintf('To be appear....')), ...
            arg({'transBandWidthGamma','TransBandWidth'},0.5,[0 Inf],sprintf('To be appear')), ...
            arg({'transBandCenter','TransBandCenter'},3,[1 Inf],sprintf('To be appear')), ...
            arg({'leakyAvgDelta','LeakyAvg'},0.05,[0 Inf],sprintf('Parameter estimation window size is proportional to 1/delta. Update rule: Rn = (1-delta)*Rn + delta*flow')), ...
            arg({'lambdaInitial','InitialFF'},0.1,[0 Inf],sprintf('Forgetting factor initial value.'))}, ...
        'adaptiveComponent', { ...
            arg({'decayRateAlpha','DecayRate'},0.02,[0 Inf],sprintf('Larger alpha leads to faster but unstable rate changes.')), ...
            arg({'upperBoundBeta','UpperBound'},0.005,[0 Inf],sprintf('To be appear....')), ...
            arg({'transBandWidthGamma','TransBandWidth'},0.5,[0 Inf],sprintf('To be appear')), ...
            arg({'transBandCenter','TransBandCenter'},3,[1 Inf],sprintf('To be appear')), ...
            arg({'leakyAvgDelta','LeakyAvg'},0.05,[0 Inf],sprintf('Parameter estimation window size is proportional to 1/delta. Update rule: Rn = (1-delta)*Rn + delta*flow')), ...
            arg({'lambdaInitial','InitialFF'},0.1,[0 Inf],sprintf('Forgetting factor initial value.'))}, ...
        'adapt_Cichocki', { ...
            arg({'alpha','DynamicRate'},0.1,[0 Inf],sprintf('The global learning rule. Larger alpha leads to faster but unstable rate changes.')), ...
            arg({'deltaV','MemoryFactorV'},0.02,[0 Inf],sprintf('Smaller the value, longer the memory of V update.')), ...
            arg({'deltaU','MemoryFactorLambda'},0.02,[0 Inf],sprintf('Smaller the value, longer the memory of lambda update.')), ...
            arg({'lambda_0','FFInitialValue'},0.1,[0 Inf],sprintf('Forgetting factor initial value.'))}, ...
    },'Adaptive forgetting factor rules.'), ...
    arg_sub({'perfmetrics','PerformanceMetrics'},{},...
    {...
        arg_subtoggle({'convergence','Convergence'},[],...
        { ...
            arg({'A','TrueMixingMatrix'},'signal.icawinv_true',[],'Ground truth mixing maxing. Used only for performance checks','type','expression'), ...
            arg({'spheremat','SpheringMatrix'},'signal.icasphere',[],'Sphering matrix.','type','expression'), ...
        },'Compute convergence metric from Akhtar et al, 2012 [1]. Requires true mixing matrix to be known'), ...
    },'Options for performance metrics'), ...
    arg({'computeConvergence','ComputeConvergence'},false,[],'Compute convergence properties.'), ...
    arg({'storeweights','StoreWeights'},false,[],'Store icaweights for all blocks. Results are stores in signal.icaweights_all. This can quickly use up memory. Use with caution'), ...
    arg({'nlfunc','Nonlinearity'},[],[],'Nonlinear ICA function. Optional. This should be a lambda function. For instance @(y)= -2*tanh(y) corresponds to the hyperbolic secant supergaussian prior density assumed in Infomax ICA (default)','type','expression'), ...
    arg({'verb','Verbosity'},false,[],'Verbose output'), ...
    arg_norep({'state','State'},unassigned));

%             arg({'beta','GlobalRate'},1000,[0 Inf],sprintf('Minimum level is proportional to beta. Larger beta, the algorithm may track changes.')), ...

% Note:
%   skipFactor: should be implemented in separate function.1
%   WeightBlowup: MAX_WEIGHT, MAX_BLOWUP, GAMMA_FAC, LAMBDA_HAT_FAC
%   perfmetrics: testing only
        
%% initialization of the filter
timeOpt = 1; % timing the function 
if timeOpt; tic; end;

[nChs,nPts] = size(signal.data);
numpcs = nChs;

% if internalSphere.arg_selection
%     if ~isempty(internalSphere.numpcs)
%         numpcs = internalSphere.numpcs;
%     end
% end
    
% initialize state structure (only for the first block)
if ~exist('state','var') || isempty(state)
    % learning related parameters:
    state.icaweights    = eye(numpcs);     % should use infinity for convergence
    state.lambda_k    	= zeros(1,options.blockSize);    % readout lambda
    state.lambdaComp    = [];
    if internalSphere.arg_selection
%         if numpcs > nChs, error('Number of components to retain cannot exceed number of channels (%d)',nChs); end
        state.icasphere = eye(numpcs,nChs); % FIXME
    end

    switch adaptiveFF.arg_selection
        case 'nonadapt_cooling'
            state.gamma         = adaptiveFF.gamma;         % store forgetting rate factor: adapt.gamma
            state.lambda_0      = adaptiveFF.lambda_0;      % store forgetting rate factor: adapt.lambda_0
            state.t_0       	= adaptiveFF.t_0*signal.srate;           % Initial time offset (in sec)
            state.lambda_ss    	= 1-exp(-1/(adaptiveFF.tau_ss*signal.srate));  % steady state lambda
            state.lambda_const  = 1-exp(-1/(adaptiveFF.tau_const*signal.srate)); % steady state constant lambda
        case 'nonadapt_const'
            state.lambda_const  = 1-exp(-1/(adaptiveFF.tau_const*signal.srate));  % steady state lambda
        case 'nonadapt_exp'
            state.gamma         = adaptiveFF.gamma/signal.srate;         % store forgetting rate factor: adapt.gamma
            state.lambda_0      = adaptiveFF.lambda_0;      % store forgetting rate factor: adapt.lambda_0
            state.t_0       	= adaptiveFF.t_0*signal.srate;           % Initial time offset (0~25 samples)
        case 'adapt_Murata_block'
            state.decayRateAlpha      = adaptiveFF.decayRateAlpha;
            state.upperBoundBeta      = adaptiveFF.upperBoundBeta;
            state.transBandWidthGamma = adaptiveFF.transBandWidthGamma;
            state.transBandCenter     = adaptiveFF.transBandCenter;
            state.leakyAvgDelta       = adaptiveFF.leakyAvgDelta;
            state.lambda_k            = adaptiveFF.lambdaInitial;
        case 'adaptiveComponent'
            state.lambdaComp          = repmat(adaptiveFF.lambdaInitial,nChs,1);
            state.decayRateAlpha      = adaptiveFF.decayRateAlpha;
            state.upperBoundBeta      = adaptiveFF.upperBoundBeta;
            state.transBandWidthGamma = adaptiveFF.transBandWidthGamma;
            state.transBandCenter     = adaptiveFF.transBandCenter;
            state.leakyAvgDelta       = adaptiveFF.leakyAvgDelta;
        case 'adapt_Cichocki'
            state.lambdaComp    = repmat(adaptiveFF.lambda_0,nChs,1);
            state.V             = zeros(nChs);
            state.deltaV        = adaptiveFF.deltaV;
            state.deltaU        = adaptiveFF.deltaU;
            state.alpha         = adaptiveFF.alpha;
    end
    
    if computeConvergence
        if ~isfield(state,'leakyAvgDelta'), state.leakyAvgDelta = 0.05; end
        state.Rn            = [];
        state.normRn        = [];
        state.minNormRn     = [];
        state.ratioOfNormRn = [];
        state.columnWiseNormRn = [];
        state.mincolumnWiseNormRn = [];
        state.columnWiseNormRatio = [];
        state.RnW = zeros(numpcs);
        state.normRnW = [];
        state.lambdaGrad = [];
        state.covyy = eye(numpcs);
        state.covfy = eye(numpcs);
        state.normCovyy = [];
        state.normCovfy = [];
        state.yEst = zeros(numpcs,1);
        state.fEst = zeros(numpcs,1);
    end
    
    % sign of kurtosis for each component: true(supergaussian), false(subgaussian)
    state.kurtsign      = ones(numpcs,1) > 0;      % store kurtosis sign for each channels
    if options.nsub ~= 0
        state.kurtsign(1:options.nsub) = false;
    end 
    
    % implementation parameters
    state.counter            = 0; % time index counter, used to keep track of time for computing lambda
    state.buffer             = zeros(nChs,options.blockSize-1);
    state.bufferIdx          = 0;
    
            % testing only
            state.convMat       = zeros(numpcs);      % store convergence matrix for sanity check
            state.statIdx       = 0;
            state.statIdxMax    = 0;
            state.time_offset   = 0;    % time_offset for annealFF option
%             if annealFF
%                 state.meanIdx       = 0;
%                 state.covIdx        = 0;
%             end
end

% handle special case of no data
if nPts == 0
    fprintf('No input data.\n');
    return;
end

% handle case of blockSize larger than number of samples 
if nPts < options.blockSize || state.bufferIdx > 0
    % concatenate buffer data to signal.data
    if state.bufferIdx + nPts < options.blockSize
        state.buffer(:,(state.bufferIdx+1):(state.bufferIdx+nPts)) = signal.data;
        state.bufferIdx = state.bufferIdx + nPts;
        fprintf('Concatenate data.\n');
        return;
    else
        signal.data = [state.buffer(:,1:state.bufferIdx) signal.data];
        state.bufferIdx = 0;
        [nChs,nPts] = size(signal.data);
        numpcs = nChs;
        fprintf('Concatenate and process data.\n');
    end
end

% divide chunk data into blocks for batch update
numsplits = floor(nPts/options.blockSize);

%% online RLS whitening - whitening inside the filter
if internalSphere.arg_selection
    if isempty(internalSphere.blockSize); internalSphere.blockSize = options.blockSize; end
    state = adaptiveWhitening(signal.data, state, options, internalSphere, adaptiveFF);
    Mixtures = state.icasphere * signal.data;
    signal.icasphere = state.icasphere;
    if storeweights; signal.sphere_all = signal.icasphere; end
    if timeOpt; signal.executeTimeWhite = toc; end
else
    if isempty(signal.icasphere); signal.icasphere = eye(numpcs,nChs);end
    Mixtures = signal.icasphere * signal.data;
    if storeweights; signal.sphere_all = signal.icasphere; end
end

% handle sphering matrix
        % testing only
        if perfmetrics.convergence.arg_selection
            if ischar(perfmetrics.convergence.A)
                perfmetrics.convergence.A = eval(perfmetrics.convergence.A);
            end
            if isempty(perfmetrics.convergence.A)
                fprintf('You must provide a mixing matrix for convergence performance metric\n');
                perfmetrics.convergence.arg_selection = false;
            end
            if ischar(perfmetrics.convergence.spheremat)
                perfmetrics.convergence.spheremat = eval(perfmetrics.convergence.spheremat);
            end
            if isempty(perfmetrics.convergence.spheremat)
                perfmetrics.convergence.spheremat = eye(numpcs,nChs);
            end
            if internalSphere.arg_selection
                perfmetrics.convergence.spheremat = state.icasphere;
            end
        end    
    
%% online recursive ICA 

% initialize running sum, output averaged weigths
% sumWeights = zeros(size(state.icaweights));
detection = 0; % detection of non-stationary

% multiple passes of current chunk data
for it = 1 : options.numPass
    
    % shuffling data to help eliminate temporal correlation
    if options.timeperm
        permIdx = randperm(nPts);
    else
        permIdx = 1:nPts;
    end   
               
    % apply online recursive ICA algorithm, dynamicORICA(), to each block
    for bi = 0 : numsplits-1
        
        % define dataRange for each block
        dataRange = 1 + floor(bi*nPts/numsplits) : min(nPts, floor((bi+1)*nPts/numsplits));

        % call dynamicORICA() function for block-wise ORICA
        state = adaptiveOrica(Mixtures(:, permIdx(dataRange)), state, options, dataRange, ...
                            adaptiveFF, computeConvergence);
    
%         % running sum
%         sumWeights = sumWeights + state.icaweights;
            
                % Compute and store performace index (error).
                if perfmetrics.convergence.arg_selection
                    % only for testing changing mixing matrix
                    if (state.counter+state.time_offset+dataRange(1)) > 601*300 && (state.counter+state.time_offset+dataRange(end))  < 1201*300
                        if isfield(signal,'etc')
                            perfmetrics.convergence.A = signal.etc.LFM{2};
                        end
                    elseif (state.counter+state.time_offset+dataRange(1)) > 1201*300 && (state.counter+state.time_offset+dataRange(end))  < 1801*300
                        if isfield(signal,'etc')
                            perfmetrics.convergence.A = signal.etc.LFM{3};
                        end
                    end
                    % This uses performance metric in 2012 ORICA paper.
                    H = state.icaweights * perfmetrics.convergence.spheremat * perfmetrics.convergence.A;
                    C = H.^2;
                    signal.convergence(dataRange) = (numpcs-sum(max(C,[],1)./sum(C,1))/2-sum(max(C,[],2)./sum(C,2))/2)/(numpcs-1);                    
                    % signal.convergence(state.counter+dataRange) = (numpcs-sum(max(C,[],1)./sum(C,1))/2-sum(max(C,[],2)./sum(C,2))/2)/(numpcs-1);
                end
    
                % store convergence matrix for sanity check
                if computeConvergence
                    % store stationary index
                    signal.statIdx(dataRange) = state.statIdx;
                    if state.statIdx > state.statIdxMax; state.statIdxMax = state.statIdx; end 
                    
                    % store lambda information
                    if isempty(state.lambdaComp)
                        signal.lambda_k(dataRange) = state.lambda_k;
                    else
                        signal.lambda_k = state.lambdaComp;
                    end
                    
                    signal.normRn(dataRange) = repmat(state.normRn,1,length(dataRange));
                    signal.lambdaGrad(dataRange) = repmat(state.lambdaGrad,1,length(dataRange));
                    signal.normRnW(dataRange) = repmat(state.normRnW,1,length(dataRange));
                    signal.normCovyy(dataRange) = repmat(state.normCovyy,1,length(dataRange));
                    signal.normCovfy(dataRange) = repmat(state.normCovfy,1,length(dataRange));
                    signal.covyy = state.covyy;
                    signal.covfy = state.covfy;
                    signal.RnW = state.RnW;
                    signal.minNormRn(dataRange) = repmat(state.minNormRn,1,length(dataRange));
                    signal.ratioOfNormRn(dataRange) = repmat(state.ratioOfNormRn,1,length(dataRange));
                    signal.columnWiseNormRn = state.columnWiseNormRn;
                end

    end % for each block

    % increment paased sample counter
    tInf = 5e6; % Roughly the number of data of 4 hrs run with srate = 300Hz.
    if state.counter > tInf
        state.counter = Inf;
    else
        if detection == 0
            state.counter = state.counter + nPts;
        else
            state.counter = state.counter + nPts - detection;
            detection = 0;
        end
    end

    % testing only
    if storeweights % only store the last result
        signal.icaweights_all = state.icaweights;
        signal.time_all = state.counter + state.time_offset;
    end

end % for each pass

%% output
        
% compute ica activation and store icaweights
signal.icaact = state.icaweights * Mixtures;
signal.icaweights = state.icaweights;

% should think about how to implement averaged icaweights, averaged over
% signal pass or all passes?
% state.icaweights = sumWeights / options.numPass / numsplits;

%% compute MIR
if computeConvergence
    % compute entropy for raw data and ica activation
    h0 = getent2(signal.data);
    h = getent2(signal.icaact);
    mir = sum(h0) - sum(h) + sum(log(abs(eig(state.icaweights*state.icasphere))));
    signal.mir = mir * ones(1,nPts);
end


if timeOpt
   signal.sampleSize = nPts;
   signal.executeTime = toc - signal.executeTimeWhite;
end

exp_endfun;
