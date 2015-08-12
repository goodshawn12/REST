function [signal, state] = flt_orica(varargin)
% This function implemented an online recursive rule for ICA, ORICA[1].
% The recursive rule is a fixed point solution derived from infomax 
% learning rule. In another point of view, this function is an adaptive
% non-linear recursive least-square (RLS) filter.

% Author:   Sheng-Hsiou Hsu, 2013-2014, SCCN/INC/UCSD
%           Tim Mullen, 2013-2014, SCCN/INC/UCSD
% References: 
% [1] S.-H. Hsu, T. Mullen, T.-P Jung, and G. Cauwenberghs, "Online recursive independent component analysis for real-time source separation of high-density EEG," in IEEE EMBS, 2014.
% [2] Akhtar, Mu. T., Jung, T.-P., Makeig, S., & Cauwenberghs, G. (2012). Recursive independent component analysis for online blind source separation. 2012 IEEE International Symposium on Circuits and Systems, (6), 2813?2816. doi:10.1109/ISCAS.2012.6271896

if ~exp_beginfun('filter') return; end

%% define input arguments
% has its own highpass filter, sometimes applied on re-referenced data
declare_properties('name','ORICA', 'experimental',true,'precedes',{'flt_fir','flt_iir'}, 'follows',{'flt_reref'}, 'independent_trials',false, 'independent_channels',false);

arg_define([0 1],varargin,...
    arg_norep({'signal','Signal'}), ...
    arg_sub({'options','Parameters'},{}, ...
    { ...
        arg({'blockSize','BlockLen','BlockSize'},8,[1 Inf],'Block size for weights batch update. Guideline: if signal is relatively stationary increasing blockSize will speed up runtime without sacrificing too much performance'), ...
        arg({'numPass','NumOfPasses'},1,[0 20],'Number of passes over the same chunk data. ORICA obtains a new icaweight at each pass.'), ...
        arg({'timeperm','TimePerm'},true,[],'Shuffle data order'), ...
        arg({'nsub','NumSubgaussian'},0,[0 Inf],'Number of subgaussian sources in EEG signal. EEG brain sources are usually supergaussian. Subgaussian sources are motstly artifact or noise.'),...
    },'Options for implementation features.'), ...
    arg_sub({'adaptff','AdaptForgettingRate'},{},...
    {...
        arg({'gamma','FFDecayRate','Gamma'},0.6,[0 Inf],sprintf('Forgetting factor decay rate. \nThis applies to the learning rule: \nlambda = (lambda_0 / t ^ gamma) + lambda_ss. \nlarger gamma -> faster convergence but may be unstable. Gamma=0 will disable decay of forgetting factor.')), ...
        arg({'lambda_0','FFInitialValue'},0.995,[0 Inf],'Forgetting factor initial value.'), ...
        arg({'tau_ss','TauSteadyState'},Inf,[0 Inf],sprintf('Memory e-folding time at steady state (in sec). \nThis determines the steady-state value of the forgetting factor (lambda) as per the equation\ntau_ss = -1 / log(1-lambda_ss). \nIn this way the contribution to W_t of past weights W_k where (k = t - tau_ss) will be (1/e)*W_k. \nThis value should be chosen to be the approximate time window for local stationarity of the process/weights. Good choices for EEG are usually in the data_range of 1-3 seconds.')), ...
    },'Parameters for algorithm convergence.'), ...
    arg_subtoggle({'constff','ConstForgettingRate'},[],...
    {...
        arg({'tau_ss','SteadyStateTimeConst'},3,[0 Inf],sprintf('Memory e-folding time at steady state (in sec). \nThis determines the steady-state value of the forgetting factor (lambda) as per the equation\ntau_ss = -1 / log(1-lambda_ss). \nIn this way the contribution to W_t of past weights W_k where (k = t - tau_ss) will be (1/e)*W_k. \nThis value should be chosen to be the approximate time window for local stationarity of the process/weights. Good choices for EEG are usually in the data_range of 1-3 seconds.')), ...
    },'Parameters for algorithm convergence.'), ...
    arg_sub({'perfmetrics','PerformanceMetrics'},{},...
    {...
        arg_subtoggle({'convergence','Convergence'},[],...
        { ...
            arg({'A','TrueMixingMatrix'},'signal.icawinv_true',[],'Ground truth mixing maxing. Used only for performance checks','type','expression'), ...
            arg({'spheremat','SpheringMatrix'},'signal.icasphere',[],'Sphering matrix.','type','expression'), ...
        },'Compute convergence metric from Akhtar et al, 2012 [1]. Requires true mixing matrix to be known'), ...
    },'Options for performance metrics'), ...
    arg_subtoggle({'internalSphere','InternalSphere'},[],...
    {...
        arg({'blockSize','BlockSize'},[],[1 Inf],'Block size for whitening block update'), ...
        arg({'centerdata','CenterData'},true,[],'Remove data mean'), ...
        arg({'numpcs','NumPCs'},[],[],'Number of principal components to retain. If empty, default to number of channels'), ...
    }, 'Run RLS whitening in ORICA.'), ...
    arg({'computeConvergence','ComputeConvergence'},false,[],'Compute convergence properties.'), ...
    arg({'nlfunc','Nonlinearity'},[],[],'Nonlinear ICA function. Optional. This should be a lambda function. For instance @(y)= -2*tanh(y) corresponds to the hyperbolic secant supergaussian prior density assumed in Infomax ICA (default)','type','expression'), ...
    arg({'verb','Verbosity'},false,[],'Verbose output'), ...
    arg_norep({'state','State'},unassigned));
        
%% initialization of the filter
[nChs,nPts] = size(signal.data);
numpcs = nChs;
    
% initialize state structure (only for the first block)
if ~exist('state','var') || isempty(state)
    % learning related parameters:
    state.icaweights    = eye(numpcs);     % should use infinity for convergence
    state.gamma         = adaptff.gamma;         % store forgetting rate factor: adapt.gamma
    state.lambda_0      = adaptff.lambda_0;      % store forgetting rate factor: adapt.lambda_0
    state.lambda_ss    	= 1-exp(-1/(adaptff.tau_ss*signal.srate));  % steady state lambda
    
    if constff.arg_selection
        state.constLambda = 1-exp(-1/(constff.tau_ss*signal.srate));
    else
        state.constLambda = [];
    end
    
    if internalSphere.arg_selection
        state.icasphere = eye(numpcs,nChs); % FIXME
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
    state.lambda_k    	= zeros(1,options.blockSize);    % readout lambda
    state.statIdx       = 0;
    
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
    state = dynamicWhitening(signal.data, state, internalSphere);
    Mixtures = state.icasphere * signal.data;
    signal.icasphere = state.icasphere;
else
    if isempty(signal.icasphere); signal.icasphere = eye(numpcs,nChs);end
    Mixtures = signal.icasphere * signal.data;
end

% if ground truth is provided
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
for it = 1 : options.numPass
    
    if options.timeperm    % shuffling data to help eliminate temporal correlation
        permIdx = randperm(nPts);
    else
        permIdx = 1:nPts;
    end   
               
    % apply online recursive ICA algorithm, dynamicORICA(), to each block
    for bi = 0 : numsplits-1
        
        dataRange = 1 + floor(bi*nPts/numsplits) : min(nPts, floor((bi+1)*nPts/numsplits));
        state = dynamicOrica(Mixtures(:, permIdx(dataRange)), state, dataRange, computeConvergence);
        
        if perfmetrics.convergence.arg_selection
            % This uses performance metric in 2012 ORICA paper.
            H = state.icaweights * perfmetrics.convergence.spheremat * perfmetrics.convergence.A;
            C = H.^2;
            signal.convergence(dataRange) = (numpcs-sum(max(C,[],1)./sum(C,1))/2-sum(max(C,[],2)./sum(C,2))/2)/(numpcs-1);
        end
        
        if computeConvergence
            signal.lambda_k(dataRange) = state.lambda_k;
            signal.statIdx(dataRange) = state.statIdx;
        end
        
    end % for each block

    % increment paased sample counter
    tInf = 5e6; % Roughly the number of data of 4 hrs run with srate = 300Hz.
    if state.counter > tInf
        state.counter = Inf;
    else
        state.counter = state.counter + nPts;
    end
    
end % for each pass

% output
signal.icaact = state.icaweights * Mixtures;
signal.icaweights = state.icaweights;

% compute MIR
if computeConvergence
    % compute entropy for raw data and ica activation
    h0 = getent2(signal.data);
    h = getent2(signal.icaact);
    mir = sum(h0) - sum(h) + sum(log(abs(eig(state.icaweights*state.icasphere))));
    signal.mir = mir * ones(1,nPts);
end

exp_endfun;