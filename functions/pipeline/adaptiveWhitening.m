function state = adaptiveWhitening(data, state, options, sphere, adaptiveFF)

% Usage: 
% Compute whitening (sphering) matrix in blockwise RLS fashion. Store
% whitening matrix in signal.icasphere and return whitened signal.data. Use
% adpative forgetting rate, i.e. exponential decay. 

% Reference:
% Zhu X.L., Zhang X.D, Ye J.I. "Natural gradient-based recursive
% least-sqrares algorithm for adaptive blind source separation", Science in
% China, 2004. 

%% Setup
if nargin < 5
    adaptiveFF.arg_selection = 'nonadapt_cooling'; end
Mixtures    = data;
[nChs,nPts] = size(Mixtures);

% divide chunk data into blocks for batch update
numsplits    = floor(nPts/sphere.blockSize);
timeperm = 1:nPts; % randperm(nPts);

if sphere.centerdata
    rowmeans = mean(Mixtures,2);
    Mixtures = bsxfun(@minus,Mixtures,rowmeans);
end

%% RLS whitening - block update for icasphere matrix
for bi = 0 : numsplits-1
    % define data_range for each block: evenly divide nPts to 'num_block' blocks
    dataRange = 1+floor(bi*nPts/numsplits) : min(nPts,floor((bi+1)*nPts/numsplits));
    nPtBlk = length(dataRange);

    % define adaptive forgetting rate: lambda
    switch adaptiveFF.arg_selection
        case 'nonadapt_cooling' % cooling FF: lambda_0 / (t+t_0)^gamma + lambda_ss
            lambda = nonadapt_cooling(state.counter+dataRange, state.gamma, state.lambda_0, state.t_0, state.lambda_ss);
            if lambda(1) < state.lambda_const, lambda = repmat(state.lambda_const,1,nPts); end
        case 'nonadapt_const' % constant FF
            lambda = repmat(state.lambda_const,1,nPts);
        case 'nonadapt_exp' % exponentially attenuated FF:  lambda_0 * exp(-gamma*(t-t_0)), t>t_0
            lambda = nonadapt_exp(state.counter+dataRange,state.lambda_0,state.gamma,state.t_0);
        case 'adapt_Murata_block' % Murata adaptive FF:
            lambda = repmat(state.lambda_k(end),1,nPts); % using previous adaptive lambda_k from adaptiveOrica
        case 'adaptiveComponent' % Murata adaptive FF:
            lambda = repmat(mean(state.lambdaComp),1,nPts); % using previous adaptive lambda_k from adaptiveOrica
        case 'adapt_Cichocki'
            lambda = repmat(mean(state.lambdaComp),1,nPts);
        otherwise
            disp('Method not specified.')
    end
        
    % blockwise whitening
    % ------------------------------------------------------------------------
    xPreWhite = state.icasphere * Mixtures(:, timeperm(dataRange)); % pre-whitened data 
    lambda_avg = 1 - lambda(ceil(end/2));    % median lambda

    % adaptive RLS whitening algorithm
    QWhite = lambda_avg/(1-lambda_avg) + trace(xPreWhite' * xPreWhite) / nPtBlk;
    state.icasphere = 1/lambda_avg * (state.icasphere - xPreWhite * xPreWhite' / nPtBlk / QWhite * state.icasphere);

end

end

% non-adaptive cooling forgetting factor 
% lambda = lambda_0 / sample^gamma
function lambda = nonadapt_cooling(t,gamma,lambda_0,t_0,lambda_ss)
% t: time
% gamma: decay rate
% lambda_0: initial value
% t_0: initial time offset (usually 0)
% lambda_ss: asymptotic value of the function
% Tim Mullen, SCCN/INC UCSD 2013
    lambda = lambda_0 ./ ((t + t_0) .^ gamma) + lambda_ss; % forgetting rate
end

% exponentially attenuated FF
function lambda = nonadapt_exp(t,lambda_0,Td,T0)
    if t(1)<=T0
        lambda = repmat(lambda_0,1,length(t));
    else
        lambda = lambda_0 * exp(-Td*(t-T0));
    end
end

