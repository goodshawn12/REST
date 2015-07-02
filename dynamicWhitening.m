function state = dynamicWhitening(data, state, sphere)

% Usage: 
% Compute whitening (sphering) matrix in blockwise RLS fashion. Store
% whitening matrix in signal.icasphere and return whitened signal.data. Use
% adpative forgetting rate, i.e. exponential decay. 

% Reference:
% Zhu X.L., Zhang X.D, Ye J.I. "Natural gradient-based recursive
% least-sqrares algorithm for adaptive blind source separation", Science in
% China, 2004. 

%% Setup
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

    % define adaptive forgetting rate
    lambda = state.lambda_ss + state.lambda_0 ./ ((state.counter+dataRange) .^ state.gamma);
    if ~isempty(state.constLambda)
        if lambda < state.constLambda, lambda = state.constLambda; end
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