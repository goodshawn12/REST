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

% cooling forgetting factor: lambda = lambda_0 / sample^gamma
function lambda = genCoolingFF(t,gamma,lambda_0)
    lambda = lambda_0 ./ (t .^ gamma);
end

