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
            + sum(evalConvergence.leakyAvgDeltaVar*bsxfun(@times,variance,(1-evalConvergence.leakyAvgDeltaVar).^(nPts-1:-1:0)),2);
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
