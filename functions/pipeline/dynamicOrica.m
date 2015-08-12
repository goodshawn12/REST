function state = dynamicOrica(blockdata, state, dataRange, computeConvergence,nlfunc)

if nargin < 5
    computeConvergence = false; end
if nargin < 6
    nlfunc = []; end

% initialize
[nChs, nPts] = size(blockdata);
weights      = state.icaweights;
f            = zeros(nChs, nPts);

% block update using online recursive ICA rule
y = weights * blockdata;

% choose different nonlinear functions for super- and sub-gaussian 
if isempty(nlfunc)
    f(state.kurtsign,:)  = -2 * tanh(y(state.kurtsign,:));                        % Supergaussian
    f(~state.kurtsign,:) = tanh(y(~state.kurtsign,:)) - y(~state.kurtsign,:);     % Subgaussian
else
    f = nlfunc(y);
end

% define adaptive forgetting rate: lambda
lambda_k = state.lambda_ss + state.lambda_0 ./ ((state.counter+dataRange) .^ state.gamma);
if ~isempty(state.constLambda)
    if lambda_k < state.constLambda, lambda_k = state.constLambda; end
end
state.lambda_k = lambda_k;

% batch update learning rule
lambda_prod = prod(1./(1-lambda_k));
Q = 1 + lambda_k .* (dot(f,y,1)-1);
weights = lambda_prod * (weights - y * diag(lambda_k./Q) * f' * weights);

% orthogonalize weight matrix - avoid inverse of signular matrix
% [Q,R] = qr(weights);
% weights = Q*diag(sign(diag(R))); % ~6 times faster than svd method
[V,D] = eig(weights * weights');
weights = V/sqrt(D)*V' * weights; 

% output
state.icaweights = weights;

if computeConvergence
    yEval = weights * blockdata;
    if isempty(nlfunc)
        fEval(state.kurtsign,:)  = -2 * tanh(yEval(state.kurtsign,:));                           % Supergaussian
        fEval(~state.kurtsign,:) = tanh(yEval(~state.kurtsign,:)) - yEval(~state.kurtsign,:);    % Subgaussian
    else
        fEval = nlfunc(yEval);
    end
    nlcov = yEval * fEval' / nPts;
    state.statIdx = norm(nlcov - eye(nChs), 'fro');       
    
end

