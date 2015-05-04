function state = dynamicOrica(blockdata, state, options, dataRange, annealFF, computeConvergence,nlfunc)
%% Description
% To be completed...
% assume weight matrix is stationary in blockSize samples
% Since W is associated with components localization, it's relatively
% stationary (at least in a few hundred milliseconds)
% 1. Learning rule for weights matrix in block-based fasion:
%   sample-based learning rule: W = 1/(1-lambda) * (W - lambda/(1+lambda*(f^T*y-1))*y*f^T*W)
%   block-based learning rule: W = cumprod(1/(1-lambda)) * (W - sum(lambda/(1+lambda*(f^T*y-1)*y*f^T*W))

% handle input
if nargin < 5
    annealFF = false; end
if nargin < 6
    computeConvergence = false; end
if nargin < 7
    nlfunc = []; end
% ------------------------------------------------------------------------
% initialize
[nChs, nPts] = size(blockdata);
weights      = state.icaweights;
f            = zeros(nChs, nPts);

% ------------------------------------------------------------------------
% block update using online recursive ICA rule
y = weights * blockdata;

% choose different nonlinear functions for super- and sub-gaussian 
if isempty(nlfunc)
    f(state.kurtsign,:)  = -2 * tanh(y(state.kurtsign,:));                        % Supergaussian
    f(~state.kurtsign,:) = tanh(y(~state.kurtsign,:)) - y(~state.kurtsign,:);     % Subgaussian
else
    f = nlfunc(y);
end

% define forgetting rate: lambda
lambda_k = exp_decay(state.counter+dataRange, state.gamma, state.lambda_0, state.t_0, state.lambda_ss);
state.lambda_k = lambda_k;

% version 1:
% batch update learning rule
lambda_prod = prod(1./(1-lambda_k));
Q = 1 + lambda_k .* (dot(f,y,1)-1);
weights = lambda_prod * (weights - y * diag(lambda_k./Q) * f' * weights);

% version 2: doesn't work
% From online RLS whitening matrix
% lambda = 1 - lambda_k(ceil(end/2));
% Q = lambda / (1-lambda) + trace(f'*y) / nPts;
% weights = 1/lambda * (weights - y * f' / nPts / Q * weights);

% key point:
%   version one doesn't average when computing inner and outer product!!

% orthogonalize weight matrix - avoid inverse of signular matrix
[V,D] = eig(weights * weights');
weights = V/sqrt(D)*V' * weights;    % runs 70% faster
% weights = V/sqrt(D)/V * weights;    % run 30% faster
% weights = real(inv(V * sqrt(D) / V)) * weights;
% if all eigenvalues is non-zero, A=WW' is non-singular and 
% A^(-1/2) = V*D^(-1/2)*V^-1; W = A^(-1/2) * W.

% -------------------------------------------------------------------------
% output
state.icaweights = weights;

% compute non-stationary index (criterion for ICA convergence)
% |E{yf'}-I|_F

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

%         % check convergence criteria
%         if computeConvergence
%             % lambdaConv = (lambda_k(1) + lambda_k(end))/2;     % average lambda
%             Ka = 3;                     % assume 3 sec window
%             alpha = 1 - 1/(Ka*nPts);    %  
%             yEval = weights * blockdata;
%             if isempty(nlfunc)
%                 fEval(state.kurtsign,:)  = -2 * tanh(yEval(state.kurtsign,:));                           % Supergaussian
%                 fEval(~state.kurtsign,:) = tanh(yEval(~state.kurtsign,:)) - yEval(~state.kurtsign,:);    % Subgaussian
%             else
%                 fEval = nlfunc(yEval);
%             end
%             % compute outer product f*y^T for convergnece test
%             state.convMat = alpha .* state.convMat + (1-alpha) .* fEval * yEval' ./ nPts;
%             % computer inner product f^T*y for convergence test: as nonstationary index
%             state.statIdx = alpha * state.statIdx + (1-alpha) * sum(dot(fEval,yEval,1)) / nPts;
%         end    

