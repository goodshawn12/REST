function [J,alpha,beta,T,history] = dynamicLoreta(Ut,Y,s2,iLV,L,options,alpha,beta)
%[J,varargout] = dynamicLoreta(V,varargin)
%
% Computes the posterior distribution of the parameters J given some data V.
% The program solves levels of inference: 1) optimization of parameters J, and
% 2) optimization of hyperparameters alpha and beta. See Trujillo-Barreto
% et. al. (2004) for details.
%
% Ut,s2, and iLV are defined as follows:
% Y: Nsensors x time points data matrix
% K: N x P predictor matrix
% L: sparse P x P square root of the precision matrix
% [U,s,V] = svd( K*inv(L) )
% iLV = inv(L)*V
% s2 = s.^2
%
% alpha, beta: hyperparameters
% J: estimated parapeters
%
% P(V|J,alpha)*P(J|beta)
% P(J|V,alpha,beta) = ----------------------
% P(V|alpha,beta)
%
% Author: Alejandro Ojeda, SCCN/INC/UCSD, Jan-2013
%
% References:
% Trujillo-Barreto, N., Aubert-Vazquez, E., Valdes-Sosa, P.A., 2004.
% Bayesian model averaging in EEG/MEG imaging. NeuroImage 21, 1300???1319
if nargin < 5, error('Not enough input arguments.');end
if nargin < 6
    options.maxTol = 1e-3;
    options.maxIter = 100;
    options.gridSize = 100;
    options.verbose = true;
    options.history = true;
    options.useGPU = false;
    options.initNoiseFactor = 0.001;
end
if options.history
    [history.alpha, history.beta, history.err] = deal(nan(1,options.maxIter));
else
    history = [];
end
s = s2.^(0.5);
n = length(s);
p = length(s);
alpha_vec = zeros(options.maxIter,1);
beta_vec = zeros(options.maxIter,1);
gcv_vec = zeros(options.maxIter,1);
err_win = 3;
% Initialize hyperparameters
if nargin < 7 || isempty(alpha)
    UtY = Ut*Y;
    tol = max([n p])*eps(max(s));
    lambda2 = logspace(log10(tol),log10(max(s)),options.gridSize);
    gcv = zeros(options.gridSize,1);
    for k=1:options.gridSize
        d = lambda2(k)./(s2+lambda2(k));
        f = mean(diag(d)*UtY,2);
        gcv(k) = dot(f,f,1)/sum(d)^2;
    end
    loc = getMinima(gcv);
    if isempty(loc), loc = 1;end
    loc = loc(end);
    lambda2 = lambda2(loc);
    alpha = options.initNoiseFactor*(Y(:)'*Y(:))/n;
    beta = alpha*lambda2;
    alpha_vec(1) = alpha;
    beta_vec(1) = beta;
    gcv_vec(1) = gcv(loc);
end
for it=2:options.maxIter
    % computing hat matrix, mse, and ||L*J_hat||^2
    H = Ut'*diag(alpha*s2./(alpha*s2+beta))*Ut;
    mse = mean(sum((Y - H*Y).^2));
    norm_JtJ = trace(Y'*Ut*diag(((alpha*s)./(alpha*s2+beta)).^2)*Ut'*Y);
    % computing gamma
    gamma = p-beta*sum(1./(alpha*s2+beta));
    % computing GCV
    gcv = mse/(1-trace(H)/n)^2;
    gcv_vec(it) = gcv;
    % updating hyperparameters
    alpha = n-gamma;
    beta = gamma/(norm_JtJ+eps); % adding eps for numerical stability
    alpha_vec(it) = alpha;
    beta_vec(it) = beta;
    if it-err_win < 1
        err = 0.5*std(alpha_vec(1:it)) + 0.5*std(beta_vec(1:it));
    else
        err = 0.5*std(alpha_vec(it-err_win:it)) + 0.5*std(beta_vec(it-err_win:it));
    end
    if options.history
        history.alpha(it) = alpha;
        history.beta(it) = beta;
        history.err(it) = err;
    end
    if options.verbose
        disp([num2str(it-1) ' => alpha: ' num2str(alpha) ' beta: ' num2str(beta) ' df: ' num2str(gamma) ' hyperp. error: ' num2str(err) ' gcv: ' num2str(gcv)]);
    end
    if err < options.maxTol, break;end
end
if it == options.maxIter, warning('Maximum iteration reached. Failed to converge.');end
% parameters's estimation
T = iLV*diag(alpha.*s./(alpha.*s2+beta))*Ut;
J = T*Y;
% standardized Loreta
E = sum(Y-H*Y,2);
sigma = E'*E/(n-trace(H));
dT = 1./(sqrt(dot(T,T,2))+eps); % adding eps for numerical stability
S = 1./sigma*dT;
S = S./std(eps+S);
T = bsxfun(@times,T,S);
J = bsxfun(@times,J,S);
end
%---
function indmin = getMinima(x)
fminor = diff(x)>=0;
fminor = ~fminor(1:end-1, :) & fminor(2:end, :);
fminor = [0; fminor; 0];
indmin = find(fminor);
end
function T = standT(T,H,Y,n)
E = Y-H*Y;
df = (n-trace(H));
sigma = E'*E/df;
dT = 1./sqrt(dot(T,T,2));
S = 1./sigma*dT;
T = bsxfun(@times,S,T);
end