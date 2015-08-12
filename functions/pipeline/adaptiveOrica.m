function state = adaptiveOrica(blockdata, state, options, dataRange, adaptiveFF, computeConvergence,nlfunc)
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
    adaptiveFF.arg_selection = 'nonadapt_cooling'; end
if nargin < 6
    computeConvergence = false; end
if nargin < 7
    nlfunc = []; end
% ------------------------------------------------------------------------
% initialize
[nChs, nPts] = size(blockdata);
weights      = state.icaweights;
f            = zeros(nChs, nPts);
if ~isempty(state.lambdaComp) && size(state.lambdaComp,2)~=nPts, state.lambdaComp = repmat(state.lambdaComp(:,end),1,nPts); end

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

if computeConvergence
%     A = pinv(weights);
%     flow = A - blockdata*f'/nPts;   
    flow = eye(nChs)+y*f'/nPts;
    variance = blockdata.*blockdata;
    if isempty(state.Rn)
        state.Rn = flow;
        state.Var = sum(variance,2)/(nPts-1);
    else
        state.Rn = (1-state.leakyAvgDelta)*state.Rn + state.leakyAvgDelta*flow; % !!! this does not account for block update!
        state.Var = (1-state.leakyAvgDeltaVar)^nPts*state.Var ...
            + sum(state.leakyAvgDeltaVar*bsxfun(@times,variance,(1-state.leakyAvgDeltaVar).^(nPts-1:-1:0)),2);
    end
    state.normRn = norm(state.Rn,'fro');
    if isempty(state.minNormRn), state.minNormRn = state.normRn; end;
    state.minNormRn = max(min(state.minNormRn, state.normRn),1);
    state.ratioOfNormRn = state.normRn/state.minNormRn;

    state.columnWiseNormRn = sum(abs(state.Rn).^2,2).^0.5;
    if isempty(state.mincolumnWiseNormRn), state.mincolumnWiseNormRn = state.columnWiseNormRn; end;
    state.mincolumnWiseNormRn = max(min(state.mincolumnWiseNormRn, state.columnWiseNormRn),ones(nChs,1));
    state.columnWiseNormRatio = state.columnWiseNormRn ./ state.mincolumnWiseNormRn;
end

% define adaptive forgetting rate: lambda
switch adaptiveFF.arg_selection
    case 'nonadapt_cooling' % cooling FF: lambda_0 / (t+t_0)^gamma + lambda_ss
        state.lambda_k = nonadapt_cooling(state.counter+dataRange, state.gamma, state.lambda_0, state.t_0, state.lambda_ss);
        if state.lambda_k(1) < state.lambda_const, state.lambda_k = repmat(state.lambda_const,1,nPts); end
    case 'nonadapt_const' % constant FF
        state.lambda_k = repmat(state.lambda_const,1,nPts);
    case 'nonadapt_exp' % exponentially attenuated FF:  lambda_0 * exp(-gamma*(t-t_0)), t>t_0
        state.lambda_k = nonadapt_exp(state.counter+dataRange,state.lambda_0,state.gamma,state.t_0);
    case 'adapt_Murata_block' % Murata adaptive FF: 
        state.lambda_k = adapt_Murata_block(dataRange,state.lambda_k,state.decayRateAlpha,state.upperBoundBeta,state.transBandWidthGamma,state.transBandCenter,state.ratioOfNormRn);
    case 'adaptiveComponent'
        state.lambdaComp = adaptiveComponent(dataRange,state.lambdaComp,state.decayRateAlpha,state.upperBoundBeta,state.transBandWidthGamma,state.transBandCenter,state.columnWiseNormRatio);
    case 'adapt_Cichocki'
        state.lambdaComp = adapt_Cichocki('component',state.lambdaComp,state.V,flow*weights,state.deltaV,state.deltaU,state.alpha);
    otherwise
        disp('Method not specified.')
end

if isempty(state.lambdaComp)
    % % version 1: block update for global lambda: doesn't average when computing inner and outer product!!
    lambda_prod = prod(1./(1-state.lambda_k));
    Q = 1 + state.lambda_k .* (dot(f,y,1)-1);
    weights = lambda_prod * (weights - y * diag(state.lambda_k./Q) * f' * weights);
else
    % version 2: block update for component-wise lambda
    % the update rule needs validation
    lambdaMat = state.lambdaComp; % repmat(state.lambdaComp,1,nPts);
    lambdaMatProd = cumprod(1./(1-lambdaMat),2);
    M = y(:,1)*f(:,1)'*diag(lambdaMat(:,1)./(1-lambdaMat(:,1)));
    Q = M./(1+trace(M));
    for it = 2:nPts
        M = y(:,it)*f(:,it)'*diag(lambdaMat(:,it)./(1-lambdaMat(:,it)));
        Q = Q + diag(1./lambdaMatProd(:,it))*M*diag(lambdaMatProd(:,it))./(1+trace(M));
    end
    weights = diag(lambdaMatProd(:,end)) * (weights - Q * weights); 
end

% orthogonalize weight matrix - avoid inverse of signular matrix
[V,D] = eig(weights * weights');
weights = V/sqrt(D)/V * weights;    % run 30% faster
% weights = real(inv(V * sqrt(D) / V)) * weights;
% if all eigenvalues is non-zero, A=WW' is non-singular and 
% A^(-1/2) = V*D^(-1/2)*V^-1; W = A^(-1/2) * W.

% -------------------------------------------------------------------------
% output
state.icaweights = weights;

% compute non-stationary index (criterion for ICA convergence)
% |E{yf'}-I|_F

if computeConvergence
    % TODO: implement some signal dependence measurement
    
    yNew = weights * blockdata;
    if isempty(nlfunc)
        fNew(state.kurtsign,:)  = -2 * tanh(yNew(state.kurtsign,:));                           % Supergaussian
        fNew(~state.kurtsign,:) = tanh(yNew(~state.kurtsign,:)) - yNew(~state.kurtsign,:);    % Subgaussian
    else
        fNew = nlfunc(yNew);
    end

    % Gradient step size: Douglas 1998
    state.lambdaGrad = -(nChs-1)/(1+state.lambda_k(end)) - (1-trace(y'*f)/nPts)/(1+state.lambda_k(end)*(1-trace(y'*f)/nPts)) ...
                       + trace(yNew'*fNew)/nPts - trace(fNew'*f*y'*yNew)/nPts^2;
             
    % Nonlinear correlation: motivated from Murata 2002
    delta = 0.02; % forgetting factor of signal dependence measurement.
    state.RnW = (1-delta) * state.RnW + delta * (eye(nChs) + yNew * fNew' / nPts);
    state.normRnW = norm(state.RnW, 'fro');       
    
    % Cross correlation: Zhang 2003
    delta = 0.02; % forgetting factor of signal dependence measurement.    
    yEst = (1-delta) * state.yEst + delta * mean(yNew,2);
    fEst = (1-delta) * state.fEst + delta * mean(fNew,2);
    yError = yEst - state.yEst;
    fError = fEst - state.fEst;
    state.covyy = (1-delta)*(state.covyy+yError*yError') + delta*bsxfun(@minus,yNew,yEst)*bsxfun(@minus,yNew,yEst)';
    state.covfy = (1-delta)*(state.covfy+fError*fError') + delta*bsxfun(@minus,fNew,fEst)*bsxfun(@minus,yNew,yEst)';
    % can calculate component-wise correlation
    state.normCovyy = norm(state.covyy, 'fro');  
    state.normCovfy = norm(state.covfy, 'fro');  
    
    state.yEst = yEst;
    state.fEst = fEst;
    
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

% Murata adaptive FF
function lambda = adapt_Murata_block(dataRange,lambda,decayRateAlpha,upperBoundBeta,transBandWidthGamma,transBandCenter,ratioOfNormRn)
    
    gainForErrors = upperBoundBeta*0.5*(1+tanh((ratioOfNormRn-transBandCenter)/transBandWidthGamma));
    f = @(n) (1+gainForErrors).^n * lambda(end) - decayRateAlpha*((1+gainForErrors).^(2*n-1)-(1+gainForErrors).^(n-1))/gainForErrors*lambda(end)^2;
    lambda = f(1:length(dataRange));
    
end

function lambda = adaptiveComponent(dataRange,lambda,decayRateAlpha,upperBoundBeta,transBandWidthGamma,transBandCenter,columnWiseNormRatio)

gainForErrors = upperBoundBeta*0.5*(1+tanh((columnWiseNormRatio-transBandCenter)/transBandWidthGamma));
    f = @(ch,n) (1+gainForErrors(ch)).^n * lambda(ch,end) - decayRateAlpha*((1+gainForErrors(ch)).^(2*n-1)-(1+gainForErrors(ch)).^(n-1))/gainForErrors(ch)*lambda(ch,end)^2;
    for ic = 1:size(lambda,1)
        lambda(ic,:) = f(ic,1:length(dataRange));
    end

end

% Cichocki adaptive FF
function U = adapt_Cichocki(opt,U,V,G,deltaV,deltaU,alpha)
    switch opt
        case 'element'
            V = V + deltaV*(G-V);
            U = U + deltaU*(alpha*abs(V)-U); % U is matrix
        case 'component'
            V = V + deltaV*(G-V);
            U = U + deltaU*(alpha*U.*sum(tanh(abs(V)),2) - U.^2);
    end

end


