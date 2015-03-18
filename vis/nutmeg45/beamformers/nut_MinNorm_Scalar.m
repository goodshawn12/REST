function [weight,eta]=nut_MinNorm_Scalar(Lp,data, flags) %---------------------------------------------------------
% NUT_SCALAR_MINNORM   scalarizes leadfield and calcs MinNorm Estimate.
%
% [weight,eta]=nut_MinNorm_Scalar(Lp,data,flags)
%
%
% Lp           lead field (required)
% data.Ryy     sample covariance (required)
% flags        (optional) 
%  flags.gamma=[]      data-driven regularization, or
%  flags.gamma=user defined regularization constant, or 
%  flags.gamma='auto'  for leadfield-based regularization, or
%  flags.snr=estimated data SNR   for noise based regularization
% 
% weight        spatial filter matrix
% eta           optimal dipole orientations

if nargin<3, flags.gamma=[]; end

L = reshape(Lp,size(Lp,1),size(Lp,2)*size(Lp,3));
G = L*L'; clear L

if isfield(flags,'snr') 
    if ~isfield(data,'C') || isempty(data.C)
        data.C = eye(size(data.Ryy));
    end    
    gamma = trace(G)/(trace(data.C)*flags.snr^2)   
    InvG = inv(G+gamma*data.C);
elseif ~isfield(flags,'gamma') || isempty(flags.gamma)
    gamma = 1e0*max(eig(data.Ryy))
    InvG = inv(G+gamma*eye(size(G)));
elseif isnumeric(flags.gamma)
    gamma = flags.gamma   %* max(eig(data.Ryy))
    InvG = inv(G+gamma*eye(size(G)));
else
    % leadfield-based normalization
    doplot=false;
    x = [-20:20];
    gamma=10.^x; %.* max(eig(data.Ryy));
    numg = length(gamma);
    InvG = zeros([size(G) numg]);
    meanuptr = zeros(numg,1);
    warning('off','MATLAB:nearlySingularMatrix')
    for k=1:numg
        InvG(:,:,k) = inv(G+gamma(k)*eye(size(G)));
        uptr = abs(triu(InvG(:,:,k),1));  
        meanuptr(k) = mean(uptr( find(uptr(:)) ));  % mean of matrix without diagonal
    end
    warning('on','MATLAB:nearlySingularMatrix')
    if doplot, figure; plot(x(1:end-2),-diff(diff(log10(meanuptr)))); end
    [dum,idx]=findpeaks(-diff(diff(log10(meanuptr))),'MINPEAKHEIGHT',.1);  % max of second derivation of InvG across gamma magnitudes
    idx = max(idx);
    fprintf('Optimal gamma: 1e%d\n',x(idx));
    InvG = InvG(:,:,idx);
end

if size(Lp,2)<3
    disp('Leadfield is already scalar, skipping orientation calculation.'),
    L = permute(Lp,[1 3 2]);
else
    GRG = InvG * data.Ryy * InvG;
    L = zeros(size(Lp,1),size(Lp,3));
    eta = zeros(size(Lp,3),3);
    for i=1:size(Lp,3);
        [vvv,ddd] = eig(Lp(:,:,i)' * GRG * Lp(:,:,i));
        [jnk,maxeig]=max(diag(ddd));
        eta(i,:)=vvv(:,maxeig);
        L(:,i) = Lp(:,:,i) * eta(i,:)';
    end
end

w1 = zeros(size(L));

for i=1:size(Lp,3)
    w1(:,i)=InvG*L(:,i);
%     J = inv(sqrt(Lp1(:,i)'*InvGLp));
%     w1(:,i) = InvGLp * J;
end

weight(:,1,:) = w1;

disp('done');
% end
