function [Gamma,s,w,cost,k,dGamma] = champagne_plain(y,F,Sigma_e,iters,Gamma_init,nd)

% This function implements the Champagne algorithm in a non-optimized, but
% readable way. Champagne is an iterative algorithm that infers the sources
% in the model y = F*s + e, where s ~ N[0,diag(Gamma)], e ~ N[0,Sigma_e] .
%
% *** USAGE ***
% [Gamma,s,w,cost,k,dGamma] = champagne_plain(y,F,Sigma_e,iters,Gamma_init,nd)
%
% *** INPUTS ***
% y             = nk X nt data matrix
% F             = nk X nd*nv lead-field (interleaved columns)
% Sigma_e       = nk X nk interference matrix
% iters         = maximum number of iterations
% Gamma_init    = nd X nd X nv array of initial gamma values (each nd x nd
%               block should be symmetric psd)
%
% *** OUTPUTS ***
% Gamma         = nd X nd X nv array of hyperparameters
% s             = nd X nv X nt arrary of source estimates
% w             = weights used for x = w*y (or s_bar = w*y)
% cost          = -log likelihood + constant terms
% k             = number of iterations used
% dGamma        = change in Gamma after k iterations

% Unless stated elsewise, all mentioned equation numbers in the code refer 
% to D. P. Wipf et al. (2010), NeuroImage.

% May 23rd, 2010 by Thomas Jochmann

%% Set constants

MIN_DELTA_GAMMA = 1e-10;% Minimum change in norm of gamma needed for early termination
MIN_GAMMA = 1e-8;          % Minimum value for pruning hyperparameters (to save time)
MIN_Z_EIGENVALUE = 1e-8;% Minimum eigenvalue of Zi before matrix inversion

IF_PLOT = 0;            % if true, cost, Gamma etc. are plotted each loop
SAVE_OUTPUT = 0;        % export_fig is required for this option. Make sure 
                        % to get it from MatlabCentral before activating it
                        
%% Initialize variables

nv = size(F,2)/nd;      % number of voxels
nk = size(F,1);         % number of sensors
nt = size(y,2);         % number of (post-stim) data vectors

F = reshape(F,nk,nd,nv);% reshape F 2d -> 3d. One slice per voxel.

cost = zeros(iters+1);  % one more value than iterations in order to consider the initial condition
X = zeros(nd,nt,nv);    % = sqrt(Q) in Attias' notes, mu in Wipf's code
s = zeros(nv*nd,nt);    % sources
Zi = zeros(nd,nd);      % = R in Attias' notes
w = zeros(nd*nv,nk);
activeVoxelList = 1:nv;
inactiveVoxelList = [];

Gamma = Gamma_init;
Gamma_old = zeros(size(Gamma)); % Gamma of the previous iteration step

y = y/sqrt(nt);

%% Learning loop
for k=1:iters+1
    % k=iters+1 in order to enter loop one more time, than the update rules
    % are applied, to calculate Sigma_b from the latest Gamma. exiting the
    % loop by using "break;" after Sigma_b is calculated.
    
    %% Prune inactive voxels
    % calculate a vector (length=nv) with the trace of each Gamma_i. The
    % trace is used as pruning criteria and for visualization.
    tr = zeros(nv,1); % trace vector of all Gammas
    for i=1:nd 
        tr = tr + squeeze(Gamma(i,i,:));
    end
    if (min(tr) < MIN_GAMMA ) % assume, small trace means inactive -> = 0
        % find active voxels
        activeVoxelList = find(tr > MIN_GAMMA)';
        inactiveVoxelList = setdiff(1:nv,activeVoxelList);
        % set Gamma of inactive voxels to zero
        for i=inactiveVoxelList
            Gamma(:,:,i)=zeros(nd,nd);
        end
        fprintf(' (pruned inactive voxels; %i of %i left (%i)\n',size(activeVoxelList,2),nv,size(activeVoxelList,2)-nv);
    end
    
    %% calculate Sigma_b and its inverse
    SumLGammaLtrans = zeros(nk,nk);
    for i=activeVoxelList
        SumLGammaLtrans = SumLGammaLtrans+F(:,:,i)*Gamma(:,:,i)*F(:,:,i)';
    end
    Sigma_b = Sigma_e + SumLGammaLtrans; % eq. (5)
    %invSigma_b = inv(Sigma_b); <- replaced by faster method below
    [V D] = eig(Sigma_b);                 % eig decomp: Sigma_b == V*D*V'
    invSigma_b = V*diag((1./diag(D)))*V'; % inv(Sigma_b) == V*inv(D)*V'
    
    %% calculate cost function and deltaGamma
    cost(k) = sum(log(svd(Sigma_b))) + sum(sum( y.*(Sigma_b\y) ));
    dGamma = sum(sum(sum(abs(Gamma-Gamma_old))));
    
    %% plot cost function and voxel activity
    if IF_PLOT
        plotChampagneStatus(cost,iters,k,nv,nt,tr,activeVoxelList,X)
        if SAVE_OUTPUT == 1; 
            export_fig('champagne_convergence','-pdf');
        end
    end
    
    %% exit conditions
    % for explanation, why this is in the middle of the loop, see comment 
    % after beginning of learning loop.
    if k>iters || dGamma < MIN_DELTA_GAMMA
        for i=inactiveVoxelList
            w(nd*i-nd+1:nd*i,:) = zeros(size(w(nd*i-nd+1:nd*i,:)));
        end
        break;
    end
    
    fprintf(' Loop %i/%i; deltaGamma = %f, cost = %f\n',k,iters,dGamma,cost(k));
    
    % sava Gamma of previous step
    Gamma_old=Gamma; 
    
    %% update aux. sources X, resolution matrix Z and therewith Gamma
    for i=activeVoxelList
        w(nd*i-nd+1:nd*i,:) = Gamma(:,:,i)*F(:,:,i)'*invSigma_b;
        X(:,:,i) = w(nd*i-nd+1:nd*i,:)*y; % eq. (10)
        Zi = F(:,:,i)'*invSigma_b*F(:,:,i); % eq. (11)
        
        %sqrtZi = real(sqrtm(Zi)); <- replaced by faster method below
        %invsqrtZi = inv(sqrtZi); <- replaced by faster method below
        [V D] = eig(Zi); % eig decomp: Zi == V*D*V' for fast computation of inv and invsqrt, exploiting Zi beeing real and symmetric.
        sqrtZi = V*diag( sqrt((sign(diag(D))+1)/2.*diag(D)) )*V'; % sqrtm(Zi) == V*sqrtm(D)*V'
        
        D = diag(max(diag(D),MIN_Z_EIGENVALUE)); % set eigen values at least to 10^-8 before inversion
        % this is important if a voxel's lead field isn't full rank like it
        % happens with 3-dimensional MEG lead fields.
        
        invsqrtZi = V*diag( 1./(sqrt((sign(diag(D))+1)/2.*diag(D))) )*V'; % inv(sqrtm(Zi)) == V*inv(sqrtm(D))*V'
        
        % Wipf update rule: update Gamma, following eq. (16)
        Gamma(:,:,i) = invsqrtZi*sqrtmRealSym(sqrtZi*X(:,:,i)*X(:,:,i)'*sqrtZi)*invsqrtZi;
    end
end

%% Calculate final values and return

k = k-1; % subtract last step of the loop, because no update occurs there
cost = cost(1:k+1); % cut off placeholders in the cost vector for iterations
                    % never made

for i=1:nv % put all voxel time courses one below the other
    s((i*nd-nd+1):(i*nd),:) = X(:,:,i);   
end

% rescale sources and weights, because they where calculated not by y but
% by y/sqrt(nt).

w = w*sqrt(nt);
s = s*sqrt(nt);

fprintf('\n\nPruned %i voxels.\n\n',(nv-size(activeVoxelList,2)));
return;




function A = sqrtmRealSym(A)
% This function performs an eigen decomposition and calculates the matrix
% square root through the scare root of the eigen values.
% If eigen values become negative (happens due to numerical reasons), they
% are set to zero.
[V D] = eig(A);  % eigen decomposition: A == V*D*V'
A = V*diag( sqrt((sign(diag(D))+1)/2.*diag(D)) )*V'; % A = sqrt(A)
% sign(diag(D))+1)/2 equals 1, if eigenvalue positive, 0 if negative and 1/2
% if 0, but in that case, the eigenvalue is zero anyway.
return;

function plotChampagneStatus(cost,iters,k,nv,nt,tr,activeVoxelList,X)
myColorOrder = sortrows(repmat([1 0 0;0 0 1; 1 1 0; 0 0 0],[3 1 1]));
subplot(3,2,1), plot(0:k-1,cost(1:k));
title(['Cost function; ' int2str(k-1) ' / ' int2str(iters) ]);
xlabel('Number of iterations');
ylabel('Cost');
set(gca(),'XLim',[0 k]);

subplot(3,2,2), plot(1:nv,tr);
title(['trace(Gamma); active voxels: ' int2str(size(activeVoxelList,2)) ' / ' int2str(nv)]);
xlabel('Voxels i');
ylabel('trace(Gamma_i)');

% find and plot most active voxels
mostActiveVoxels = find(tr>0.1*max(tr));
[neverused,index] = sort(tr(mostActiveVoxels),'descend');
for i=1:min(5,size(mostActiveVoxels,1))
    v = mostActiveVoxels(index(i));
    subplot(3,2,3:6), plot(1:nt,X(:,:,v)');
    hold on;
end
hold off;

title('source time courses)');
xlabel('time (samples)');
ylabel('source strength');

drawnow;
%maximize; % maximizes the figure window after plotting. requires the 
          % maximize function from matlabcentral. no effect on the results.
return;