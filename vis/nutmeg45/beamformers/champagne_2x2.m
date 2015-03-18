function [gamma,s,w,cost,k,dgamma] = champagne_2x2(y,F,Sigma_e,iters,gamma_init);
% *************************************************************************
% 
% *** Purpose ***
% Iterative algorithm for inferring sources s in the model
% y = F*s + e, where s ~ N[0,diag(gamma)], e ~ N[0,Sigma_e]
%
%
% *** USAGE ***
% [gamma,s_bar,like,k,dgamma] = localize_sources(y,F,Sigma_e,iters,gamma_init);
%
%
% *** INPUTS ***
% y             = nk X nt data matrix
% F             = 2 X nk X nv lead-field for two components
% Sigma_e       = nk X nk interference matrix
% iters         = maximum number of iterations
% gamma_init    = 2 X 2 X nv array of initial gamma values (each 2 x 2
%               block should be symmetric psd)
%
%
% *** OUTPUTS ***
% gamma         = 2 X 2 X nv array of hyperparameters
% s             = 2 X nv X nt arrary of source estimates
% w             = weights used for x = w*y
% cost          = -log likelihood + constant terms
% k             = number of iterations used
% dgamma        = change in gamma after k iterations
%
%
% *************************************************************************

MIN_DELTA_GAMMA = 1e-10;     % Minimum change in norm of gamma needed for early termination
MIN_GAMMA = 1e-5;           % Minimum value for pruning hyperparameters

% Check conditioning of initial gamma and exit if bad
ok = check_gamma_cond(gamma_init,MIN_GAMMA);
if (~ok)
    disp('Initial gamma is ill-conditioned');
    return;
end;

% Initializations
[dummy nk nv] = size(F);	% nk = num sensors, nl = (num voxels * dip comps.)
[nk nt] = size(y);          % nt = num time points 
nl = 2*nv;                  % nl = num dipole comps.: assume 2 dipole components per voxel

gamma_sp = sparse2x2(gamma_init);
cost = zeros(iters,1);
keep_list = [1:nv]';

% Fuse leadfield into single nk X nl matrix with related dipoles components
% in adjacent columns
temp = zeros(nk,nl);
temp(:,[1:2:nl]) = F(1,:,:);
temp(:,[2:2:nl]) = F(2,:,:);
F = temp;

% Create interference normalized version of F for use with svd implementation
% of s_bar computation below
[U,S,V] = svd(Sigma_e);
temp = 1./( sqrt(diag(S)) );
temp=double(temp);
inv_sqrt_S = sparse(diag(temp));
inv_sqrt_Sigma_e = (V*inv_sqrt_S*U');
eF = inv_sqrt_Sigma_e*F;


%gamma1x1 = ones(nl,1);
%gamma1x1_old = 0;

% Learning loop
k = 1;
gamma_old_sp = 0;
m = nv;

while (1)
    
	% *** Compute -log likelihood and display ***
    
	Sigma_y = F*gamma_sp*F' + Sigma_e;
	cost(k) = sum(log(svd(Sigma_y))) + sum(sum( y.*(inv(Sigma_y)*y) ))/nt;
    
	%d = gamma_sp - gamma_old_sp;
    %dgamma = full(sum(sum(d'*d)));
    
    d = full( diag(gamma_sp) - diag(gamma_old_sp) );
    dgamma = max(abs(d));
   
    
	%Sigma_y = F*sparse(diag(gamma1x1))*F' + Sigma_e;
	%cost1x1(k) = sum(log(svd(Sigma_y))) + sum(sum( y.*(inv(Sigma_y)*y) ))/nt;
    
    % disp([k cost(k) cost1x1(k) dgamma max(abs(gamma1x1_old - gamma1x1))]);
    disp([k cost(k) dgamma m ]);
    
	% *** Check stopping criteria and update k ***

    if ( dgamma < MIN_DELTA_GAMMA )
        break;
    end;
    if (k == iters)
        break;
    end;
    k = k + 1;
    
    %gamma1x1_old = gamma1x1;
    %gamma1x1 = update_gamma1x1(F,y,1,gamma1x1);
     
    % *** Prune things as hyperparameters go to zero ***
    
	gamma = full2x2(gamma_sp);
	tr = squeeze(gamma(1,1,:)+gamma(2,2,:));
   
    if (min(tr) < MIN_GAMMA )
		index = find(tr > MIN_GAMMA);
		gamma = gamma(:,:,index);
        gamma_sp = sparse2x2(gamma);
        F = F(:,sort([(2*index-1); 2*index]));
        eF = eF(:,sort([(2*index-1); 2*index]));

		keep_list = keep_list(index);
        m = length(keep_list);
        nl = 2*m;
     
        if (m == 0)   break;  end;
    end;
    

	% *** Compute mean of s given gamma ***
    
    gamma = full2x2(gamma_sp);
	[U_gamma S_gamma] = eig2x2(gamma);
	U_gamma_sp = sparse2x2(U_gamma);
	S_gamma_sp = sparse2x2(S_gamma);
    
	sqrt_gamma = U_gamma_sp*sqrt(S_gamma_sp)*U_gamma_sp';
	eFg = eF*sqrt_gamma;

	[U S V] = svd(eFg,'econ');
    S = diag(S);
	temp = V*sparse(diag(S./(S.^2 + 1)))*U' ;
	w = sqrt_gamma*temp*inv_sqrt_Sigma_e;

	s_bar=w*y;


    
    % *** Update gamma using Attias decomposition method ***

    % Compute 2x2 block diagonal elements of Q from Attias notes
	b0 = sum(s_bar.^2,2)/nt;
	bm = sum(s_bar.*[s_bar(2:nl,:);zeros(1,nt)],2)/nt;
	bm(2:2:nl) = 0;
	bp = [0;bm(1:nl-1)];
	mu2_sp = spdiags([bm b0 bp],[-1 0 1],nl,nl);  % This is called Q

    
    % Compute 2x2 block diagonal elements of S from Attias notes
    % (here it is called R b/c it is similar to the resolution matrix
	b0 = sum(w'.*F(:,1:nl),1)';
    
	bm = sum(w'.*[zeros(nk,1) F(:,1:nl-1)],1)'; 
    bm=[bm(2:nl);0];
    bm(2:2:nl)=0;
    
	bp=sum(w'.*[F(:,2:nl) zeros(nk,1)],1)';
    bp=[0;bp(1:nl-1)];
    bp(1:2:nl)=0;
    
	temp = spdiags([bm b0 bp],[-1 0 1],nl,nl);
	R_sp = inv(gamma_sp)*temp;

   
    % Now compute inv sqrt of R and sqrt of C from Attias notes    
    sqrt_R = sqrt2x2(R_sp);
    inv_sqrt_R = inv(sqrt_R);

    C_sp = sqrt_R*mu2_sp*sqrt_R;
    sqrt_C = sqrt2x2(C_sp);
   

    % Finally update gamma
    gamma_old_sp = gamma_sp;
    gamma_sp = inv_sqrt_R*sqrt_C*inv_sqrt_R;;

end;


% *** Expand hyperparameters and sources ***
gamma = full2x2(gamma_sp);
temp = zeros(2,2,nv);
if (m > 0) temp(:,:,keep_list) = gamma;  end;
gamma = temp;

temp = zeros(2,nv,nt);
if (m > 0) 
    temp(1,keep_list,:) = s_bar(1:2:nl,:);
    temp(2,keep_list,:) = s_bar(2:2:nl,:);
end;
s = temp;

cost = cost(1:k);

return;



function gamma = full2x2(gamma_sp);
%-----------------------------------------------------------------------
% Convert sparse matrix to full representation
%-----------------------------------------------------------------------

nl=size(gamma_sp,1);
nv=nl/2;
gamma = zeros(2,2,nv);

b0=diag(gamma_sp,0);
bp=diag(gamma_sp,1);
bm=diag(gamma_sp,-1);

gamma(1,1,:)=b0(1:2:nl);
gamma(2,2,:)=b0(2:2:nl);
gamma(1,2,:)=bp(1:2:nl);
gamma(2,1,:)=bm(1:2:nl);

return;



function gamma_sp=sparse2x2(gamma)
%-----------------------------------------------------------------------
% Converts gamma to a sparse matrix structure for efficient computations
%-----------------------------------------------------------------------

nl = 2*size(gamma,3);

b0=zeros(nl,1);bm=zeros(nl,1);bp=zeros(nl,1);

b0(1:2:nl)=squeeze(gamma(1,1,:));
b0(2:2:nl)=squeeze(gamma(2,2,:));
bm(1:2:nl)=squeeze(gamma(2,1,:));
bp(2:2:nl)=squeeze(gamma(1,2,:));

gamma_sp = spdiags([bm b0 bp],(-1:1),nl,nl);

return;



function [U,S]=eig2x2(gamma);
%-----------------------------------------------------------------------
% Compute block diagonal eigendecomposition of each 2x2 gamma block
% U is composed of 2x2 orthogonal blocks of eigenvectors
% S is diagonal with associated eigenvalues
%-----------------------------------------------------------------------

nv=size(gamma,3);
U=zeros(2,2,nv);
S=zeros(2,2,nv);

U(1,1,:)=1;U(2,2,:)=1;
S(1,1,:)=gamma(1,1,:);S(2,2,:)=gamma(2,2,:);

% Only compute decomposition when significant off-diagonal elements are
% present
ff=find(abs(squeeze(gamma(1,2,:)))>1e-16);

g11=squeeze(gamma(1,1,ff));
g12=squeeze(gamma(1,2,ff));
g22=squeeze(gamma(2,2,ff));

lam1 = (g11 + g22 + sqrt((g11-g22).^2 + 4*(g12.^2)))/2;
lam2 = (g11 + g22 - sqrt((g11-g22).^2 + 4*(g12.^2)))/2;

s1 = sqrt((g11-lam1).^2 + g12.^2);
U(1,1,ff) = g12./s1;
U(2,1,ff) = (lam1-g11)./s1;

s2 = sqrt((g11-lam2).^2 + g12.^2);
U(1,2,ff) = g12./s2;
U(2,2,ff) =(lam2-g11)./s2;

S(1,1,ff)=lam1;
S(2,2,ff)=lam2;

return;



function [sqrt_X] = sqrt2x2(X_sp)
%-----------------------------------------------------------------------
% Compute the sqrt of a sparse, 2x2 block diagonal matrix X_sp
%-----------------------------------------------------------------------

X = full2x2(X_sp);
[U_X S_X] = eig2x2(X);
U_X_sp = sparse2x2(U_X);
S_X_sp = sparse2x2(S_X);
    
temp = full(diag(S_X_sp));
temp(find(temp<1e-16)) = 1e-16;
%S_X_sp = sparse(diag(temp)); 
ind1=[1:length(temp)];
ind2=[1:length(temp)];
S_X_sp = sparse(ind1,ind2,temp);
sqrt_X = U_X_sp*sqrt(S_X_sp)*U_X_sp';

return;



function [gamma,cost] = update_gamma1x1(Phi,T,lambda,gamma)

[N,M] = size(Phi);
[N,L] = size(T);

G = repmat(sqrt(gamma)',N,1);
PhiG = Phi.*G; 
[U,S,V] = svd(PhiG,'econ');
    
[d1,d2] = size(S);
if (d1 > 1)     diag_S = diag(S);  
else            diag_S = S(1);      end;
    
U_scaled = U(:,1:N).*repmat((diag_S./(diag_S.^2 + lambda + 1e-16))',N,1);       
Xi = G'.*(V*U_scaled'); 
        
mu = Xi*T; 
    
    
% *** Update hyperparameters ***
mu2_bar = sum(abs(mu).^2,2);
R_diag = real( (sum(Xi.'.*Phi)).' );
gamma = mu2_bar./(L*R_diag);  


return;



function [ok] = check_gamma_cond(gamma,MIN_GAMMA);
%-----------------------------------------------------------------------
% Check gamma for ill-conditioned blocks
%-----------------------------------------------------------------------

nv = size(gamma,3);
tr = squeeze(gamma(1,1,:)+gamma(2,2,:));

ok = 1;

for i = 1:nv
    g = gamma(:,:,i);
    
    if ( (tr(i) > MIN_GAMMA) & (cond(g) > 1e6) )
        ok = 0;
        break;
    end;    
end;

return;
