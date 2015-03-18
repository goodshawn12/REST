function rnk = nut_rank(A,tol)
% rank = nut_rank(A,tol)
% rank = nut_rank(A)  % default tol=1e-5*sum(svd(A))
% replaces matlab's rank, which sucks with real data (esp. lead fields)
% tolerance is a *ratio* of each component to sum of all component
% (i.e., percent of total energy)

if(nargin==1)
    tol = 1e-5;
end
eigvals=svd(A);
ratios=eigvals./sum(eigvals);
rnk = max(find(eigvals > tol*sum(eigvals)));
