function [InvES,q,u,v]=nut_eiginv(Rzz,signalspace)
% signalspace should be vector of which eigenvalues to include
[u,q,v]=svd(Rzz,'econ');
q=diag(q);

% compute qinv by truncating and taking reciprocal of q;
% we want to keep original q for future reference and troubleshooting
qinv=zeros(size(q));
qinv(signalspace)=1./q(signalspace);

InvES=v*diag(qinv)*u';
