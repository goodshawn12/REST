function [a,lam,alp,xbar,psi,p,d]=nut_vbfa(y,nl,nem,a_init,lam_init,plotflag);

% vb-em algorithm for inferring the factor analysis model   y = a*x + v
%
% y(nk,nt) = data
% a(nk,nl) = mixing matrix
% lam(nk,nk) = diagonal noise precision matrix 
% alp(nl,nl) = diagonal hyperparmaeter matrix 
% xbar(nl,nt) = posterior means of the factors 
%
% nk = number of data points
% nt = number of time points
% nl = number of factors
% nem = number of em iterations

nk=size(y,1);
nt=size(y,2);

% initialize by svd

ryy=y*y';
if a_init==0
   [p d]=svd(ryy/nt);
   d=diag(d);
   a=p*diag(sqrt(d));
   a=a(:,1:nl);
   lam=diag(nt./diag(ryy));
else
   a=a_init;
   lam=lam_init;
end

alp=diag(1./diag(a'*lam*a/nk));
alp=min(diag(alp))*diag(ones(nl,1));
psi=eye(nl)/nt;

% em iteration

like=zeros(nem,1);
alapsi=a'*lam*a+nk*psi;

for iem=1:nem
   gam=alapsi+eye(nl);
   igam=pinv(gam);  % for large inputs pinv is better than nut_pinv.
   xbar=igam*a'*lam*y;

   ldlam=sum(log(diag(lam/(2*pi))));
   ldgam=sum(log(svd(gam)));
   ldalp=sum(log(diag(alp)));
   ldpsi=sum(log(svd(psi)));
   if (plotflag | iem==nem)
   like(iem)=.5*nt*(ldlam-ldgam)-.5*sum(sum(y.*(lam*y)))+.5*sum(sum(xbar.*(gam*xbar)))+.5*nk*(ldalp+ldpsi);
   end
if plotflag
    figure(10);
   subplot(3,3,7);plot((1:iem)',like(1:iem));title('vbfa');ylabel('vbfa');
   subplot(3,3,8);plot([mean(a.^2,1)' 1./diag(alp)]);
   subplot(3,3,9);plot(1./diag(lam));
   drawnow;
end
   ryx=y*xbar';
   rxx=xbar*xbar'+nt*igam;
   psi=inv(rxx+alp);

   a=ryx*psi;
   lam=diag(nt./diag(ryy-a*ryx'));
   alapsi=a'*lam*a+nk*psi;
   alp=diag(nk./diag(alapsi));
end

return



%-----------------------------------------------------------------------




