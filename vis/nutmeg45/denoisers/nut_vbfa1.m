function [a,nu,alp,xbar,ubar1,cy,like,www,psiaa,outalp,igam]=nut_vbfa1(y,nl,b,lam,nem,a_init,infer_nu,plotflag);

% vb-em algorithm for inferring the evoked mixing matrix in the factor 
% analysis model   y = a*x + b*u + v
% learn a,nu
% (source estimates) = www * y     (where y = observations)
% (cleaned observations) = a * (source estimates) = a * www * y
%
% y(nk,nt) = data
% a(nk,nl) = evoked mixing matrix
% b(nk,nm) = interference mixing matrix
% lam(nk,nk) = diagonal noise precision matrix 
% alp(nl,nl) = diagonal hyperparmaeter matrix 
% xbar(nl,nt) = posterior means of the evoked factors 
% nu(nm,1) = precision of interference factors
% infer_nu = 1 means infer nu from data, otherwise nu = ones(nm,1)
%
% nk = number of data points
% nt = number of time points
% nl = number of evoked factors
% nm = number of interference factors
% nem = number of em iterations

nk=size(y,1);
nt=size(y,2);
nm=size(b,2);

% initialize by svd

ryy=y*y';
if a_init==0 
%   [p d q]=svd(ryy/nt);d=diag(d);
%   a=p*diag(sqrt(d));
%   a=a(:,1:nl);
   sig0=b*b'+diag(1./diag(lam));
   [p0 d0]=svd(sig0);
   d0=diag(d0);
   s=p0*diag(sqrt(d0))*p0';
   invs=p0*diag(1./sqrt(d0))*p0';
   [p d]=svd(invs*ryy*invs/nt);
   d=diag(d);
%   a=s*p(:,1:nl)*diag(sqrt(max(d(1:nl)-1,0)));
   a=s*p(:,1:nl)*diag(sqrt(abs(d(1:nl)-1)));
else
   a=a_init;
end

alp=diag(1./diag(a'*lam*a/nk));
alp=min(diag(alp))*diag(ones(nl,1));
psi=eye(nl)/nt;
nu=ones(nm,1);

% em iteration

like=zeros(nem,1);

for iem=1:nem
   ab=[a b];
%   gam=ab'*lam*ab+eye(nl+nm);
   gam=ab'*lam*ab+diag([ones(nl,1);nu]);
   gam(1:nl,1:nl)=gam(1:nl,1:nl)+nk*psi;
   igam=inv(gam);
   xubar=igam*ab'*lam*y;
   xbar=xubar(1:nl,:);
   ubar=xubar(nl+1:nl+nm,:);
   www = igam(1:nl,:)*ab'*lam;

   ldlam=sum(log(diag(lam/(2*pi))));
   ldgam=sum(log(svd(gam)));
   ldalp=sum(log(diag(alp)));
   ldpsi=sum(log(svd(psi)));
   if (plotflag | iem==nem)
   like(iem)=.5*nt*(ldlam-ldgam)-.5*sum(sum(y.*(lam*y)))+.5*sum(sum(xubar.*(gam*xubar)))+.5*nk*(ldalp+ldpsi)+.5*nt*sum(log(nu));
   end
   if plotflag
figure(3)
   subplot(3,3,2);plot((1:iem)',like(1:iem));title('vbfa1');
   subplot(3,3,5);plot([mean(a.^2,1)' 1./diag(alp)]);
   drawnow;
   end
   
   ryx=y*xbar';
   rux=ubar*xbar'+nt*igam(nl+1:nl+nm,1:nl);
   rxx=xbar*xbar'+nt*igam(1:nl,1:nl);
   psi=inv(rxx+alp);

   a=(ryx-b*rux)*psi;
   alapsi=a'*lam*a+nk*psi;
   alp=diag(nk./diag(alapsi));

   if infer_nu==1
      ruu=ubar*ubar'+nt*igam(nl+1:nl+nm,nl+1:nl+nm);
      nu=nt./diag(ruu);
   end
   if plotflag
   subplot(3,3,8);plot(1./nu);
   end
end
cy0=a*rxx*a';
cy_reg=trace(rxx*psi)*diag(1./diag(lam));
cy=cy0+cy_reg;
ubar1=ubar;
psiaa=psi;
outalp=diag(alp);

return


