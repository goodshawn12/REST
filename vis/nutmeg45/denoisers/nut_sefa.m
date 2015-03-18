function [a,b,lam,nu,alp,bet,xbar,ybar,ubar1,cy,like,www,psiaa,psibb,outalp,igam]=nut_sefa(ypre,ypost,nl,nm,nem,nem1,dosefa0);

% function [a,b,lam,nu,alp,bet,xbar,ybar,cy]=sefa(ypre,ypost,nl,nm,nem0,nem,nem1);

% xbar = estimates of sources in source space (post stim period)
% ybar = estimates of sources in sensor space (post stim period)
% ubar = estimates of interferences in source space (pre stim period)
% (source estimates) = www * ypost     (where ypost = observations) (same as xbar)
% (cleaned observations) = a * (source estimates) = a * www * ypost (same as ybar)

%[a,b,lam,alp,bet,xbar]=sefa0(ypost,ypre,nl,nm,nem0);

%b_init=b;
%lam_init=lam;
b_init=0;
lam_init=0;
[b,lam,bet,ubar,psibb]=nut_vbfa(ypre,nm,nem,b_init,lam_init,0);

%a_init=a;
a_init=0;
infer_nu=0;
[a,nu,alp,xbar,ubar1,cy,like,www,psiaa,outalp,igam]=nut_vbfa1(ypost,nl,b,lam,nem1,a_init,infer_nu,0);

ybar=a*xbar;

%% this sefa0 uses vbfa/vbfa1 as initialization
if dosefa0 & ~infer_nu
    [a,b,lam,alp,bet,xbar,like]=sefa0_sriinit(y,y0,nl,nm,nem,a_init,b_init,lam_init,alp_init,bet_init);
end

return



% %-----------------------------------------------------------------------
% 
% 
% 
% function [a,lam,alp,xbar]=vbfa(y,nl,nem,a_init,lam_init);
% 
% % vb-em algorithm for inferring the factor analysis model   y = a*x + v
% %
% % y(nk,nt) = data
% % a(nk,nl) = mixing matrix
% % lam(nk,nk) = diagonal noise precision matrix 
% % alp(nl,nl) = diagonal hyperparmaeter matrix 
% % xbar(nl,nt) = posterior means of the factors 
% %
% % nk = number of data points
% % nt = number of time points
% % nl = number of factors
% % nem = number of em iterations
% 
% nk=size(y,1);
% nt=size(y,2);
% 
% % initialize by svd
% 
% ryy=y*y';
% if a_init==0
%    [p d q]=svd(ryy/nt);d=diag(d);
%    a=p*diag(sqrt(d));
%    a=a(:,1:nl);
%    lam=diag(nt./diag(ryy));
% else
%    a=a_init;
%    lam=lam_init;
% end
% 
% alp=diag(1./diag(a'*lam*a/nk));
% alp=min(diag(alp))*diag(ones(nl,1));
% psi=eye(nl)/nt;
% 
% % em iteration
% 
% like=zeros(nem,1);
% alapsi=a'*lam*a+nk*psi;
% 
% for iem=1:nem
%    gam=alapsi+eye(nl);
%    igam=inv(gam);
%    xbar=igam*a'*lam*y;
% 
%    ldlam=sum(log(diag(lam/(2*pi))));
%    ldgam=sum(log(svd(gam)));
%    ldalp=sum(log(diag(alp)));
%    ldpsi=sum(log(svd(psi)));
%    like(iem)=.5*nt*(ldlam-ldgam)-.5*sum(sum(y.*(lam*y)))+.5*sum(sum(xbar.*(gam*xbar)))+.5*nk*(ldalp+ldpsi);
% figure(3)
%    subplot(3,3,1);plot((1:iem)',like(1:iem));title('vbfa');
%    subplot(3,3,4);plot([mean(a.^2,1)' 1./diag(alp)]);
%    subplot(3,3,7);plot(1./diag(lam));
%    drawnow;
% 
%    ryx=y*xbar';
%    rxx=xbar*xbar'+nt*igam;
%    psi=inv(rxx+alp);
% 
%    a=ryx*psi;
%    lam=diag(nt./diag(ryy-a*ryx'));
%    alapsi=a'*lam*a+nk*psi;
%    alp=diag(nk./diag(alapsi));
% end
% 
% return
% 


%-----------------------------------------------------------------------




%-----------------------------------------------------------------------



function [a,b,lam,alp,bet,xbar]=sefa0(y,y0,nl,nm,nem);

% vb-em algorithm for inferring the sefa analysis model   y = a*x + b*u + v
% learn a,b,lam
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
nt0=size(y0,2);
			    
% initialize by svd

ryy=y*y';
[p d q]=svd(ryy/nt);d=diag(d);
a=p*diag(sqrt(d));
a=a(:,1:nl);

ryy0=y0*y0';
[p d q]=svd(ryy0/nt0);d=diag(d);
b=p*diag(sqrt(d));
b=b(:,1:nm);
lam=diag(nt0./diag(ryy0));

alp=diag(1./diag(a'*lam*a/nk));
alp=min(diag(alp))*diag(ones(nl,1));
bet=diag(1./diag(b'*lam*b/nk));
bet=min(diag(bet))*diag(ones(nm,1));

ab=[a b];
alpbet=diag([diag(alp);diag(bet)]);
nlm=nl+nm;
psi=eye(nlm)/(nt+nt0);

% em iteration

like=zeros(nem,1);
alapsi=ab'*lam*ab+nk*psi;

for iem=1:nem
   gam=alapsi+eye(nlm);
   igam=inv(gam);
   xubar=igam*ab'*lam*y;

   b=ab(:,nl+1:nlm);
   psib=psi(nl+1:nlm,nl+1:nlm);
   gam0=b'*lam*b+nk*psib+eye(nm);
   igam0=inv(gam0);
   ubar0=igam0*b'*lam*y0;

   ldlam=sum(log(diag(lam/(2*pi))));
   ldgam=sum(log(svd(gam)));
   ldgam0=sum(log(svd(gam0)));
   ldalpbet=sum(log(diag(alpbet)));
   ldpsi=sum(log(svd(psi)));
   like0=.5*nt0*(ldlam-ldgam0)-.5*sum(sum(y0.*(lam*y0)))+.5*sum(sum(ubar0.*(gam0*ubar0)));
   like(iem)=.5*nt*(ldlam-ldgam)-.5*sum(sum(y.*(lam*y)))+.5*sum(sum(xubar.*(gam*xubar)))+.5*nk*(ldalpbet+ldpsi)+like0;

   subplot(3,3,3);plot((1:iem)',like(1:iem));title('vbfa2');
   subplot(3,3,6);plot([mean(ab.^2,1)' 1./diag(alpbet)]);
   subplot(3,3,9);plot(1./diag(lam));
   drawnow;

   rxuxu=xubar*xubar'+nt*igam;
   ruu0=ubar0*ubar0'+nt0*igam0;
   rxuxu(nl+1:nlm,nl+1:nlm)=rxuxu(nl+1:nlm,nl+1:nlm)+ruu0;
   psi=inv(rxuxu+alpbet);

   ryxu=y*xubar';
   ryu0=y0*ubar0';
   ryxu(:,nl+1:nlm)=ryxu(:,nl+1:nlm)+ryu0;
   ab=ryxu*psi;
   lam=diag((nt+nt0)./diag(ryy+ryy0-ab*ryxu'));
   alapsi=ab'*lam*ab+nk*psi;
   alpbet=diag(nk./diag(alapsi));
end

a=ab(:,1:nl);
b=ab(:,nl+1:nlm);
alp=alpbet(1:nl,1:nl);
bet=alpbet(nl+1:nlm,nl+1:nlm);
xbar=xubar(1:nl,:);

return



%-----------------------------------------------------------------------




