% SEFAICA
% HT Attias, Golden Metallic 11/8/04


% Denoise and separate MEG signals
%
% Step 1 -- initialization: estimate initial noise model from pre-stimulus 
%           data by FA
% Step 2 -- denoising: estimate noise and signal models and estimate 
%           post-stimulus factors from pre+post-stimulus data by SEFA
% Step 3 -- separation: estimate mixing model and post-stimulus sources 
%           by ICAMOG from factors estimated in step 2



function [yproj,zbar,a,cy,cysep,BBTL,perf,b,lam,alp,g]=sefaica(ypre,ypost,xpost,flg,nl_x,dflag);

if ~exist('flg','var')
   flg = 0; % flg=0 uses step3, flg=1 uses Jade/MRMI-SIG
end

flg2 = 1;

if flg2
   mpre = mean(ypre');
   mpost = mean(ypost');
   ypost = ypost - repmat(mpost',1,size(ypost,2));
   ypre = ypre - repmat(mpre',1,size(ypre,2));
end

x0=xpost;y0=ypost;
u0=ypre;

sigy=sqrt(mean(mean(y0.^2)));
y0=y0/sigy;x0=x0/sigy;u0=u0/sigy;

nk=size(y0,1);nt=size(y0,2);


%nem_u = number of EM iterations for step 1
%nl_u  = number of factors in the noise model
%nem_x = number of EM iterations for step 2
%nl_x  = number of factors in the signal model
%nem_z = number of EM iterations for step 3
%eta   = learning rate for mixing model
 
nem_u=200;nl_u=75;
nem_x=200;%nl_x=2;
nem_z=1500;eta=.005;

%ns=2;w=[.5 .5]';nu=[.27 3.8]';mu=[0 0]';% this is randn*randn
%ns=3; w=[.31 .52 .17]';nu=[.36 3.38 134.2]';
%ns=3;w=[1 1 1]'/3;nu=[4 2 4]';mu=[-3.5 0 3.5]';
if strcmp(dflag,'peaky')
   ns=3;w=[0.4 0.2 0.4]';nu=[1 50 1000]';mu=[0 0 0]';
elseif strcmp(dflag,'sin')
   ns=3;w=[.2 .6 .2]';nu=[80 2 80]';mu=[-1.3 0 1.3]';% this is sin(50*t)
end


% estimate noise model (FA)

figure(3); clf
[b,lam,like_u,bet]=step1(u0,nl_u,nem_u);

hsub=subplot(3,3,4);plot((1:nl_u)',sqrt([mean(b.^2,1)' 1./bet]));
title('1/beta');
sig=b*b'+diag(1./lam);
prec_u=inv(sig);


% estimate signal model (SEFA)

a_init=randn(nk,nl_x);
[a,b,lam,xbar,like_x,alp,bet,rxx,raa]=step2(u0,y0,nl_x,nem_x,a_init,lam,b,bet,nl_u);

perf = [like_u(end) like_x(end)];% a_init(:)'];

BBTL = b*b' + diag(1./lam);

hsub=subplot(3,3,6);plot((1:nl_x)',sqrt([mean(a.^2,1)' 1./alp]));
title('1/alpha');
hsub=subplot(3,3,5);plot((1:nl_u)',sqrt([mean(b.^2,1)' 1./bet]));
title('1/beta');


% estimate mixing matrix (ICAMOG)

figure(4); clf
if flg == 0
   [g,zbar,like_z]=step3(xbar,nem_z,eta,ns,w,nu,mu);
else
   %[g,zbar] = Jade(xbar,eye(size(xbar,1))); like_z = NaN; g = g';
   [g,zbar,wmtrx]=mrmi_inst_sig(xbar,0.1,size(xbar,2),500,[0.25 1],randn(size(xbar,1)));
   clf, plot(wmtrx')
   keyboard
end


% project on sensors

ig=inv(g);
h=a*ig;
rzz=g*rxx*g';
rhh=ig'*raa*ig;

%ybar=a*xbar;
ybar=h*zbar;

%cy0=a*rxx*a';
cy0=h*rzz*h';
%cy_reg=trace(rxx*raa)*diag(1./lam);
cy_reg=trace(rzz*rhh)*diag(1./lam);
cy=cy0+cy_reg;

cysep=zeros(nk,nk,nl_x);
for il=1:nl_x
%   cysep(:,:,il)=rxx(il,il)*a(:,il)*a(:,il)'+rxx(il,il)*raa(il,il)*diag(1./lam);
   cysep(:,:,il)=rzz(il,il)*h(:,il)*h(:,il)'+rzz(il,il)*rhh(il,il)*diag(1./lam);
end



% project sources on sensor space

figure(5); clf
yproj=zeros(nk,nt,nl_x);
for il=1:nl_x
   yproj(:,:,il)=h(:,il)*zbar(il,:);
   if il <= 4
      hsub=subplot(2,2,il);plot(yproj(1:5:end,:,il)');
      title(['projected source ' num2str(il)]);
      %axis([1 400 -12 12]);
      axis('tight')
   end
end

hsub=subplot(2,2,3);
plot([zbar(1,:)'/std(zbar(1,:))]);
axis('tight')
title('source 1')
if nl_x > 1
   hsub=subplot(2,2,4);
   plot([zbar(2,:)'/std(zbar(2,:))]);
   axis('tight')
   title('source 2')
end


% compute SNR
if ~isempty(xpost)
   snr1=10*log10(mean(mean(x0.^2))/mean(mean((y0-x0).^2)));
   snr2=10*log10(mean(mean(x0.^2))/mean(mean((ybar-x0).^2)));
end

tit0='Clean signal';
%tit1=['Noisy signal: SNR=' num2str(snr1)];
%tit2=['SEFA: SNR=' num2str(snr2)];
tit1=['Original observations'];
tit2=['Observations after Bayesian denoising'];

figure(6);

%hsub=subplot(3,2,1);plot(x0(1:5:nk,:)');title(tit0);axis([1 400 -12 12]);
%hsub=subplot(3,2,2);plot(y0(1:5:nk,:)');title(tit1);axis([1 400 -12 12]);
%hsub=subplot(3,2,3);plot(ybar(1:5:nk,:)');title(tit2);axis([1 400 -12 12]);
hsub=subplot(2,1,1);set(hsub,'FontSize',12);plot(y0(1:5:nk,:)');axis tight;title(tit1,'FontSize',12);
hsub=subplot(2,1,2);set(hsub,'FontSize',12);plot(ybar(1:5:nk,:)');axis tight;title(tit2,'FontSize',12);xlabel('Time (ms)','FontSize',12)

yproj = yproj*sigy;
cy = cy*sigy^2/size(ypost,2); %%%%% + (mpost-mpre)'*(mpost-mpre); % TRY THIS OUT WITH REAL DATA (assumes interference plus sensor noise is zero mean so that mpre is sensor offsets)
cysep = cysep*sigy^2/size(ypost,2); % not sure just yet how to apply the idea above to cysep
%xbar = xbar*sigy;
BBTL = BBTL*sigy^2;
cy = cy*sigy^2;


% if flg2
%    for i=1:size(yproj,3)
%       yproj(:,:,i) = yproj(:,:,i) + repmat(mpost',1,size(ypost,2));
%    end
% end

return





%===============================


% y=Bu+v   estimate B and Lam (=precision of v)  



function [b,lam,mlike,bet]=step1(y,nl,nem);

nk=size(y,1);
nt=size(y,2);
%disp([nk nt]);
 
ryy=y*y';
[p d q]=svd(ryy/nt);
d=diag(d);
b=p*diag(sqrt(d));
b=b(:,1:nl);
lam=1./diag(ryy/nt);

bet=ones(nl,1);

like=zeros(nem,1);
rbb=eye(nl)/nt;

for iem=1:nem
   dbet=diag(bet);ldbet=sum(log(bet));
   dlam=diag(lam);ldlam=sum(log(lam));

   gam=b'*dlam*b+eye(nl)+nk*rbb*0;
   igam=inv(gam);   
   ubar=igam*b'*dlam*y;
   ryu=y*ubar';
   ruu=ubar*ubar'+nt*igam;

   [p d q]=svd(gam);ldgam=sum(log(diag(d)));
   temp1=-.5*ldgam*ones(1,nt)+.5*sum(ubar.*(gam*ubar),1);
   temp2=.5*ldlam*ones(1,nt)-.5*lam'*(y.^2);
   f=temp1+temp2;	
   f3=.5*nl*ldlam+.5*nk*ldbet-.5*trace(b'*dlam*b*dbet);
   like(iem)=mean(f)+f3/nt;
							   
   betbar=ruu+dbet;
   ibetbar=inv(betbar);
   b=ryu*ibetbar;

   ilam=diag(ryy-b*ryu')/(nt+nl);
   lam=1./ilam;
   dlam=diag(lam);
       					
   bet=1./(diag(b'*dlam*b)/nk+diag(ibetbar));

   if rem(iem,50) == 0
      hsub=subplot(3,3,1);plot((1:iem)',like(1:iem));title('Step1   likelihood')
      hsub=subplot(3,3,4);plot((1:nl)',sqrt([mean(b.^2,1)' 1./bet]));title('1/beta');
      hsub=subplot(3,3,7);plot(1./lam);title('1/lam');
      drawnow;
   end
   
   rbb=ibetbar;
end

mlike=like(iem);


return





%===============================


% y=Ax+Bu+v   estimate A,B,Lam



function [a,b,lam,xbar,mlike,alp,bet,rxx,raa]=step2(y0,y,nl,nem,a_init,lam,b,bet,nl0);

nk=size(y,1);
nt=size(y,2);
nt0=size(y0,2);
%disp([nk nt]);
   
nlz=nl0+nl;
ntz=nt0+nt;

ryy=y*y';
[p d q]=svd(ryy/nt);
d=diag(d);
a=p*diag(sqrt(d));
a=a(:,1:nl);
%a=a_init;

ryy0=y0*y0';


%%% alp?
sig_u=b*b'+diag(1./lam);
temp=sum(diag(ryy/nt-sig_u))/sum(1./lam);
alp0=1/(temp/nl);
if alp0<=0;alp0=1;end
%disp(['alp0=' num2str(alp0)]);
%%%%%%

alp=ones(nl,1)*alp0;
betz=[bet;alp];

like=zeros(nem,1);
rbbz=eye(nlz);

for iem=1:nem
   bz=[b a];
   dbetz=diag(betz);ldbetz=sum(log(betz));
   dalp=diag(alp);dbet=diag(bet);
   dlam=diag(lam);ldlam=sum(log(lam));

   rbb=rbbz(1:nl0,1:nl0);
   gam0=b'*dlam*b+eye(nl0)+nk*rbb*0;
   igam0=inv(gam0);   
   ubar0=igam0*b'*dlam*y0;
   ryu0=y0*ubar0';
   ruu0=ubar0*ubar0'+nt0*igam0;

   [p d q]=svd(gam0);ldgam0=sum(log(diag(d)));
   temp1=-.5*ldgam0*ones(1,nt0)+.5*sum(ubar0.*(gam0*ubar0),1);
   temp2=.5*ldlam*ones(1,nt0)-.5*lam'*(y0.^2);
   f0=temp1+temp2;

   gamz=bz'*dlam*bz+eye(nlz)+nk*rbbz*0;
   igamz=inv(gamz);
   ubarz=igamz*bz'*dlam*y;
   ryuz=y*ubarz';
   ruuz=ubarz*ubarz'+nt*igamz;

   [p d q]=svd(gamz);ldgamz=sum(log(diag(d)));
   temp1=-.5*ldgamz*ones(1,nt)+.5*sum(ubarz.*(gamz*ubarz),1);
   temp2=.5*ldlam*ones(1,nt)-.5*lam'*(y.^2);
   f=temp1+temp2;	
   f3=.5*nlz*ldlam+.5*nk*ldbetz-.5*trace(bz'*dlam*bz*dbetz);
   like(iem)=(nt0*mean(f0)+nt*mean(f)+f3)/ntz;
   %like(iem)=(nt0*mean(f0)+nt*mean(f))/ntz;

   ruuz(1:nl0,1:nl0)=ruuz(1:nl0,1:nl0)+ruu0;
   ryuz(:,1:nl0)=ryuz(:,1:nl0)+ryu0;				   
   betbarz=ruuz+dbetz;
   ibetbarz=inv(betbarz);
   %bz=ryuz*ibetbarz;
   
   %>>>
   ryx=ryuz(:,nl0+1:nlz);
   rxx=ruuz(nl0+1:nlz,nl0+1:nlz);
   aa=ryx*inv(rxx+dalp);
   bz=[b aa];
   %>>>

   ilam=diag(ryy0+ryy-bz*ryuz')/(ntz+nlz);
   %lam=1./ilam;
   dlam=diag(lam);

   betz=1./(diag(bz'*dlam*bz)/nk+diag(ibetbarz));
   b=bz(:,1:nl0);
   a=bz(:,nl0+1:nlz);
   bet=betz(1:nl0);
   alp=betz(nl0+1:nlz);
   
   if rem(iem,50) == 0
      hsub=subplot(3,3,2);plot((1:iem)',like(1:iem));title('Step2   likelihood')
      hsub=subplot(3,3,5);plot((1:nl0)',sqrt([mean(b.^2,1)' 1./bet]));title('1/beta');
      hsub=subplot(3,3,8);plot(1./lam);title('1/lam');
      hsub=subplot(3,3,6);plot((1:nl)',sqrt([mean(a.^2,1)' 1./alp]));title('1/alp');
      drawnow;
   end

   rbbz=ibetbarz;
end

mlike=like(iem);
xbar=ubarz(nl0+1:nlz,:);
rxx=ruuz(nl0+1:nlz,nl0+1:nlz);
raa=rbbz(nl0+1:nlz,nl0+1:nlz);


return





%===============================


% x=Gy (ICA)



function [g,x,like]=step3(y,nem,eta,ns,w,nu,mu);

nl=size(y,1);
nt=size(y,2);
%disp([nl nt]);

c=y*y'/nt;
[p d q]=svd(c);
g=diag(1./sqrt(diag(d)))*p';
x=g*y;

like=zeros(nem,1);logp=zeros(nl,nt);phi=zeros(nl,nt);
logw0=(.5*log(nu/(2*pi))+log(w))*ones(1,nt);

for it=1:nem
   g1=g;

   for il=1:nl
      x1=ones(ns,1)*x(il,:)-mu*ones(1,nt);
      logw1=logw0-.5*(nu*ones(1,nt)).*(x1.^2);
      w1=exp(logw1-ones(ns,1)*max(logw1,[],1));z1=sum(w1,1);
      w1=w1./(ones(ns,1)*z1);
      temp=zeros(ns,nt);i0=find(w1>0);temp(i0)=log(w1(i0));
      logp(il,:)=squeeze(sum(w1.*(logw1-temp),1));
      phi(il,:)=(w1'*nu)'.*x(il,:)-(w1'*(nu.*mu))';
   end

   like(it)=log(abs(det(g1)))+sum(mean(logp,2));

   dg=eye(nl)-phi*x'/nt;
   g=g1+eta*dg*g1;      
   x=g*y;

   if rem(it,50) == 0
      hsub=subplot(2,2,1);plot((1:it)',like(1:it));title('ICAMOG likelihood');
      hsub=subplot(2,2,2);plot(reshape(g,nl^2,1));title('g elements');
      drawnow
   end
end

return





%===============================
