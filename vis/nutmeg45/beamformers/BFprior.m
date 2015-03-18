function [out]=nut_BFprior(Lp,data, flags) %---------------------------------------------------------
% note from jmz:  this was an idea I was playing around with in 2008.  it
% didn't work as-is, but I have several reasons now why that was the case,
% and ideas how to improve.  Haven't had time to try them out, but might as
% well includes this in svn version for now

% Lp : lead field ( channel X 3 )

% Rprior is covariance, not precision
bfprior=flags.bfprior;
nk=size(Lp,1);
nv=size(Lp,3);
nf=size(Lp,2);
t0=dsearchn(data.latency,0);
% yall=bolts.meg'*1e15;
yall=data.y';

% yall=yall-mean(yall(:,1:t0-1),2)*ones(1,size(yall,2));
% ypre=yall(:,1:t0-1);
% y=yall(:,t0:bfprior.timept(2));

nt=size(yall,2);
ypre=yall(:,1:t0-1);
ypre=ypre-mean(ypre,2)*ones(1,t0-1);
y=yall(:,t0:end);
y=y-mean(y,2)*ones(1,nt-t0+1);
yall=yall-mean(yall,2)*ones(1,nt);

% y=yall;
% alpha=1;  %% should this be scaled in some way to match the data??
Rprior=bfprior.alpha*bfprior.Rprior;  
iRpri=nut_inv(diag(Rprior),'tikhonov');
Rpriorb=repmat(Rprior,1,nf);
Rpriora=reshape(Rpriorb',nv*nf,1);
clear Rpriorb
iRpria=nut_inv(diag(Rpriora),'tikhonov');
% nt=size(y,2);
Ryy=y*y';
iRyy=nut_inv(Ryy,'tikhonov');
reg=mean(Rprior)*ones(size(Rprior));% as if prior was scalar*identity for source cov
ireg=nut_inv(diag(reg),'tikhonov');

Ryall=yall*yall';
iRyall=nut_inv(Ryall,'tikhonov');
[u,q,v]=svd(Ryall);
q=diag(q);
qinv=zeros(size(q));
qinv(flags.signalspace)=1./q(flags.signalspace);
iESRyall=v*diag(qinv)*u';

nem=bfprior.nem;
nm=bfprior.nm_pre;
lam_init=0;
b_init=0;
[b,lam,bet,ubar,psibb]=nut_vbfa(ypre,nm,nem,b_init,lam_init,0);
clear bet ubar psibb
Rnn=(b*b'+inv(lam));
iRnn=nut_inv(Rnn,'tikhonov');
Rypre=ypre*ypre';
iRypre=nut_inv(Rypre,'tikhonov');

out.W1=zeros(nk,nf,nv,47);
out.lqv=zeros(nv,1);

ldlam=sum(log(diag(lam/(2*pi))));
ldrnn=sum(log(diag(iRnn/(2*pi))));
Lp1=squeeze(Lp(:,1,:));
Lp2=squeeze(Lp(:,2,:));
Lpa = reshape(Lp,nk,nf*nv);

% Wv11(:,1,:)=nut_inv(Lp1'*iRnn*Lp1,'tikhonov')*Lp1'*iRnn;
% Wv11(:,2,:)=nut_inv(Lp2'*iRnn*Lp2,'tikhonov')*Lp2'*iRnn;
Wv1(:,1,:)=nut_inv(Lp1'*iRnn*Lp1+iRpri,'tikhonov')*Lp1'*iRnn;
Wv1(:,2,:)=nut_inv(Lp2'*iRnn*Lp2+iRpri,'tikhonov')*Lp2'*iRnn;
Wv2=reshape(nut_inv(Lpa'*iRnn*Lpa+iRpria,'tikhonov')*Lpa'*iRnn,nf,nv,nk);
Wv3(:,1,:)=diag(Rprior)*Lp1'*nut_inv(Rnn+Lp1*diag(Rprior)*Lp1','tikhonov');
Wv3(:,2,:)=diag(Rprior)*Lp2'*nut_inv(Rnn+Lp2*diag(Rprior)*Lp2','tikhonov');
Wv4=reshape(diag(Rpriora)*Lpa'*nut_inv(Rnn+Lpa*diag(Rpriora)*Lpa','tikhonov'),nf,nv,nk);

Wv14(:,1,:)=nut_inv(Lp1'*iESRyall*Lp1+iRpri,'tikhonov')*Lp1'*iESRyall;
Wv14(:,2,:)=nut_inv(Lp2'*iESRyall*Lp2+iRpri,'tikhonov')*Lp2'*iESRyall;
Wv15=reshape(nut_inv(Lpa'*iESRyall*Lpa+iRpria,'tikhonov')*Lpa'*iESRyall,nf,nv,nk);
% Wv16(:,1,:)=nut_inv(Lp1'*iRnn*Lp1+iRpri,'tikhonov')*Lp1'*iESRyall;
% Wv16(:,2,:)=nut_inv(Lp2'*iRnn*Lp2+iRpri,'tikhonov')*Lp2'*iESRyall;
% Wv17=reshape(nut_inv(Lpa'*iRnn*Lpa+iRpria,'tikhonov')*Lpa'*iESRyall,nf,nv,nk);
% Wv18(:,1,:)=diag(Rprior)*Lp1'*nut_inv(Ryall+Lp1*diag(Rprior)*Lp1','tikhonov');
% Wv18(:,2,:)=diag(Rprior)*Lp2'*nut_inv(Ryall+Lp2*diag(Rprior)*Lp2','tikhonov');
% Wv19=reshape(diag(Rpriora)*Lpa'*nut_inv(Ryall+Lpa*diag(Rpriora)*Lpa','tikhonov'),nf,nv,nk);

% Wv31(:,1,:)=(nut_pinv(Lp1'*iRyy*Lp1)-diag(Rprior))*Lp1'*iRyy; 
% Wv31(:,2,:)=(nut_pinv(Lp2'*iRyy*Lp2)-diag(Rprior))*Lp2'*iRyy; 
% Wv32=reshape((nut_pinv(Lpa'*iRyy*Lpa)-diag(Rpriora))*Lpa'*iRyy,nf,nv,nk); 
% Wv41(:,1,:)=nut_pinv(Lp1'*iRyy*Lp1)*Lp1'*data.InvES - diag(Rprior)*Lp1'*iRyy;
% Wv41(:,2,:)=nut_pinv(Lp2'*iRyy*Lp2)*Lp2'*data.InvES - diag(Rprior)*Lp2'*iRyy;
% Wv42=reshape(nut_pinv(Lpa'*iRyy*Lpa)*Lpa'*data.InvES - diag(Rpriora)*Lpa'*iRyy,nf,nv,nk);
Wv43(:,1,:)=(nut_pinv(Lp1'*iRyy*Lp1)-diag(Rprior))*Lp1'*iESRyall; 
Wv43(:,2,:)=(nut_pinv(Lp2'*iRyy*Lp2)-diag(Rprior))*Lp2'*iESRyall; 
Wv44=reshape((nut_pinv(Lpa'*iRyy*Lpa)-diag(Rpriora))*Lpa'*iESRyall,nf,nv,nk); 


out.W1(:,:,:,1)=permute(Wv1,[3 2 1]);
out.W1(:,:,:,3)=permute(Wv3,[3 2 1]);
out.W1(:,:,:,14)=permute(Wv14,[3 2 1]);
% out.W1(:,:,:,16)=permute(Wv16,[3 2 1]);
% out.W1(:,:,:,18)=permute(Wv18,[3 2 1]);
% out.W1(:,:,:,31)=permute(Wv31,[3 2 1]);
% out.W1(:,:,:,41)=permute(Wv41,[3 2 1]);
out.W1(:,:,:,43)=permute(Wv43,[3 2 1]);

out.W1(:,:,:,2)=permute(Wv2,[3 1 2]);
out.W1(:,:,:,4)=permute(Wv4,[3 1 2]);
out.W1(:,:,:,15)=permute(Wv15,[3 1 2]);
% out.W1(:,:,:,17)=permute(Wv17,[3 1 2]);
% out.W1(:,:,:,19)=permute(Wv19,[3 1 2]);
% out.W1(:,:,:,32)=permute(Wv32,[3 1 2]);
% out.W1(:,:,:,42)=permute(Wv42,[3 1 2]);
out.W1(:,:,:,44)=permute(Wv44,[3 1 2]);

% s_th = transpose(squeeze(out.W1(:,1,:,1)))*y;
% s_ph = transpose(squeeze(out.W1(:,2,:,1)))*y;
% Rsth=s_th*s_th';
% Rsph=s_ph*s_ph';
% sa=[s_th; s_ph];
% Rsa=sa*sa';
% FRssF=Lp1*Rsth*Lp1'+Lp2*Rsph*Lp2';
FRsspriF=Lp1*diag(Rprior)*Lp1'+Lp2*diag(Rprior)*Lp2';

for iv=1:nv  %for now, nasty for-loop over voxels
    Lpv=squeeze(Lp(:,:,iv));
    sbar=squeeze(out.W1(:,:,iv,1))'*y;  % how choose which weight to initialize?  iterate??
    Rsy=sbar*y';
    Rss=sbar*sbar';  % plus gam?
%     Rnnv=Ryy-2*Lpv*Rsy+Lpv*Rss*Lpv';
%     iRnnv=nut_inv(Rnnv,'tikhonov');
    Rsspri=Rprior(iv)*eye(nf);
    regpri=reg(iv)*eye(nf);
    iRsspri=inv(Rsspri);
    iregpri=inv(regpri);
    Ryyr=Ryy-2*Lpv*Rsy;
    iRyyr=nut_inv(Ryyr,'tikhonov');

%     ldphi=sum(log(svd(iRsspri)));
    %initialize sbar
    Wv5=nut_pinv(Lpv'*iRnn*Lpv+iRsspri)*Lpv'*iRnn; %min norm + sefa (Rnn)
    Wv6=nut_pinv(Lpv'*iRypre*Lpv+iRsspri)*Lpv'*iRypre; %min norm 
    Wv7=Rsspri*Lpv'*nut_inv(Rnn+Lpv*Rsspri*Lpv','tikhonov');  % dale;  is Rnnv init correctly?
    Wv8=Rsspri*Lpv'*nut_inv(Rypre+Lpv*Rsspri*Lpv','tikhonov');  % dale;  is Rnnv init correctly?
    Wv9=(nut_pinv(Lpv'*iRyy*Lpv)-Rsspri)*Lpv'*iRnn; 
    Wv10=Rsspri*Lpv'*inv(Rnn+FRsspriF);
    Wv11=(nut_pinv(Lpv'*iRyy*Lpv)-Rsspri)*Lpv'*iRyy; 
    
    Wv12=(nut_pinv(Lpv'*iRyy*Lpv)-Rsspri)*Lpv'*data.InvES;
%     Wv13=nut_pinv(Lpv'*iRnn*Lpv+iRsspri)*Lpv'*data.InvES; %min norm + sefa (Rnn)
    Wv20=(nut_pinv(Lpv'*iRyall*Lpv)-Rsspri)*Lpv'*data.InvES;
    Wv21=nut_pinv(Lpv'*iRyall*Lpv+iRsspri)*Lpv'*data.InvES; %min norm + sefa (Rnn)
%     Wv22=nut_pinv(Lpv'*iRyall*Lpv+iRsspri)*Lpv'*iRyall; %min norm + sefa (Rnn)
    Wv23=nut_pinv(Lpv'*iRyall*Lpv+iRsspri)*Lpv'*iESRyall; %min norm + sefa (Rnn)
%     Wv24=nut_pinv(Lpv'*iRnn*Lpv+iRsspri)*Lpv'*iESRyall; %min norm + sefa (Rnn)
%     Wv25=Rsspri*Lpv'*nut_inv(Ryall+Lpv*Rsspri*Lpv','tikhonov');  % dale;  is Rnnv init correctly?
%     Wv26=Rsspri*Lpv'*nut_inv(Ryall+FRsspriF);
    
%     Wv27=(nut_pinv(Lpv'*iRyy*Lpv)-Rsspri)*Lpv'*iRyall; 
    Wv28=(nut_pinv(Lpv'*iRyy*Lpv)-Rsspri)*Lpv'*iESRyall; 
%     Wv29=(nut_pinv(Lpv'*iRyall*Lpv)-Rsspri)*Lpv'*iRyy; 
%     Wv30=(nut_pinv(Lpv'*iRnn*Lpv)-Rsspri)*Lpv'*iRyy; 
%     Wv33=nut_pinv(Lpv'*iRyy*Lpv)*Lpv'*data.InvES - Rsspri*Lpv'*iRyy;
%     Wv34=nut_pinv(Lpv'*iRyy*Lpv)*Lpv'*iRyall - Rsspri*Lpv'*iRyy;
    Wv35=nut_pinv(Lpv'*iRyy*Lpv)*Lpv'*data.InvES - Rsspri*Lpv'*iRnn;
    Wv36=nut_pinv(Lpv'*iRyy*Lpv)*Lpv'*data.InvES - Rsspri*Lpv'*iRypre;
%     Wv37=nut_pinv(Lpv'*iRyall*Lpv)*Lpv'*iRyall - Rsspri*Lpv'*iRyy;
%     Wv38=nut_pinv(Lpv'*iRyall*Lpv)*Lpv'*iRyall - Rsspri*Lpv'*iRyall;
%     Wv39=nut_pinv(Lpv'*iRyall*Lpv)*Lpv'*iESRyall - Rsspri*Lpv'*iRyy;
    Wv40=nut_pinv(Lpv'*iRyall*Lpv)*Lpv'*iESRyall - Rsspri*Lpv'*iRyall;
    Wv45=(nut_pinv(Lpv'*iRyall*Lpv)-Rsspri)*Lpv'*iESRyall; 
    Wv46=(nut_pinv(Lpv'*iRyall*Lpv)-Rsspri)*Lpv'*iRnn; 
    Wv47=(nut_pinv(Lpv'*iRyall*Lpv)-Rsspri)*Lpv'*iRypre; 
 

    out.W1(:,:,iv,5)=Wv5';
    out.W1(:,:,iv,6)=Wv6';
    out.W1(:,:,iv,7)=Wv7';
    out.W1(:,:,iv,9)=Wv9';
    out.W1(:,:,iv,10)=Wv10';
    out.W1(:,:,iv,11)=Wv11';
    out.W1(:,:,iv,8)=Wv8';
    out.W1(:,:,iv,12)=Wv12';
%     out.W1(:,:,iv,13)=Wv13';
    out.W1(:,:,iv,20)=Wv20';
    out.W1(:,:,iv,21)=Wv21';
%     out.W1(:,:,iv,22)=Wv22';
    out.W1(:,:,iv,23)=Wv23';
%     out.W1(:,:,iv,24)=Wv24';
%     out.W1(:,:,iv,25)=Wv25';
%     out.W1(:,:,iv,26)=Wv26';
%     out.W1(:,:,iv,27)=Wv27';
    out.W1(:,:,iv,28)=Wv28';
%     out.W1(:,:,iv,29)=Wv29';
%     out.W1(:,:,iv,30)=Wv30';
%     out.W1(:,:,iv,33)=Wv33';
%     out.W1(:,:,iv,34)=Wv34';
    out.W1(:,:,iv,35)=Wv35';
    out.W1(:,:,iv,36)=Wv36';
%     out.W1(:,:,iv,37)=Wv37';
%     out.W1(:,:,iv,38)=Wv38';
%     out.W1(:,:,iv,39)=Wv39';
    out.W1(:,:,iv,40)=Wv40';
    out.W1(:,:,iv,45)=Wv45';
    out.W1(:,:,iv,46)=Wv46';
    out.W1(:,:,iv,47)=Wv47';

    %         sbar=Wv2*y;
    %         Rss=sbar*sbar'+nt*igam;
    %         Rss=sbar*sbar';
    %         like(iter)=.5*nt*(ldlam+ldphi-ldgam)-.5*sum(sum(y.*(lam*y)))+.5*sum(sum(sbar.*(gam*sbar)));
    %         like(iter)=.5*nt*(ldrnn+ldphi-ldgam)-.5*sum(sum(y.*(iRnn*y)))+.5*sum(sum(sbar.*(gam*sbar)));
    %         if bfprior.plot==1
    %             figure(5);plot(like)
    %         end
    %     end
    %     out.lqv(iv)=like(end);
    if(~mod(iv,20))  % update progress every 20th voxel
        disp(['Please wait. BFprior has finished ' num2str(iv) ' out of ' num2str(nv) ' voxels.']);
    end
end

out.flags=bfprior;
