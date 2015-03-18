function tb1aout=nut_nsefaloc_sub(Lp,data,flags);
% function [lqv,g]=algo1a(y,f,sig0,lam0,nem,phi,nm_x,plotflag);
% model y=F*G*phi+A*x+v

% global bolts 
global ndefaults
tbflags=flags.tb;

% first do regular sefa to get temporal basis functions phi
% t0=dsearchn(bolts.latency,0);

% ypre=bolts.meg(1:t0-1,:)';
% ypre=ypre-mean(ypre,2)*ones(1,size(ypre,2));
% ypost=bolts.meg(t0:end,:)';
% ypost=ypost-mean(ypost,2)*ones(1,size(ypost,2));
% yall=bolts.meg'*1e15; %can't assume data is in T.

%changed so now 'active' and 'control' can be any time, pre or post t0
yall=double(data.y');  
yall=yall-mean(yall(:,tbflags.timeptc(1):tbflags.timeptc(2)),2)*ones(1,size(yall,2));
yall=yall/std(yall(:));
% yall=yall-mean(yall(:,1:t0-1),2)*ones(1,size(yall,2));
% if tbflags.timept(1)>=t0
%     error('you must select a time interval that includes some pre-stim');
% end
% if tbflags.timept(2)<=t0
%     error('you must select a time interval that includes some post-stim');
% end
ypre=yall(:,tbflags.timeptc(1):tbflags.timeptc(2));
ypost=yall(:,tbflags.timepta(1):tbflags.timepta(2));

nk=size(yall,1);
if ~isfield(tbflags,'nl')
%         tbflags.nl=max(bolts.params.signalspace);
        error('need to specify tbflags.nl');
end
if ~isfield(tbflags,'nm')
    tbflags.nm=ceil(nk/3);
end
if ~isfield(tbflags,'nm_x')
    tbflags.nm_x=ceil(nk*.1);
end

nem_sefa=50;
nem1_sefa=80;
nem=tbflags.nem;
[a,b,lam,nu,alp,bet,xbar,ybar,ubar1,cy,like,Ws,psiaa,psibb]=nut_sefa(ypre,ypost,tbflags.nl,tbflags.nm,nem_sefa,nem1_sefa,0);
sig0=b*b'+diag(1./diag(lam));
lam0=lam;
if tbflags.usephiall
%     phi=www*yall;
    y=yall;
else
    y=ypost;
%     phi=xbar;
end
phi=Ws*y;
tb1aout.Ws=Ws;
clear a b lam nu alp bet xbar ybar ubar1 cy like Ws psiaa psibb


f=permute(Lp,[1 3 2]);
clear Lp;
nt=size(y,2);
nv=size(f,2);
nf=size(f,3);
nl=tbflags.nl;
nu=1;


% compute correlation matrices and hyperparameters
ryy=y*y';
ryp=y*phi';
rpp=phi*phi';
irpp=inv(rpp);

% mg1=inv((ryy+sig0-ryp*irpp*ryp')/(nt+nu));
noisecov=ryy-ryp*irpp*ryp';
Sinv=inv((ryy+sig0-ryp*irpp*ryp')/(nt+nu));
% mg2=Sinv*ryp*irpp;

% compute likelihood for all voxels by em

tb1aout.lmap=zeros(nv,2);
tb1aout.g=zeros(nf,nl,nv,2);
tb1aout.W1=zeros(nk,nf,nv,2);
% lam0=lam;

tic
for iv=1:nv
%     if(~mod(iv,10))  % update progress every 10th voxel
%          if(ndefaults.bf.viewwaitbar)
%             waitbar(iv/nv,bolts.barhandle)
%          end
%     end

    f1=squeeze(f(:,iv,:));

    %-------

    % algo 1
    g1=nut_pinv(f1'*Sinv*f1)*f1'*Sinv*ryp*irpp;
    fg1=f1*g1;
    tt=ryp*fg1';
    sig=(ryy+sig0-tt-tt'+fg1*rpp*fg1')/(nt+nu);
    d=svd(sig);
    ldisig=-sum(log(d/(2*pi)));
    tb1aout.lmap(iv,1)=.5*(nt+nu)*ldisig-.5*(nt+nu)*nk;
    tb1aout.g(:,:,iv,1)=g1;
    tb1aout.W1(:,:,iv,1)=(g1*tb1aout.Ws)';

    %-------

    % algo 1a

    lam=lam0;
    flf=nut_pinv(f1'*lam*f1)*f1'*lam;
    g1=flf*ryp*inv(rpp);

    ryy=y*y';
    y1=y-f1*g1*phi;
    ryy1=y1*y1';
    [p d]=svd(ryy1/nt);d=diag(d);
    a=p*diag(sqrt(d));
    a=a(:,1:tbflags.nm_x);
    alp=diag(1./diag(a'*lam*a/nk));
    alp=min(diag(alp))*diag(ones(tbflags.nm_x,1));
    psi=eye(tbflags.nm_x)/nt;

    % em iteration

    % nem=5;
    like=zeros(nem,1);
    alapsi=a'*lam*a+nk*psi;

    for iem=1:nem
        gam=alapsi+eye(tbflags.nm_x);
        igam=inv(gam);
        xbar=igam*a'*lam*y1;

        ldlam=sum(log(diag(lam/(2*pi))));
        ldgam=sum(log(svd(gam)));
        ldalp=sum(log(diag(alp)));
        ldpsi=sum(log(svd(psi)));
        if (tbflags.plot | iem==nem)
            like(iem)=.5*nt*(ldlam-ldgam)-.5*sum(sum(y1.*(lam*y1)))+.5*sum(sum(xbar.*(gam*xbar)))+.5*nk*(ldalp+ldpsi);
        end
        if (tbflags.plot)
            subplot(3,3,8);plot((1:iem)',like(1:iem));title('algo1a');
            subplot(3,3,9);plot([1./mean(a.^2,1)' diag(alp)]);
            %subplot(3,3,9);plot(1./diag(lam));
            drawnow;
        end

        ryx=y*xbar';
        rxx=xbar*xbar'+nt*igam;
        %ryx=(ryy-ryp*g1'*f1')*lam*a*igam;
        %ryy1=ryy-f1*g1*ryp'-ryp*g1'*f1'+f1*g1*rpp*g1'*f1';
        %rxx=igam*a'*lam*ryy1*lam*a*igam+nt*igam;
        psi=inv(rxx+alp);

        rpx=phi*xbar';
        g1=flf*(ryp-ryx*psi*rpx')*inv(rpp-rpx*psi*rpx');
        a=(ryx-f1*g1*rpx)*psi;

        y1=y-f1*g1*phi;
        ryy1=y1*y1';
        ryx1=y1*xbar';
        lam=diag(nt./diag(ryy1-a*ryx1'));
        flf=nut_pinv(f1'*lam*f1)*f1'*lam;
        alapsi=a'*lam*a+nk*psi;
        alp=diag(nk./diag(alapsi));
    end
    tb1aout.lmap(iv,2)=like(iem);
    tb1aout.g(:,:,iv,2)=g1;
    tb1aout.phi=phi;
    tb1aout.W1(:,:,iv,2)=(g1*tb1aout.Ws)';  %weight to apply to y to get s (y=Fs + stuff)
    tb1aout.tbflags=tbflags;

    %-------
    if tbflags.plot
        subplot(3,3,4);plot((1:iv)',lqv(1:iv,1));ylabel('log q(voxel)');
        title(['algo1   voxel=' num2str(iv)]);
        eg1=squeeze(mean(mean(abs(g(:,:,1:iv,1)),1),2));
        eg2=squeeze(mean(mean(abs(g(:,:,1:iv,2)),1),2));
        subplot(3,3,5);plot((1:iv)',[eg1 eg2]);ylabel('g(voxel)');xlabel('voxel');

        subplot(3,3,7);plot((1:iv)',lqv(1:iv,2));ylabel('log q(voxel)');
        title(['algo1a   voxel=' num2str(iv)]);
        drawnow
    end
%     toc,disp(iv)
    if(~mod(iv,20))  % update progress every 20th voxel
        disp(['Please wait. NSEFALoc has finished ' num2str(iv) ' out of ' num2str(nv) ' voxels.']);
    end
end



