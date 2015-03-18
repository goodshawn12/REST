function out=nut_BF_LFerror(Lp,data,flags);
%
% dipole source localization for the model y = F*s + A*x + B*u + v
% s~N(0,phi) , x~N(0,I), v~N(0,lam)
%
% INPUTS:
% Lp(nk,nf,nv) = lead field
% data.y and data.latency
% flags.sake - list of options for update rules and initialization
%     flags.sake.nem = number of EM iterations, recommended at least 20
%     flags.sake.sbar_init = either ['bf_nut', 'bf_ES', 'bf_reg', 'algo1'], meaning beamformer by nutmeg, bf by eigenspace (recommended), regular beamformer, or algo1 (aka sefaloc1)
%     flags.sake.thisrun = simulation index/number for your own records
%     flags.sake.iryy_init = either ['vbfa','sefa_aa','sefa_cy'], meaning vbfa to poststim (recommended), vbfa1 to poststim using AA', or vbfa1 to poststim using cleaned cov.
%     flags.sake.phi_init = either ['ssT','magic_y'], meaning s*s', or magic (recommended, see code)
%     flags.sake.a_init = either ['vbfa','vbfa1','magic_y','magic_y1'], meaning vbfa to y1, vbfa1 to y1, magic using y, magic using y1 (recommended); few lines down for y1 def.
%     flags.sake.xubar = either ['reg','eig'], meaning regular weight (recommended), or eigenspace corrected weight
%     flags.sake.updatelam = binary flag, 1 = yes, update lam in EM (recommended = 0)
%     flags.sake.plot = binary flag, 1= yes, plot some intermediate results, useful for debugging, not for batching
% flags.sake.nl = number of x-factors (optional, but recommended to be related to how many significant eigenvalues)
% flags.sake.nm = number of u interference factors (optional, but recommended at least 20 up to 100, depending on number of sensors)
% flags.sake.eigs = number of eigenvalues to choose if use eig-space bf for initialization (optional, but recommended to be related to how many significant eigenvalues)
%
% INTERMEDIATE stuff
% nf = number of lead field components
% s(nf,nt) = source factors
% nk = number of sensors
% nv = number of voxels
%
% number of em iterations:
% nem0 = for vbfa used to compute ryy
% nem1 = for vbfa used to initialize em at the given voxel
% nem  = for em at the given voxel
%
% OUTPUTS:
% sakeout.lmap0(nv,1) = initial likelihood map
% sakeout.lmap(nv,1) = final likelihood map  (<<-- MOST IMPORTANT OUTPUT for localization)
% sakeout.pmap0(nf,nf,nv) = initial power map (source-component covariance matrix)
% sakeout.pmap(nf,nf,nv) = final power map
% sakeout.pmapt(nv,1) = trace of source cov power per voxel
% sakeout.pmapt0(nv,1) = initial trace of source cov power per voxel
% sakeout.sbar_out(nf,nt,nv) = final source estimates
% sakeout.sbar0(nv,nt,nf) = initial source estimates
% sakeout.W1(nk,nf,nv) = weight to apply to data to get source estimate
% sakeout.lratmap = attempt at likelihood ratio (doesn't work)


% global bolts ndefaults

flags=flags.LFerror;

t0=dsearchn(data.latency,0);

yall=data.y'*1e15;
yall=yall-mean(yall(:,1:t0-1),2)*ones(1,size(yall,2));
ypre=yall(:,1:t0-1);
y=yall(:,t0:flags.timept(2));


f=permute(Lp,[1 3 2]);
clear Lp;

nk=size(y,1);
nt=size(y,2);
nv=size(f,2);
nf=size(f,3);

% if ~isfield(sakeflags,'nm')
%     sakeflags.nm=ceil(nk/3);    
% end
% if ~isfield(sakeflags,'eigs')
%     sakeflags.eigs=max(bolts.params.signalspace);
% end
% if ~isfield(sakeflags,'nl')
%     sakeflags.nl=sakeflags.eigs;  % used to be eigs + 1
% end

% nem0=40;nem1=10;
% % nem=20;
nem=flags.nem;
nm_pre=flags.nm_pre;
% nl=sakeflags.nl;
nl=flags.nl;

like=zeros(nem,1);
tphi=zeros(nem,1);
out.sbar0=zeros(nv,nt,nf);
out.sbar_out=zeros(nf,nt,nv);
out.pmap=zeros(nf,nf,nv);
out.pmap0=zeros(nf,nf,nv);
out.lmap0=zeros(nv,1);
out.lmap=zeros(nv,1);
out.lratmap=zeros(nv,1);
out.pmapt0=zeros(nv,1);
out.pmapt=zeros(nv,1);
out.W1=zeros(nk,nf,nv);
% 
a_init=0;lam_init=0;b_init=0;
switch flags.sbar_init  %test this right away before start computing
    case 'bf_nut'
        if ischar(flags.thisrun) %beamname
            load(flags.thisrun);
        else %hack for simulations
            beam=run_bf_onsim2(flags.thisrun);
        end
        sbar=zeros(nf,nt);
        beam.ts{1}=beam.s_th(:,end-nt+1:end);
        beam.ts{2}=beam.s_ph(:,end-nt+1:end);
        if nf==3
            beam.ts{3}=beam.s_z(:,end-nt+1:end);
        end
end
[b,lam,bet,ubar,psibb]=nut_vbfa(ypre,nm_pre,nem,b_init,lam_init,0);
clear bet ubar psibb

% switch sakeflags.iryy_init
%     case 'vbfa'
%         [a1,lam1,alp1,xbar1,psiaa1]=nut_vbfa(y,nl,nem0,a_init,lam_init,0);
%         ryy=a1*a1'+inv(lam1);
%         iryy=inv(ryy);
%         clear a1 lam1 alp1 psiaa1
%     case {'sefa_cy','sefa_aa'}
%         nem0=30;nem=30;nem1=30;
%         %      [a,b,lam,nu,alp,bet,xbar,ybar,ubar,cy,psi]=sefa_oct2005_3_kensuke(ypre,y,nl,nm,nem0,nem,nem1);
%         [a1,nu1,alp1,xbar1,ubar1,cy1,like1,www1,psiaa1]=nut_vbfa1(y,nl,b,lam,nem1,a_init,infer_nu,plotflag);
%         %      [a,b,lam,nu,alp,bet,xbar,ybar,ubar1,cy,like,www,psiaa,psibb]=nut_sefa(ypre,y,nl,nm,nem,nem1,0);
%         %    std_noise=sqrt(1 ./diag(lam));
%         %    y=zeros(size(ybar));
%         %    for ink=1:nk
%         %        sen_noise=std_noise(ink)*randn(1,nt);
%         %        y(ink,:)=ybar(ink,:)+sen_noise;
%         %    end
%         y_ori=y;
%         %   y=y-b*ubar;
%         switch sakeflags.iryy_init
%             case 'sefa_aa'
%                 ryy=a1*a1'+inv(lam);
%             case 'sefa_cy'
%                 ryy=cy1;
%         end
%         iryy=inv(ryy);
%         nem0=40;nem1=10;nem=10;
%         clear a1 nu1 alp1 ubar1 www1 psiaa1
% end
ryy=y*y';
[U,S]=eig(ryy);
eig_els=abs(diag(S));
[eig_el, iorder]=sort(-eig_els); %, eig_el=-eig_el;
U(:,:)=U(:,iorder);
ES = U(:,1:flags.nl);  % used to be .eigs when nl = eigs +1;
ESES=ES*ES';

% %% this gets called here in Bu version, assuming iryy_init='vbfa'
% [a,b,lam,nu,alp,bet,xbar,ybar,ubar1,cy,like,www,psiaa,psibb]=nut_sefa(ypre,y,nl,nm,nem,nem1,0);
% clear a nu alp bet xbar ybar ubar1 www psiaa psibb

% switch sakeflags.sbar_init  % needs some xbar1 estimate above
%     case 'algo1'
%         ryy=y*y';
%         ryx=y*xbar1';
%         rxx=xbar1*xbar1';
%         irxx=nut_pinv(rxx);
% 
%         sig0=inv(lam);
%         nu=1;
%         mg1=inv((ryy+sig0-ryx*irxx*ryx')/(nt+nu));
%         mg2=mg1*ryx*irxx;
% end
switch flags.phi_init
    case 'magic_y'
        ryy=y*y';

        sig0=b*b'+diag(1./diag(lam));
        [p0 d0]=svd(sig0);
        d0=diag(d0);
        s=p0*diag(sqrt(d0))*p0';
        invs=p0*diag(1./sqrt(d0))*p0';

        %         invs=diag(sqrt(diag(lam)));
        [p d]=svd(invs*ryy*invs/nt);
        d=diag(d);
%         a=inv(sqrtm(lam))*p(:,nf+1:nf+nl)*diag(sqrt(max(d(nf+1:nf+nl)-1,0)));
        phi_premult=inv(diag(sqrt(max(d(1:nf)-1,0))))*p(:,1:nf)'*sqrtm(lam);
end


%%%%%%%%%%%%%%% begin looping through voxels  %%%%%%%%%%%%%%%%

% tic
for iv=1:nv
    if(~mod(iv,20))  % update progress every 20th voxel
%         if(ndefaults.bf.viewwaitbar)
%             waitbar(iv/nv,h);
%         end
        disp(['Please wait. nut_BF_LFerror has finished ' num2str(iv) ' out of ' num2str(nv) ' voxels.']);
    end

%     iv,toc
    %     iv
    f1=squeeze(f(:,iv,:));

    % initialize
    switch flags.sbar_init
        case 'bf_reg'
            phi=f1'*iryy*f1;
            w=nut_pinv(phi)*f1'*iryy;
            sbar=w*y;
        case 'bf_nut'
            for ii=1:nf
                sbar(ii,:)=beam.ts{ii}(iv,:);
            end
            %             sourcecov=sbar*sbar';
            %             sourcesvd=svd(sourcecov);
            %             if ~sourcesvd(end)
            %                 phi=sourcesvd;
            %             end
        case 'bf_ES'
            phi=f1'*iryy*f1;
            w=nut_pinv(phi)*f1'*iryy;
            sbar=w*ESES*y;
        case 'algo1'
            g=nut_pinv(f1'*mg1*f1)*f1'*mg2;
            phi=nut_pinv(g*rxx*g'/nt);
            sbar=g*xbar1;
        otherwise
            warning('you didn''t give a proper sbar_init');
    end
    out.sbar0(iv,:,:)=sbar';
    switch flags.phi_init
        case 'ssT'
            phi=sbar*sbar';
        case 'magic_y'
            phi=(phi_premult*f1)^2;
        case 'fmri'
            phi=betamap;  %define from fmri
    end
%     switch sakeflags.a_init
%         case 'vbfa'
%             y1=y-f1*sbar;
%             [a,lam0,alp,xbar,psiaa]=nut_vbfa(y1,nl,nem1,a_init,lam_init,0);
%             clear lam0 xbar
%         case 'vbfa1'
%             y1=y_ori-f1*sbar;
%             infer_nu=0;
%             psiaa=[];
%             [a,nu,alp,xbar0,ubar0,cy,like0,www0,psiaa]=nut_vbfa1(y1,nl,b,lam,nem,a_init,infer_nu,0);
%             clear nu xbar0 ubar0 cy like0 www0
%         case 'magic_y1'  %use with fixed non-voxel-specific lam (test if from vbfa or sefa?)
%             y1=y-f1*sbar;
%             ryy1=y1*y1';
%             sig0=b*b'+diag(1./diag(lam));
%             [p0 d0]=svd(sig0);
%             d0=diag(d0);
%             s=p0*diag(sqrt(d0))*p0';
%             invs=p0*diag(1./sqrt(d0))*p0';
% 
%             %             invs=diag(sqrt(diag(lam)));
%             [p d]=svd(invs*ryy1*invs/nt);
%             d=diag(d);
%             a=inv(sqrtm(lam))*p(:,1:nl)*diag(sqrt(max(d(1:nl)-1,0)));
%             %             % only to get bet and psi for now
%             %             [tmp1,tmp2,bet,tmp3,psi]=nut_vbfa(y1,nl,nem1,a_init,lam_init,0);
%             %             clear tmp1 tmp2 tmp3
%             alp=diag(1./diag(a'*lam*a/nk));
%             alp=min(diag(alp))*diag(ones(nl,1));
%             psiaa=eye(nl)/nt;
%         case 'magic_y'
%             %             y1=y-f1*sbar;
%             %             % only to get bet and psi for now
%             %             [tmp1,tmp2,bet,tmp3,psi]=nut_vbfa(y1,nl,nem1,a_init,lam_init,0);
%             %             clear tmp1 tmp2 tmp3
%             alp=diag(1./diag(a'*lam*a/nk));
%             alp=min(diag(alp))*diag(ones(nl,1));
%             psiaa=eye(nl)/nt;
%     end

psi=eye(nf)/nt;
alp=eye(nf)/nk;

    out.pmap0(:,:,iv)=nut_pinv(phi);

    ryy=y*y';
%     nu=ones(nm,1);

    % em iteration

    fg=f1;

%     fab=[fg a b];
%     alapsi=fab'*lam*fab;
    alapsi=fg'*lam*fg + nk*psi + phi;
%     alapsi(nf+1:nf+nl,nf+1:nf+nl)=alapsi(nf+1:nf+nl,nf+1:nf+nl)+nk*psiaa;
    ldlam=sum(log(diag(lam/(2*pi))));

    for iem=1:nem
        gam=alapsi;
%         gam(1:nf,1:nf)=gam(1:nf,1:nf)+phi;
%         gam(nf+1:end,nf+1:end)=gam(nf+1:end,nf+1:end)+diag([ones(nl,1);nu]);
        igam=nut_pinv(gam);
        W=igam*fg'*lam;
        switch flags.sbar_up
            case 'reg'
                sbar=W*y;
            case 'eig'
                sbar=W*ESES*y;
                %             case 'ybar'
                %                 sxubar=W*ybar;
        end
        if flags.updatelam
            ldlam=sum(log(diag(lam/(2*pi))));
        end
        ldphi=sum(log(svd(phi)));
        ldgam=sum(log(svd(gam)));
        ldalp=sum(log(diag(alp)));
        ldpsi=sum(log(svd(psi)));

%         sbar=sxubar(1:nf,:);
%         xbar=sxubar(nf+1:nf+nl,:);
%         ubar=sxubar(nf+nl+1:end,:);
%         xubar=sxubar(nf+1:end,:);
%         gamxu=gam(nf+1:end,nf+1:end);
        rys=y*sbar';
%         ryx=y*xbar';
%         ryu=y*ubar';
        rss=sbar*sbar'+nt*igam;
%         rxx=xbar*xbar'+nt*igam(nf+1:nf+nl,nf+1:nf+nl);
%         ruu=ubar*ubar'+nt*igam(nf+nl+1:nf+nl+nm,nf+nl+1:nf+nl+nm);
%         rsx=sbar*xbar'+nt*igam(1:nf,nf+1:nf+nl);
%         rux=ubar*xbar'+nt*igam(nf+nl+1:end,nf+1:nf+nl);
%         rsu=sbar*ubar'+nt*igam(1:nf,nf+nl+1:end);
        if (flags.plot | iem==nem | iem==1)
            like(iem)=.5*nt*(ldlam+ldphi-ldgam)-.5*sum(sum(y.*(lam*y)))+.5*sum(sum(sbar.*(gam*sbar)))+.5*nk*(ldalp+ldpsi);
            %             noise(iem)=.5*nt*log(det([a*rxx*a'+ b*ruu*b' + lam]));
            noise(iem)=.5*nt*(ldlam-ldgam)-.5*sum(sum(y.*(lam*y)))+.5*nk*(ldalp+ldpsi);
            likerat(iem)=like(iem)-noise(iem);
        end

        tphi(iem)=trace(nut_pinv(phi));
        if flags.plot
            figure(5);
            subplot(3,3,3);plot((1:iem)',tphi(1:iem));title(iv);
            subplot(3,3,6);plot((1:iem)',like(1:iem)/nt);
            subplot(3,3,7);plot((1:iem)',noise(1:iem)/nt);
            drawnow;
        end


        psi=nut_pinv(rss+alp);
        phi=nut_pinv(rss/nt);
        fg=rys*psi;
%         a=(ryx-fg*rsx-b*rux)*psiaa;

        if flags.updatelam
            %             % next two lines for saketini (not Bu)
            %             ryy1=ryy-fg*rys'-rys*fg'+fg*rss*fg';
            %             lam=diag(nt./diag(ryy1-a*(ryx-fg*rsx)'));
            lam=diag(nt./diag(ryy-2*fg*rys'+fg*rss*fg'));
        end

        %         aa=[fg a];
        %         alapsi=aa'*lam*aa;
        %         alapsi(nf+1:end,nf+1:end)=alapsi(nf+1:end,nf+1:end)+nk*psi;
        %         alp=diag(nk./diag(alapsi(nf+1:end,nf+1:end)));
%         fab=[fg a b];
%         alapsi=fab'*lam*fab;
%         alapsi(nf+1:nf+nl,nf+1:nf+nl)=alapsi(nf+1:nf+nl,nf+1:nf+nl)+nk*psiaa;
%         alp=diag(nk./diag(alapsi(nf+1:nf+nl,nf+1:nf+nl)));
        alp=diag(nk./diag( fg'*lam*fg + rss + alp ));
    alapsi=fg'*lam*fg + nk*psi + phi;
%     alapsi(nf+1:nf+nl,nf+1:nf+nl)=alapsi(nf+1:nf+nl,nf+1:nf+nl)+nk*psiaa;
    ldlam=sum(log(diag(lam/(2*pi))));

%         if sakeflags.infer_nu  % leaving this off for now.
%             nu=nt./diag(ruu);
%             if sakeflags.plot
%                 subplot(3,3,8);plot(1./nu);
%             end
%         end
    end
    out.W1(:,:,iv)=W';
    %     Wout(:,:,iv)=W(1:nf,:);
    out.pmap(:,:,iv)=nut_pinv(phi);
    out.sbar_out(:,:,iv)=sbar;

    out.lmap0(iv)=like(1)/nt;
    out.lmap(iv)=like(iem)/nt;
%     out.bic(iv)=-2*out.lmap(iv)+nl*log(nt);
%     out.aic(iv)=-2*out.lmap(iv)+nl*2;
    out.lratmap(iv)=likerat(iem)/nt;
    out.pmapt0(iv)=trace(out.pmap0(:,:,iv));
    out.pmapt(iv)=trace(out.pmap(:,:,iv));
    out.Fg(:,:,iv)=fg;
    [uu,dd,vv]=svd(psi);
    out.lampsi(:,:,iv)=diag(lam)*diag(dd)';
    if flags.plot
        figure(5);
        subplot(3,3,1);plot((1:iv)',out.pmapt0(1:iv));title('initial');ylabel('pmap')
        subplot(3,3,2);plot((1:iv)',out.pmapt(1:iv));title('final');
        subplot(3,3,4);plot((1:iv)',out.lmap0(1:iv));ylabel('lmap');xlabel('voxel');
        subplot(3,3,5);plot((1:iv)',out.lmap(1:iv));xlabel('voxel');
        subplot(3,3,9);plot((1:iv)',out.lratmap(1:iv));xlabel('voxel');
    end
end
% out.W1=permute(Wout,[2 1 3]);
out.flags=flags;
pack;
return




