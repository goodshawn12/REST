function champ_out = nut_Champagne(Lp,data,flags)

% function to run Champagne:  Y_post = Fs + Bu + v
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%inputs:
%Lp - leadfield matrix: (sensorsxdirxvoxels)
% data - uses data.y for sensor data: (timexsensors)
%flags: uses the control and active time markers
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%output
%will output hyperparameter and the source time-courses 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

col_norm=flags.cn;
AX=flags.champ.ax;
multf = flags.champ.multf;

if(AX==1||AX==3)
    display(['using ' num2str(multf) ' as scalar in noise estimate'])
end

%find num of time, sensors, and voxels
nt=size(data.y(flags.tb.timeptc(1):flags.tb.timepta(2)),1);
ns=size(data.y,2);
nv=size(Lp,3);

if flags.champ.nem_ch<50
    warning('That is not enough iterations.  We recommend more than 50.')
end

if size(data.y,3)>1
    display('Champagne only runs on averaged data')
    return
end

%column normalize lead field
[Lp, foo] = norm_leadf_col(Lp);

%make data on scale 1-10
m=max(max(max(abs(data.y))));

data.y=data.y*(1/m);
data.y=double(data.y);

%extract pre- and post-stim data
prestim=data.y(flags.tb.timeptc(1):flags.tb.timeptc(2),:)';
prestim=prestim-mean(prestim,2)*ones(1,size(prestim,2));
poststim=data.y(flags.tb.timepta(1):flags.tb.timepta(2),:)';
poststim=poststim-mean(poststim,2)*ones(1,size(poststim,2));

timepts=data.latency(flags.tb.timepta(1):flags.tb.timepta(2));

%initialize ialp (hyper parameter matrix)
    for i=1:nv
        for j=1:size(Lp,2)
            for k=1:size(Lp,2)
                if j==k
                    ialp(j,k,i)=1;
                else
                    ialp(j,k,i)=0;
                end
            end
        end
    end
    ialp0=zeros(size(Lp,2),size(Lp,2),nv);
    %puts random # in iapl0 with off diags the same
    for iv=1:nv
        a0=randn(size(Lp,2),size(Lp,2));ialp0(:,:,iv)=a0*a0';
    end
    ialp1=max(max(max(ialp)));
    ialp_init=ialp+ialp0*ialp1/10000; %moves values sligtly off

%leadfield should be [num dir (2) X ns X nv] lead-field for two components
Lp = permute(Lp,[2 1 3]);

display('computing b and lam from pre-stim data using SEFA')
nl=20; %num of factors
nem_sefa=100;
[a1,b,lam,alp,bet,xbar1]=sefa0(poststim,prestim,nl,25,nem_sefa,0);

% could use vbfa
%b_init=0;
%lam_init=0;
%[b,lam,bet,ubar]=nut_vbfa(prestim,nl,nem_vbfa,b_init,lam_init,0);

%computer noise covariance
if(AX==1)
    Sigma_e =diag(multf*ones(ns,1)); %3 works most of the time
    poststim=a1*xbar1;
elseif(AX==0)
    multf=0;
    Sigma_e = b*b' + inv(lam);
elseif(AX==3)
    Sigma_e =diag(multf*ones(ns,1));
elseif(AX==2)
    Sigma_e =inv(lam);
else
    error('you somehow didn''t choose a valid AX')
end

%number of em iterations
nem_ch=flags.champ.nem_ch;

%interleave columns for champagne_plain
F=zeros(size(Lp,2), size(Lp,1)*size(Lp,3));
for i=1:size(Lp,1)
    F(:,i:size(Lp,1):size(F,2))=permute(Lp(i,:,:), [2 3 1]);
end
used_plain=false;
display('running Champagne')
if size(Lp,1)==2
    if nv<7000
        [ialp_c3,s_c3,W1,cost,k,d_ialp] = champagne_2x2(poststim,Lp,Sigma_e,nem_ch,ialp_init);
    else
        used_plain=true;
        [ialp_c3,s_c3,W1,cost,k,d_ialp] = champagne_plain(poststim,F,Sigma_e,nem_ch,ialp_init,2);
    end
elseif size(Lp,1)==3
    used_plain=true;
    [ialp_c3,s_c3,W1,cost,k,d_ialp] = champagne_plain(poststim,F,Sigma_e,nem_ch,ialp_init,3);
elseif size(Lp,1)==1
    used_plain=true;
    [ialp_c3,s_c3,W1,cost,k,d_ialp] = champagne_plain(poststim,F,Sigma_e,nem_ch,ialp_init,1);
end

if used_plain
    %de-interleave columns from champagne_plain
    s2=zeros(size(Lp,1),size(Lp,3),length(timepts));
    for i=1:size(Lp,1)
        s2(i,:,:)=permute(s_c3(i:size(Lp,1):size(Lp,1)*size(Lp,3),:), [3 1 2]);
    end
    s_c3=s2;
    clear s2;
end


    
%output the source time courses and hyperparameters
champ_out.sources=s_c3;
champ_out.hyper=ialp_c3;

%output the post-stim time points
champ_out.timepts=timepts;

%separate out the time courses in each direction
%champ_out.s_th=squeeze(s_c3(1,:,:));
%champ_out.s_ph=squeeze(s_c3(2,:,:));
%champ_out.s=s_c3;

%get one hyperparameter value per voxel by summing
temp=ialp_c3(1,1,:);
for i=2:size(Lp,1)
    temp=temp+ialp_c3(i,i,:);
end
champ_out.hyper1=squeeze(temp);
clear temp;

%save the weights from nut_Champagne_2x2
champ_out.W1=W1;
champ_out.mn=m;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [L,colnorm] = norm_leadf_col(L)

L = permute(L,[1 3 2]);

if size(L,3)>1
    for i=1:size(L,3)
        colnorm(:,i) = sqrt(sum(L(:,:,i).^2));
        L(:,:,i) = L(:,:,i)./repmat(colnorm(:,i)',[size(L,1) 1]);
    end
else
    for i=1:size(L,2)
        colnorm(:,i) = sqrt(sum(L(:,i).^2));
        L(:,i) = L(:,i)./repmat(colnorm(:,i)',[size(L,1) 1]);
    end
end

L = permute(L,[1 3 2]);

%%
function [a,b,lam,alp,bet,xbar]=sefa0(y,y0,nl,n_inf,nem_init,ifplot);

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
% nem_init = number of em iterations

nk=size(y,1);
nt=size(y,2);
nt0=size(y0,2);
			    


disp('VBFA initialization of SEFA0');
b_init=0;
lam_init=0;
[b,lam,bet,ubar]=vbfa(y0,n_inf,nem_init,b_init,lam_init,ifplot);



a_init=0;
ryy=y*y';
ryy0=y0*y0';
if a_init==0 
%   [p d q]=svd(ryy/nt);d=diag(d);
%   a=p*diag(sqrt(d));
%   a=a(:,1:nl);
   sig0=b*b'+diag(1./diag(lam));
   [p0 d0 q0]=svd(sig0);d0=diag(d0);
   s=p0*diag(sqrt(d0))*p0';
   invs=p0*diag(1./sqrt(d0))*p0';
   [p d q]=svd(invs*ryy*invs/nt);d=diag(d);
%   a=s*p(:,1:nl)*diag(sqrt(max(d(1:nl)-1,0)));
   a=s*p(:,1:nl)*diag(sqrt(abs(d(1:nl)-1)));
else
   a=a_init;
end



% initialize by svd
if 1>2
ryy=y*y';
[p d q]=svd(ryy/nt);d=diag(d);
a=p*diag(sqrt(d));
a=a(:,1:nl);

ryy0=y0*y0';
[p d q]=svd(ryy0/nt0);d=diag(d);
b=p*diag(sqrt(d));
b=b(:,1:n_inf);
lam=diag(nt0./diag(ryy0));
end



alp=diag(1./diag(a'*lam*a/nk));
alp=min(diag(alp))*diag(ones(nl,1));
bet=diag(1./diag(b'*lam*b/nk));
bet=min(diag(bet))*diag(ones(n_inf,1));

ab=[a b];
alpbet=diag([diag(alp);diag(bet)]);
nlm=nl+n_inf;
psi=eye(nlm)/(nt+nt0);

% em iteration

like=zeros(nem_init,1);
alapsi=ab'*lam*ab+nk*psi;

for iem=1:nem_init
   gam=alapsi+eye(nlm);
   igam=inv(gam);
   xubar=igam*ab'*lam*y;

   b=ab(:,nl+1:nlm);
   psib=psi(nl+1:nlm,nl+1:nlm);
   gam0=b'*lam*b+nk*psib+eye(n_inf);
   igam0=inv(gam0);
   ubar0=igam0*b'*lam*y0;

   ldlam=sum(log(diag(lam/(2*pi))));
   ldgam=sum(log(svd(gam)));
   ldgam0=sum(log(svd(gam0)));
   ldalpbet=sum(log(diag(alpbet)));
   ldpsi=sum(log(svd(psi)));
   like0=.5*nt0*(ldlam-ldgam0)-.5*sum(sum(y0.*(lam*y0)))+.5*sum(sum(ubar0.*(gam0*ubar0)));
   like(iem)=.5*nt*(ldlam-ldgam)-.5*sum(sum(y.*(lam*y)))+.5*sum(sum(xubar.*(gam*xubar)))+.5*nk*(ldalpbet+ldpsi)+like0;
if(ifplot)
   subplot(3,3,1);plot((1:iem)',like(1:iem));title('SEFA0');
   subplot(3,3,4);plot([mean(ab.^2,1)' 1./diag(alpbet)]);
   subplot(3,3,7);plot(1./diag(lam));
   drawnow;
end
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

%%
function [a,lam,alp,xbar]=vbfa(y,nl,nem_init,a_init,lam_init,ifplot);

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
% nem_init = number of em iterations

nk=size(y,1);
nt=size(y,2);

% initialize by svd

ryy=y*y';
if a_init==0
   [p d q]=svd(ryy/nt);d=diag(d);
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

like=zeros(nem_init,1);
alapsi=a'*lam*a+nk*psi;


for iem=1:nem_init
   gam=alapsi+eye(nl);
   igam=inv(gam);
   xbar=igam*a'*lam*y;


if(ifplot)
   subplot(3,3,1);plot((1:iem)',like(1:iem)/nt);title('VBFA: like');
   subplot(3,3,4);plot([mean(a.^2,1)' 1./diag(alp)]);title('1/alp');
   subplot(3,3,7);plot(1./diag(lam));title('1/lam');
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
%%
