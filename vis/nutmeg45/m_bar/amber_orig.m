
function [yhat,shat]=amber(ypre,ypost,f,nd,nl,nm);


% Reconstruct voxel activity from multisensor stimulus evoked data
% by probability-based time-frequecy modeling
%
% --- Version 1/3/11 -----
% 
% © 2010-11 Convex Imaging
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Output: 
% yhat(nk,nt) = time-domain evoked voxel contribution to sensor data 
% shat(nv*nd,nt) = time-domain evoked voxels 
%
% Input:
% ypre(nk,nt0) = time-domain prestim sensor data 
% yost(nk,nt) = time-domain poststim sensor data 
% f(nk,nv*nd) = lead field matrix
% nd = dipole vector dimensionality
% nl = # of evoked factors
% nm = # of interference factors
%
% #s:
% nk = # of sensors
% nv = # of voxels
% nt0 = # of prestim time points
% nt = # of poststim time points
% nem = # of EM iterations
% nw = window length
% nf = FFT length
% eta = factor spectral smoothness


nem=30;
nw=12;
nf=16;
eta=.95;

%nem1=nem;
%nem2=nem;
%nem3=min(nem/2,15);

[nk nt]=size(ypost);
nt0=size(ypre,2);

ifplot=1;
[yhat xhat a b lam]=prepoststim(ypost,ypre,nl,nm,nem,ifplot,nw,nf,eta);

cb=b*b'+diag(1./diag(lam));
nem1=nem;
[gamma shat]=champ(ypost,f,cb,nem1,nd);
yhat=f*shat;


figure;
t=linspace(-nt0/nt,1,nt0+nt)';
ym=max(max(abs([ypre ypost])));

subplot(2,2,1);plot(t,[ypre ypost]');
title('Sensor data','fontsize',16);
xlabel('time','fontsize',12);ylabel('amplitude','fontsize',12);
axis([min(t) max(t) -ym ym]);

subplot(2,2,2);plot(t,[zeros(nk,nt) yhat]');
title('Voxel contribution','fontsize',16);
xlabel('time','fontsize',12);ylabel('amplitude','fontsize',12);
axis([min(t) max(t) -ym ym]);

nv=floor(size(f,2)/nd);
trg=zeros(nv,1);
for id=1:nd
    trg=trg+squeeze(gamma(id,id,:));
end
subplot(2,2,3);plot(trg);
title('Voxel power','fontsize',16);
xlabel('voxel index','fontsize',12);ylabel('power','fontsize',12);


return

