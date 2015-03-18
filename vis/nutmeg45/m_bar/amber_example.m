

% An example of using Amber
%
% --- Version 1/3/11 -----
% 
% © 2010-11 Convex Imaging
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Please read the comments in this file as well as in amber.m


% Simulate data using the underlying model

nk=30; % # of sensors
nv=100; % # of voxels
nm=10; % # of interference factors
nt0=900; % # of prestim time points
nt=900; % # of poststim time points

snr=-8; % signal-to-noise ratio

f=randn(nk,nv); % lead field matrix
b=randn(nk,nm); % interference mixing matrix
lam=ones(nk,1)*10; % noise precision
u0=randn(nm,nt0); % prestim interference factors
u=randn(nm,nt); % poststim interference factors
v0=randn(nk,nt0)./(sqrt(lam)*ones(1,nt0)); % prestim noise
v=randn(nk,nt)./(sqrt(lam)*ones(1,nt)); % poststim noise


% Simulate time course of active voxels

t=linspace(0,1,nt)';
g=exp(-25*((t-.5).^2));
sactive=[cos(2*pi*4*t) cos(2*pi*5*t+pi/5) cos(2*pi*6*t+pi/3)];
sactive=(ones(3,1)*g').*sactive';
s=zeros(nv,nt);
s([10 41 83],:)=sactive; % voxel activity


% Scale according to SNR

yposts=f*s; % sensor data from voxels
ypostn=b*u+v; % sensor data from noise and interference

scl=sqrt(mean(mean(yposts.^2))*(10^(-snr/10))/mean(mean(ypostn.^2)));
u0=scl*u0;
v0=scl*v0;
u=scl*u;
v=scl*v;

ypre=b*u0+v0; % prestim data
ypost=f*s+b*u+v; % poststim data


% Call Amber

nd=1;
nl=3;
nm=10;

[yhat,shat]=amber(ypre,ypost,f,nd,nl,nm);


% Compute sensor SNR gain

es=mean(mean((f*s).^2));
en=mean(mean((b*u+v).^2));
ey=mean(mean((ypost-f*s).^2));
eyhat=mean(mean((yhat-f*s).^2));
%disp(['snr = ',num2str(10*log10(es/en))]);
disp(['Sensor SNR in = ',num2str(10*log10(es/ey))]);
disp(['Sensor SNR out = ',num2str(10*log10(es/eyhat))]);

% compute voxel reconstruction SNR

ev=mean(mean(s.^2))/mean(mean((s-shat).^2));
disp(['Voxel reconstruction SNR = ',num2str(10*log10(ev))]);


% Regarding nd = dipole vector dimensionality
%
% The above example uses nd=1. For nd>=1, the lead field matrix f has 
% nv*nd rows and columns. Amber assumes that the rows are ordered by 
% voxel index - component index. For, e.g., nd=3, the row order is:
%
% row 1: voxel 1, component 1
% row 2: voxel 1, component 2
% row 3: voxel 1, component 3
%     .....
% row nv*3-2: voxel nv, component 1
% row nv*3-1: voxel nv, component 2
% row nv*3: voxel nv, component 3


% Regarding nm = interference factor dimensionality
%
% The pre-stimulus activity is modeled using a factor based approach,
% y = B * u + v
% where u = interference factors, B = interference mixing matrix, v= noise.
% The maximum number of factors nm is input by the user. The optimal number
% of factors is chosen by the code. It is somewhat similar to the number of
% above-threshold principal components of the pre-stimulus data. Too large
% nm can slow down processing. To ensure nm is not too small, run the code 
% and look at the subplot titled 'Prestim eigenvalues'. A good choice of nm 
% is in the regime where the curve approaches zero. However, results are
% sometime sensitive to nm, and some experimentation may be required.


% Regarding nl = evoked factor dimensionality
%
% The post-stimulus activity -- at an intermediate processing step -- is also
% modeled using factors, including both interference and evoked factors. The 
% maximum number of evoked factors nl is input by the user. The optimal 
% number is chosen by the code, and is related to the number of evoked dipole 
% components. To ensure nl is not too small, run the code and
% look at the subplot titled 'Poststim eigenvalues'. A good choice of nl is
% in the regime where the curve approaches zero. However, results are
% sometime sensitive to nl, and some experimentation may be required.


% Regarding nt0 = length of pre-stimulus period
%
% To obtain an accurate estimate of the interenference and noise model, the
% pre-stimulus period should be sufficiently long. If it's too short, there
% are not enough data; if it's too long, non-stationarity might apply; in
% eithe case, the resulting estimate may be inaccurate. Best results are
% obtained when the pre-stimulus period length is similar to the post-
% stimulus period length. If it's much shorted, concatenate 2 or more
% copies of the pre-stimulus data.







