%**************************************************************************
% 1D Interpolation
%**************************************************************************

x = 0:1.25:10; 
y = sin(x); 
xi = 0:.1:10; 

%Matlab
yi = interp1(x,y,xi); 
subplot(2,1,1); plot(x,y,'o',xi,yi, xi, sin(xi),'r'); title('Interpolation using Matlab function interp1');

%RBF
%op=rbfcreate(x, y,'RBFFunction', 'thinplate'); rbfcheck(op);
%op=rbfcreate(x, y,'RBFFunction', 'linear'); rbfcheck(op);
%op=rbfcreate(x, y,'RBFFunction', 'cubic'); rbfcheck(op);
%op=rbfcreate(x, y,'RBFFunction', 'gaussian'); rbfcheck(op);
op=rbfcreate(x, y,'RBFFunction', 'multiquadric', 'RBFConstant', 2); rbfcheck(op);
op=rbfcreate(x, y); rbfcheck(op);
%op=rbfcreate(x, y,'RBFFunction', 'gaussian');
%op=rbfcreate(x, y);
fi = rbfinterp(xi, op);
subplot(2,1,2); plot(x, y,'o', xi, fi,xi, sin(xi),'r'); title('RBF interpolation');

yi=rbfinterp_ez(x,y,xi);
plot(x, y,'o', xi, fi,xi, sin(xi),'r'); title('RBF interpolation');

clear all; figure;
%**************************************************************************
% 2D Interpolation
%**************************************************************************
%Matlab standard interpolation using griddata 
rand('seed',0)
x = rand(50,1)*4-2; y = rand(50,1)*4-2;
z = x.*exp(-x.^2-y.^2);

ti = -2:.05:2; 
[XI,YI] = meshgrid(ti,ti);
ZI = griddata(x,y,z,XI,YI,'cubic');

subplot(2,2,1); mesh(XI,YI,ZI), hold, axis([-2 2 -2 2 -0.5 0.5]); 
plot3(x,y,z,'.r'), hold off; title('Interpolation using Matlab function griddata(method=cubic)');

subplot(2,2,3); pcolor(abs(ZI - XI.*exp(-XI.^2-YI.^2))); colorbar; title('griddata(method=cubic) interpolation error');

%RBF interpolation
%op=rbfcreate([X(:)'; Y(:)'], Z(:)','RBFFunction', 'thinplate'); rbfcheck(op);
%op=rbfcreate([X(:)'; Y(:)'], Z(:)','RBFFunction', 'linear'); rbfcheck(op);
%op=rbfcreate([X(:)'; Y(:)'], Z(:)','RBFFunction', 'cubic'); rbfcheck(op);
%op=rbfcreate([X(:)'; Y(:)'], Z(:)','RBFFunction', 'gaussian'); rbfcheck(op);
op=rbfcreate([x'; y'], z','RBFFunction', 'multiquadric', 'RBFConstant', 2);
rbfcheck(op);
%op=rbfcreate(x, y,'RBFFunction', 'gaussian');
%op=rbfcreate(x, y);
ZI = rbfinterp([XI(:)'; YI(:)'], op);
ZI = reshape(ZI, size(XI));
subplot(2,2,2); mesh(XI,YI,ZI), hold
plot3(x,y,z,'.r'), hold off; title('RBF interpolation'); axis([-2 2 -2 2 -0.5 0.5]);
subplot(2,2,4); pcolor(abs(ZI - XI.*exp(-XI.^2-YI.^2))); colorbar; title('RBF interpolation error');

clear all;
%**************************************************************************
% 1D smooth interpolation 
%**************************************************************************

x = 0:1.25:10; 
y = sin(x)+rand(size(x))-0.5; 
xi = 0:.1:10; 

%Matlab
yi = interp1(x,y,xi); 
subplot(2,1,1); plot(x,y,'o',xi,yi, xi, sin(xi),'r'); title('Interpolation using Matlab function interp1');

%RBF
op=rbfcreate(x, y,'RBFFunction', 'multiquadric', 'RBFConstant', 2, 'RBFSmooth', 0.1); rbfcheck(op);
fi = rbfinterp(xi, op);
subplot(2,1,2); plot(x, y,'o', xi, fi,xi, sin(xi),'r'); title('RBF interpolation');

clear all;
%**************************************************************************
% 1D  interpolation of random data
%**************************************************************************
N = 30;
smooth = 1e-2;
x = rand(1, N); 
y = rand(1, N); 
xi = [x rand(1, 10*N)]; 

%RBF
%op=rbfcreate(x, y,'RBFFunction', 'multiquadric', 'RBFConstant', 2, 'RBFSmooth', 1e-6); rbfcheck(op);
op=rbfcreate(x, y,'RBFFunction', 'cubic', 'RBFSmooth', smooth); 
yi = rbfinterp(xi, op);
[xi,index] = sort(xi);
subplot(2,2,1); plot(x, y,'.r', xi, yi(index)); title('RBF cubic interpolation'); 

op=rbfcreate(x, y,'RBFFunction', 'linear', 'RBFSmooth', smooth); 
yi = rbfinterp(xi, op);
[xi,index] = sort(xi);
subplot(2,2,2); plot(x, y,'.r', xi, yi(index)); title('RBF linear interpolation');

op=rbfcreate(x, y,'RBFFunction', 'multiquadric', 'RBFSmooth', smooth); 
yi = rbfinterp(xi, op);
[xi,index] = sort(xi);
subplot(2,2,3); plot(x, y,'.r', xi, yi(index)); title('RBF multiquadric interpolation');

op=rbfcreate(x, y,'RBFFunction', 'gaussian', 'RBFSmooth', smooth); 
yi = rbfinterp(xi, op);
[xi,index] = sort(xi);
subplot(2,2,4); plot(x, y,'.r', xi, yi(index)); title('RBF gaussian interpolation');
