clear all; figure;
%**************************************************************************
% 3D Interpolation (takes some time)
%**************************************************************************
%Matlab
[x,y,z,v] = flow(10); 
[xi,yi,zi] = meshgrid(.1:.25:10, -3:.25:3, -3:.25:3);
vi = interp3(x,y,z,v,xi,yi,zi); 
subplot(2,2,1); slice(xi,yi,zi,vi,[6 9.5],2,[-2 .2]), shading flat; ; title('Interpolation using Matlab function interp3');
subplot(2,2,2); slice(x,y,z,v,[6 9.5],2,[-2 .2]), shading flat; title('data used for interpolation');

%RBF
op=rbfcreate([x(:)'; y(:)'; z(:)'], v(:)','RBFFunction', 'multiquadric', 'Stats', 'on');
rbfcheck(op);
rbfvi = rbfinterp([xi(:)'; yi(:)'; zi(:)'], op);
rbfvi = reshape(rbfvi, size(xi));

subplot(2,2,3); slice(xi,yi,zi,rbfvi,[6 9.5],2,[-2 .2]), shading flat; ; title('RBF Interpolation');

[xr,yr,zr,vr] = flow(100); 
subplot(2,2,4); slice(xr,yr,zr,vr,[6 9.5],2,[-2 .2]), shading flat; title('Original 3D data');

clear all; figure;
%**************************************************************************
% 3D Interpolation using isosurface vizualization
%**************************************************************************
%Matlab
[x,y,z,v] = flow(10); 
[xi,yi,zi] = meshgrid(.1:.55:10, -3:.55:3, -3:.55:3);

vi = interp3(x,y,z,v,xi,yi,zi); 
subplot(2,2,1); p = patch(isosurface(xi,yi,zi,vi,-3));  title('Interpolation using Matlab function interp3');
isonormals(xi,yi,zi,vi,p)
set(p,'FaceColor','red','EdgeColor','none');
daspect([1 1 1])
view(3); axis tight
camlight 
lighting gouraud
% show data used for interpolation
subplot(2,2,2); p = patch(isosurface(x,y,z,v,-3));  title('data used for interpolation');
isonormals(x,y,z,v,p)
set(p,'FaceColor','red','EdgeColor','none');
daspect([1 1 1])
view(3); axis tight
camlight 
lighting gouraud

%RBF
op=rbfcreate([x(:)'; y(:)'; z(:)'], v(:)','RBFFunction', 'multiquadric', 'Stats', 'on');
rbfcheck(op);
rbfvi = rbfinterp([xi(:)'; yi(:)'; zi(:)'], op);
rbfvi = reshape(rbfvi, size(xi));

subplot(2,2,3); p = patch(isosurface(xi,yi,zi,rbfvi,-3)); title('RBF Interpolation');
isonormals(xi,yi,zi,rbfvi,p)
set(p,'FaceColor','red','EdgeColor','none');
daspect([1 1 1])
view(3); axis tight
camlight 
lighting gouraud

[x,y,z,v] = flow(100);
subplot(2,2,4); p = patch(isosurface(x,y,z,v,-3)); title('Original 3D data');
isonormals(x,y,z,v,p)
set(p,'FaceColor','red','EdgeColor','none');
daspect([1 1 1])
view(3); axis tight
camlight 
lighting gouraud
