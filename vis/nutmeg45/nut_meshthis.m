function [meshx,meshy,Z]=nut_meshthis(coords,Zin,gridsize,method);
% [meshx,meshy,Z]=nut_meshthis(coords,Zin,gridsize)
% [meshx,meshy,Z]=nut_meshthis(coords,Zin,gridsize,method)
% function previously named fu*kmatlab
%
% creates MATLAB's prized meshgrid and remaps your corresponding Z values.
% lets you use contour, surf, mesh without getting pissed off.
% e.g.:
% contour(meshx,meshy,Z);
%
% optional "method" is passed on to griddata (cubic, linear, or nearest)

if(~exist('method','var'))
    method=[];
end

[meshx,meshy] = meshgrid(min(coords(:,1)):gridsize:max(coords(:,1)),min(coords(:,2)):gridsize:max(coords(:,2)));
if(0)
    [k,d] = dsearchn(coords,[meshx(:) meshy(:)]);
    Z = Zin(k);
    % Z(d>.1) = NaN;
    Z=reshape(Z,size(meshx,1),size(meshx,2));
else
    Z = griddata(coords(:,1),coords(:,2),Zin(:),meshx,meshy,method);
end
