function [r,p]=nut_phs(x,y)
% [R,P]=NUT_PHS(X,Y)
%
% AK: Some way to convert cartesian coordinates into polar. I can't
% catch whether there are any differences, but there may be some.

r=sqrt(x.^2+y.^2);
in=find(r~=0);
p=0*r;
p(in)=acos(x(in)./r(in)).*nut_mysign(y(in));
return
