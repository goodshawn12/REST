function rr=nut_car2pol(xx)
% function rr=nut_car2pol(xx)
%  conversion from cartesian to polar
%  xx=[x, y, z]        : cartesian coordinate
%  rr=[r, theta, fai]  : polar coordinate
%
rr=zeros(1,3);
r=norm(xx);
if(r==0), 
   rr=[0 0 0];
else
  rrxy=sqrt(xx(1)^2 + xx(2)^2);
  if(rrxy~=0), theta=pi/2-atan(xx(3)/rrxy);,end
  if(rrxy==0)
     if(xx(3)>0)
         theta=0;
     else
         theta=pi;
     end
  end
  fai=atan2(xx(2),xx(1));
  if(xx(2)<0), fai=2*pi + fai;, end
   rr=[r theta fai];
end

