function [w,r0,err]=nut_cnv3(r1,r2,st)
%[ROTATION,TRANSLATION,ERROR]=NUT_CNV3(DATA1,DATA2,1)
% This function is used to calculate the transformation
% matrix. Used by nut_adjust_coregistry.

n=size(r1,1);
err=ones(n,1); 
Wp=spdiags(err/sum(err),0,n,n);
%find the centroid for the two incoming data
r10=sum(Wp*r1,1);
r20=sum(Wp*r2,1);
%normalize the two incoming data
p1=r1-nut_spread(r1,r10);
p2=(r2-nut_spread(r2,r20))*st;
%find the rotation matrix
[u,q,v]=svd(p2'*Wp*Wp*p1,0);
w=u*v';
%find error in rotation
err=sqrt(sum((p2*w-p1).^2,2));
err=1./(err+0.01*max(abs(err)));
w=st*w;
%find the translation
r0=-r20*w+r10;
%find error in transformation
err=r2*w+nut_spread(r2,r0)-r1;
return;

% %%-----------------------------
% function p=nut_spread(head,r0)
% p=r0(:)*ones(1,size(head,1),1);
% p=p';
