function r=nut_mysign(t)
% R=NUT_MYSIGN(T)
% 
% Works much like SIGN(), but while SIGN(0)==0, MYSIGN(0)==1.

r=t;
r(find(t>=0))=1;
r(find(t<0))=-1;
return;
