function [tout,df]=nut_ttest_meanStd(m1,m2,s1,s2,n1,n2)
% function [tout,df]=nut_ttest_meanStd(m1,m2,s1,s2,n1,n2)

num=m1-m2;
den=sqrt((s1.^2)/n1+(s2.^2)/n2);
tout=num./den;

df=(((s1.^2)/n1+(s2.^2)/n2).^2)./((((s1.^2)/n1).^2)/(n1-1)+(((s2.^2)/n2).^2)/(n2-1));




