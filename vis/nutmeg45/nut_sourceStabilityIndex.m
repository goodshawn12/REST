function [ssimean,ssistd]=nut_sourceStabilityIndex(W,data,N)
% function nut_sourceStabilityIndex(W,data,N)
%
% W = inverse weight matrix: size sensors x voxels
% data = sensor data: size time x sensors x trials
% N = number of permutations 
%
% ssi: mean of half of trials correlated with mean of other half of trials,
% permuted over N times.
%
% Based on Source Stability Index paper, see paper Gary Green last author
%
% jmz, 13 Aug, 2010

T=size(data,3);
for ii=1:N
    trialperm=randperm(T);
    s1=W'*mean(data(:,:,1:floor(T/2)),3)';
    s2=W'*mean(data(:,:,(1+floor(T/2)):(2*floor(T/2))),3)';
    for jj=1:size(s1,1)
        scorr(jj,ii)=corr(s1(jj,:)',s2(jj,:)');
    end
    
end
ssimean=mean(scorr,2);
ssistd=std(scorr,[],2);