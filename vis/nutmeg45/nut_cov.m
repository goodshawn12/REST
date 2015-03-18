function [Rzz1,mineigUnave]=nut_cov(data,avecovflag);
% [Rzz1,mineigUnave]=nut_cov(data,avecovflag);
% 
% data:      time by sensors by trials (num trials can equal 1 if input averaged data)
% avecovflag: =1 means average over trials then cov, =0 cov each trial then ave
%
% Rzz1:     output covariance matrix
% mineigUnave:  minimum eigenvalue of the output of avecovflag=0

global ndefaults

if ndefaults.rzzcov
    Rzz1=zeros(size(data,2),size(data,2));
    for jj=1:size(data,3)
        Rzz1=Rzz1+cov(data(:,:,jj));
    end
    Rzz1=Rzz1/size(data,3);
else % this is for not-mean-subtracting the data prior to covariance
    dataC = zeros([size(data,2) size(data,1) size(data,3)]);
    for ii = 1:size(data,3)
        dataC(:,:,ii) = data(:,:,ii)';
    end
    dataC = reshape(dataC,[size(data,2) size(data,1)*size(data,3)]);
    Rzz1 = dataC*dataC';
%    Rzz1=Rzz1/size(data,3); % FIXME: JZ asks: why does this line not exist here?
end
mineigUnave=min(eig(Rzz1));  % saved for 'minEigUnAve' way of regularising

if avecovflag
    clear Rzz1
    if ndefaults.rzzcov
        Rzz1=cov(mean(data,3));
    else
        Rzz1=mean(data,3)'*mean(data,3);
    end
    % else
    %     Rzz1=zeros(size(data,2),size(data,2));
    %     for jj=1:size(data,3)
    %         Rzz1=Rzz1+cov(data(:,:,jj));
    %     end
    %     Rzz1=Rzz1/size(data,3);
end
