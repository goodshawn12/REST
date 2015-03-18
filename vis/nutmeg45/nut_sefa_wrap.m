function [y_clean,cleaning,www,cy,sig,cyall,cytr,lam,scaleout] = nut_sefa_wrap(datain,nl,nm,actind,conind,avecov,sefaave,numave)

if isempty('actind') | isempty('conind')
    global nuts
    t0=dsearchn(nuts.meg.latency,0);
    actind=t0:length(nuts.meg.latency);
    conind=1:(t0-1);
end

if sefaave
    data=double(mean(datain(:,:,:),3))';
    datam=data-repmat(mean(data(:,conind),2),1,size(data,2));
    stddm=std(datam(:));
    if stddm>1e4; %mainly to catch data in T
        datam=datam/stddm;
        scaleout=stddm;
    else 
        scaleout=0;
    end
    ypre=datam(:,conind);
    ypost=datam(:,actind);
    y1cov0=[];
else
    [datcatc,datcata]=ave_sometrials(datain,numave,conind,actind);
    data=[datcatc datcata];
    datam=data-repmat(mean(data(:,1:size(datcatc,2)),2),1,size(data,2));
    stddm=std(datam(:));
    if stddm>1e4; %mainly to catch data in T
        datam=datam/stddm;
        scaleout=stddm;
    else
        scaleout=0;
    end
    ypre=datam(:,1:size(datcatc,2));
    ypost=datam(:,(size(datcatc,2)+1):end);
    y1cov0=cov(datcata');
end



% ypre=datam(:,1:t0-1);
% ypost=datam(:,t0:end);

% ypre=double(nuts.meg.data(1:t0-1,:))';
% ypost=double(nuts.meg.data(t0:end,:))';
% ypre=ypre-repmat(mean(ypre,2),1,size(ypre,2));
% ypost=ypost-repmat(mean(ypost,2),1,size(ypost,2));

nem=50;nem1=50;dosefa0=0;

[a,b,lam,nu,alp,bet,xbar,ybar,ubar1,cy,like,www,psiaa,psibb,outalp,igam]=nut_sefa(ypre,ypost,nl,nm,nem,nem1,dosefa0);
cleaning=a*www;
y_clean=a*www*datam;

sig=a*a'+diag(1./diag(lam));

if avecov
    cyall=[];
else
    cyall=zeros(size(datain,2),size(datain,2));
    cytr=zeros(size(datain,2),size(datain,2),size(datain,3));
    for ii=1:size(datain,3)
%         xtrial=www*datain(:,:,ii)';
        xtrial=www*datain(actind,:,ii)';
        rxx=xtrial*xtrial'+length(actind)*igam(1:nl,1:nl);
        cy0=a*rxx*a';
%         cy0=a*rxx*a'*stddm^2;
        cy_reg=trace(rxx*psiaa)*diag(1./diag(lam));
%         cy_reg=trace(rxx*psiaa)*diag(1./diag(lam))/stddm^2;
        cytr(:,:,ii)=cy0+cy_reg;
%         cyall=cyall+cy0+cy_reg;
    end
%     cyall=cyall/size(datain,3);
    cyall=mean(cytr,3);
end




