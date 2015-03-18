function fcm_plotroimatrix(M,frhz,tims,selroi)
% fcm_plotroimatrix(M,freq_in_Hz,time_in_ms,selroi)

if nargin>1 && ~isempty(frhz)
    switch length(frhz)
        case 1
            f = dsearchn(mean(M.frq,2),frhz);
        case 2
            f = find(M.frq(:,1)>=frhz(1) & M.frq(:,2)<=frhz(2));
        otherwise
            error('Invalid frequency limits.')
    end
else
    f = 1;
end

if nargin>2 && ~isempty(tims)
    switch length(tims)
        case 1
            t = dsearchn(mean(M.time,2),tims);
        case 2
            t = find(M.time(:,1)>=tims(1) & M.time(:,2)<=tims(2));
        otherwise
            error('Invalid time limits.')
    end
else
    t = 1;
end

if isfield(M,'r')
    X=M.r;
else
    X=M.coh;
end
if isfield(M,'p_uncorr')
    X(M.p_uncorr>.05)=0;
end

if nargin>3
    X=X(selroi,selroi,:,:);
    M.roilabel = M.roilabel(selroi);
end

figure('units','normalized','position',[0.3327    0.2381    0.4768    0.6638]);

nL=length(M.roilabel);
clim = max(max(abs(X(:,:,t,f))));
imagesc(X(:,:,t,f),[-clim clim])
set(gca,'xtick',1:nL,'ytick',1:nL,'xticklabel',[],'yticklabel',[])
for k=1:nL
    text(0,k,M.roilabel{k},'fontsize',6,'horizontalalignment','right','Interpreter','none')
end
for k=1:nL
    text(k,nL+1,M.roilabel{k},'fontsize',6,'rotation',90,'horizontalalignment','right','Interpreter','none')
end

hc=colorbar; 
if isfield(M,'info')
    axes(hc); ylabel(M.info,'fontweight','bold')
end