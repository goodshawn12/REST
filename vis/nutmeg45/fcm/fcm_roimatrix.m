function M=fcm_roimatrix(CC,typ,doz,flag,interpfrq)
% M=fcm_roimatrix(CC,typ,do_z_tf)
% typ       coherence type (has no effect for other functional connectivity measures)
%           1 = abs imag coh (default), 2 = mag squared coh, 3 = abs real coh
% do_z_tf   transform connectivity values to z-scores (default is false).

if nargin<2 || isempty(typ), typ=1; end
if nargin<3 || isempty(doz), doz=false; end
if nargin<4, flag=true; end  % output as structure

load(CC.roifile);


numroi=length(R.goodroi);
[nconn,nt,nf] = size(CC.coh);

if strcmp(CC.method,'ccohere')
    switch typ
        case 1
            X = abs(imag(CC.coh));
            zo = 0;
            M.info = 'Abs Imag Coh';
        case 2
            X = abs(CC.coh).^2;
            zo = 1;
            M.info = 'Mag Squared Coh';
        case 3
            X = abs(real(CC.coh));
            zo = 1;
            M.info = 'Abs Real Coh';            
        otherwise
            error('Invalid type.')
    end
else
    X = abs(CC.coh);
    zo = 1;
end

if doz
    M = nanmean(X,1);
    X = X-repmat(M,[nconn 1 1]); clear M
    S = nanstd(X,[],1);
    X = X./repmat(S,[nconn 1 1]); clear S
end

if nargin>4 && ~isempty(interpfrq)
    interpfrq=interpfrq(:);
    X2=interp1(mean(CC.frq,2),permute(X,[3 1 2]),interpfrq,'linear');
    X = permute(X2,[2 3 1]); clear X2
    d = (interpfrq(2)-interpfrq(1))/2;
    CC.frq = [interpfrq-d interpfrq+d]; clear d
    nf = length(interpfrq);
end

tmp = nan(numroi);
tmp(logical(eye(numroi))) = zo;
M.coh = repmat(tmp,[1 1 nt nf]); clear tmp

for k=1:size(CC.comps,1)
    M.coh(CC.comps(k,1),CC.comps(k,2),:,:) = X(k,:,:);
    M.coh(CC.comps(k,2),CC.comps(k,1),:,:) = X(k,:,:);
end

if flag
    M.method = CC.method;
    M.frq = CC.frq;
    M.time = CC.time;
    M.N = CC.N;
    M.goodroi  = R.goodroi;
    M.roilabel = R.roilabel(R.goodroi);
else
    M=M.coh;
end