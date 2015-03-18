function beam=fcm_comcohSC2Lbeam(comcohfile)
% FCM_COMCOHSC2LBEAM  calculates mean imaginary coherence of each voxel and
%       compares to homologous contralateral voxel (i.e., creates L-images).
%       SC --> for configurations with selected seed voxels and homologous
%              contralateral voxels.
%
% beam = fcm_comcohSC2Lbeam(comcohfile)
%
%  COMCOHFILE  name of file containing complex coherence results.


%if (nargin<2 || isempty(remspur)), remspur = false;  end

global nuts
load(comcohfile);
if ~isfield(CC,'method'), CC.method='ccohere'; end
[nc,nt,nf]=size(CC.coh);

if nc==size(CC.frq,1)     % Legacy compatibility
    CC.coh=permute(CC.coh,[2 3 1]); 
    [nc,nt,nf]=size(CC.coh);
end    

switch CC.method
    case {'ccohere' 'nccohere'}
        CC.coh = CC.coh./abs(CC.coh) .* atanh(abs(CC.coh));     % Do Z-transform 
        CC.coh = abs(imag(CC.coh));                             % get abs of imaginary component
    case 'ampcorr'
        CC.coh = abs(atanh(CC.coh));
    case 'pli'
        CC.coh = abs(CC.coh);
end                    

% if remspur
%     varreal = 1/(2*CC.N);                % Variance of real part of coherence
%     varimag = (1-varreal^2)/(2*CC.N);    % Variance of imaginary part, for large N this approaches varreal
%     uci = norminv(0.95,0,sqrt(varimag)); % Upper bound of one-tailed 95% confidence interval
%     spurious = ( CC.coh < uci );
%     
%     CC.coh(spurious)=0;         % Set spurious to 0
%     clear spurious varreal varimag
% end

nipsi = length(nuts.selvox.ipsi);
ncontra = length(nuts.selvox.contra);
vvv=unique(CC.comps(:,2));
nconn = length(vvv);

IM=nan(nipsi,nconn,nt,nf);
for v=1:nipsi
    f = ( CC.comps(:,1)==nuts.selvox.ipsi(v) );
    c = ismember(vvv,CC.comps(f,2));
    IM(v,c,:,:)=CC.coh(f,:,:);
end
CM=nan(size(IM));
for v=1:ncontra
    f = (CC.comps(:,1)==nuts.selvox.contra(v));
    c = ismember(vvv,CC.comps(f,2));    
    CM(v,c,:,:)=CC.coh(f,:,:);
end

MM = IM - repmat(nanmean(CM,1),[nipsi 1 1 1]);
clear IM CM

T=zeros(nipsi,nt,nf);
p=ones(nipsi,nt,nf);
for v=1:nipsi
    %f=find(TC(:,1)==tuv(v));
    numconn = sum(isfinite(MM(v,:,1,1)));
    T(v,:,:) = nanmean(MM(v,:,:,:),2) ./ ( nanstd(MM(v,:,:,:),[],2)./sqrt(numconn) );
    p(v,:,:) = 2*tcdf(-abs(T(v,:,:)), numconn-1);
end

% Output structure
isspec = ( nt==1 & nf>1 );
if isspec   % if we have a spectrogram, put freq data in time dimension.
    T = permute(T,[1 3 2]);
    p = permute(p,[1 3 2]);
    timepts=mean(CC.frq,2);     
    timewin=CC.frq;    
    nt=nf; nf=1;
    bands=[CC.frq(1,1) CC.frq(end,2)];
else
    timepts=mean(CC.time,2);
    timewin=CC.time;
    bands=CC.frq;
end
if length(timepts)>1
    srate = 1/(timepts(2)-timepts(1));
else
    srate = 1;  % arbitrary
end
beam=struct('s',{{T}},'timewindow',timewin, ...
    'timepts',timepts,'bands',bands,'voxels',nuts.voxels, ...
    'voxelsize',nuts.voxelsize,'srate',srate,'coreg',nuts.coreg);    

beam.sinfo={'T'};
beam.ttest.tail='both';
beam.ttest.T=T;
beam.ttest.p_uncorr = p;
beam.ttest.FDR = 0.01;    

if isspec
    beam.labels.xaxis = 'Frequency (Hz)';
    beam.labels.yaxis = 'T';
else
    beam.labels.xaxis = 'Time (ms)';
    if nf>1
        beam.labels.yaxis = 'Frequency (Hz)';
        beam.labels.colorbar = 'T';
    else
        beam.labels.yaxis = 'T';
    end
end
