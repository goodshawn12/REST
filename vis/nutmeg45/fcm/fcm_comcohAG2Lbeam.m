function beam=fcm_comcohAG2Lbeam(comcohfile,radius)
% FCM_COMCOHAG2LBEAM  calculates mean imaginary coherence of each voxel and
%       compares to homologous contralateral voxel (i.e., creates L-images).
%       AG --> for connections between 'All' seed voxels and a test 'Grid'.
%
% beam = fcm_comcohAG2Lbeam(comcohfile,¦remspur¦,¦radius¦)
%   ( parameters in ¦¦ are optional )
%
%  COMCOHFILE  name of file containing complex coherence results.
%  RADIUS      Radius of contralateral region in mm to be used for
%              t-test of each voxel. Default is 20 mm.

%if (nargin<3 || isempty(remspur)), remspur = false;  end
if nargin<2, radius=20; end

global nuts
% nuts = load(sessionfile,'coreg','voxels','voxelsize');

load(comcohfile);
if ~isfield(CC,'method'), CC.method='ccohere'; end
[nc,nt,nf]=size(CC.coh);
nv=size(nuts.voxels,1);

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

vvv = unique(CC.comps(:,2));    % grid voxels
T=zeros(nv,nt,nf);
p=ones(nv,nt,nf);

% Calculate distances
voxma=[abs(nuts.voxels(:,1)) nuts.voxels(:,2:3)];
voxma=repmat(reshape(voxma,[nv 1 3]),[1 nv 1]);
dist=zeros(nv,nv,3);
for k=1:3
    dist(:,:,k)=voxma(:,:,k)-voxma(:,:,k)';
end
clear voxma
dist=sqrt(sum(dist.^2,3));
dist=dist+diag(Inf(nv,1));

% Bring to matrix form for faster performance
CM=NaN(nv,nv,nt,nf);
for k=1:size(CC.comps,1)
    CM(CC.comps(k,1),CC.comps(k,2),:,:)=CC.coh(k,:,:);
    %CM(CC.comps(k,2),CC.comps(k,1),:)=CC.coh(:,k);
end

voxsign=sign(nuts.voxels(:,1));
for vv=1:nv
    e = find( voxsign~=sign(nuts.voxels(vv,1)) );    
    f = find( dist(vv,e) < radius  );
    if ~isempty(f)
        currconn = isfinite(CM(vv,:,1,1));
        numconn = sum(currconn);
        I = shiftdim(CM(vv,currconn,:,:),1);   % we cannot use squeeze here, because it will but the frequency dim in the beginning if there is only 1 freq bin
        C = shiftdim(nanmean(CM(e(f),currconn,:,:),1),1);
        test = I - C;
        T(vv,:,:)=mean(test,1) ./ (std(test,[],1) ./ sqrt(numconn));
        p(vv,:,:)=2*tcdf(-abs(T(vv,:,:)),numconn-1);
    end                
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
