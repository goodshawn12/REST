function beam=fcm_comcohA2seedbeam(comcohfile,idx)
% beam=fcm_comcohA2seedbeam(comcohfile,idx)


global nuts
idx=idx(:)';    % make sure idx is row vector

% nuts = load(sessionfile,'coreg','voxels','voxelsize');
load(comcohfile);
if ~isfield(CC,'method'), CC.method='ccohere'; end
[nc,nt,nf]=size(CC.coh);
nv=size(nuts.voxels,1);
ni = length(idx);

if nc==size(CC.frq,1)     % Legacy compatibility
    CC.coh=permute(CC.coh,[2 3 1]); 
    [nc,nt,nf]=size(CC.coh);
end    

isspec = ( nt==1 & nf>1 );
if isspec   % if we have a spectrogram, put freq data in time dimension.
    CC.coh = permute(CC.coh,[1 3 2]);
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

beam=struct('s',{{[]}},'timewindow',timewin, ...
    'timepts',timepts,'bands',bands,'voxels',nuts.voxels, ...
    'voxelsize',nuts.voxelsize,'srate',srate,'coreg',nuts.coreg);    

% Do Z-transform 
CC.coh = fcm_fisherZtf(CC.coh,CC.method,1);    

temp=nan(nv,nt,nf,ni);

% Bring comps to matrix form for faster performance
CM=zeros(nv,nv);
for k=1:size(CC.comps,1)
    CM(CC.comps(k,:),CC.comps(k,:))=k;
end
CM=CM+diag(NaN(nv,1));

for vv=1:ni
    %g=find(any(CC.comps==vv,2));
    g = CM(idx(vv),:);
    g(idx(vv))=[];
    temp(setdiff(1:nv,idx(vv)),:,:,vv) = CC.coh(g,:,:);
end
temp = nanmean(temp,4);

% Inverse Fisher
temp = fcm_fisherZtf(temp,CC.method,-1);

switch CC.method
    case 'ampcorr'
        beam.s{1} = temp;  clear temp
        ylab = 'Amplitude Correlation';
        clab = 'AEC';
        beam.sinfo = {ylab};
    case {'ccohere' 'nccohere'}
        beam.s{1} = imag( temp );
        beam.s{2} = abs(temp).^2;
        beam.s{3} = real( temp ); clear temp
        beam.sinfo = {'Abs Imaginary Coherence' 'Magnitude Squared Coherence' 'Abs Real Coherence'};
        beam.labels.contrasts = {'Imaginary Coherence';'Magnitude Squared Coherence';'Real Coherence'};
        ylab = 'Coherence Magnitude';
        clab = 'Coh';
    case 'pli'
        beam.s{1} = temp; clear temp
        ylab = 'PLI';
        clab = ylab;
        beam.sinfo={ylab};
    case 'glcohere'
        beam.s{1} = temp;  clear temp
        ylab = 'General Lagged Coherence';
        clab = 'Coh';
        beam.sinfo = {ylab};        
end

if isspec
    beam.labels.xaxis = 'Frequency (Hz)';       
    beam.labels.yaxis = ylab;
else
    beam.labels.xaxis = 'Time (ms)';
    if nf>1
        beam.labels.yaxis = 'Frequency (Hz)';
        beam.labels.colorbar = clab;
    else
        beam.labels.yaxis = ylab;
    end
end   
% beam.connectionfile = comcohfile;

