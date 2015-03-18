function beam=fcm_comcohE2meanbeam(comcohfile)
% FCM_COMCOHE2MEANBEAM calculates mean coherence of each voxel and
%      creates beam structure. For connections with external channels (eg, EMG).
%      E --> Extracerebral channels are seed.
%
%  beam = fcm_comcohE2meanbeam(comcohfile)
%
%  COMCOHFILE  name of file containing complex coherence results.


global nuts
% nuts = load(sessionfile,'coreg','voxels','voxelsize');

load(comcohfile);
if ~isfield(CC,'method'), CC.method='ccohere'; end
[nc,nt,nf]=size(CC.coh);
numext = length(unique(CC.comps(:,2)));
if numext>1, disp('Creating separate s_beam file for each extra channel.'); end 

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

switch CC.method
    case {'ccohere' 'nccohere'}
        ylab = 'Coherence Magnitude';
        clab = 'Coh';
        sinfo = {'Imaginary Coherence' 'Magnitude Squared Coherence' 'Real Coherence'};
        labels.contrasts = {'Imaginary Coherence';'Abs Imaginary Coherence';'Magnitude Squared Coherence';'Real Coherence'};        
    case 'ampcorr'
        ylab = 'Amplitude Correlation';
        clab = 'AEC';
        sinfo = {ylab};
    case 'pli'
        ylab = 'PLI';
        clab = ylab;
        sinfo={ylab};  
    case 'glcohere'
        ylab = 'Generalized Lagged Coherence';
        clab = ylab;
        sinfo={ylab};  
end
        
if isspec
    labels.xaxis = 'Frequency (Hz)';       
    labels.yaxis = ylab;
else
    labels.xaxis = 'Time (ms)';
    if nf>1
        labels.yaxis = 'Frequency (Hz)';
        labels.colorbar = clab;
    else
        labels.yaxis = ylab;
    end
end  

beam=struct('s',cell(numext,1),'timewindow',timewin, ...
    'timepts',timepts,'bands',bands,'voxels',nuts.voxels, ...
    'voxelsize',nuts.voxelsize,'srate',srate,'coreg',nuts.coreg, ...
    'sinfo',{sinfo},'labels',labels);  

for bb=1:numext
    cc = ( CC.comps(:,2)==bb );

    switch CC.method
        case {'ccohere' 'nccohere'}
            beam(bb).s = cell(1,3);
            beam(bb).s{1}=imag(CC.coh(cc,:,:));
            beam(bb).s{2}=abs(CC.coh(cc,:,:)).^2;
            beam(bb).s{3}=real(CC.coh(cc,:,:));     
            
        otherwise
            beam(bb).s = {abs(CC.coh(cc,:,:))};
    end
end
