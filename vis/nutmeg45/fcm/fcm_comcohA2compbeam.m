function beam=fcm_comcohA2compbeam(comcohfile1,comcohfile2,sessionfile)
% calculates functional connectivity changes from one recording to another.
% beam=fcm_comcohA2compbeam(comcohfile1,comcohfile2,sessionfile)

CC1=load(comcohfile1);
CC1=CC1.CC;
if ~isfield(CC1,'method'), CC1.method='ccohere'; end
switch CC1.method
    case 'ccohere'
        CC1.coh = CC1.coh./abs(CC1.coh) .* atanh(abs(CC1.coh));     % Do Z-transform 
        CC1.coh = abs(imag(CC1.coh));                               % get abs of imaginary component
    case 'ampcorr'
        CC1.coh = abs(atanh(CC1.coh));
    case 'pli'
        CC1.coh = abs(CC1.coh);
end
CC2=load(comcohfile2);
CC2=CC2.CC;
if ~isfield(CC2,'method'), CC2.method='ccohere'; end
switch CC2.method
    case 'ccohere'
        CC2.coh = CC2.coh./abs(CC2.coh) .* atanh(abs(CC2.coh));     % Do Z-transform 
        CC2.coh = abs(imag(CC2.coh));                               % get abs of imaginary component
    case 'ampcorr'
        CC2.coh = abs(atanh(CC2.coh));
    case 'pli'
        CC2.coh = abs(CC2.coh);
end                        

if( ~isequal(CC1.comps,CC2.comps) || ~isequal(CC1.frq,CC2.frq) )
    error('The 2 comcohfiles are not similar.')
end

nuts = load(sessionfile,'coreg','voxels','voxelsize');

[nc,nt,nf]=size(CC1.coh);
nv=size(nuts.voxels,1);

if nc==size(CC1.frq,1)     % Legacy compatibility
    CC1.coh=permute(CC1.coh,[2 3 1]); 
    if nc==size(CC2.frq,1)     
        CC2.coh=permute(CC2.coh,[2 3 1]);
    end
    [nc,nt,nf]=size(CC1.coh);
end    

isspec = ( nt==1 & nf>1 );
if isspec   % if we have a spectrogram, put freq data in time bins.
    CC1.coh = permute(CC1.coh,[1 3 2]);
    CC2.coh = permute(CC2.coh,[1 3 2]);
    timepts=mean(CC1.frq,2);     
    timewin=CC1.frq;
    nt=nf; nf=1;    
    bands=[CC1.frq(1,1) CC1.frq(end,2)];
else
    timepts=mean(CC1.time,2);
    timewin=CC1.time;
    bands=CC1.frq;
end
if length(timepts)>1
    srate = 1/(timepts(2)-timepts(1));
else
    srate = 1;  % arbitrary
end

beam=struct('s',{{zeros(nv,nt,nf)}},'timewindow',timewin, ...
    'timepts',timepts,'bands',bands,'voxels',nuts.voxels, ...
    'voxelsize',nuts.voxelsize,'srate',srate,'coreg',nuts.coreg);
beam.ttestFDR=struct('T',zeros(nv,nt,nf),'tail','both','p_uncorr',ones(nv,nt,nf));

for vv=1:nv
    g=find(any(CC1.comps==vv,2));
    numconn=length(g);
    temp= ( CC2.coh(g,:,:)-CC1.coh(g,:,:) );
    beam.s{1}(vv,:,:) = mean(temp,1);
    beam.ttestFDR.T(vv,:,:) = beam.s{1}(vv,:,:) ./ ( std(temp,[],1) ./ sqrt(numconn) );
    beam.ttestFDR.p_uncorr(vv,:,:) = 2.*tcdf(-abs(beam.ttestFDR.T(vv,:,:)), numconn-1);
end
beam.ttestFDR.FDR=0.01;
[dum,beam.ttestFDR.cutoff]=nut_FDR(beam.ttestFDR.p_uncorr,[],beam.ttestFDR.FDR);

switch CC1.method
    case 'ccohere'
        lab='Imaginary Coherence';
        clab='ImCoh';
    case 'ampcorr'
        lab='Amplitude Correlation';
        clab='AEC';
    case 'pli'
        lab='PLI';
        clab=lab;
end
beam.sinfo = {lab};
beam.labels.contrasts = {lab};
if isspec
    beam.labels.xaxis = 'Frequency (Hz)';
    beam.labels.yaxis = ['\Delta ' lab];
else
    beam.labels.xaxis = 'Time (ms)';
    if nf>1
        beam.labels.yaxis = 'Frequency (Hz)';
        beam.labels.colorbar = ['\Delta ' clab];
    else
        beam.labels.yaxis = ['\Delta ' lab];
    end
end

