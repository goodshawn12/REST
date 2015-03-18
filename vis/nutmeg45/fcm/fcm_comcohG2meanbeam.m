function beam=fcm_comcohG2meanbeam(comcohfile)
% FCM_COMCOHG2MEANBEAM calculates mean coherence of each voxel and
%   creates beam structure. 
%   G --> for configurations with connections to test grid.
%
% beam = fcm_comcohG2meanbeam(comcohfile,¦remspur¦)
%     ( parameters in ¦¦ are optional )
%
%  COMCOHFILE  name of file containing complex coherence results.

%if nargin<2, dostats = false; end

global nuts
% nuts = load(sessionfile,'coreg','voxels','voxelsize');

load(comcohfile);
if ~isfield(CC,'method'), CC.method='ccohere'; end
[nc,nt,nf]=size(CC.coh);

if nc==size(CC.frq,1)     % Legacy compatibility
    CC.coh=permute(CC.coh,[2 3 1]); 
    [nc,nt,nf]=size(CC.coh);
end    

% Do Z-transform 
CC.coh = fcm_fisherZtf(CC.coh,CC.method,1);

% if dostats
%     uci = norminv(0.95,0,1/sqrt(2*CC.N));  % Upper bound of one-tailed 95% confidence interval
%     spurious = ( abs(CC.coh) < uci );
%     
%     CC.coh(spurious)=0+0i;         % Set spurious to 0
%     clear spurious
%     
% %     if strcmp(comcohfile(end-3:end),'.mat'), comcohfile=comcohfile(1:end-4); end
% %     comcohfile = [comcohfile 'r'];
% %     save(comcohfile,'CC');   % save new comcoh file 
% end   

V1 = unique(CC.comps(:,1));
n1 = length(V1);

temp=zeros(n1,nt,nf);
for vv=1:n1
    g = (CC.comps(:,1)==V1(vv));
    temp(vv,:,:) = mean(CC.coh(g,:,:),1);
end
    
% Inverse Fisher
temp = fcm_fisherZtf(temp,CC.method,-1);

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
    'timepts',timepts,'bands',bands,'voxels',nuts.voxels(V1,:), ...
    'voxelsize',nuts.voxelsize,'srate',srate,'coreg',nuts.coreg);   

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
        beam.s{1} = temp; clear temp
        ylab = 'Generalized Lagged Coherence';
        clab = ylab;
        beam.sinfo={ylab};          
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
