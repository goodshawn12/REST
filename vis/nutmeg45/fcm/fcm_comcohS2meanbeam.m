function beam=fcm_comcohS2meanbeam(comcohfile,typ)
% FCM_COMCOHS2MEANBEAM calculates mean coherence of each voxel and
%      creates beam structure. 
%      S --> for configurations with connections between selected seed voxels
%            and the rest of the brain.
%
% beam = fcm_comcohS2meanbeam(comcohfile,¦avgtyp¦)
%     ( parameters in ¦¦ are optional )
%
%  COMCOHFILE  name of file containing complex coherence results.
%  AVGTYP      'all' to average across all connections of each voxel (default).
%              'interhomo' to average across homologous interhemispheric
%              connections only.
%              'interhetero' to average across heterologous interhemispheric 
%              connections only.
%              'intra' to average across intrahemispheric connections only.
%              'none' to average across seed voxels.
%              If there is only 1 seed voxel, no averaging will be done and
%              the beam output structure will contain functional connectivity of
%              each brain voxel with the seed voxel.


if(nargin<2 || isempty(typ)), typ='all'; end
%if nargin<3, dostats = false; end

global nuts
% nuts = load(sessionfile,'coreg','voxels','voxelsize');

load(comcohfile);
if ~isfield(CC,'method'), CC.method='ccohere'; end
[nc,nt,nf]=size(CC.coh);
nv=size(nuts.voxels,1);
seedvox = unique(CC.comps(:,1));
n1 = length(seedvox);

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
    'timepts',timepts,'bands',bands, ...
    'voxels',nuts.voxels(seedvox,:), 'voxelsize',nuts.voxelsize, ...
    'srate',srate,'coreg',nuts.coreg);    

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

if length(seedvox)==1   % in case of a single seed voxel
    temp = zeros(nv,nt,nf);
    temp(CC.comps(:,2),:,:) = CC.coh;
    beam.voxels=nuts.voxels;
else                    % in case of seed region
    switch lower(typ)
        case 'all'
            % Bring comps to matrix form for faster performance
            CM=nan(nv,nv);
            for k=1:size(CC.comps,1)
                CM(CC.comps(k,1),CC.comps(k,2))=k;
                CM(CC.comps(k,2),CC.comps(k,1))=k;
            end
            
            temp=zeros(n1,nt,nf);
            for vv=1:n1
                %g=find(any(CC.comps==vv,2));
                g = CM(seedvox(vv),:);
                g(seedvox(vv))=[];
                temp(vv,:,:)=mean(CC.coh(g,:,:),1);
            end
        case 'interhomo'
            % Calculate distances
            voxma=[abs(beam.voxels(:,1)) beam.voxels(:,2:3)];
            voxma=repmat(reshape(voxma,[nv 1 3]),[1 nv 1]);
            dist=zeros(nv,nv,3);
            for k=1:3
                dist(:,:,k)=voxma(:,:,k)-voxma(:,:,k)';
            end
            clear voxma
            dist=sqrt(sum(dist.^2,3));
            dist=dist+diag(Inf(nv,1));

            % Bring comps to matrix form for faster performance
            CM=nan(nv,nv);
            for k=1:size(CC.comps,1)
                CM(CC.comps(k,1),CC.comps(k,2))=k;
                CM(CC.comps(k,2),CC.comps(k,1))=k;
            end
            %CM=CM+diag(NaN(nv,1));
            
            temp=zeros(n1,nt,nf);
            voxsign=sign(beam.voxels(:,1));
            for vv=1:n1
                e = find( voxsign~=sign(beam.voxels(seedvox(vv),1)) );
                %distances = nut_rownorm(nut_coord_diff(nuts.voxels(e,:),cc));
                f = find( dist(seedvox(vv),e)<20 );
                %cons=sort([vv*ones(length(f),1) e(f)],2);
                %g=find(ismember(CC.comps,cons,'rows'));
                if ~isempty(f)
                    temp(vv,:,:)=mean(CC.coh(CM(seedvox(vv),e(f)),:,:),1);
                else
                    temp(vv,:,:)=0;
                end
            end
        case 'interhetero'
            % Calculate distances
            voxma=[abs(beam.voxels(:,1)) beam.voxels(:,2:3)];
            voxma=repmat(reshape(voxma,[nv 1 3]),[1 nv 1]);
            dist=zeros(nv,nv,3);
            for k=1:3
                dist(:,:,k)=voxma(:,:,k)-voxma(:,:,k)';
            end
            clear voxma
            dist=sqrt(sum(dist.^2,3));
            dist=dist+diag(Inf(nv,1));

            % Bring comps to matrix form for faster performance
            CM=nan(nv,nv);
            for k=1:size(CC.comps,1)
                CM(CC.comps(k,1),CC.comps(k,2))=k;
                CM(CC.comps(k,2),CC.comps(k,1))=k;                
            end
            %CM=CM+diag(NaN(nv,1));
            
            temp=zeros(n1,nt,nf);
            voxsign=sign(beam.voxels(:,1));
            for vv=1:nv
                e = find( voxsign~=sign(beam.voxels(seedvox(vv),1)) );
                f = ( dist(seedvox(vv),e)>20 );
                temp(vv,:,:)=mean(CC.coh(CM(seedvox(vv),e(f)),:,:),1);
            end
        case 'intra'          
            % Bring comps to matrix form for faster performance
            CM=nan(nv,nv);
            for k=1:size(CC.comps,1)
                CM(CC.comps(k,1),CC.comps(k,2))=k;
                CM(CC.comps(k,2),CC.comps(k,1))=k;                
            end
            %CM=CM+diag(NaN(nv,1));

            temp=zeros(n1,nt,nf);
            voxsign=sign(beam.voxels(:,1));
            for vv=1:nv
                e = ( voxsign==sign(beam.voxels(seedvox(vv),1)) );
                g = CM(seedvox(vv),e);
                g = g(isfinite(g));
                temp(vv,:,:)=mean(CC.coh(g,:,:),1);
            end
            
        case 'seed'
            % Bring comps to matrix form for faster performance
            CM=nan(nv,nv);
            for k=1:size(CC.comps,1)
                CM(CC.comps(k,1),CC.comps(k,2))=k;
                CM(CC.comps(k,2),CC.comps(k,1))=k;
            end
            
            temp=nan(nv,nt,nf,n1);
            for vv=1:n1
                %g=find(any(CC.comps==vv,2));
                g = CM(seedvox(vv),:);
                g(seedvox(vv))=[];
                temp(setdiff(1:nv,vv),:,:,vv)=CC.coh(g,:,:);
            end            
            temp=nanmean(temp,4);
            beam.voxels=nuts.voxels;
            
        otherwise
            error('typ must be ''all'', ''interhomo'', ''interhetero'', ''intra'', or ''seed''')
    end
end
    
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
if length(seedvox)>1
    beam.connectionfile = comcohfile;
end

