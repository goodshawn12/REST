function beam=fcm_comcohA2meanbeam(comcohfile,typ)
% FCM_COMCOHA2MEANBEAM calculates mean coherence of each voxel and
%      creates beam structure. A --> All voxels are seed voxels.
%
% beam = fcm_comcohA2meanbeam(comcohfile,¦typ¦)
%     ( parameters in ¦¦ are optional )
%
%  COMCOHFILE  name of file containing complex coherence results.
%  TYP         'all' to average across all connections of each voxel (default).
%              'interhomo' to average across homologous interhemispheric
%              connections only.
%              'interhetero' to average across heterologous interhemispheric 
%              connections only.
%              'intra' to average across intrahemispheric connections only.


if(nargin<2 || isempty(typ)), typ='all'; end
%if nargin<4, dostats = false; end

global nuts

% nuts = load(sessionfile,'coreg','voxels','voxelsize');
load(comcohfile);
if ~isfield(CC,'method'), CC.method='ccohere'; end
useroi = isfield(CC,'roifile');
[nc,nt,nf]=size(CC.coh);
if useroi
    load(CC.roifile);
    if ~isfield(R,'goodroi'), R.goodroi=1:82; end
    R.roilabel = R.roilabel(R.goodroi);
    nv=length(R.goodroi);
else
    nv=size(nuts.voxels,1);
end

ismrx = (nc==nt);  % this is usually false
if ismrx
    [nv2,dum,nt,nf]=size(CC.coh); clear dum
else
    nv2=length(unique(CC.comps(:)));
end
if nv~=nv2
    error('Session file does not seem to match the connectivity file.')
end
clear nv2
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
%     CC.coh(spurious)=0;         % Set spurious to 0
%     clear spurious
%     
%     %     if strcmp(comcohfile(end-3:end),'.mat'), comcohfile=comcohfile(1:end-4); end
%     %     comcohfile = [comcohfile 'r'];
%     %     save(comcohfile,'CC');   % save new comcoh file 
% end   

temp=zeros(nv,nt,nf);
switch lower(typ)
    case 'all'
        if ~ismrx
            % Bring comps to matrix form for faster performance
            CM=zeros(nv,nv);
            for k=1:size(CC.comps,1)
                CM(CC.comps(k,:),CC.comps(k,:))=k;
            end
            CM=CM+diag(nan(nv,1));

            for vv=1:nv
                %g=find(any(CC.comps==vv,2));
                g = CM(vv,:);
                g(vv)=[];
                temp(vv,:,:)=nanmean(CC.coh(g,:,:),1);
            end
        else
            % Remove diagonals
            for tt=1:nt
                for ff=1:nf
                    CC.coh(:,:,tt,ff)=CC.coh(:,:,tt,ff)+diag(nan(nv,1));
                end
            end
            temp = permute(nanmean(CC.coh,2),[1 3 4 2]);
        end
    case 'interhomo'
        if useroi
            % Bring comps to matrix form for faster performance
            CM=zeros(nv,nv);
            for k=1:size(CC.comps,1)
                CM(CC.comps(k,:),CC.comps(k,:))=k;
            end

            for rr=1:nv
                f = strmatch(R.roilabel{rr}(1:end-1),R.roilabel);
                if ~isempty(f)
                    temp(rr,:,:) =  CC.coh(CM(rr,f),:,:);
                else
                    temp(rr,:,:) = 0;
                end
            end
            
        else
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

            % Bring comps to matrix form for faster performance
            CM=zeros(nv,nv);
            for k=1:size(CC.comps,1)
                CM(CC.comps(k,:),CC.comps(k,:))=k;
            end
            CM=CM+diag(nan(nv,1));

            voxsign=sign(nuts.voxels(:,1));
            for vv=1:nv
                e = find( voxsign~=sign(nuts.voxels(vv,1)) );
                %distances = nut_rownorm(nut_coord_diff(nuts.voxels(e,:),cc));
                f = find( dist(vv,e)<20 );
                %cons=sort([vv*ones(length(f),1) e(f)],2);
                %g=find(ismember(CC.comps,cons,'rows'));
                if ~isempty(f)
                    temp(vv,:,:)=mean(CC.coh(CM(vv,e(f)),:,:),1);
                else
                    temp(vv,:,:)=0;
                end
            end
        end
    case 'interhetero'
        if useroi
            side=repmat(' ',[1 nv]);  % prealloc
            for k=1:nv
                side(k) = R.roilabel{k}(end);  % should contain L or R for each ROI
            end
            
            % Bring comps to matrix form for faster performance
            CM=zeros(nv,nv);
            for k=1:size(CC.comps,1)
                CM(CC.comps(k,:),CC.comps(k,:))=k;
            end
            CM=CM+diag(nan(nv,1));
            
            for rr=1:nv
                e = find(side(rr)~=side);  % all contralateral ROIs
                f = setdiff(e,strmatch(R.roilabel{rr}(1:end-1),R.roilabel));    % except the homologous
                temp(rr,:,:) = mean(CC.coh(CM(rr,f),:,:),1);
            end
            
        else
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

            % Bring comps to matrix form for faster performance
            CM=zeros(nv,nv);
            for k=1:size(CC.comps,1)
                CM(CC.comps(k,:),CC.comps(k,:))=k;
            end
            CM=CM+diag(nan(nv,1));

            voxsign=sign(nuts.voxels(:,1));
            for vv=1:nv
                e = find( voxsign~=sign(nuts.voxels(vv,1)) );
                f = ( dist(vv,e)>20 );
                temp(vv,:,:)=mean(CC.coh(CM(vv,e(f)),:,:),1);
            end
        end
    case 'intra'      
        if useroi
            side=repmat(' ',[1 nv]);  % prealloc
            for k=1:nv
                side(k) = R.roilabel{k}(end);  % should contain L or R for each ROI
            end
            
            % Bring comps to matrix form for faster performance
            CM=zeros(nv,nv);
            for k=1:size(CC.comps,1)
                CM(CC.comps(k,:),CC.comps(k,:))=k;
            end
            CM=CM+diag(nan(nv,1));
            
            for rr=1:nv
                e = find(side(rr)==side);   % all ipsilateral ROIs
                f = setdiff(e,rr); % except the current one itself
                temp(rr,:,:) = mean(CC.coh(CM(rr,f),:,:),1);
            end
            
        else        
            % Bring comps to matrix form for faster performance
            CM=zeros(nv,nv);
            for k=1:size(CC.comps,1)
                CM(CC.comps(k,:),CC.comps(k,:))=k;
            end
            CM=CM+diag(nan(nv,1));

            voxsign=sign(nuts.voxels(:,1));
            for vv=1:nv
                e = ( voxsign==sign(nuts.voxels(vv,1)) );
                g = CM(vv,e);
                g = g(isfinite(g));
                temp(vv,:,:)=mean(CC.coh(g,:,:),1);
            end
        end
    case 'seed'
        if useroi
            if ~iscell(nuts.selvox.ipsilab)
                nuts.selvox.ipsi={nuts.selvox.ipsilab};
            end
            ni = length(nuts.selvox.ipsilab);
            for k=1:ni
                idx(k) = find(ismember(R.roilabel,nuts.selvox.ipsilab{k}));
            end
        else
            idx = nuts.selvox.ipsi; % ROI voxel indices
            ni  = length(idx);
        end
        temp= nan(nv,nt,nf,ni);

        % Bring comps to matrix form for faster performance
        CM=zeros(nv,nv);
        for k=1:size(CC.comps,1)
            CM(CC.comps(k,:),CC.comps(k,:))=k;
        end
        CM=CM+diag(nan(nv,1));

        for vv=1:ni
            %g=find(any(CC.comps==vv,2));
            g = CM(idx(vv),:);
            g(idx(vv))=[];
            temp(setdiff(1:nv,idx(vv)),:,:,vv) = CC.coh(g,:,:);
            temp(idx(vv),:,:,vv) = atanh(.99999);
        end
        temp = nanmean(temp,4);

    otherwise
        error('typ must be ''all'', ''interhomo'', ''interhetero'', ''intra'', or ''seed''')
end

% Inverse Fisher
temp = fcm_fisherZtf(temp,CC.method,-1);

% Output structure
isspec = ( nt==1 & nf>1 );
if isspec   % if we have a spectrogram, put freq data in time dimension.
    temp = permute(temp,[1 3 2]);
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
    case 'glcohere'
        beam.s{1} = temp; clear temp
        ylab = 'General Lagged Coherence';
        clab = ylab;
        beam.sinfo={ylab};        
    case 'pli'
        beam.s{1} = temp; clear temp
        ylab = 'PLI';
        clab = ylab;
        beam.sinfo={ylab};
end
if useroi
    beam.rois=beam.s; beam.s=[];
    if isfield(R,'goodvoxels')
        beam.voxels = beam.voxels(R.goodvoxels,:);
    end
    for ss=1:length(beam.rois)
        for k3=1:size(beam.rois{ss},3)
            beam.s{ss}(:,:,k3) = R.roi2voxel_tfm * beam.rois{ss}(:,:,k3);
        end
    end
    beam.R = R;
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
beam.connectionfile = comcohfile;
if isfield(CC,'baseline'), beam.baseline=CC.baseline; end

