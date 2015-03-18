function tfrun(sessionfile, covfile, timeselect, algo, covfile_usechar)
% tfrun(sessionfile, covfile, timeselect, algo)

error('tfrun will be deprecated!  call nut_tfsrc instead!');

tic
if ( ~strcmp(upper(algo),'SAM') && ~params.dualstate )
    error('Single state analyses are only implemented for SAM so far.')
end
% print params for preservation in qsub log files
%fprintf(['to rerun this bit:\nqsub -t ' num2str(timeselect) ' ~/bin/qtfrun.csh ' sessionfile ' ' covfile ' ' algo '\n']);
% disp(['running ' num2str(timeselect) ' out of '  ' covfile ']);

if strcmp(sessionfile((end-1):end),'Lp')
    % old way with *Lp.mat files
    load(sessionfile); % must contain nuts.Lp, nuts.voxels, nuts.voxelsize, nuts.coreg
else
    nuts = load(sessionfile,'Lp','voxels','voxelsize','coreg','goodvoxels');
end
load(covfile);     % must contain R, Rcon, filtERF, params
nuts.meg.data = []; % discard raw data to save RAM

if ~exist('covfile_usechar')
%     covfile_usechar=4:length(covfile);
    covfile_usechar=1:(length(covfile)-4);
end

%saveweights = false; 

if(ischar(timeselect))
    timeselect = str2num(timeselect);
end

beam.params.cn = params.cn;
beam.params.dualstate = params.dualstate;
if isfield(params,'minEigUnAve')
    beam.params.minEigUnAve=params.minEigUnAve;
elseif isfield(params,'mineig')
    beam.params.minEigUnAve=params.mineig; % set in nut_calltfcov
end
beam.bands = params.band;
beam.voxels = nuts.voxels;
if isfield(nuts,'goodvoxels')
    beam.voxels=beam.voxels(nuts.goodvoxels,:);
end
beam.voxelsize = nuts.voxelsize;
beam.coreg = nuts.coreg;

Ract = R(:,:,timeselect); clear R
if size(Rcon,3)>1, Rcon = Rcon(:,:,timeselect); end
beam.timewindow = params.active(timeselect,:);
beam.timepts = mean(beam.timewindow);

if params.dualstate
    Rall = (Ract + Rcon)/2;
else
    Rall = Ract;
end

if isfield(nuts.meg,'goodchannels')
    nuts.Lp = nuts.Lp(nuts.meg.goodchannels,:,:);
end

if isfield(nuts,'goodvoxels')
    nuts.Lp = nuts.Lp(:,:,nuts.goodvoxels);
end

switch(algo)
    case 'SAM'
%         beam.params.cn = 0 % REMEMBER TO REMOVE THIS*************
        condRall = cond(Rall)
        condRact = cond(Ract)
        condRcon = cond(Rcon)
        [W,eta] = nut_SAM(nuts.Lp,Rall,beam.params);
        Sact = sum(W.*(Ract*W))';
        if params.dualstate
            Scon = sum(W.*(Rcon*W))';
        else
            Scon = ones(size(Sact));     % fill with 1's for single state beamformer
        end
        
        [u,s,v]=svd(Rall);
        sig=s(end,end);
        noise = sum(W.*(sig*W))';
        
        beam.params.beamformertype='SAM-timef';
    case 'SAMfd'
%         beam.params.cn = 0 % REMEMBER TO REMOVE THIS*************
        load(['W_SAMfd_' covfile]);
        Sact = sum(W.*(Ract*W))';
        Scon = sum(W.*(Rcon*W))';
        
        [u,s,v]=svd(Rall);
        sig=s(end,end);
        noise = sum(W.*(sig*W))';
        
        beam.params.beamformertype='SAM-timef';
    case 'SAMbb'
%         beam.params.cn = 0 % REMEMBER TO REMOVE THIS*************
        load(['W_SAMbb_' sessionfile]);
        Sact = sum(W.*(Ract*W))';
        Scon = sum(W.*(Rcon*W))';
        
        [u,s,v]=svd(Rall);
        sig=s(end,end);
        noise = sum(W.*(sig*W))';
        
        beam.params.beamformertype='SAM-timef';
%     case 'sLORETA'
%         [W] = nut_sLORETA(nuts.Lp,Rall);
%         Sact = sum(W.*(Ract*W))';
%         Scon = sum(W.*(Rcon*W))';
%
%         [u,s,v]=svd(Rall);
%         sig=s(end,end);
%         noise = sum(W.*(sig*W))';
%
%         beam.params.beamformertype='sLORETA-timef';
    case 'sLORETA'
        load(['W_sLORETA_' sessionfile]);
        W1 = squeeze(W(:,1,:));
        W2 = squeeze(W(:,2,:));
        W3 = squeeze(W(:,3,:));
        
        Sact = sum(W1.*(Ract*W1))' + sum(W2.*(Ract*W2))' + sum(W3.*(Ract*W3))';
        Scon = sum(W1.*(Rcon*W1))' + sum(W2.*(Rcon*W2))' + sum(W3.*(Rcon*W3))';
        
        % does noise calculation apply to sLORETA???
        [u,s,v]=svd(Rall);
        sig=s(end,end);
        noise = sum(W1.*(sig*W1))' + sum(W2.*(sig*W2))' + sum(W3.*(sig*W3))';
        
        beam.params.beamformertype='sLORETA-timef';
    case 'LCMV'
        [Wact,Wcon] = nut_TF_LCMV_Beamformer(nuts.Lp,Ract,Rcon,beam.params.cn);
        Sact = sum(Wact.*(Ract*Wact))';
        Scon = sum(Wcon.*(Rcon*Wcon))';
        beam.params.beamformertype='LCMV-timef';
        
        [u,s,v]=svd(Ract);
        sig=s(end,end);
        noise = sum(Wact.*(sig*Wact))';
    case 'ES'
        [Wact,Wcon] = nut_TF_ES_Beamformer(nuts.Lp,Ract,Rcon,beam.params.cn);
        Sact = sum(Wact.*(Ract*Wact))';
        Scon = sum(Wcon.*(Rcon*Wcon))';
        beam.params.beamformertype='ES-timef';
        
        [u,s,v]=svd(Ract);
        sig=s(end,end);
        noise = sum(Wact.*(sig*Wact))';
    case 'PW10'
        beam.params.numeigs = 10
        [Sact,Wact] = nut_Prewhitened_Beamformer(nuts.Lp,Ract,Rcon,beam.params.cn,beam.params.numeigs);
        [Scon,Wcon] = nut_Prewhitened_Beamformer(nuts.Lp,Rcon,Ract,beam.params.cn,beam.params.numeigs);
        clear crap
        Wact=[];
        Wcon = [];
        beam.params.beamformertype='PW-timef';
    case 'PW20'
        beam.params.numeigs = 20
        [Sact,Wact] = nut_Prewhitened_Beamformer(nuts.Lp,Ract,Rcon,beam.params.cn,beam.params.numeigs);
        [Scon,Wcon] = nut_Prewhitened_Beamformer(nuts.Lp,Rcon,Ract,beam.params.cn,beam.params.numeigs);
        clear crap
        
        [u,s,v]=svd(Rall);
        sig=s(end,end);
        noise = sum(Wact.*(sig*Wact))';
        noisecon = sum(Wcon.*(sig*Wcon))';
        
        beam.params.beamformertype='PW-timef';
    case 'PW30'
        beam.params.numeigs = 30
        [Sact,crap] = nut_Prewhitened_Beamformer(nuts.Lp,Ract,Rcon,beam.params.cn,beam.params.numeigs);
        [Scon,crap] = nut_Prewhitened_Beamformer(nuts.Lp,Rcon,Ract,beam.params.cn,beam.params.numeigs);
        clear crap
        Wact=[];
        Wcon = [];
        beam.params.beamformertype='PW-timef';
    case 'PW40'
        beam.params.numeigs = 40
        [Sact,crap] = nut_Prewhitened_Beamformer(nuts.Lp,Ract,Rcon,beam.params.cn,beam.params.numeigs);
        [Scon,crap] = nut_Prewhitened_Beamformer(nuts.Lp,Rcon,Ract,beam.params.cn,beam.params.numeigs);
        clear crap
        Wact=[];
        Wcon = [];
        beam.params.beamformertype='PW-timef';
end

beam.s{1} = Sact;
beam.s{2} = Scon;
beam.s{3} = noise;


if(isfield(params,'z'))
    
%     for trial=goodactivetrials
%         curractive  = active  + meg.markers.active(trial);   % take care of timewins that are relative to markers
%         for ii = 1:size(active,1)
%             timepts = dsearchn(meg.latency,curractive(ii,1)):dsearchn(meg.latency,curractive(ii,2));
%             Ra(:,:,ii,trial) = cov(meg.data(timepts,:,trial));
%             %             R(:,:,ii) = R(:,:,ii) + meg.data(timepts(ii,:),:,trial)'*meg.data(timepts(ii,:),:,trial);
%         end
%     end
%
%     for trial=goodcontroltrials
%         currcontrol = control + meg.markers.control(trial);
%         for jj = 1:size(control,1)
%             timepts = dsearchn(meg.latency,currcontrol(jj,1)):dsearchn(meg.latency,currcontrol(jj,2));
%             Rc(:,:,jj,trial) = cov(meg.data(timepts,:,trial));
%             %             R(:,:,ii) = R(:,:,ii) + meg.data(timepts(ii,:),:,trial)'*meg.data(timepts(ii,:),:,trial);
%         end
%     end


% or load filtered file from params.savefiltereddata???
if isfield(params,'filtereddata')
    if params.dualstate
        load(params.filtereddata);
        atimepts=dsearchn(meg.latency,params.active(timeselect,1)):dsearchn(meg.latency,params.active(timeselect,2));
        ctimepts=dsearchn(meg.latency,params.control(1,1)):dsearchn(meg.latency,params.control(1,2));
    for trial = 1:size(meg.data,3) % loop over trials
        %%%% rather than having precomputed R across trials, we can refilter
        %%%% original data and recompute covariance per trial here
        Ra = cov(meg.data(atimepts,:,trial));
        Pa(:,trial)  = sum(W.*(Ra*W))';
        Rc = cov(meg.data(ctimepts,:,trial));
        Pc(:,trial) = sum(W.*(Rc*W))';
    end
        clear Ra Rc
        for ii = 1:size(Pa,1) % loop over voxels
            [beam.p(ii),h(ii),stats(ii)] = nut_signrank(Pa(ii,:),Pc(ii,:));
            beam.z(ii) = stats(ii).zval;
        end
    end
else
load(params.filtereddata);


    for trial = 1:size(R,3) % loop over trials
        %%%% rather than having precomputed R across trials, we can refilter
        %%%% original data and recompute covariance per trial here
        Ra = cov(meg.data(timepts,:,trial));
        Rc = cov(meg.data(timepts,:,trial));
        Pa(:,trial)  = sum(W.*(Ra*W))';
        Pc(:,trial) = sum(W.*(Rc*W))';
    end

    clear R Rc

    for ii = 1:size(Pa,1) % loop over voxels
        [beam.p(ii),h(ii),stats(ii)] = nut_signrank(Pa(ii,:),Pc(ii,:));
        beam.z(ii) = stats(ii).zval;
    end
        
end

end





beam.params.active = params.active(timeselect,:);
beam.params.control = params.control;

if params.savepower
    save(['s_beamtf_' covfile(covfile_usechar) '_' num2str(params.active(timeselect,1)) 'to' num2str(params.active(timeselect,2)) 'ms_' algo '.mat'],'beam');
end
if params.saveweights
    if beam.params.cn
        W=single(W);
    end
    if isfield(params,'weightsvoxels') % only save some voxels if already determined which ones for ROI
        W=W(:,params.weightsvoxels);
    end
    save(['weights_' covfile(covfile_usechar) '_' num2str(params.active(timeselect,1)) 'to' num2str(params.active(timeselect,2)) 'ms_' algo '.mat'],'W*');
end
%if(saveweights)
%     save(['s_beamtf_' covfile '_' num2str(params.active(timeselect,1)) 'to' num2str(params.active(timeselect,2)) 'ms_' algo '.mat'],'beam','W');
%    save(['s_beamtf_' covfile '_' num2str(params.active(timeselect,1)) 'to' num2str(params.active(timeselect,2)) 'ms_' algo '.mat'],'beam','Wact','Wcon');
%else
%    save(['s_beamtf_' covfile '_' num2str(params.active(timeselect,1)) 'to' num2str(params.active(timeselect,2)) 'ms_' algo '.mat'],'beam');
%end

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Wact,Wcon] = nut_TF_LCMV_Beamformer(Lp, Ract, Rcon, cn)

InvRact = inv(Ract);
InvRcon = inv(Rcon);

% InvRact = nut_inv(Ract,'tikhonov');
% InvRcon = nut_inv(Rcon,'tikhonov');

flags.LCMVcn = cn;
flags.wn = false;
flags.dualstate = true;

[Wact, etaact] = nut_LCMV_Scalar_Beamformer(Lp,InvRact, flags);
[Wcon, etacon] = nut_LCMV_Scalar_Beamformer(Lp,InvRcon, flags);


function [W,eta] = nut_SAM(Lp, Rall, params)

global ndefaults
% InvRall = inv(Rall);
if isfield(params,'minEigUnAve')
    data.InvRyy=nut_inv(Rall,'minEigUnAve',[],[],params.minEigUnAve);
else
    data.InvRyy=nut_inv(Rall,'vanilla');
end


flags.LCMVcn = params.cn;
flags.wn = false;
flags.dualstate = params.dualstate;
flags

[W, eta] = nut_LCMV_Scalar_Beamformer(Lp,data, flags);


function [Wact,Wcon] = nut_TF_ES_Beamformer(Lp, Ract, Rcon, Rall, cn)

flags.cn = cn;

InvRact = inv(Ract);
InvRcon = inv(Rcon);
InvRall = inv(Rall);

% condRact = cond(Ract);
% disp(['I just dropped in to see what condition my condition was in: ' num2str((condRact),'%0.5g')]);

[u,q,v]=svd(Ract);
q=diag(q);

signalspace = 1:4;

% compute qinv by truncating and taking reciprocal of q;
% we want to keep original q for future reference and troubleshooting
qinv=zeros(size(q));
qinv(signalspace)=1./q(signalspace);

InvES=v*diag(qinv)*u';

[Wact, etaact] = nut_Eigenspace_Beamformer(Lp,InvRact, InvES, 0, flags);
[Wcon, etacon] = nut_Eigenspace_Beamformer(Lp,InvRcon, InvES, 0, flags);




function [Sact,Wact]=nut_Prewhitened_Beamformer(Lp, Ract, Rcon,cn,numeigs)
%---------------------------------------------------------
% Lp : lead field ( channel X 3 )
% InvRact : Inverse data covariance matrix
% InvES : Inverse signal-component-data covariance matrix

%%%%%% Prewhitened Beamformer NOT YET READY FOR NUTMEG GUI!!! Need to pass
%%%%%% Rzz and ES, as opposed to InvRzz and InvES...

sqrtinvRcon = Rcon^(-1/2);
sqrtRcon = Rcon^(1/2);

Rpw = sqrtinvRcon * Ract * sqrtinvRcon;

[u,s,v]=svd(Rpw);
q = diag(s);

% clear numeigs

% find significant eigenvectors
% defined as cumulative sum = 400

sigeigs = find(cumsum(q) < 400);
if(isempty(sigeigs))
    sigeigs=1;
end
% sigeigs = find((cumsum(q/sum(q)*100)) < 50);
numeigs = sigeigs(end)

if(~exist('numeigs','var'))
    qmax = length(q);
    % qmax = min(50,length(q));
    eigfig=figure;subplot(1,2,1); plot(q(1:qmax)/sum(q)*100,'og-');
    subplot(1,2,2); plot(cumsum(q(1:qmax))/sum(q)*100,'og-');
    xlabel('Eigenspectrum'); ylabel('% of Total');
    % saveas(eigfig,'eigsup.fig');
    
    
    prompt   = 'What''s the signalspace of your prewhitened covariance?';
    title    = 'Input signalspace';
    lines = 1;
    def{1}   = num2str(20);
    answer   = inputdlg(prompt,title,lines,def);
    if (isempty(answer)) msgbox('What''s wrong with you?!');error('user is whack.'); end;
    numeigs = str2num(answer{1});
    close(eigfig);
end


Us = u(:,1:numeigs);

Hs = sqrtRcon * Us * Us' * sqrtinvRcon;

switch(1)
    case 1
        HRH = Hs * Ract * Hs';
    case 2
        sumvectors = 0;
        for ii=1:numeigs
            sumvectors = sumvectors + (q(ii)-1) * u(:,ii)*u(:,ii)';
        end
        HRH = sqrtRcon * sumvectors * sqrtRcon;
    case 3
        sumvectors = 0;
        for ii=1:numeigs
            sumvectors = sumvectors + (q(ii)) * u(:,ii)*u(:,ii)';
        end
        HRH = sqrtRcon * sumvectors * sqrtRcon;
    case 4
        HRH=sqrtRcon*Us*Us'*sqrtinvRcon*Ract*sqrtinvRcon*sqrtRcon;
end


Lp2 = reshape(Lp,size(Lp,1),size(Lp,2)*size(Lp,3));
InvRactLp = reshape(inv(Ract)*Lp2,size(Lp,1),size(Lp,2),size(Lp,3));

muI = 1e-10*norm(Rcon)*eye(size(Hs,1));
JJ = inv(HRH + muI);  % JJ = inv(R_(s+n))

Sact = zeros(size(Lp,3),1);
% Snoise = Sact;
Wact = zeros(size(Lp,1),size(Lp,3));

% when performing SVD, what is the minimum nonzero eigenvalue?
% generally, rank = 2 for single sphere, rank = 3 for more complex models
rnkLp = nut_rank(Lp(:,:,1));
switch(rnkLp)
    case 2
        disp('Lead field of rank 2 detected: single sphere model assumed');
    case 3
        disp('Lead field of rank 3 detected: fancy head model assumed');
    otherwise
        error(['Doooood, your lead field is cracked out! It''s rank is ' num2str(rnkLp) ' but needs to be either 2 or 3.']);
end

invRact = inv(Ract);

% orientation optimization
for i=1:size(Lp,3)
    % equivalent to: Lp' * InvRact* Lp
    invJ = Lp(:,:,i)'*JJ*Lp(:,:,i);
%    invJ = Lp(:,:,i)'*InvRactLp(:,:,i);
%     invJall = Lp(:,:,i)'*InvRallLp(:,:,i);
    
    % We use pinv here to accomodate single-sphere head models (which
    % will be of rank 2 and fail with a true inverse).
    % For more sophisticated models, the rank will generally be 3, and
    % then pinv yields the same results as inv.
    % see nut_pinv for modifications to pinv.
    
    J = nut_pinv(invJ);
%     Jall = nut_pinv(invJall);
    
    [v,d]=svd(invJ);
    eta(i,:) = v(:,rnkLp);
    L = Lp(:,:,i)*v(:,rnkLp);
    
    if(cn)
        L = L/norm(L);
    end
    
    JJL = JJ*L;
    
    Sact(i) = inv(L' * JJL);
    
%     [u,s,v]=svd(Ract);
%     invsig=eye(size(Ract))./s(end,end);
%
%     Snoise(i)=inv(L'*invsig*L);
    
    Wact(:,i) = Hs'*JJL*Sact(i);
    
% weight normalized  -- doesn't make sense with muI regularizer...
% Sact(i) = JJ*inv(L'*JJ*JJL)*JJ;
    
    
end


% Scon = ones(size(Sact));

subtractnoise = false;
if(subtractnoise)
    [u,s,v]=svd(Rall);
    sig=s(end,end)*eye(size(Rall));
    % noiseact = [sum(Wact1.*(sig*Wact1))' sum(Wact2.*(sig*Wact2))' sum(Wact3.*(sig*Wact3))'];
    noiseact = sum(Wact.*(sig*Wact))';
    noisecon = sum(Wcon.*(sig*Wcon))';
    
    Sact = Sact - noiseact;
%     Scon = Scon - noisecon;
end

% [u,s,v]=svd(Ract);
% sig=s(end,end)*eye(size(Ract));
% noiseact = sum(W.*(sig*W))';
% Scon = noiseact;

% Scon =ones(size(Sact));
