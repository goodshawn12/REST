function nut_tfsrc(sessionfile, covfile, timeselect, algo)
% nut_tfsrc(sessionfile, covfile, timeselect, algo)

if nargin<4, algo='SAM'; end

% print params for preservation in qsub log files
%fprintf(['to rerun this bit:\nqsub -t ' num2str(timeselect) ' ~/bin/qtfrun.csh ' sessionfile ' ' covfile ' ' algo '\n']);
% disp(['running ' num2str(timeselect) ' out of '  ' covfile ']);

% if ( ~strcmp(upper(algo),'SAM') && ~params.dualstate )
%     error('Single state analyses are only implemented for SAM so far.')
% end

if strcmpi(covfile(end-3:end),'.mat'), covfile=covfile(1:end-4); end
load([covfile '.mat']);     % must contain R, Rcon, params
if ~isfield(params,'do_avref_leadpotential'), params.do_avref_leadpotential=false; end

warning('off','MATLAB:load:variableNotFound');
nuts=load(sessionfile,'Lp','voxels','voxelsize','coreg','goodvoxels','voxor'); 
warning('on','MATLAB:load:variableNotFound');
if ~isfield(nuts,'Lp')      % in case of old saving style
    load(sessionfile);      % then this should be the ...Lp file
    if isfield(nuts.meg,'data')   % but if something went wrong
        nuts.meg.data = [];       % then discard raw data to save RAM
    end
end

if(ischar(timeselect))
    timeselect = str2num(timeselect);
end

beam.params.cn = params.cn;
beam.params.dualstate = params.dualstate;
beam.bands = params.band;
beam.voxels = nuts.voxels;
if isfield(nuts,'goodvoxels')
    beam.voxels=beam.voxels(nuts.goodvoxels,:);
end
beam.voxelsize = nuts.voxelsize;
beam.coreg = nuts.coreg;
if isfield(params,'minEigUnAve')
    beam.params.minEigUnAve=params.minEigUnAve; % should be only 1 value
    params.mineig=params.minEigUnAve; % use mineig later on    
elseif isfield(params,'mineig')
    beam.params.minEigUnAve=params.mineig(timeselect); % set in tfbf.m
else
    params.mineig=[]; 
end

Ract = R(:,:,timeselect); clear R
if size(Rcon,3)>1, Rcon = Rcon(:,:,timeselect); end
beam.timewindow = params.active(timeselect,:);
beam.timepts = mean(beam.timewindow);

if (isfield(params,'eegflag') && params.eegflag)
    if ( length(params.referenceidx)>1 && params.do_avref_leadpotential )
        disp('Average referencing lead potential...');
        meanLp = mean(nuts.Lp(params.referenceidx,:,:),1);
        nuts.Lp = nuts.Lp - repmat(meanLp,[size(nuts.Lp,1) 1 1]);
    elseif isscalar(params.referenceidx)
        disp('Referencing lead potential...');
        nuts.Lp = nuts.Lp - repmat(nuts.Lp(params.referenceidx,:,:),[size(nuts.Lp,1) 1 1]);
        goodchannels(goodchannels==params.referenceidx)=[]; 
    end
end

if size(nuts.Lp,1)>size(Ract,1)
    if exist('goodchannels','var')
        nuts.Lp = nuts.Lp(goodchannels,:,:);
    elseif( isfield(nuts,'meg') && isfield(nuts.meg,'goodchannels') )
        nuts.Lp = nuts.Lp(nuts.meg.goodchannels,:,:);
    else
        error('Leadfield channel size does not match the data.')
    end
end

if isfield(nuts,'goodvoxels')
    num_voxels=length(nuts.goodvoxels)
    if length(nuts.goodvoxels) < size(nuts.Lp,3)
        nuts.Lp = nuts.Lp(:,:,nuts.goodvoxels);
    end
end

if isfield(nuts,'voxor')        % Surface normal voxel orientations
    disp('Using surface normal voxel orientations.')
    L = zeros(size(nuts.Lp,1),1,size(nuts.Lp,3));
    for k=1:size(nuts.Lp,1)
        L(k,1,:) = dot(squeeze(nuts.Lp(k,:,:)),nuts.voxor');
    end
    nuts.Lp=L; clear L
end

switch strtok(algo,'_')
    case 'SAM'
%         beam.params.cn = 0 % REMEMBER TO REMOVE THIS*************
        if params.dualstate
            Rall = (Ract + Rcon)/2;
        else
            Rall = Ract;
        end
        condRall = cond(Rall)
        condRact = cond(Ract)
        condRcon = cond(Rcon)
        W = nut_TF_invsol(nuts.Lp,Rall,[],params,algo);
               
    case { 'LCMV' 'MinNorm' 'sLORETA' 'dSPM' }
        [Wact,Wcon] = nut_TF_invsol(nuts.Lp,Ract,Rcon,params,algo); 
%         if params.dualstate
%             Wcon = nut_TF_invsol(nuts.Lp,Rcon,params,algo);  
%         end
            
    case 'ES'
        [Wact,Wcon] = nut_TF_ES_Beamformer(nuts.Lp,Ract,Rcon,beam.params.cn);
        if params.savepower
            Sact = sum(Wact.*(Ract*Wact))';
            Scon = sum(Wcon.*(Rcon*Wcon))';
            beam.params.beamformertype='ES-timef';

            [u,s,v]=svd(Ract);
            sig=s(end,end);
            noise = sum(Wact.*(sig*Wact))';
        end
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
    otherwise
        error('NUT_TFSRC: Unknown localization algorithm.')
end

if params.savepower
   switch strtok(algo,'_')
       case 'SAM'
            Sact = sum(W.*(Ract*W))';
            if params.dualstate
                Scon = sum(W.*(Rcon*W))';
            else
                Scon = ones(size(Sact));     % fill with 1's for single state beamformer
            end

            if params.dualstate
                Rall = (Ract + Rcon)/2;
            else
                Rall = Ract;
            end
            [u,s,v]=svd(Rall);
            sig=s(end,end);
            noise = sum(W.*(sig*W))';

            beam.params.beamformertype='SAM-timef';
            
%        case 'LCMV'
%             Sact = sum(Wact.*(Ract*Wact))';
%             if params.dualstate
%                 Scon = sum(Wcon.*(Rcon*Wcon))';
%             else
%                 Scon = ones(size(Sact));     % fill with 1's for single state beamformer
%             end            
% 
%             beam.params.beamformertype='LCMV-timef';
% 
%             [u,s,v]=svd(Ract);
%             sig=s(end,end);
%             noise = sum(Wact.*(sig*Wact))';            
              
       case { 'LCMV' 'MinNorm' 'sLORETA' 'dSPM' }
            [nel,nor,nvox] = size(Wact);
            if nor>3 && nvox==1, nvox=nor; nor=1; Wact=reshape(Wact,[nel 1 nvox]); end
            Sact = zeros(nvox,1);
            
            noise = zeros(nvox,1);        % does noise calculation apply to MinNorm???
            [u,s,v]=svd(Ract);
            sig=s(end,end);
            
            for oo = 1:nor;
                Wtmp = squeeze(Wact(:,oo,:));
                Sact = Sact + sum(Wtmp .* (Ract * Wtmp))';
                noise = noise + sum(Wtmp .* (sig * Wtmp))';
            end
            
            if params.dualstate
                Scon = zeros(nvox,1);
                for oo = 1:nor;
                    Wtmp = squeeze(Wcon(:,oo,:));
                    Scon = Scon + sum(Wtmp .* (Rcon * Wtmp))';
                end              
            else
                Scon = ones(nvox,1);     % fill with 1's for single state beamformer
            end
            clear Wtmp

            beam.params.beamformertype=[algo '-timef'];         
    end

    beam.s{1} = Sact;
    beam.s{2} = Scon;
    beam.s{3} = noise;

    beam.params.active = params.active(timeselect,:);
    beam.params.control = params.control;
    save(['s_beamtf_' covfile(5:end) '_' num2str(params.active(timeselect,1)) 'to' num2str(params.active(timeselect,2)) 'ms_' algo '.mat'],'beam');
end
if params.saveweights
    save(['W_' covfile(5:end) '_' num2str(params.active(timeselect,1)) 'to' num2str(params.active(timeselect,2)) 'ms_' algo '.mat'],'W*');
end
%if(saveweights)
%     save(['s_beamtf_' covfile '_' num2str(params.active(timeselect,1)) 'to' num2str(params.active(timeselect,2)) 'ms_' algo '.mat'],'beam','W');
%    save(['s_beamtf_' covfile '_' num2str(params.active(timeselect,1)) 'to' num2str(params.active(timeselect,2)) 'ms_' algo '.mat'],'beam','Wact','Wcon');
%else
%    save(['s_beamtf_' covfile '_' num2str(params.active(timeselect,1)) 'to' num2str(params.active(timeselect,2)) 'ms_' algo '.mat'],'beam');
%end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [W,Wcon] = nut_TF_invsol(Lp,Ract,Rcon,params,algo)

% if ~isfield(params,'calceta'), params.calceta=false; end

switch strtok(algo,'_')
    case 'SAM'
        
        data.InvRyy = nut_inv(Ract,params.regularization,params.regulthres,[],params.mineig);
        
        flags.LCMVcn = params.cn;
        flags.wn = params.wn;
        %flags.dualstate = params.dualstate;
        flags.progressbar=false;
        flags
        
        if size(Lp,2)==1, algo='LCMV_Vector';   % if we already have a scalar leadfield, we do not need the scalar beamformer
        else algo='LCMV_Scalar'; 
        end

        W = feval(['nut_' algo '_Beamformer'],Lp,data, flags);
        
    case 'LCMV'
        condRact = cond(Ract)
        condRcon = cond(Rcon)
        
        flags.LCMVcn = params.cn;
        flags.wn = params.wn;
        %flags.dualstate = params.dualstate;
        flags.progressbar=false;
        flags
        
        if size(Lp,2)==1, algo='LCMV_Vector'; end   % if we already have a scalar leadfield, we do not need the scalar beamformer
        
        data.InvRyy = nut_inv(Ract,params.regularization,params.regulthres,[],params.mineig);
        W = feval(['nut_' algo '_Beamformer'],Lp,data, flags);
        W = squeeze(W);
        
        if params.dualstate
            data.InvRyy = nut_inv(Rcon,params.regularization,params.regulthres,[],params.mineig);
            Wcon = feval(['nut_' algo '_Beamformer'],Lp,data, flags);
            Wcon = squeeze(Wcon);
        else
            Wcon=[];
        end
        
    case {'sLORETA' 'MinNorm' 'dSPM'}
        data.Ryy = Ract;
        data.C   = Rcon;
%         if ( size(Lp,2)>1 && params.calceta )
%             data.InvRyy = nut_inv(R,params.regularization,params.regulthres);
%             Lp = nut_dipoleorientation(Lp,data);
%         end
        
        W = feval(['nut_' algo],Lp,data,params);
        if params.dualstate
            data.Ryy = Rcon;
            data.C   = eye(size(Rcon));
            Wcon = feval(['nut_' algo],Lp,data,params);
        else
            Wcon=[];
        end
end

%-------------------------------
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
