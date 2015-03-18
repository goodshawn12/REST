function [beam,outname] = nut_generate_beamforming_activations(outname,timept1,timept2,algorithm)
% NUT_GENERATE_BEAMFORMING_ACTIVATIONS
%
% Call appropriate beamforming algorithm
%
% Uses global nuts bolts.

savefullvolume = true;
saveweights = false;
numoutbeams=1;  %usually only one output, but algo1/1a has two

if(~exist('algorithm','var'))
    algorithm = 'LCMV Scalar Beamformer';
end


% Beamcore selection:
% true = the choice of a new generation -- partially vectorized and otherwise optimized for matlab
%        ~45x faster than original
% false = original for geezers
% fastway = true;

restrict_interval = false;  % true = output time series only contains chosen interval
% false = output t.s. contains whole interval,
% but weights computed using only covariance of chosen interval


global nuts bolts ndefaults;
% BOLTS contains temporary variables needed for communication between the
% different beamformer m-files: nut_generate_beamforming_activations.m,
% nut_cov_eigen.m, and nut_beamforming_gui.m/fig
% These variables are not needed for subsequent steps, so bolts is cleared
% upon successful beamformer generation.
% Might be a good idea to replace with appropriate function arguments/outputs.

% wn=0; % flag for weight normalization
% cn=1; % flag for column normalization ( aka lead field normalization )

Lp = nuts.Lp(nuts.meg.goodchannels,:,:);

% pack;

% oh yes, you guessed it... it's that time again when we come up with an
% ugly hack to get matlab to do what we want. what now? well, because
% matlab didn't design their cancel button properly, we have to do a
% little jig, cough twice, lucidly dream about making love to the
% matlab deities, and then offer a libation of premium malt liquor
% (Olde English tested, your results may vary with other brands) while
% partaking ourselves. simple, huh?
%
% bolts.stop=false; % keep track of cancel button hit, start off false
% if(ndefaults.bf.viewwaitbar)
%     bolts.barhandle = waitbar(0,'Please wait... computing inverse...','CreateCancelBtn','global bolts;bolts.stop=true;delete(gcbf)');
% end
% 
%%%%%%% end jig, begin substance abuse

switch algorithm
    case {'LCMV Scalar Beamformer','Eigenspace Scalar Beamformer'}
        if bolts.flags.cn
            display('leadfield will be normalized')
        else
            display('leadfield will not be normalized')
        end
        bolts.flags.LCMVcn=bolts.flags.cn;
        bolts.flags.cn=0;
    case {'Region Suppression'}
        if bolts.flags.cn
            display('leadfield will be normalized')
        else
            display('leadfield will not be normalized')
        end
        bolts.flags.RScn=bolts.flags.cn;
        bolts.flags.cn=0;
%     case {'Yu Yeh GEIB'} % deprecated
%         bolts.flags.YYcn=bolts.flags.cn;
%         bolts.flags.cn=0;
    case {'Thresholded Lead Field'}
        if bolts.flags.cn
            warning('column normalizing lead field not allowed for Thresholded lead field algo')
        end
        bolts.flags.cn=0;
%     case {'sLORETA','sLORETA xyz'}
    case {'sLORETA','swLORETA','dSPM','MinNorm','MinNorm_Scalar','TRACS Beamformer'}
        bolts.flags.cn=1; warning('enforcing lead-field normalization for this inverse method');
        bolts.flags.wn=0; warning('enforcing no weights normalization for this inverse method');
    case {'Beamspace noES','LCMV Vector Beamformer','Eigenspace Vector Beamformer','Beamspace','Point Suppression','Saketini','NSEFALoc','Correlate Columns','BF LFerror','BF LFerror vector','BFprior','Champagne','SAM2 Beamformer','Amber'}
        % enforce no lead-field column normalization
        bolts.flags.cn = 0; warning('enforcing no lead-field normalization for this inverse method');
    otherwise
        warning('this algo too cool to care about cn/wn issues?');
end


switch algorithm
    case {'LCMV Scalar Beamformer','LCMV Vector Beamformer','SAM2 Beamformer'}
        data.InvRyy=bolts.params.InvRzz1;
    case {'agmnrug'}
        data.meg = bolts.meg;
    case {'Eigenspace Scalar Beamformer','Eigenspace Vector Beamformer','Region Suppression','Point Suppression'}
        data.InvRyy=bolts.params.InvRzz1;
        data.InvES=bolts.params.InvES1;
    case {'Thresholded Lead Field'}
        data.InvRyy=bolts.params.InvRzz1;
        data.Ryy=bolts.params.Rzz1;
    case {'Beamspace'}
        data.InvRyy=bolts.params.InvRzz1;
        data.Ryy=bolts.params.Rzz1;
        data.InvES=bolts.params.InvES1;
    case {'Beamspace noES','sLORETA','swLORETA','dSPM','MinNorm','MinNorm_Scalar'}
        data.Ryy=bolts.params.Rzz1;
    case {'Saketini','NSEFALoc','Correlate Columns','BF LFerror','BF LFerror Vector'}
        data.y=mean(bolts.meg,3);
        data.latency=bolts.latency;
    case {'Champagne','Amber'}
        m=max(max(max(abs(Lp))));ms=num2str(m);mm=['1' ms(end-3:end)];mn=str2num(mm);
        Lp=Lp*(1/mn); % mn is generic rescaling of order/magnitude of Lp
%         Lp=Lp*10^round(-log10(max(abs(nuts.Lp(:))))); % bring to order 1
%         Lp=Lp*10^floor(-log10(max(abs(nuts.Lp(:))))); % bring to order 1
        data.y=mean(bolts.meg,3);
        data.latency=bolts.latency;
    case {'BFprior'}
        data.InvES=bolts.params.InvES1;
        data.y=mean(bolts.meg,3);
        data.latency=bolts.latency;        
    case {'TRACS Beamformer'}
        data.evoked = bolts.evoked;
        data.meg = bolts.meg;
        data.InvRzz = bolts.params.InvRzz1;
        data.latency = bolts.latency;
        data.prestim = dsearchn(nuts.meg.latency,nuts.preprocessing.bf_ptimeinterval');
        data.goodchannels = nuts.meg.goodchannels;
        data.voxels = nuts.voxels;
        data.eegflag = isfield(nuts.meg,'eegflag') & nuts.meg.eegflag == 1;
    case {'TMP Beamformer'}
        if bolts.flags.cn
            display('leadfield will be normalized')
        else
            display('leadfield will not be normalized')
        end
        bolts.flags.LCMVcn=bolts.flags.cn;
        bolts.flags.cn=0;
    case {'FrequencyVector Beamformer2', 'FrequencyVector Beamformer3'}
        data = bolts;
    otherwise 
        error('****  this algo not added to data input list  ****');
end
if bolts.flags.doW2
    Na=(timept2-timept1);
    Nc=bolts.flags.ptime(1)-bolts.flags.ptime(2);
    [Czz,CzzMineigUnave]=nut_cov(double(bolts.meg(bolts.flags.ptime(1):bolts.flags.ptime(2),:,:)),bolts.flags.avecov);
    bolts.params.Czz=Czz;
    [TotC,TotCMineigUnave]=nut_cov(double(bolts.meg([bolts.flags.ptime(1):bolts.flags.ptime(2) bolts.flags.time(1):bolts.flags.time(2)],:,:)),bolts.flags.avecov);
%     bolts.params.TotC=(bolts.params.Czz+bolts.params.Rzz1)/2;
    bolts.params.TotC=TotC;
    bolts.params.InvC=nut_inv(bolts.params.Czz,nuts.preprocessing.invtype,[],[],CzzMineigUnave);
    bolts.params.InvT=nut_inv(bolts.params.TotC,nuts.preprocessing.invtype,[],[],TotCMineigUnave);
    if ndefaults.bf.signalspace
        if(length(nuts.preprocessing.signalspace)==1)
            signalspace = 1:nuts.preprocessing.signalspace;
        end
    else
        signalspace=nuts.preprocessing.signalspace;
    end
    bolts.params.InvCES=nut_eiginv(bolts.params.Czz,signalspace);
    bolts.params.InvTES=nut_eiginv(bolts.params.TotC,signalspace); 
   switch algorithm
        case {'LCMV Scalar Beamformer','LCMV Vector Beamformer'}
            dataC.InvRyy=bolts.params.InvC;
            dataT.InvRyy=bolts.params.InvT;
        case {'Eigenspace Scalar Beamformer','Eigenspace Vector Beamformer','Region Suppression','Point Suppression'}
            dataC.InvRyy=bolts.params.InvC;
            dataC.InvES=bolts.params.InvCES;
            dataT.InvRyy=bolts.params.InvT;
            dataT.InvES=bolts.params.InvTES;
        otherwise
            disp('this algo not coded for control period comparison yet')
    end
end

condwarning = true;
if(isunix) % suppress condition warning for knowledgeable users :)
    [crap,whoami]=unix('whoami');
    if(any(strcmp(whoami(1:(end-1)),{'sarang','zumer'})))
        condwarning = false;
    end
end
switch algorithm
    case {'Saketini','NSEFALoc'}
        condwarning=false;
end
% cond_meaningless must match logcond_meaningless found in nut_beamforming_gui.m
cond_meaningless = 1e19;
if(condwarning && (bolts.params.cond >= cond_meaningless))
    warndlg('Your results will be meaningless. Get that condition number down! (Choose a bigger time window, remove lowpass filters, and/or discard channels.)')
    %return
end


bolts.flags.signalspace=bolts.params.signalspace;
bolts.flags.progressbar=ndefaults.bf.viewwaitbar;

if isfield(nuts.meg,'chanmixMtx') && length(nuts.meg.chanmixMtx)>1
    rotateMtx=nuts.meg.chanmixMtx;
    rotateMtx(1)=[]; % first part (coil/sensor arrangement) already applied during compute_lead_field
    nut_rotatearray(Lp,rotateMtx);
end

if(bolts.flags.cn) % column normalization ( lead field normalization )
    Lp = nut_lfnorm(Lp);
end

% NUTEEG mod - Ajust for average reference
if isfield(nuts.meg,'eegflag') && nuts.meg.eegflag
    % Data is thought to be avg referenced if not otherwise specified in nuts.meg.referenceidx.
    if bolts.flags.do_avref_leadpotential && ( ~isfield(nuts.meg,'referenceidx') || (length(nuts.meg.referenceidx)>1) )
        disp('Average referencing lead potentials...');
        meanLp = mean(Lp,1);
        Lp = Lp - repmat(meanLp,[length(nuts.meg.goodchannels) 1 1]);
    elseif isscalar(nuts.meg.referenceidx) % UNTESTED!
        disp('Referencing lead potential...');
        Lp = Lp - repmat(Lp(nuts.meg.referenceidx,:,:),[length(nuts.meg.goodchannels) 1 1]);
        % We should probably remove the reference channel (which is all 0) from data and Lp?
        %goodchannels(goodchannels==nuts.meg.referenceidx)=[]; 
    end
end

% reconstruct corresponding m-file name for chosen algorithm:
% prepend "nut_" and replace spaces with underscores
algorithm = strrep(['nut_' algorithm],' ','_');


% retain the essentials for placement in s_beam
beam.params.filename = nuts.meg.filename;
%beam.meg.pathname = nuts.meg.pathname;
beam.params.goodchannels = nuts.meg.goodchannels;
beam.coreg = nuts.coreg;
beam.srate = nuts.meg.srate;
% beam.VOIvoxels = nuts.VOIvoxels;
beam.voxels = nuts.voxels;
beam.voxelsize = nuts.voxelsize;

% retain parameters for troubleshooting and for people who can't keep track of
% the analyses they've done
beam.params = bolts.params;

% if isfield(bolts,'eigs'), beam.params.eigs = bolts.params.eigs; end
% if isfield(bolts,'signalspace'), beam.params.signalspace = bolts.params.signalspace; end
% beam.params.cond = bolts.params.cond;
% if isfield(bolts,'v') & isfield(bolts,'signalspace'), beam.params.eig_timeplot=bolts.meg(time_ndx,:)*bolts.params.v(:,bolts.params.signalspace); end

beam.params.algorithm = algorithm;
beam.params.beamformertype = algorithm;
beam.params.preprocessing = nuts.preprocessing;
beam.params.invtype = nuts.preprocessing.invtype;
beam.timepts=bolts.latency; %nuts.meg.latency;
beam.params.flags=bolts.flags;
beam.bands = [0 inf];

tic
% W1 = feval(algorithm,Lp(nuts.meg.goodchannels,:,:),bolts.params.InvRzz1,bolts.params.InvES1,wn,cn);
% W1 = feval(algorithm,Lp(nuts.meg.goodchannels,:,:),bolts.params.InvRzz1,bolts.params.InvES1,bolts.flags);
%W1 = feval(algorithm,Lp(nuts.meg.goodchannels,:,:),data,bolts.flags);
W1 = feval(algorithm,Lp,data,bolts.flags);
if isempty(W1); return; end;    % Cleanly exit if beamformer was aborted
if bolts.flags.doW2
%     W2 = feval(algorithm,Lp(nuts.meg.goodchannels,:,:),bolts.params.InvC,bolts.params.InvCES,bolts.flags);
%     Wall = feval(algorithm,Lp(nuts.meg.goodchannels,:,:),bolts.params.InvT,bolts.params.InvTES,bolts.flags);
%   W2 = feval(algorithm,Lp(nuts.meg.goodchannels,:,:),dataC,bolts.flags);
%    Wall = feval(algorithm,Lp(nuts.meg.goodchannels,:,:),dataT,bolts.flags);
    W2 = feval(algorithm,Lp,dataC,bolts.flags);
    Wall = feval(algorithm,Lp,dataT,bolts.flags);
    beam.params.Czz=bolts.params.Czz;
    beam.params.TotC=bolts.params.TotC;
end

% if(bolts.stop) % kill waitbar and bust out of function if user hit cancel
%     return
% end;
toc

% if(ndefaults.bf.viewwaitbar)
%     delete(bolts.barhandle);  % close progress bar
% end
disp('please wait.... saving output')

if isstruct(W1)  %happens in saketini, tbf1/1a, etc.
    beamout=W1;
    clear W1
    W1=beamout.W1;
    switch algorithm
        case {'nut_Saketini'}
            beam.params.sakeflags=beamout.sakeflags;
            nut_save_likelihood(beamout.lmap,algorithm,outname,beam);
        case {'nut_BF_LFerror'}
            beam.params.flags=beamout.flags;
            beam.lampsi=beamout.lampsi;
            beam.Fg=beamout.Fg;
            nut_save_likelihood(beamout.lmap,algorithm,outname,beam);            
        case {'nut_BFprior'}
            numoutbeams=size(W1,4);
            beam.params.flags=beamout.flags;
            nut_save_likelihood(beamout.lqv,algorithm,outname,beam);
%             algname{1}='1';algname{2}='2';algname{3}='3';algname{4}='4';
            tmpoutname=outname;
            clear outname;
            for ii=1:numoutbeams
                outname{ii}=[tmpoutname algorithm(4:end) '_' num2str(ii) '.mat'];
            end
        case {'nut_NSEFALoc'}
            %             applyweights=false;
            numoutbeams=2;
            beam.params.tbflags=beamout.tbflags;
            beam.params.phi=beamout.phi;
            nut_save_likelihood(squeeze(beamout.lmap(:,1)),[algorithm '_1'],outname,beam);
            nut_save_likelihood(squeeze(beamout.lmap(:,2)),[algorithm '_2'],outname,beam);
            %             if ~bolts.flags.tb.applytoalldata
            %                 beam.timewindow=beam.timewindow(dsearchn(beam.timewindow,0):end);
            %             end
            algname{1}='1';algname{2}='2';
            tmpoutname=outname;
            clear outname;
            for ii=1:2
                %                 g_th=squeeze(beamout.g(1,:,:,ii));
                %                 g_ph=squeeze(beamout.g(2,:,:,ii));
                %                 if size(beamout.g,1)==3
                %                     g_z=squeeze(beamout.g(3,:,:,ii));
                %                 else
                %                     g_z=0;
                %                 end
                %                 if beamout.tbflags.nl==1
                %                     g_th=g_th';
                %                     g_ph=g_ph';
                %                     g_z=g_z';
                %                 end
                %                 beam.s_th=g_th'*beamout.phi;
                %                 beam.s_ph=g_ph'*beamout.phi;
                %                 if size(beamout.g,1)==3
                %                     beam.s_z=g_z'*beamout.phi;
                %                 else
                %                     beam.s_z=0;
                %                 end
                %                 pack;
                %                 save([outname(1:end-4) algorithm(4:end) '_' algname{ii} '.mat'],'beam');
                outname{ii}=[tmpoutname(1:end-4) algorithm(4:end) '_' algname{ii} '.mat'];
            end
        case {'nut_Correlate_Columns'}
            bolts.meg=beamout.meg;
        case {'nut_Champagne'}
            nut_save_hyperparam(beamout, outname, beam)
        otherwise
            error('why is your weight a structure?')
    end
end


% if (applyweights)      %%%%%%%%%% Solve by Beamforming : W*meg' (6.1)
if(0) % try post-normalization
    Nvoxels=size(Lp,3);
%     [Nchannel1,Ncomponent,Nvoxels]=size(W1); % W=[Nchannels X 3 X Nvoxels]
%     Nsamples=length(time_ndx); %size(bolts.meg,1);

    for i=1:Nvoxels
        solution1(1,i,:)=solution1(1,i,:)/norm(Lp(beam.meg.goodchannels,1,i));
        solution1(2,i,:)=solution1(2,i,:)/norm(Lp(beam.meg.goodchannels,2,i));
        %solution1(3,i,:)=solution1(3,i,:)/norm(Lp(beam.meg.goodchannels,3,i));
    end;
end;

if(restrict_interval)
    time_ndx = timept1:timept2;
    beam.timepts = beam.timepts(time_ndx);
else  % TEST USING WINDOW FOR Rzz, BUT BEAMFORMING OVER WHOLE INTERVAL
    time_ndx = 1:size(bolts.meg,1);
end


if(savefullvolume)
    if isfield(bolts,'yc_full')
        bolts.meg.yc_full=bolts.params.yc_full;
        finaldata=mean(bolts.params.yc_full,3);
    else
        finaldata=bolts.evoked(time_ndx,:)';
%        finaldata=mean(bolts.meg(time_ndx,:,:),3)';
    end
    if strcmp(algorithm,'nut_Champagne')
       %already saved with nut_save_hyperparam earlier
    elseif(ndims(W1)==2) % for compatibility with orientation-optimized weights
%         beam.sa{1} = W1'*finaldata;

        if(saveweights)
            beam.W{1} = W1;
            if bolts.flags.doW2
                beam.W{2}=W2;
                beam.W{3}=Wall;
            end
            % writing out the weights first in case applying them to
            % data fails for memory reasons.
            hooz=whos('beam');
            if hooz.bytes>2e9
                save(outname,'-v7.3','-struct','beam');
            else
                save(outname,'-struct','beam');
            end
        end
        if bolts.flags.doW2
            beam.sa{1}=Wall'*finaldata;
%             beam.sc{1}=Wall'*finaldata;
        else
            beam.sa{1} = W1'*finaldata;
        end
        hooz=whos('beam');
        if hooz.bytes>2e9
            save(outname,'-v7.3','-struct','beam');
        else
            save(outname,'-struct','beam');
        end
    else
        for ii=1:numoutbeams
            if any(reshape(squeeze(W1(:,1,:,ii)),size(W1,1)*size(W1,3),1))
                % writing out the weights first in case applying them to
                % data fails for memory reasons.
                if numoutbeams==1
                    beam.W{1} = W1;
                    if bolts.flags.doW2
                        beam.W{2}=W2;
                        beam.W{3}=Wall;
                    end
                    save(outname,'-struct','beam');
                else
                    %                 beam.W = squeeze(W1(:,:,:,ii));
                    save(outname{ii},'-struct','beam');
                    if ii==numoutbeams
                        tmpoutname=outname{ii};
                        clear outname
                        outname=tmpoutname;
                    end
                end
                for jj=1:size(W1,2)
                    beam.sa{1}(:,:,1,jj)=transpose(squeeze(W1(:,jj,:,ii)))*finaldata;
                end
                if any(~isreal(beam.sa{1})) % If we have complex data
                    beam.sa{1} = abs(beam.sa{1});   % Convert to magnitude
                end
                if numoutbeams==1
                    if(saveweights)
                        beam.W{1} = W1;
                        if bolts.flags.doW2
                            beam.W{2}=W2;
                            beam.W{3}=Wall;
                        end
                    end
                    hooz=whos('beam');
                    if hooz.bytes>2e9
                        save(outname,'-v7.3','-struct','beam');
                    else
                        save(outname,'-struct','beam');
                    end
                else
                    %                 beam.W = squeeze(W1(:,:,:,ii));
                    hooz=whos('beam');
                    if hooz.bytes>2e9
                        save(outname{ii},'-v7.3','-struct','beam');
                    else
                        save(outname{ii},'-struct','beam');
                    end
                    %                     save(outname{ii},'-struct','beam');
                    if ii==numoutbeams
                        tmpoutname=outname{ii};
                        clear outname
                        outname=tmpoutname;
                    end
                end
            end
        end
    end
else
    save(outname,'-struct','beam');
end
disp('done saving')
toc
    % close(gcbf); % close remaining figures
% end


