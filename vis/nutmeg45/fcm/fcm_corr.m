function beam=fcm_corr(imagingdatafile,clindata,clindata2)
% FCM_CORR calculates the linear (partial) correlation between imaging data
%          and clinical parameters.
%
%   beam = fcm_corr(imagingdatafile, clindata [,covar])
%
% imagingdatafile   name of file containing pop structure created with
%                   nut_beampopulation.
% clindata          vector with clinical parameters. The number of entries must
%                   match the number of subjects in the imagingdata file.
% covar             optional confounding covariate(s). If given, partial
%                   correlations are calculated.

% Settings
FDR=0.1;
MinClusterSize = 180; % value is based on Monte Carlo simulation, but may be different for other datasets. You can use SnPM to do find Cluster Threshold of current dataset!

beam=[];
if nargin<2, help fcm_corr, return, end
dopartcorr = (nargin>2 && ~isempty(clindata2));

if ~exist(imagingdatafile,'file') && ~exist([imagingdatafile '.mat'],'file')
    error('Imagingdata file does not exist.')
end
load(imagingdatafile);

useroi=isfield(pop,'R');
if useroi
    pop.s=pop.rois;  % for now
end

[nsubj,nv,nt,nf] = size(pop.s);
if( nsubj ~= length(clindata) )
    error('Behavioral data and imaging data do not have the same number of subjects.')
end

clindata=clindata(:);

% if ~useroi
%     goodtim = find(squeeze(all(all(all(isfinite(pop.s),1),2),4)))';
%     goodfrq = find(squeeze(all(all(all(isfinite(pop.s(:,:,goodtim,:)),1),2),3)))';
% else
    % TODO: allow for bad time-freq bins in ROI data
    goodtim = 1:nt; goodfrq=1:nf;
% end

R = zeros(nv,nt,nf);
P = ones(nv,nt,nf);

%warning('off','MATLAB:divideByZero');
for ff = goodfrq
    for tt = goodtim
        goodvox = ( sum( ( pop.s(:,:,tt,ff)==0 | isnan(pop.s(:,:,tt,ff))) , 1) < 3 );
        if any(goodvox)
            if dopartcorr
                [R(goodvox,tt,ff),P(goodvox,tt,ff)] = partialcorr(pop.s(:,goodvox,tt,ff),clindata,clindata2,'rows','complete'); %,'type','Spearman'
            else
                [R(goodvox,tt,ff),P(goodvox,tt,ff)] = corr(pop.s(:,goodvox,tt,ff),clindata,'rows','complete'); %,'type','Spearman'
            end
        end
    end
end
%warning('on','MATLAB:divideByZero');

if strncmp( spm('ver'),'SPM8',4 ) 
    spmmri = [fileparts(which('spm')) filesep 'canonical' filesep 'avg152T1.nii'];
elseif strcmp( spm('ver'),'SPM2' )
    spmmri = [fileparts(which('nutmeg')) filesep 'templates' filesep 'wNormT1neuroax.img'];
else
    error('This version of SPM is currently not supported.')
end

if useroi
    Rr = R; Pr=P;
    bad = isnan(Rr);
    Rr(bad)=0; Pr(bad)=1;
    nv = size(pop.R.voxel2roi_tfm,1);
    R  = zeros(nv,nt,nf); 
    P  = ones(nv,nt,nf);
    for k3=1:nf
        R(:,:,k3) = pop.R.roi2voxel_tfm * Rr(:,:,k3);
        P(:,:,k3) = pop.R.roi2voxel_tfm * Pr(:,:,k3);
    end
end

beam = struct('s',{{R}}, 'timepts',pop.timepts, 'timewindow',pop.timewindow, ...
    'bands',pop.bands,'voxels',pop.voxels,'voxelsize',pop.voxelsize,'srate',pop.srate);
beam.coreg = struct('mripath', spmmri, ...
     'meg2mri_tfm',eye(4),'orientation',1, ...
     'norm_mripath',spmmri, ...
     'brainrender_path',[fileparts(which('spm.m')) filesep 'rend' filesep 'render_single_subj.mat']);
beam.corr = struct('behavdata',clindata,'imagingdata',imagingdatafile,'nr',pop.subjnr, ...
    'type','Pearson','tail','both','p_uncorr',P,'FDR',FDR);
if dopartcorr
    beam.corr.type = 'Partial';
    beam.corr.confoundingcovar = clindata2;
end

if useroi
    beam.rois={Rr};
    beam.corr.roi_p_uncorr=Pr;
    beam.corr.roilabel=pop.R.roilabel;
    beam.R = pop.R;
else
    beam.corr.FDR = FDR;
    voxelsblob = nut_voxels2blob(beam);
    if max(max(abs(round(voxelsblob)-voxelsblob)))<.05
        beam.corr.clusterthreshold= MinClusterSize;
        beam.corr.p_cluster_corr  = nut_clusterstats(R,P,beam,MinClusterSize); % Apply cluster size threshold
    end
end
        
% add labels
beam.sinfo = {'Correlation'};
if nt>1
    if nf>1     % TF plot
        beam.labels.xaxis = 'Time (ms)';
        beam.labels.yaxis = 'Frequency (Hz)';
        beam.labels.colorbar = 'Corr';
    else        % Freq spectrogram for now!
        beam.labels.xaxis = 'Frequency (Hz)';
        beam.labels.yaxis = 'Correlation (r)';
    end
elseif nt==1
    beam.labels.xaxis = 'Time (ms)';
    beam.labels.yaxis = 'Frequency (Hz)';
    beam.labels.colorbar = 'Corr';
end

