function nut_applyTFBFweights2(nuts, R, Rcon, params, weightsfile, timeselect, algo, covfile,covfile_usechar, voxind)
% function nut_applyTFBFweights2(sessionfile, covfile, weightsfile, timeselect, algo, covfile_usechar, voxind)
% tfrun(sessionfile, covfile, timeselect, algo)

tic
% print params for preservation in qsub log files
%fprintf(['to rerun this bit:\nqsub -t ' num2str(timeselect) ' ~/bin/qtfrun.csh ' sessionfile ' ' covfile ' ' algo '\n']);
% disp(['running ' num2str(timeselect) ' out of '  ' covfile ']);

% load(sessionfile); % must contain nuts.Lp, nuts.voxels, nuts.voxelsize, nuts.coreg
% nuts.meg.data = []; % discard raw data to save RAM
% load(covfile);     % must contain R, Rcon, filtERF, params

% tmp, jz add, remove soon
% params.dualstate=1;
% params.savepower=1;
% params.saveweights=0;

if ( ~strcmp(upper(algo),'SAM') && ~params.dualstate )
    error('Single state analyses are only implemented for SAM so far.')
end

if ~exist('covfile_usechar')
    covfile_usechar=5:length(covfile);
end

%saveweights = false; 

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

Ract = R(:,:,timeselect); clear R
if size(Rcon,3)>1, Rcon = Rcon(:,:,timeselect); end
beam.timewindow = params.active(timeselect,:);
beam.timepts = mean(beam.timewindow);

if params.dualstate
    Rall = (Ract + Rcon)/2;
else
    Rall = Ract;
end

% if isfield(nuts.meg,'goodchannels')
%     nuts.Lp = nuts.Lp(nuts.meg.goodchannels,:,:);
% end

if isfield(nuts,'goodvoxels')
    nuts.Lp = nuts.Lp(:,:,nuts.goodvoxels);
end

load(weightsfile)

Sact = sum(W.*(Ract*W))';
if params.dualstate
    Scon = sum(W.*(Rcon*W))';
else
    Scon = ones(size(Sact));     % fill with 1's for single state beamformer
end

[u,s,v]=svd(Rall);
sig=s(end,end);
noise = sum(W.*(sig*W))';

beam.s{1} = Sact;
beam.s{2} = Scon;
beam.s{3} = noise;

beam.params.active = params.active(timeselect,:);
beam.params.control = params.control;

if ~exist('voxind')
    if params.savepower
        save(['s_beamtf_' covfile(covfile_usechar) '_' num2str(params.active(timeselect,1)) 'to' num2str(params.active(timeselect,2)) 'ms_' algo '.mat'],'beam');
    end
    params.saveweights=0;
    if params.saveweights
        save(['weights_' covfile(covfile_usechar) '_' num2str(params.active(timeselect,1)) 'to' num2str(params.active(timeselect,2)) 'ms_' algo '.mat'],'W*');
    end
else
    if length(voxind)~=length(beam.voxels)
        beam.voxels=beam.voxels(voxind,:);
        if length(voxind)~=length(beam.s{1})
            beam.s{1}=beam.s{1}(voxind);
            beam.s{2}=beam.s{2}(voxind);
            beam.s{3}=beam.s{3}(voxind);
        end
    end
    save(['s_beamtf_' covfile(covfile_usechar) '_' num2str(params.active(timeselect,1)) 'to' num2str(params.active(timeselect,2)) 'ms_' algo '.mat'],'beam');
end

%if(saveweights)
%     save(['s_beamtf_' covfile '_' num2str(params.active(timeselect,1)) 'to' num2str(params.active(timeselect,2)) 'ms_' algo '.mat'],'beam','W');
%    save(['s_beamtf_' covfile '_' num2str(params.active(timeselect,1)) 'to' num2str(params.active(timeselect,2)) 'ms_' algo '.mat'],'beam','Wact','Wcon');
%else
%    save(['s_beamtf_' covfile '_' num2str(params.active(timeselect,1)) 'to' num2str(params.active(timeselect,2)) 'ms_' algo '.mat'],'beam');
%end

toc