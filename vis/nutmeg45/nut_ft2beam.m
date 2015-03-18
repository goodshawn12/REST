function beam = nut_ft2beam(ft,mripath,mristruct,fid,flags)
% beam = nut_ft2beam(ft,mripath,mristruct,fid,flags)
%
% converts FieldTrip analysis structure to s_beam for NUTMEG visualization
%
% ft:  structure from fieldtrip output from ft_sourceanalysis and/or
%      ft_sourcedescriptives
% mripath, mristruct, fid: (see nut_ftmriN2coreg) (may be empty)
% flags:    .avetime
%           .mmvoxflag (see nut_ftmriN2coreg)

global coreg
if(exist('mripath','var')) & size(mripath,1)~=0
   if ~exist('fid','var')
       error('you must also give fiducial file along with MRI file')
   end
   if ~exist('mristruct','var')
       error('you must also give mristruct along with MRI file; may be empty!')
   end
   if ~exist('mmvoxflag','var')
       error('you must also give mmvoxflag along with MRI file')
   end
   nut_ftmriN2coreg(mripath,mristruct,fid,mmvoxflag)
elseif (exist('mristruct','var'))
   nut_ftmriN2coreg([],mristruct)
else
   nut_ftmriN2coreg([]);
end
beam.coreg=coreg;
clear global coreg
beam.coreg.orientation = 1;


if(isfield(ft,'time'))
    if flags.avetime % ave time window edges
        beam.timepts=mean(1000*ft.time);
        beam.timewindow=1000*[ft.time(1) ft.time(end)];
        beam.srate=1;
    else % use actually .time field as real time points
        beam.timepts = 1000*ft.time';
        beam.timewindow=[ft.time-.5*ft.cfg.previous.t_ftimwin ft.time+.5*ft.cfg.previous.t_ftimwin];
        if numel(ft.time)>1
            beam.srate = 1/(ft.time(2)-ft.time(1));
        end
    end
elseif(isfield(ft,'latency'))
    beam.timepts=1000*ft.latency';
    if length(ft.latency)>1
        beam.srate=1/(ft.time(2)-ft.time(1)); % rate of time window spacing, not of data acquisition
    else
        beam.srate=1000; % fake
    end
    beam.timewindow=[beam.timepts-.5 beam.timepts+.5]; % FIXME
else  % no time if this is stats, so for now we need to pretend
    beam.timepts=[1];
    beam.timewindow = [0.5 1.5]; % fake!
    beam.srate = 1000; % fake!
end
if(isfield(ft.cfg,'frequency'))
    if strcmp(ft.cfg.method,'dics')
        if isfield(ft.cfg,'foi')
            cfg=ft.cfg;
        elseif isfield(ft.cfg.previous,'foi')
            cfg=ft.cfg.previous;
        elseif isfield(ft.cfg.previous.previous,'foi')
            cfg=ft.cfg.prevoius.previous;
        else
            warning('how deep into your previous cfg is your foi?')
        end
        foiind=dsearchn(cfg.foi',ft.freq);
        if isfield(cfg,'tapsmofrq')
            beam.bands=[cfg.foi(foiind)-.5*cfg.tapsmofrq(foiind) cfg.foi(foiind)-.5*cfg.tapsmofrq(foiind)];
        else
            try 
                lowerf=cfg.foi(foiind-1);
            catch
                lowerf=cfg.foi(foiind)-2; % FIXME: how else get lower bound on band?
            end
            try 
                higherf=cfg.foi(foiind+1);
            catch
                higherf=cfg.foi(foiind)+2; % FIXME: how else get upper bound on band?
            end
            beam.bands=[cfg.foi(foiind)-.5*(cfg.foi(foiind)-lowerf) cfg.foi(foiind)+.5*(higherf-cfg.foi(foiind))];
        end
    elseif strcmp(ft.cfg.method,'lcmv')
        if isfield(ft.cfg,'bpfreq')
            beam.bands=ft.cfg.bpfreq;
        end
    end
end
    
if(isfield(ft,'pos'))
    % FT provides cm, we want mm
    beam.voxels = ft.pos(ft.inside,:)*10;
else
    % standard MNI voxel grid
    warning('no voxel grid found. using standard MNI 2mm spaced grid');
    beam.voxels=nut_coordgrid(-89:2:91,-125:2:91,-71:2:109);
end
    
% beam.voxels = voxels(ft.inside,:);
% voxelsize = voxels(2,1)-voxels(1,1);
voxelsize = max(nut_coord_diff(beam.voxels(2,:),beam.voxels(1,:)));
beam.voxelsize = [voxelsize voxelsize voxelsize];

beam.params.algorithm = 'FieldTrip';

if isfield(ft,'avg')
    avtr=1;
elseif isfield(ft,'trial')
    avtr=2;
else
    avtr=0;
end

if isfield(ft,'mom')
    beam.params.beamformertype = 'FieldTrip-mom';
    s = cat(3,ft.mom{ft.inside});
    for jj=1:size(s,1)
        beam.s{1}(:,:,:,jj) = squeeze(s(jj,:,:))';
    end
elseif [avtr==1 && isfield(ft.avg,'mom')] 
    beam.params.beamformertype = 'FieldTrip-avg-mom';
    s = cat(3,ft.avg.mom{ft.inside});
    for jj=1:size(s,1)
        beam.s{1}(:,:,:,jj) = squeeze(s(jj,:,:))';
    end
elseif(isfield(ft,'stat'))
    beam.params.beamformertype = 'FieldTrip-stats';
    beam.s{1} = ft.stat(ft.inside);
    beam.p = ft.prob(ft.inside);
elseif avtr==1 && isfield(ft.avg,'pow')
    beam.s{1}=ft.avg.pow(ft.inside);
    if size(beam.s{1},1)==1
        beam.s{1}=beam.s{1}';
    end
elseif avtr==2 && isfield(ft.trial,'pow')
    beam.params.beamformertype = 'FieldTrip-trial-pow';
    for ii=1:length(ft.trial)
        beam.s{1}(:,ii)=ft.trial(ii).pow(ft.inside);
    end
    beam.dimord='pos_rpt';
else
    error('I only know about converting "mom" and "stat" and "avg" and "trial" -- what''s this that you''ve given me?!');
end
% if(isfield(ft,'avg'))
%     beam.noise = ft.avg.noise(ft.inside);
%     beam.power = ft.avg.pow(ft.inside);
% end
if isfield(ft,'sem')
    beam.s{length(beam.s)+1}=ft.sem.pow(ft.inside);
end
if isfield(ft,'var')
    beam.s{length(beam.s)+1}=ft.var.pow(ft.inside);
end
beam.params.ft_cfg=ft.cfg;

if isfield(ft,'filter') 
    beam.W=cell2mat({ft.filter{ft.inside}}');
elseif [avtr==1 && isfield(ft.avg,'filter')]
    beam.W=cell2mat({ft.avg.filter{ft.inside}}');
end
