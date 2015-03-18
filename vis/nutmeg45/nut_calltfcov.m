function nut_calltfcov(sessionfile,paramfile,lowfreq,hifreq,varargin)
% tfcov(sessionfile,paramfile,lowfreq,highfreq)
% wrapper script for nut_tfcov, intended for compilation into standalone
% varargin{1} is .mat file with filt. and params. stuff set in it

error('This is to be deprecated. Please use tfbf.m instead')

tic

if(ischar(lowfreq))
    lowfreq = str2num(lowfreq);
end
if(ischar(hifreq))
    hifreq = str2num(hifreq);
end

if(nargin<5)
    filt.class='butter';
    filt.order=4;
    filt.type='bp';
    filtname='butterbp4';
%     algo={'SAM'};
    params.cn=1
    params.qsub=0;
else
    load(varargin{1}); % load additional params (like filter settings)
    [filtpath,filtname,ext] = fileparts(varargin{1});
%     algo=cell(length(varargin)-1,1);
%     [algo{1:size(algo,1)}]=varargin{2:end};
% algo={'SAM'};
% params.qsub=0;
end

[filepath,filename,ext] = fileparts(sessionfile);
if isempty(filepath)
    filtsessionname = [filename '_' filtname '_' num2str(lowfreq) 'to' num2str(hifreq) 'Hz.mat'];
else
    filtsessionname = [filepath '/' filename '_' filtname '_' num2str(lowfreq) 'to' num2str(hifreq) 'Hz.mat'];
end
load(filtsessionname);

% load(sessionfile); % must contain nuts.Lp, nuts.voxels, nuts.voxelsize, nuts.coreg, nuts.meg
load(paramfile);  % must contain time windows ("active" and "control")

% print run parameters (for preservation in log files)
disp(filtsessionname)
disp(paramfile)



if isfield(meg,'numtrial')
    if size(meg.numtrial.stim)==1
        meg.numtrial.stim=1:meg.numtrial.stim;
        meg.numtrial.rest=meg.numtrial.stim(end)+1:meg.numtrial.stim(end)+meg.numtrial.rest;
    end        
    if isfield(params,'trialact')
        trialselect{1}=meg.numtrial.(params.trialact);
        trialselect{2}=meg.numtrial.(params.trialcon);
    end
end

oldway = true;
if(oldway)
    if exist('trialselect','var')
        [R,Rcon]=nut_tfcov(meg,active,control,trialselect);
    else
        [R,Rcon]=nut_tfcov(meg,active,control);
    end

    % ~isempty needed because matlab thinks ([] | 1) = []
    if(~isempty(strmatch('SAMbbW',varargin,'exact')) | ~isempty(strmatch('SAMfd',varargin,'exact')))
        % note this is not for general use!!!! ******
        active2=[-250 250];
        control2=[-950 -450];
%         Rtmp=nut_tfcov(nuts.meg,[active2;control2]);
        Rtmp=nut_tfcov(meg,[active2;control2]);
        Rall = mean(Rtmp,3);
        clear Rtmp
    end
    % ************************************************
else
    % else we can get a covariance per time sample, and time window after
    % e.g., size(R) = 275 x 275 x timesamples
    error('not implemented yet.');
end
toc

filtERF = mean(meg.data,3);
params.megds = meg.filename;
% clear nuts
params.active = active;
params.control = control;
params.band = [lowfreq hifreq];
params.session = sessionfile;
for kk=1:size(R,3)
    if ~isempty(Rcon)
        params.mineig(kk)=min(eig([R(:,:,kk)+Rcon]/2));
    else
        params.mineig(kk)=min(eig([R(:,:,kk)]/2));
    end
end
params.filtereddata=filtsessionname;
% params.cn=1;


[filepath,filename,ext] = fileparts(sessionfile);
% if ~isempty(control)
%     outname = [filename '_' filtname '_' num2str(lowfreq) 'to' num2str(hifreq) 'Hz_con' num2str(control(1,1)) 'to' num2str(control(1,2)) 'ms.mat' ];
% else
if strcmp(paramfile((end-3):end),'.mat');
    outname = [filename '_' filtname '_' num2str(lowfreq) 'to' num2str(hifreq) 'Hz_' paramfile(1:(end-4)) '.mat' ];
else
    outname = [filename '_' filtname '_' num2str(lowfreq) 'to' num2str(hifreq) 'Hz_' paramfile(1:end) '.mat' ];
end
% end
% outname = [filename '_' filtname '_' num2str(lowfreq) 'to' num2str(hifreq) 'Hz_con'  num2str(control(1,2)) 'ms' ];
% save(outname,'R','Rcon','filtERF','params'); % don't need filtERF if got filtered raw data
save(outname,'R','Rcon','params');

