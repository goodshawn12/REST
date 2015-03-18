function nut_calltfrun(sessionfile,paramfile,lowfreq,hifreq,varargin)
% tfcov(sessionfile,paramfile,lowfreq,highfreq)
% wrapper script for nut_tfcov, intended for compilation into standalone
% varargin{1} is .mat file with filt. and params. stuff set in it
% varargin{2} is covfile_usechar (which characters to use from the covfile name)
% varargin{3} is covwindowname

tic
global ndefaults
nut_defaults;

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
    algo={'SAM'};
    params.cn=1
    params.qsub=0;
else
    load(varargin{1}); % load additional params (like filter settings)
    [filtpath,filtname,ext] = fileparts(varargin{1});
%     algo=cell(length(varargin)-1,1);
%     [algo{1:size(algo,1)}]=varargin{2:end};
algo={'SAM'};
end
if nargin>5
    covfile_usechar=varargin{2};
    if isempty(covfile_usechar)
        clear covfile_usechar
    end
end
if nargin>6
    windowfile=varargin{3};
end

% load(sessionfile); % must contain nuts.Lp, nuts.voxels, nuts.voxelsize, nuts.coreg, nuts.meg
load(paramfile);  % must contain time windows ("active" and "control")

% print run parameters (for preservation in log files)
% disp(sessionfile)
% disp(paramfile)

[filepath,filename,ext] = fileparts(sessionfile);
% be sure you called nut_liposession(filename) already prior to this step
sessionLp = [filename 'Lp.mat'];  % load stripped session file
if exist('windowfile')
    rcovname=[filename '_' filtname '_' num2str(lowfreq) 'to' num2str(hifreq) 'Hz_' windowfile '.mat' ];
else
    rcovname=[filename '_' filtname '_' num2str(lowfreq) 'to' num2str(hifreq) 'Hz_con' num2str(control(1,1)) 'to' num2str(control(1,2)) 'ms.mat' ];
end
paramstmp=params;
load(rcovname);
params.cn=paramstmp.cn;
params.qsub=paramstmp.qsub;
params.dualstate=paramstmp.dualstate;
params.algo=paramstmp.algo;
params.savepower=paramstmp.savepower;
params.saveweights=paramstmp.saveweights;
if params.saveweights
    if isfield(params,'weightsvoxels');
        if isempty(params.weightsvoxels)
            params=rmfield(params,'weightsvoxels');
        end
    end
end
save(rcovname,'R','Rcon','params');
if exist('windowfile')
    outname=[filename '_' filtname '_' num2str(lowfreq) 'to' num2str(hifreq) 'Hz_' windowfile];
else
    outname=[filename '_' filtname '_' num2str(lowfreq) 'to' num2str(hifreq) 'Hz_con' num2str(control(1,2)) 'ms' ];
end

 if(strmatch('SAMfd',varargin{1},'exact'))
    [W,eta] = nut_SAM(nuts.Lp,Rall,0);
    save(['W_SAMfd_' outname],'W');
elseif(strmatch('SAMbbW',varargin{1},'exact'))
    [W,eta] = nut_SAM(nuts.Lp,Rall,0);
    save(['W_SAMbb_' sessionLp],'W');
    return
end

if params.qsub
    % set priority for local cluster, leave empty everywhere else
    [crap,whereami] = unix('id -gn');
    if(strcmp(whereami(1:3),'bil'))
        pristring = '-p -500 ';
    else
        pristring = '';
    end


    for jj=1:size(algo,1)
        unix(['qsub ' pristring '-t 1-' num2str(size(R,3)) ' /data/research_meg/tfbf/bin/qtfrun.csh ' sessionfile ' ' outname ' ' algo{jj}]);
    end
else
    tic
    for ii=1:size(params.active,1)
        if exist('covfile_usechar')
%             tfrun(sessionLp,rcovname,ii,algo{1},covfile_usechar);
            tfrun(filename,rcovname,ii,algo{1},covfile_usechar);
        elseif exist('windowfile')
%             tfrun(sessionLp,rcovname,ii,algo{1});
            tfrun(filename,rcovname,ii,algo{1});
        else
%             tfrun(sessionLp,outname,ii,algo{1});
            tfrun(filename,outname,ii,algo{1});
        end
        disp(['finished running ' num2str(ii) ' out of ' num2str(length(params.active)) ' on ' outname ]);
        toc
    end
end

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [W,eta] = nut_SAM(Lp, Rall, cn)

InvRall = inv(Rall);

flags.LCMVcn = cn;
flags.wn = false;
flags.dualstate = true;
flags

[W, eta] = nut_LCMV_Beamformer(Lp,InvRall, flags);
