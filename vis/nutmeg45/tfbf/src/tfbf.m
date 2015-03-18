function tfbf(sessionfile,timewinfile,lowfreq,hifreq,varargin)
% TFBF   runs the time-frequency beamformer
%
%   tfbf(sessionfile,timewinfile,lowfreq,highfreq,�filterconfig�,�invsolconfig�)
%
% Input arguments in �� are optional
%
% sessionfile       path to file containing NUTMEG session (nuts structure)
% timewinfile       path to file containing timewindow settings (active and control)
% lowfreq           bandpass filter low freq cutoff in Hz
% highfreq          bandpass filter high freq cutoff in Hz
% filterconfig      path to file containing "filt" structure with the fields:
%                   - class     'firls'/'butter'/'firpm'. Default is 'butter'.
%                   - order     filter order. Default is 4.
%                   - type      'bp' for bandpass filter/'low'/'high'. Default is 'bp'
%                   - reference EEG data can be rereferenced to the average
%                               ('AVG') or to a single ref electrode ('eleclabel'). 
%                               Default is no rereferencing.
%                   Alternatively, you can call your own filter:
%                   - filt.callback = 'nuts.meg.data = yourfilter(nuts.meg.data, ...
%                     lowfreq,hifreq,nuts.meg.srate);'
%                   - Output must be nuts.meg.data, available inputs are
%                     lowfreq, hifreq, nuts.meg.srate
%                   - This will not work if you run tfbf as a compiled function.
% invsolconfig      path to file containing "params" structure with the fields:
%                   - algo: cell with source loc algorithms. Default is {'SAM'}.
%                     Options include, among others, 'sLORETA' and 'MinNorm'
%                   - cn: power normalization to avoid overestimating
%                     the power of deep sources. 1 (default) or 0.
%                   - regularization: algorithm used for matrix inversion.
%                     'none'=matlab inv function (default); options include
%                     'tikhonov', 'pinv', 'bayesian'.
%                   - regulthres: threshold for regularization algorithm.
%                     'always' (default) or 'auto'=inversion error >1/10 of absolute 
%                     data minimum.
%                   - qsub: send to qsub cluster. 0 (default) or 1.
%                   - savepower: 1 (default) or 0.
%                   - saveweights: 0 (default) or 1. Set to 1 if you would
%                     like to calculate virtual channels.
%                   - savefiltereddata: 0 (default) or 1. Set to 1 if you
%                     would like to calculate virtual channels.

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
else
    load(varargin{1}); % load filter settings
    [filtpath,filtname,ext] = fileparts(varargin{1});    
end
if ~exist('params','var')       % if filt and params were not saved in same single file
    if nargin<6
        params=[];
        paramname='dsSAMcn';
    else
        load(varargin{2}); % load parameter file
        [parampath,paramname,ext] = fileparts(varargin{2});
    %     algo=cell(length(varargin)-1,1);
    %     [algo{1:size(algo,1)}]=varargin{2:end};
    end
end

nut_defaults;

% Load everything
[filepath,filename,ext] = fileparts(sessionfile);

nuts=load(sessionfile);
isoldsavestyle = isfield(nuts,'nuts');
if isoldsavestyle
    nuts=nuts.nuts;
end
nuts.Lp = [];  % kill it to save memory... don't need it here

load(timewinfile);  % must contain time windows ("active" and "control")
if isfield(params,'dualstate') && ~params.dualstate && ~isempty(control)
    warning('NUTMEG:ControlWindowsDefinedForSingleStateBeamformer','Ignoring control windows for single state beamformer!')
    control=[]; 
end

% Default settings
if ~isfield(params,'algo'), params.algo={'SAM'}; end
if ischar(params.algo), params.algo={params.algo}; end
params.dualstate = ~isempty(control);
params.eegflag = ( isfield(nuts.meg,'eegflag') && nuts.meg.eegflag );
if ~isfield(params,'cn'), params.cn = ~params.eegflag; end  % weight normalization is default for EEG data, leadfield normalization for MEG data
if ~isfield(params,'wn'), params.wn = (~params.cn && params.eegflag); end
if ~isfield(params,'qsub'), params.qsub=false; end
if ~isfield(params,'savepower'), params.savepower=true; end
if ~isfield(params,'saveweights'), params.saveweights=false; end
if ~isfield(params,'savefiltereddata'), params.savefiltereddata=false; end
if ~isfield(params,'regularization')
    if params.eegflag, params.regularization='pinv'; 
    else params.regularization='none'; 
    end
end
if ~isfield(params,'regulthres'), params.regulthres='always'; end
if ~exist('paramname','var')
    if params.dualstate && params.cn
        paramname=['ds' params.algo{1} 'cn'];
    elseif params.dualstate && ~params.cn
        paramname=['ds' params.algo{1}];
    elseif ~params.dualstate && params.cn
        paramname=[params.algo{1} 'cn'];
    else
        paramname=[params.algo{1}];
    end 
end

% Remove bad channels if necessary
if isfield(nuts.meg,'goodchannels')
    if size(nuts.meg.data,2)>length(nuts.meg.goodchannels)
        nuts.meg.data=nuts.meg.data(:,nuts.meg.goodchannels,:);
    end
    goodchannels=nuts.meg.goodchannels;
else
    goodchannels=1:size(nuts.meg.data,2);
end
numchannels=size(nuts.meg.data,2)

% Deal with EEG specific problems
if params.eegflag
    params.eegflag=true;
    if isfield(nuts.meg,'referenceidx')
        params.referenceidx=nuts.meg.referenceidx;
    else  % if reference not defined we assume average reference
        params.referenceidx=goodchannels;
    end
	if isscalar(params.referenceidx)  % If single electrode reference
        nuts.meg.data(:,goodchannels==params.referenceidx,:)=[];  % Remove reference channel (which is all 0)
    else        % If average reference
        if ~isfield(params,'do_avref_leadpotential')        % if not defined in beamformer config file
            params.do_avref_leadpotential = false;      % In my experience, average referencing the leadpotential can lead to problems and does not bring any advantage (Adrian). But this may be different in other settings.
        end
    end
end

% print run parameters (for preservation in log files)
disp(sessionfile)
disp(timewinfile)
filt
params

% Filter and preprocessing
outname = ['filtdata_' filename '_' filtname '_' num2str(lowfreq) 'to' num2str(hifreq) 'Hz.mat'];
if ~exist(outname,'file')
    if isfield(filt,'reference')        % EEG data can be rereferenced to AVG or given electrode here.
        nuts.meg = nut_eegref(nuts.meg,filt.reference);
    end
    if isfield(filt,'callback')
        if ~strcmp(filt.callback(1:13),'nuts.meg.data'), error('Invalid filter callback specified.'), end
        eval(filt.callback);
        %e.g., filt.callback = 'nuts.meg.data = yourfilt(nuts.meg.data,nuts.meg.srate,lowfreq,hifreq);'
    else
        nuts.meg.data = nut_filter2(nuts.meg.data,filt.class,filt.type,filt.order,lowfreq,hifreq,nuts.meg.srate,1);
    end
    if params.savefiltereddata
        fprintf('Saving filtered data as %s...\n',outname)
        meg=nuts.meg;
        save(outname,'meg');
        clear meg
    end
else
    warning(['Loading filtered data ' outname ' saved earlier!'])
    load(outname);
    if exist('F','var') % old way
        activetimewins = ( (nuts.meg.latency > active(1,1)-1) & (nuts.meg.latency < active(end,end)+1) );
        nuts.meg.data(activetimewins,:,:) = F.data;
        clear F
    elseif exist('meg','var')
        nuts.meg=meg;
        clear meg
    end
end

if isfield(nuts.meg,'numtrial') % probably only affects Johanna
    error('nuts.meg.numtrial is deprecated. please use meg.markers')
end

% Covariance
oldway = true;
if(oldway)
    [R,Rcon]=nut_tfcov(nuts.meg,active,control);
    %Rcon = R(:,:,end);
    %R(:,:,end) = [];
else
    % else we can get a covariance per time sample, and time window after
    % e.g., size(R) = 275 x 275 x timesamples
    error('not implemented yet.');
end

%filtERF = mean(nuts.meg.data,3);
params.megds = nuts.meg.filename;
% clear nuts
params.active = active;
params.control = control;
params.band = [lowfreq hifreq];
params.session = sessionfile;
for kk=1:size(R,3) % Johanna uses this mineig for inverse regularization
    % FIXME: in future, use 2nd output of nut_cov rather than recompute mineig here
    if ~isempty(Rcon)
        params.mineig(kk)=min(eig([R(:,:,kk)+Rcon]/2));
    else
        params.mineig(kk)=min(eig([R(:,:,kk)]/2));
    end
end

outname = ['cov_' filename '_' filtname '_' paramname '_' num2str(lowfreq) 'to' num2str(hifreq) 'Hz'];
save([outname '.mat'],'R','Rcon','params','goodchannels');   % ,'filtERF'

if isoldsavestyle
% be sure you called nut_liposession(filename) already prior to this step
    sessionfile = [filename 'Lp.mat'];  % load stripped session file
end

% Inverse Solution

% switch(params.algo{1}) % might need to change this if multiple algos specified...
%     case 'sLORETA'
%         data.y = nuts.meg.data;
%         [W] = nut_sLORETA(nuts.Lp,data);
%         save(['W_sLORETA_' outname],'W');
%     otherwise
%         disp('onk.')
% end

if(params.qsub)
    if(length(getenv('BQSCLUSTER')) > 0)
        clustertype = 'bqs'
    elseif(length(getenv('SGE_TASK_ID')) > 0)
        clustertype = 'sge'
    elseif(length(getenv('XGRID')) > 0)
        clustertype = 'xgrid'
    else
        warning('unknown cluster type... running TFBF serially...');
        clustertype = '1cpu'
    end
else
    clustertype = '1cpu';
end

switch(clustertype)
    case 'sge' % e.g., BIL and QB3 clusters
        [crap,whereami] = unix('id -gn');
        if(strcmp(whereami(1:3),'bil'))
            pristring = '-p -500 ';
        else
            pristring = '';
        end

        for jj=1:length(params.algo)
            unix(['qsub ' pristring '-t 1-' int2str(size(R,3)) ' ~/bin/qtfrun.csh ' sessionfile ' ' outname ' ' params.algo{jj}]);
        end
    case 'bqs' % e.g., French IN2P3 cluster, some TeraGrid sites
        for jj=1:length(params.algo)
            for ii=1:size(R,3)
                unix(['setenv QSUB_PWD ' pwd '; setenv SGE_TASK_ID ' int2str(ii) '; setenv ARGV "' sessionfile ' ' outname ' ' params.algo{jj} '";  qsub -e $QSUB_PWD/qlogs -o $QSUB_PWD/qlogs -v QSUB_PWD,SGE_TASK_ID,ARGV ~/bin/qtfrun.csh ']);
            end
        end
    case 'xgrid' % Mac cluster
        % NOT READY YET!!!!
        for jj=1:length(params.algo)
            unix(['xgrid ' pristring '-t 1-' int2str(size(R,3)) ' ~/bin/qtfrun.csh ' sessionfile ' ' outname ' ' params.algo{jj}]);
        end
    case '1cpu' % no cluster, just your lonely computer
        tic
        for jj=1:length(params.algo)
            for ii=1:size(params.active,1)
                nut_tfsrc(sessionfile,outname,ii,params.algo{jj});
                disp(['Finished running ' num2str(ii) ' out of ' num2str(size(params.active,1)) ' on ' outname ' using ' params.algo{jj} '.' ]);
                toc
            end
        end
end

