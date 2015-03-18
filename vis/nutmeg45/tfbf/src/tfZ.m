function tfZ(sessionfile,timewinfile,lowfreq,hifreq,varargin)
% TFZ   post-processes the time-frequency beamformer to give Wilcoxon Z scores
%
%   tfZ(sessionfile,timewinfile,lowfreq,highfreq,filterconfig,beamformerconfig)
%
% Input arguments in  are optional
%
% sessionfile       path to file containing NUTMEG session (nuts structure)
%                   OR (johanna's option) the file of save filtered data
% timewinfile       path to file containing timewindow settings (active and control)
% lowfreq           bandpass filter low freq cutoff in Hz
% highfreq          bandpass filter high freq cutoff in Hz
% filterconfig      varargin{1}: path to file containing "filt" structure with the fields:
%                   - class     'firls'/'butter'/'firpm'. Default is 'butter'.
%                   - order     filter order. Default is 4.
%                   - type      'bp' for bandpass filter/'low'/'high'. Default is 'bp'
%                   Alternatively, you can call your own filter:
%                   - filt.callback = 'nuts.meg.data = yourfilter(nuts.meg.data, ...
%                     lowfreq,hifreq,nuts.meg.srate);'
%                   - Output must be nuts.meg.data, available inputs are
%                     lowfreq, hifreq, nuts.meg.srate
%                   - This will not work if you run tfbf as a compiled function.
% beamformerconfig  varargin{2}: path to file containing "params" structure with the fields:
%                   - algo: cell with source loc algorithms. Default is {'SAM'}.
%                     Options include 'sLORETA' and 'LCMV'
%                   - dualstate: 1 (default) or 0. Set to 0 for single
%                     state beamformers (= no control/baseline time window)
%                   - cn: power normalization to avoid overestimating
%                     the power of deep sources. 1 (default) or 0.
%                   - qsub: send to qsub cluster. 0 (default) or 1.
%                   - savepower: 1 (default) or 0.
%                   - saveweights: 0 (default) or 1. Set to 1 if you would
%                     like to calculate virtual channels.
%                   - savefiltereddata: 0 (default) or 1. Set to 1 if you
%                     would like to calculate virtual channels.
% varargin{3}       ask Johanna


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
if ~isfield(params,'algo'), params.algo={'SAM'}; end
if ischar(params.algo), params.algo={params.algo}; end
if ~isfield(params,'dualstate'), params.dualstate=true; end
if ~isfield(params,'cn'), params.cn=true; end
if ~isfield(params,'qsub'), params.qsub=false; end
if ~isfield(params,'savepower'), params.savepower=true; end
if ~isfield(params,'saveweights'), params.saveweights=false; end
if ~isfield(params,'savefiltereddata'), params.savefiltereddata=false; end
if ~isfield(params,'doSumSquares'), params.doSumSquares=false; end
if ~isfield(params,'noiseproject'), params.noiseproject=false; end
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

[filepath,filename,ext] = fileparts(sessionfile);

load(timewinfile);  % must contain time windows ("active" and "control")

% print run parameters (for preservation in log files)
disp(sessionfile)
disp(timewinfile)
filt
params

if ~params.dualstate && ~isempty(control)
    warning('NUTMEG:ControlWindowsDefinedForSingleStateBeamformer','Ignoring control windows for single state beamformer!')
    control=[]; 
end

load(sessionfile); % must contain nuts.voxels, nuts.voxelsize, nuts.coreg, nuts.meg
if exist('nuts')
    nuts.Lp = [];  % kill it to save memory... don't need it here
    
    
    if isfield(nuts.meg,'goodchannels')
        nuts.meg.data=nuts.meg.data(:,nuts.meg.goodchannels,:);
    end
    numchannels=size(nuts.meg.data,2)
    
    if isfield(filt,'callback')
        if ~strcmp(filt.callback(1:13),'nuts.meg.data'), error('Invalid filter callback specified.'), end
        eval(filt.callback);
        %e.g., filt.callback = 'nuts.meg.data = yourfilt(nuts.meg.data,nuts.meg.srate,lowfreq,hifreq);'
    else
        nuts.meg.data = nut_filter2(nuts.meg.data,filt.class,filt.type,filt.order,lowfreq,hifreq,nuts.meg.srate,1);
    end
    meg=nuts.meg;
elseif exist('meg')
    % just a sanity check
    disp('you loaded previously filtered MEG data');
else
    error('your session file is wrong');    
end

% if params.savefiltereddata
%     activetimewins =  find( (nuts.meg.latency > active(1,1)-1) & (nuts.meg.latency < active(end,end)+1) );
%     F.data = nuts.meg.data(activetimewins,:,:);
%     F.latency= nuts.meg.latency;
%     F.srate=nuts.meg.srate;
%     outname = ['filtdata_' filename '_' filtname '_' num2str(lowfreq) 'to' num2str(hifreq) 'Hz'];
%     fprintf('Saving filtered data as %s...\n',outname)
%     save(outname,'F');
%     clear F activetimewins outname
% end


%% new Wilcoxon stuff
% params.megds = nuts.meg.filename;
% params.active = active;
% params.control = control;
% params.band = [lowfreq hifreq];
% params.session = sessionfile;

if nargin>6
    % this option works for johanna (let me know if it messes you up!)
    [twfilepath,twfilename,twext] = fileparts(timewinfile);
    covfile = ['cov_' filename '_' twfilename];
    covfile_usechar=[5:24 29:length(covfile)];
    load(varargin{3});
else
    covfile = ['cov_' filename '_' filtname '_' paramname '_' num2str(lowfreq) 'to' num2str(hifreq) 'Hz'];
    covfile_usechar=[5:length(covfile)];
end


for timeselect = 1:size(active,1);
    load(['W_' covfile(covfile_usechar) '_' num2str(active(timeselect,1)) 'to' num2str(active(timeselect,2)) 'ms_' params.algo{1} '.mat'],'W*');
    
    if exist('at') %from markers file
        [Pa,Pc] = deal(zeros(size(W,2),length(find(at))));
        fat=find(at);
        fct=find(ct);
        for trial = 1:min(length(fat),length(fct)) % loop over trials
            [R] = nut_tfcov(meg,active(timeselect,:),[],fat(trial));
            [Rc] = nut_tfcov(meg,active(timeselect,:),[],fct(trial));
            Pa(:,trial)  = sum(W.*(R*W))';
            Pc(:,trial) = sum(W.*(Rc*W))';
        end
    else %this is default for almost everyone
        [Pa,Pc] = deal(zeros(size(W,2),size(meg.data,3)));
        for trial = 1:size(meg.data,3) % loop over trials
            [R,Rc] = nut_tfcov(meg,active(timeselect,:),control,trial);
            Pa(:,trial)  = sum(W.*(R*W))';
            if ~isempty(Rc)
                Pc(:,trial) = sum(W.*(Rc*W))';
            elseif params.noiseproject % no by default
                [u,s,v]=svd(R);
                sig=s(end,end);
                Pc(:,trial) = sum(W.*(sig*W))'; % poor man's control
            end
        end
    end
    
    beam=load(['s_beamtf_' covfile(covfile_usechar) '_' num2str(active(timeselect,1)) 'to' num2str(active(timeselect,2)) 'ms_' params.algo{1} '.mat']);
    if isfield(beam,'beam')
        beam=beam.beam;
    end
    if isfield(beam,'z')
        beam=rmfield(beam,'z');
    end
    if isfield(beam,'sav')
        beam=rmfield(beam,'sav');
        beam=rmfield(beam,'scv');
    end
    [beam.z{1},beam.p{1}] = deal(zeros(size(W,2),1));
    
    
    for ii = 1:size(Pa,1) % loop over voxels
        [beam.p{1}(ii),h,stats] = nut_signrank(Pa(ii,:),Pc(ii,:));
        beam.z{1}(ii) = stats.zval;
    end
    
    if isfield(params,'keeptrials')
        if params.keeptrials
            beam.trial{1}=Pa;
            if max(Pc(:))>0
                beam.trial{2}=Pc;
            end
        end
        if isfield(params,'keepzp') % yes, I know the point of this function is Z, but for memory reasons maybe not want
            if ~params.keepzp
                beam=rmfield(beam,'z');
                beam=rmfield(beam,'p');
            end
        end
    end
        
    if params.doStdDev
        if ~isfield(beam,'sa')
            beam=nut_oldbeam2newbeam(beam);
        end
        beam.sasd{1}=zeros(size(beam.sa{1}));
        beam.scsd{1}=zeros(size(beam.sc{1}));
        for jj=1:size(Pa,2) %loop over trials
            beam.sasd{1}=beam.sasd{1}+(beam.sa{1}-Pa(:,jj)).^2;
            beam.scsd{1}=beam.scsd{1}+(beam.sc{1}-Pc(:,jj)).^2;
        end
        beam.sasd{1}=sqrt(beam.sasd{1}/size(Pa,2));
        beam.scsd{1}=sqrt(beam.scsd{1}/size(Pa,2));
        beam.numtrial=size(Pa,2);
    end
%    timeselect
%     save(['s_beamtf_' covfile(covfile_usechar) '_' num2str(active(timeselect,1)) 'to' num2str(active(timeselect,2)) 'ms_' params.algo{1} '.mat'],'-append','beam');
    save(['s_beamtf_' covfile(covfile_usechar) '_' num2str(active(timeselect,1)) 'to' num2str(active(timeselect,2)) 'ms_' params.algo{1} '.mat'],'-struct','beam');
end
%%%%% save nut_timef2tfbf??????

clear R Rc


%%


% 
% if(params.qsub)
%     if(length(getenv('BQSCLUSTER')) > 0)
%         clustertype = 'bqs';
%     elseif(length(getenv('SGE_TASK_ID')) > 0)
%         clustertype = 'sge';
%     elseif(length(getenv('XGRID')) > 0)
%         clustertype = 'xgrid';
%     else
%         warning('unknown cluster type... running TFBF serially...');
%         clustertype = '1cpu';
%     end
% else
%     clustertype = '1cpu';
% end
% 
% switch(clustertype)
%     case 'sge' % e.g., BIL and QB3 clusters
%         [crap,whereami] = unix('id -gn');
%         if(strcmp(whereami(1:3),'bil'))
%             pristring = '-p -500 ';
%         else
%             pristring = '';
%         end
% 
%         for jj=1:size(params.algo,1)
%             unix(['qsub ' pristring '-t 1-' int2str(size(R,3)) ' ~/bin/qtfrun.csh ' sessionfile ' ' outname ' ' params.algo{jj}]);
%         end
%     case 'bqs' % e.g., French IN2P3 cluster, some TeraGrid sites
%         for jj=1:size(params.algo,1)
%             for ii=1:size(R,3)
%                 unix(['setenv QSUB_PWD ' pwd '; setenv SGE_TASK_ID ' int2str(ii) '; setenv ARGV "' sessionfile ' ' outname ' ' params.algo{jj} '";  qsub -v QSUB_PWD,SGE_TASK_ID,ARGV ~/bin/qtfrun.csh ']);
%             end
%         end
%     case 'xgrid' % Mac cluster
%         % NOT READY YET!!!!
%         for jj=1:size(params.algo,1)
%             unix(['xgrid ' pristring '-t 1-' int2str(size(R,3)) ' ~/bin/qtfrun.csh ' sessionfile ' ' outname ' ' params.algo{jj}]);
%         end
%     case '1cpu' % no cluster, just your lonely computer
%         tic
%         for jj=1:size(params.algo,1)
%             for ii=1:size(params.active,1)
%                 nut_tfsrc(sessionfile,outname,ii,params.algo{jj});
%                 disp(['Finished running ' num2str(ii) ' out of ' num2str(size(params.active,1)) ' on ' outname ' using ' params.algo{jj} '.' ]);
%                 toc
%             end
%         end
% end

