%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REAL-TIME BCILAB-LSL DEMONSTRATION SCRIPT
%
% This script demonstrates how to set up a default BCILAB-LSL pipeline for
% real-time data extraction and processing from an EEG device or file
%
% This script requires the following tools to be installed:
% BCILAB: https://www.sccn.ucsd.edu/wiki/BCILAB
% LSL:    https://code.google.com/p/labstreaminglayer/
%
% Optional (if using SourceLocalization)
% MobiLab: http://sccn.ucsd.edu/wiki/MoBILAB
%
% Author: Tim Mullen,  October 14, 2013, SCCN/INC/UCSD
% Email:  tmullen@ucsd.edu
%
% .........................................................................
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
% .........................................................................
%
% If you use this script as a basis for your published work, kindly
% consider acknowledging us in your manuscript. A possible example:
% "Thanks to Tim Mullen and Christian Kothe (SCCN, UCSD) for providing
%  scripting support"
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Force EEGLAB to use double precision
pop_editoptions( 'option_single', false);

%% Define config options (USER MUST DEFINE THIS)
% -------------------------------------------------------------------------
opts.runlsl = false;

% Set up the name of the LSL stream
opts.lsl.StreamName          = 'EEGDATA';
opts.lsl.SelectionProperty   = 'name';  % can also be 'name'
opts.lsl.SelectionValue      = 'Emotiv EPOC';

% Determine a fixed 'sliding' window length (seconds) for pulling data from LSL
% If winlen=0, pull the most recent (variable length) chunk of unprocessed data
opts.winlen = 0;

% suppress console output in processing loop?
% if false, this can reusult in a LOT of text printed to console
opts.silence = true;

% display benchmarking
opts.dispBenchmark = true;

% select calibration data segment
opts.calibEpoch = [0 1];

% Establish file paths
% ..........................................................................
% User Tip:
%
% All paths are relative to opts.datapath. This is platform-independent and
% itself can be relative to bcilab root folder (i.e. data:/ is the userdata
% folder in the bcilab root dir)
% ..........................................................................
opts.datapath = 'D:\Matlab Coding\OnlineICA\P3_realtime\data\';
% (required) calibration data set
opts.TrainingDataFile = 'EmotivTrain_EyeClose_icainfo.set';%'Colin27_Biosemi_1010_0_norot.set'; %'Flanker_SH_cleaned_calib_testicainfo.set';%'sim_64ch_superlaplaceblink1sub7_800s.set';
% (optional) file to "playback" if opts.runlsl = false
opts.PlaybackDataFile = 'EmotivTrain_EyeClose_icainfo.set';%'EmotivTrain_EyeClose_icainfo.set'; %'Flanker_SH_cleaned_test_icainfo.set';%'sim_64ch_superlaplaceblink1sub7_800s.set';
% (optional) pipeline config file
opts.BCILAB_PipelineConfigFile = 'ORICA_pipeline_config_realtime.mat'; % make sure this file doesn't have 'signal' entry

[filename, pathname] = uigetfile('*.set');
if filename ~= 0
    opts.TrainingDataFile = filename;
    opts.PlaybackDataFile = filename;
    opts.datapath         = pathname;
end


%% Load head model object
% -------------------------------------------------------------------------
% opts.HeadModelDataFile= ''; %'NameOfHeadModelFile.mat';   % Look in <SIFT_RootFolder>\resources\headmodels\standard-Colin27-385ch.mat
% if ~isempty(opts.HeadModelDataFile)
%     hmObj = headModel.loadFromFile(env_translatepath([opts.datapath opts.HeadModelDataFile]));
% end

%% load calibration recording
% -------------------------------------------------------------------------
calibData = exp_eval_optimized(io_loadset([opts.datapath opts.TrainingDataFile], ...
    'markerchannel',{'remove_eventchns',false}));
% select data
% calibData = set_selinterval(calibData, opts.calibEpoch);
calibData = pop_select(calibData,'time',opts.calibEpoch);
flt_pipeline('update');

%% Set up the pipeline
% -------------------------------------------------------------------------
% load existing config file(s) (if it exists)
try    fltPipCfg = exp_eval(io_load([opts.datapath opts.BCILAB_PipelineConfigFile]));
catch, disp('-- no existing pipeline --'); fltPipCfg = {}; end

% enforce usage of user-defined head model
% if ~isempty(opts.HeadModelDataFile)
%     fltPipCfg.psourceLocalize.hmObj = [opts.datapath opts.HeadModelDataFile]; end

% display GUI to allow user to configure the BCILAB pipeline...
% if an error appears here, delete config file and re-run the script
fltPipCfg = arg_guipanel('Function',@flt_pipeline, ...
    'Parameters',[{'signal',calibData} fltPipCfg], ...
    'PanelOnly',false);

% ..........................................................................
% User Tip:
%
% You might want to save the pipeline config (fltPipCfg) at this point for
% later reload (from 'opts.BCILAB_PipelineConfigFile')
if ~isempty(fltPipCfg)
    if isfield(fltPipCfg,'signal'); fltPipCfg = rmfield(fltPipCfg,'signal'); end
    save(env_translatepath([opts.datapath opts.BCILAB_PipelineConfigFile]),...
        '-struct','fltPipCfg');
end
% ..........................................................................

% ...apply the pipeline to calibration data
disp('-- Applying pipeline to calibration data (please wait) --');
cleaned_data = exp_eval(flt_pipeline('signal',calibData,fltPipCfg));

% ..........................................................................
% User Tip:
%
% Instead of using the GUI, you can also easily define pipelines in the
% script itself. Here are some examples:
%
% % remove flatline channels:
% noflats      = flt_clean_flatlines(flt_seltypes(calibData,'EEG'));
% % apply a high-pass filter:
% highpassed   = flt_fir(noflats,[0.5 1],'highpass');
% % discard noisy channels:
% goodchannels = flt_clean_channels('signal',highpassed,'min_corr',0.35,'ignored_quantile',0.3);
%
% % you always need to make sure to apply the pipeline
% cleaned_data = exp_eval(goodchannels);
% ..........................................................................
if isfield(calibData,'icawinv_true')
    ground_truth_ica = exp_eval(flt_project(calibData,inv(calibData.icawinv_true)));
end

% ..........................................................................
% User Tip:
%
% After this step, cleaned_data will contain the results of having applied
% the entire pipeline to your calibration data. Thus, the steps up to this
% point essentially define an offline data processing pipeline which you
% can use for your publications. It is also helpful at this stage to
% examine the contents of cleaned_data for deficiencies, which can help you
% fine-tune your pipeline for later online use.
% ..........................................................................

%% Start LSL streaming
% -------------------------------------------------------------------------
if opts.runlsl == true
    run_readlsl('MatlabStream',opts.lsl.StreamName,                 ...
        'SelectionProperty',opts.lsl.SelectionProperty,     ...
        'SelectionValue',opts.lsl.SelectionValue);
else
    % ... OR read from datafile
    run_readdataset('MatlabStream',opts.lsl.StreamName,                         ...
        'Dataset',io_loadset([opts.datapath opts.PlaybackDataFile], ...
        'markerchannel',{'remove_eventchns',false}));
end

% initialize the pipeline for streaming data
pipeline     = onl_newpipeline(cleaned_data,opts.lsl.StreamName);

% render a panel for viewing streams
gui_vis_filtered;

%% Main loop
% -------------------------------------------------------------------------
disp('-- Running pipeline --');
chunk_len = round(opts.winlen*cleaned_data.srate);
% initialize parameters to store
% testPI = []; ConnData = []; convEst = []; statIdx = []; meanIdx = []; covIdx = []; alarm = [];
% icaweights_all = []; sphere_all = []; time_all = []; lambda_k = [];
% sampleSize = []; executeTime = []; sampleSizeWhite = []; executeTimeWhite = [];
updateTime = 0; 
plotcomponent = 1; plotconvergence = 0; ploticaact = 0; plotstatIdx = 0; plotmir = 0;
flag = 1;

while true
    
    tic;
    % grab a chunk of data from the stream and send through the filter pipeline
    [eeg_chunk,pipeline] = onl_filtered(pipeline, chunk_len, opts.silence);
    benchmarking.pipeline = toc;
    
    % .....................................................................
    % User Tip:
    %
    % The structure 'eeg_chunk' now contains the results of applying the
    % pipeline to the chunk of raw data acquired from the LSL stream.
    % This is an augmented EEGLAB data structure.
    %
    % - You can add code here to apply your own additional processing steps
    % - You'll also want to add code to store the results you'd like keep
    %
    % Hint:
    %  eeg_chunk.data contains processed/cleaned data
    %
    %  (if using source localization)
    %  eeg_chunk.srcpot contains current density for ROIs
    %
    %  (if using source localization and keepFullCsd is enabled)
    %  eeg_chunk.srcpot_all contains current density for all mesh vertices
    %
    %  (if using ICA)
    %  eeg_chunk.icaact contains ica activations
    %
    % (if using SIFT)
    % eeg_chunk.CAT  contains a SIFT object/datastructure
    %  CAT
    %   .Conn       Connectivity object
    %       .winCenterTimes
    %       .freqs
    %       .dims   dimensions for connectivity matrix
    %               e.g. (chans_to, chans_from, freqs, time)
    %       .<measure_name>    Connectivity matrix (conforms to dims above)
    %   .MODEL      Model structure
    %
    % .....................................................................
    
    
    if opts.dispBenchmark
        % TODO: implement this as a figure
        vis_benchmark(benchmarking);
        
        if plotconvergence
            % plot ORICA convergence
            if fltPipCfg.porica.arg_selection && isfield(eeg_chunk,'convergence')
                if ~exist('bpfig','var') || isempty(bpfig) || ~ishandle(bpfig)
                    CONVERGE_BUF = 60*eeg_chunk.srate; % size of convergence buffer
                    bpfig = figure;
                    bpaxes= axes('parent',bpfig);
                    convergence = nan(1,CONVERGE_BUF);
                    convergence(end-length(eeg_chunk.convergence)+1:end) = 10*log(eeg_chunk.convergence);
                    h=plot(linspace(0,CONVERGE_BUF/eeg_chunk.srate,CONVERGE_BUF),convergence);
                    set(bpfig,'userdata',h);
                    ylabel('10*ln(Convergence)');
                    xlabel('Time');
                else
                    blkpnts = length(eeg_chunk.convergence); %options.blockSize*(skipCount);
                    convergence(:,1:end-blkpnts) = convergence(:,blkpnts+1:end);
                    % insert new samples ...
                    convergence(:,end-blkpnts+1:end) = 10*log(eeg_chunk.convergence);
                    set(h,'ydata',convergence);
                end
            end
        end
        
        if plotstatIdx
            % plot ORICA convergence
            if fltPipCfg.porica.arg_selection && isfield(eeg_chunk,'statIdx')
                if ~exist('sifig','var') || isempty(sifig) || ~ishandle(sifig)
                    CONVERGE_BUF_si = 60*eeg_chunk.srate; % size of convergence buffer
                    sifig = figure;
                    siaxes= axes('parent',sifig);
                    statIdx = nan(1,CONVERGE_BUF_si);
                    statIdx(end-length(eeg_chunk.statIdx)+1:end) = eeg_chunk.statIdx;
                    hsi=plot(linspace(0,CONVERGE_BUF_si/eeg_chunk.srate,CONVERGE_BUF_si),statIdx);
                    set(sifig,'userdata',hsi);
                    ylabel('Stationary Idx');
                    xlabel('Time');
                else
                    blkpnts = length(eeg_chunk.statIdx); %options.blockSize*(skipCount);
                    statIdx(:,1:end-blkpnts) = statIdx(:,blkpnts+1:end);
                    % insert new samples ...
                    statIdx(:,end-blkpnts+1:end) = eeg_chunk.statIdx;
                    set(hsi,'ydata',statIdx);
                end
            end
        end
        
        
        if plotmir
            % plot ORICA convergence
            if fltPipCfg.porica.arg_selection && isfield(eeg_chunk,'mir')
                if ~exist('mirfig','var') || isempty(mirfig) || ~ishandle(mirfig)
                    CONVERGE_BUF_mir = 60*eeg_chunk.srate; % size of convergence buffer
                    mirfig = figure;
                    miraxes= axes('parent',mirfig);
                    mir = nan(1,CONVERGE_BUF_mir);
                    mir(end-length(eeg_chunk.mir)+1:end) = eeg_chunk.mir;
                    hmir=plot(linspace(0,CONVERGE_BUF_mir/eeg_chunk.srate,CONVERGE_BUF_mir),mir);
                    set(mirfig,'userdata',hmir);
                    ylabel('MIR');
                    xlabel('Time');
                else
                    blkpnts = length(eeg_chunk.mir); %options.blockSize*(skipCount);
                    mir(:,1:end-blkpnts) = mir(:,blkpnts+1:end);
                    % insert new samples ...
                    mir(:,end-blkpnts+1:end) = eeg_chunk.mir;
                    set(hmir,'ydata',mir);
                end
            end
        end
        
        
        if ploticaact
            % plot ORICA ica activation
            if fltPipCfg.porica.arg_selection && isfield(eeg_chunk,'icaact')
                if ~isempty(eeg_chunk.icaact)
                    if ~exist('srcfig','var') || isempty(srcfig) || ~ishandle(srcfig)
                        SRCPOT_BUF = 10*eeg_chunk.srate; % size of convergence buffer
                        srcfig = figure;
                        srcaxes= axes('parent',srcfig);
                        plotdata = nan(eeg_chunk.nbchan,SRCPOT_BUF);
                        plotdata(:,end-size(eeg_chunk.icaact,2)+1:end) = eeg_chunk.icaact;
                        keyboard;
                        % zero-mean
                        plotdata = bsxfun(@minus, plotdata, mean(plotdata,2));
                        
                        % arrange for plotting
                        plotoffsets = (0:size(plotdata,1)-1); %*stream.opts.datascale;
                        plotdata = bsxfun(@plus, plotdata, plotoffsets');
                        plottime = linspace(0,SRCPOT_BUF/eeg_chunk.srate,size(plotdata,2));
                        % === actual drawing ===
                        
                        % draw the block contents...
                        if ~isempty(plotdata)
                            if ~exist('lines','var') || isempty(lines)
                                lines = plot(srcaxes,plottime,plotdata);
                                title(srcaxes,'EEG icaact');
                                xlabel(srcaxes,'Time (sec)','FontSize',12);
                                ylabel(srcaxes,'Activation','FontSize',12);
                            else
                                for k=1:length(lines)
                                    set(lines(k),'Ydata',plotdata(:,k));
                                    set(lines(k),'Xdata',plottime);
                                end
                            end
                            
                            % update the axis limit and tickmarks
                            %                             axis(srcaxes,[stream.xmin stream.xmax -stream.opts.datascale size(plotdata,2)*stream.opts.datascale + stream.opts.datascale]);
                            set(srcaxes, 'YTick',plotoffsets, 'YTickLabel',{eeg_chunk.chanlocs(1:eeg_chunk.nbchan).labels});
                        end
                        
                        drawnow;
                        
                        %                         hsrc=plot(linspace(0,SRCPOT_BUF/eeg_chunk.srate,SRCPOT_BUF),srcpot);
                        %                         set(srcfig,'userdata',hsrc);
                        %                         ylabel('10*ln(Convergence)');
                        %                         xlabel('Time');
                    else
                        blkpnts = size(eeg_chunk.icaact,2); %options.blockSize*(skipCount);
                        plotdata(:,1:end-blkpnts) = plotdata(:,blkpnts+1:end);
                        % insert new samples ...
                        plotdata(:,end-blkpnts+1:end) = eeg_chunk.icaact;
                        %                         set(lines,'ydata',plotdata);
                        for k=1:length(lines)
                            set(lines(k),'Ydata',plotdata(:,k));
                            set(lines(k),'Xdata',plottime);
                        end
                        
                    end
                end
            end
        end
        
        
        updateTime = updateTime + size(eeg_chunk.data,2) / eeg_chunk.srate; % sec
        
        if plotcomponent && updateTime > 0.2
            if fltPipCfg.porica.arg_selection && isfield(eeg_chunk,'icaweights') && isfield(eeg_chunk,'icasphere')
                if ~isempty(eeg_chunk.icaweights)
                    eeg_chunk.icawinv = inv(eeg_chunk.icaweights*eeg_chunk.icasphere);
                    if ~exist('compfig','var') % || isempty(bpfig) || ~ishandle(bpfig)
                        compfig = figure;
                    end
                    
                    compID = 1:14;
                    for plotIdx = 1:length(compID)
                        subplot(4,4,plotIdx);
                        if flag
                            topoplot( eeg_chunk.icawinv(:, compID(plotIdx)), eeg_chunk.chanlocs, ...
                                'style', 'map', 'electrodes', 'off', 'verbose', 'off','gridscale',32);
                            if plotIdx == length(compID);flag = 0; end
                        else
                            handleComp = get(gca,'children');
                            handleComp = handleComp(end);
                            [Zi,amax,amin] = topoplot_easy( eeg_chunk.icawinv(:, compID(plotIdx)), eeg_chunk.chanlocs, ...
                                'style', 'map', 'electrodes', 'off', 'verbose', 'off', 'gridscale',32);
                            set(handleComp,'cdata',Zi);
                            caxis([amin,amax]);
                        end
                    end
                end
            end
            updateTime = 0;
        end
        
    end
    
    % (optionally) give Matlab a little break...
    pause(0.01);
    
end


%%
%         plotconvergence = 1;
%         if plotconvergence
%             % plot ORICA convergence
%             if fltPipCfg.porica.arg_selection && isfield(eeg_chunk,'convergence')
%                 if ~exist('bpfig','var') || isempty(bpfig) || ~ishandle(bpfig)
%                     CONVERGE_BUF = 60*eeg_chunk.srate; % size of convergence buffer
%                     bpfig = figure;
%                     bpaxes= axes('parent',bpfig);
%                     convergence = nan(1,CONVERGE_BUF);
%                     convergence(end-length(eeg_chunk.convergence)+1:end) = 10*log(eeg_chunk.convergence);
%                     h=plot(linspace(0,CONVERGE_BUF/eeg_chunk.srate,CONVERGE_BUF),convergence);
%                     set(bpfig,'userdata',h);
%                     ylabel('10*ln(Convergence)');
%                     xlabel('Time');
%                 else
%                     blkpnts = length(eeg_chunk.convergence); %options.blockSize*(skipCount);
%                     convergence(:,1:end-blkpnts) = convergence(:,blkpnts+1:end);
%                     % insert new samples ...
%                     convergence(:,end-blkpnts+1:end) = 10*log(eeg_chunk.convergence);
%                     set(h,'ydata',convergence);
%                 end
%             end
%         end
%
%         plotnonStatIdx = 0;
%         if plotnonStatIdx
%             if fltPipCfg.porica.arg_selection && isfield(eeg_chunk,'statIdx')
%                 if ~exist('bpfig','var') || isempty(bpfig) || ~ishandle(bpfig)
%                     STAT_BUF = 60*eeg_chunk.srate; % size of convergence buffer
%                     bpfig = figure;
%                     bpaxes= axes('parent',bpfig);
%                     nonStatIdx = nan(1,STAT_BUF);
%                     nonStatIdx(end-length(eeg_chunk.statIdx)+1:end) = eeg_chunk.statIdx;
%                     h=plot(linspace(0,STAT_BUF/eeg_chunk.srate,STAT_BUF),nonStatIdx);
%                     set(bpfig,'userdata',h);
%                     ylabel('Non Stationary Index');
%                     xlabel('Time');
%                 else
%                     blkpnts = length(eeg_chunk.statIdx); %options.blockSize*(skipCount);
%                     nonStatIdx(:,1:end-blkpnts) = nonStatIdx(:,blkpnts+1:end);
%                     % insert new samples ...
%                     nonStatIdx(:,end-blkpnts+1:end) = eeg_chunk.statIdx;
%                     set(h,'ydata',nonStatIdx);
%                 end
%             end
%         end
%
%         plotcomponent = 0;
%         timer = timer + toc;
%         if plotcomponent && timer > 2
%             if fltPipCfg.porica.arg_selection && isfield(eeg_chunk,'icaweights') && isfield(eeg_chunk,'icasphere')
%                 compID = 1:16;
%                 eeg_chunk.icawinv = pinv(eeg_chunk.icaweights*eeg_chunk.icasphere);
%                 if ~exist('compfig','var') % || isempty(bpfig) || ~ishandle(bpfig)
%                     compfig = figure;
% %                     compfig = cell(1,length(compID));
%                 end
%                 if refreshFig > 10;
%                     clf(compfig,'reset');
%                     refreshFig = 0;
%                 end
%                 refreshFig = refreshFig + 1;
%
%                 for plotIdx = 1:length(compID)
%                     subplot(4,4,plotIdx);
%                     topoplot( eeg_chunk.icawinv(:, compID(plotIdx)), eeg_chunk.chanlocs, ...
%                         'style', 'map', 'electrodes', 'off', 'verbose', 'off');
%                 end
% %                 tmpobj = pop_topoplot(eeg_chunk,0,3:4,{},[],'colorbar','off','electrodes', 'off');%, ...
%                                    %'style', 'map', 'electrodes', 'off', 'verbose', 'off');
% %                         refreshdata(compfig{plotIdx},'base')
% %                         drawnow;
% %                 end
%             end
%             timer = 0;
%         end
%
%         storeData = 1;
%         if storeData
%             % store convergence, convEst, statIdx / computeConvergence on
%             if isfield(eeg_chunk,'statIdx') && isfield(eeg_chunk,'alarm')
%                 if isempty(testPI) && isempty(statIdx) && isempty(alarm)
%                     statIdx = cleaned_data.statIdx;
%                     alarm = cleaned_data.alarm;
%                 end
%                 statIdx = [statIdx eeg_chunk.statIdx];
%                 alarm = [alarm eeg_chunk.alarm];
%             end
%
%             if isfield(eeg_chunk,'convergence')
%                 if isempty(testPI)
%                     testPI = cleaned_data.convergence;
%                 end
%                 testPI = [testPI eeg_chunk.convergence];
%             end
%             if isfield(eeg_chunk,'lambda_k')
%                 if isempty(lambda_k)
%                     lambda_k = cleaned_data.lambda_k;
%                 end
%                 lambda_k = [lambda_k eeg_chunk.lambda_k];
%             end
%
%
%             % store meanIdx covIdx alarm / annealFF on
%             if isfield(eeg_chunk,'meanIdx') && isfield(eeg_chunk,'covIdx') && isfield(eeg_chunk,'alarm')
%                 if isempty(meanIdx) && isempty(covIdx) && isempty(alarm)
%                     meanIdx = cleaned_data.meanIdx;
%                     covIdx = cleaned_data.covIdx;
%                     alarm = cleaned_data.alarm;
%                 end
%                 meanIdx = [meanIdx eeg_chunk.meanIdx];
%                 covIdx = [covIdx eeg_chunk.covIdx];
%                 alarm = [alarm eeg_chunk.alarm];
%             end
%
%             % store time and icaweights / storeweights on
%             if isfield(eeg_chunk,'icaweights_all') && isfield(eeg_chunk,'time_all') && isfield(eeg_chunk,'sphere_all')
%                 if isempty(icaweights_all) && isempty(time_all) && isempty(sphere_all)
%                     icaweights_all = cleaned_data.icaweights_all;
%                     sphere_all = cleaned_data.sphere_all;
%                     time_all = cleaned_data.time_all;
%                 end
%                 icaweights_all = cat(3,icaweights_all,eeg_chunk.icaweights_all);
%                 sphere_all = cat(3,sphere_all,eeg_chunk.sphere_all);
%                 time_all = [time_all eeg_chunk.time_all];
%             end
%
%             % store sampleSize and executeTime for timing the filter
%             if isfield(eeg_chunk,'sampleSize') && isfield(eeg_chunk,'executeTime')
%                 if isempty(sampleSize) && isempty(executeTime)
%                     sampleSize = cleaned_data.sampleSize;
%                     executeTime = cleaned_data.executeTime;
%                 end
%                 sampleSize = [sampleSize eeg_chunk.sampleSize];
%                 executeTime = [executeTime eeg_chunk.executeTime];
%             end
%
%             % store sampleSize and executeTime for timing the filter
%             if isfield(eeg_chunk,'executeTimeWhite')
%                 if isempty(executeTimeWhite)
%                     executeTimeWhite = cleaned_data.executeTimeWhite;
%                 end
%                 executeTimeWhite = [executeTimeWhite eeg_chunk.executeTimeWhite];
%             end
%         end
