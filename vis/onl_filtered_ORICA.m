function onl_filtered_ORICA(~,~,stream_name)
% Obtain processed data from a filter pipeline online.
%
% This function returns a chunk of most recent filtered output from a filter pipeline.
% 
% A filter pipeline is a recursive data structure (like a tree) whose nodes are the filter stages and 
% whose edges represent the output of one stage going into the input of another stage. The leaf
% nodes refer to raw online streams (structs in the workspace) and are queried via onl_peek, all other
% nodes are evaluated by calling the filter function on its input data (recursively).
%
% In:
%   Pipeline : previous filter pipeline struct
%
%   DesiredLength : number of samples to get (or 0 to get all new samples) (default: 0)
%
%   SuppressOutput : suppress console output (default: true)
%
%   SetOnlineScope : set the regular online-processing scope (can be turned off for efficiency if
%                    that scope is already set for some reason) (default: true)
%
% Out:
%   Chunk : EEGLAB dataset struct representing the desired data chunk
%           Can be shorter than desired length if not enough data is available at the moment; if the
%           chunk is epoched, the desired length is ignored.
%
%   Pipeline : updated filter pipeline struct
%
%
% Example:
%   % load calibration set
%   raw = io_loadset('calib.set')
%
%   % apply a series of filter to it (the processed set now has a filter expression and initial state)
%   processed = exp_eval(flt_iir(flt_resample(raw,128),[0.5 1],'highpass'));
%
%   % start streaming some data
%   run_readdataset('mystream','action.set');
%   % and put a pipeline on top of it that replicates the processing applied to processed and continues it on new data
%   pip = onl_newpipeline(processed,{'mystream'});
%
%   while 1
%      % generate a 200-sample view into the processed stream
%      [EEG,pip] = onl_filtered(pip,200);
%   end
%
% See also:
%   onl_newpipeline, onl_newstream, onl_append, onl_peek
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2012-05-13

p = evalin('base','pipeline');

% run update_pipeline() with appropriate options
[console_output,chunk,p] = evalc('hlp_scope({''disable_expressions'',1,''is_online'',1},@update_pipeline,p)'); %#ok<ASGLU>

if isempty(chunk.data)
    return
end

% save pipeline in base workspace
assignin('base','pipeline',p);

% save convergence info in base workspace % !!! clean for output
conv_statIdx = evalin('base','conv_statIdx');
conv_mir = evalin('base','conv_mir');
learning_rate = evalin('base','learning_rate');
if isfield(chunk,'statIdx')
    len = length(chunk.statIdx);
    assignin('base','conv_statIdx',[conv_statIdx(len+1:end) chunk.statIdx]);
    assignin('base','conv_mir',[conv_mir(len+1:end) chunk.mir]);
    assignin('base','learning_rate',[learning_rate(len+1:end) db(chunk.lambda_k)]);
end

% raw data: p.parts{2}.parts{2}.out
% iir filtered data: p.out
% whitened data: chunk.icasphere
% ica activations: chunk.icaact

% update vis_stream_ORICA buffer
buffername = ['lsl_' stream_name '_stream'];
buffer = evalin('base',buffername);
try
    % save data to buffer
    strin = '] = deal(buffer.smax + size(p.out,2),';
    strout = '[buffer.smax,';
    [strin,strout] = build_command(p,0,[],strin,strout,'buffer',false);
    eval([strout(1:end-1) strin(1:end-1) ');']);
catch e
    % reformat buffer: !!! this won't work if we're adding back in channels or if channles are removed later on
    [strin,strout] = build_command(p,0,[],'] = deal(','[','buffer',true);
    eval([strout(1:end-1) strin(1:end-1) ');']);
    for it = 1:length(t)
        buffer.data{it}(t(it)+1:end,:) = []; end
    
    % save data to buffer
    strin = '] = deal(buffer.smax + size(p.out,2),';
    strout = '[buffer.smax,';
    [strin,strout] = build_command(p,0,[],strin,strout,'buffer',false);
    eval([strout(1:end-1) strin(1:end-1) ');']);
end

% save buffer
assignin('base',buffername,buffer)


function [chunk,p,n] = update_pipeline(p)
% Update the given filter pipeline and get a chunk of the newly appended output
% [Chunk,Pipeline] = update_pipeline(Pipeline)
%
% A pipeline is a recursive data structure (like a tree) whose nodes are the filter stages and 
% whose edges represent the output of one stage going into the input of another stage. The leaf
% nodes refer to raw data streams (structs in the workspace) and are queried via onl_peek, all other
% nodes are evaluated by calling the filter function on its input data (recursively).
%
% At each node we store the filter function (.head) and its arguments (.parts), some of which may be
% input filter pipelines themselves, plus some miscellaneous book-keeping data. These include:
%  * .israw : true for nodes that represent raw data
%  * .pipeline_indices : indices of those input arguments that are pipelines themselves (if not raw)
%  * .stateful : true if the node has state
%  * .state : previous filter state, if stateful
%  * .smax : number of samples seen so far (if israw)
%
% In:
%   Pipeline : previous filter pipeline struct
%
% Out:
%   Chunk : EEGLAB dataset struct representing newly appended data, filtered
%
%   Pipeline : updated filter pipeline struct

inputs = p.parts;
if p.subnodes
    % update input pipelines to the current node and store the results in inputs
    for k=p.subnodes
        [inputs{k},p.parts{k},n] = update_pipeline(p.parts{k}); end
    % process the inputs by calling the respective filter function
    if p.stateful
        [chunk,p.state] = p.head(inputs{:},'state',p.state,'arg_direct',true);
    else
        chunk = p.head(inputs{:});
    end
    if isempty(chunk.icaact)
        p.out = chunk.data(:,end-n+1:end);
    else
        p.out = chunk.icaact(:,end-n+1:end);
    end
else
    % get the most recent samples since our buffer's smax from a raw stream: 
    % inputs holds the cell array {stream_name,channel_range)
    chunk = onl_peek(inputs{1},p.smax,'index',inputs{2});
    p.smax = chunk.smax;
    p.out = chunk.data;
    n = size(chunk.data,2);
end

function [strin,strout,count] = build_command(p,count,index,strin,strout,buffername,flag_channelCount)
if p.subnodes
    % recursive call for deeper parts of pipeline
    for k=p.subnodes
            [strin,strout,count] = build_command(p.parts{k},count,[index '.parts{' num2str(k) '}'],strin,strout,buffername,flag_channelCount);
    end
end
% build the outputs
count = count + 1;
if ~flag_channelCount
    strout = [strout 'buffer.data{' num2str(count) '}(:,1+mod(buffer.smax:buffer.smax+size(p.out,2)-1,buffer.pnts)),'];
    strin = [strin 'p' index '.out,'];
else
    strout = [strout 't(' num2str(count) '),'];
    strin = [strin 'size(p' index '.out,1),'];
end





