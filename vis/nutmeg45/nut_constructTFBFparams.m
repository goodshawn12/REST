function [Active,Control] = nut_constructTFBFparams(startpt,endpt,windowsize,stepsize,controlstart,controlend,filename)
% nut_constructTFBFparams   creates time parameter files for the time-frequency beamformer.
%
%  nut_constructTFBFparams(activestart,activeend,windowsize,stepsize,controlstart,controlend,filename)
%    saves time parameter files to the current directory
%
%  [active,control] = nut_constructTFBFparams(...)
%    creates parameters in output variables (no saving)
%
% activestart   [scalar] begin of active time period in ms. 
% activeend     [scalar] end of active time period in ms. 
% windowsize    [vector if you want to create different window lengths for different 
%               frequency bands, otherwise scalar] length of time windows in ms. 
% stepsize      [optional: scalar] latency is ms to start of next time window. 
%               Determines overlap, default is stepsize=windowsize -> no overlap. 
% controlstart  [optional: scalar or vector] begin of baseline time period in ms. 
%               For single-state beamformers, omit this input or leave empty.
% controlend    [optional: scalar or vector] end of baseline time period in ms. 
%               The time range controlend-controlstart can either be of the same 
%               size as windowsize (to calculate activations), or of the same size 
%               as the time range activeend-activestart (to calculate constrasts 
%               between conditions). For single-state beamformers, omit
%               this input.
% filename      [optional] name of new timewindow file
%
% Examples:
% - nut_constructTFBFparams(0,500,[100 200 300],50,-600,-500)  is interpreted as
%   nut_constructTFBFparams(0,500,[100 200 300],50,[-600 -600 -600],[-500 -400 -300])
%   --> activations, creates 3 files with different windowsizes
% - nut_constructTFBFparams(0,500,[100 200 300],50,1000,1100)  is interpreted as
%   nut_constructTFBFparams(0,500,[100 200 300],50,[1000 900 800],[1100 1100 1100])
%   --> activations, creates 3 files with different windowsizes
% - nut_constructTFBFparams(0,500,[100 200 300],50,0,500)
%   --> contrasts, creates 3 files with different windowsizes
% - nut_constructTFBFparams(0,60000,60000,60000)
%   --> single-state beamformer, saves 1 file with 1 windowsize
% - [active,control]=nut_constructTFBFparams(0,500,100,50,-600,-500)
%   --> activations, creates output variables with 1 windowsize, no saving

if nargin<3
    help nut_constructTFBFparams
    return
end
if nargin<4 || isempty(stepsize)
    stepsize=windowsize; 
end
if nargin<6
    controlstart=[]; controlend=[];
end
if nargin<7
    filename='';
end


numwins=length(windowsize);

if nargout>0
    if numwins>1, error('You can not use output variables with more than one windowsize.'), end
    issave=false;
else
    issave=true;
end

isdual = ( ~isempty(controlstart) && ~isempty(controlend) );

if isdual
    if controlend(1)-controlstart(1)==endpt-startpt  % contrast
        iscontrast=true;
    elseif controlend(1)-controlstart(1)-windowsize(1)<10    % activation, tolerate max 10ms of difference
        iscontrast=false;
    else
        error('Your control window size must either be of the same size as windowsize or as activeend-activestart.')
    end  
end

if numwins==1   
    active=[(startpt:stepsize:(endpt-windowsize))' ((startpt+windowsize):stepsize:endpt)']
    if isempty(active)
        error('No active time window found. Make sure that windowsize is not greater than activeend-activestart.')
    end
    if isdual
        if iscontrast
            control=[(controlstart:stepsize:(controlend-windowsize))' ((controlstart+windowsize):stepsize:controlend)']
        else
            control=[controlstart controlend]
            if ((controlstart(1)<startpt && control(end,end)>active(1,1)) || (controlstart(1)>startpt && control(1,1)<active(end,end)))
                disp('WARNING: Your baseline window is partly overlapping with the active time period!!')
            end
        end
    else
        control=[];
    end
    if issave
        save_timewinfile(filename,active,control);
    else
        Active=active;      % this prevents MATLAB from sending active to the ans variable if no output variable is given.
        Control=control;
    end
    
else    % if more than one windowsize, we need to adjust starting and end points
    if isdual
        if iscontrast
            if length(controlstart)>1 || length(controlend)>1
                error('For contrasts, controlstart and controlend must be scalars.')
            end
            controlstart=controlstart.*ones(1,numwins);     % need to replicate for compatibility with activation setup
            controlend  =controlend  .*ones(1,numwins);
        else      % adjust length of baseline window to length of each windowsize, if not already done by user
            if ~isequal(numwins,length(controlstart)) || ~isequal(numwins,length(controlend)) || ~all(controlend-controlstart-windowsize<10)
                if controlstart(1)<startpt     % if baseline before active period
                    controlstart=controlstart(1).*ones(1,numwins);
                    controlend  =controlstart+windowsize;
                else                        % if baseline after active period
                    controlend  =controlend(1).*ones(1,numwins);
                    controlstart=controlend-windowsize;
                end
            end    
        end
    end
    
    for k=1:numwins
        fprintf('\n%d ms windows:\n',windowsize(k))
        active=[((startpt-windowsize(k)+stepsize):stepsize:(endpt-stepsize))' ((startpt+stepsize):stepsize:(endpt+windowsize(k)-stepsize))']
        if isdual
            if iscontrast
                control=[((controlstart(k)-windowsize(k)+stepsize):stepsize:(controlend(k)-stepsize))' ((controlstart(k)+stepsize):stepsize:(controlend(k)+windowsize(k)-stepsize))']
            else
                control=[controlstart(k) controlend(k)]
                if length(startpt)>1
                    starttest=startpt(k)
                else
                    starttest=startpt;
                end
                if ((controlstart(k)<starttest && control(end,end)>active(1,1)) || (controlstart(1)>startpt && control(1,1)<active(end,end)))
                    disp('WARNING: Your baseline window is partly overlapping with the active time period!!')
                end
            end
        else
            control=[];
        end
        save_timewinfile(sprintf('%dms',windowsize(k)),active,control);
    end
end

%-------------------
function save_timewinfile(filename,active,control)

if isempty(filename)
   [fMAT, pMAT] = uiputfile('*.mat', 'Save timewindows as...')
   if isequal(fMAT,0), return, end
   filename = fullfile(pMAT,fMAT);
end
save(filename,'active','control');
