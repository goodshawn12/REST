function h = nut_progress(iter,total,dispmode,varargin)
% h = nut_progress(iter,total,dispmode,varargin)
% called to display progress for various situations in NUTMEG
% dispmode = 0 for no display
%          = 1 to display a dot per iteration in the command window
%          = 2 to display a percentage in the command window
%          = 3 for GUI progress bar
%          = 4 for GUI progress bar with a cancel button
% if a GUI progress bar is used, h is the handle. To determine if the
% cancel button was pressed, use ishandle(h).

if(~exist('dispmode','var'))
    dispmode = 1;
end

if(~length(varargin))
    h = [];
else
    h = varargin{1};
end

switch(dispmode)
    case 0
        return
    case 1
        fprintf('.')
    case 2
        if(iter==1)
            disp(char(32*ones(1,64))); % 64 blanks
            t=0;
            tic
        end
        
%         checkpoint_vec = 0.05:0.05:1;
%         checkpoint = floor(total*checkpoint_vec);
%         [idplease,checkpointnum]=ismember(iter,checkpoint);
%         if(idplease)
%             percent = 100*checkpoint_vec(checkpointnum);

        percent = iter/total;
        if(~mod(percent,0.01))
            t=toc;
            str=['progress: ' num2str(100*percent) '% done; ETA: ' num2str(ceil(t*(1-percent)/percent)) ' seconds'];
            str((end+1):64)=' '; % pad with blanks
            strkill=char(8*ones(1,65)); % char(8) is backspace
            disp([strkill str]);
            
            if(percent == 100 & t >= 120)
                load nutdone
                sound(nutdone.snd,nutdone.Fs);
                toc
            end
        end
    case {3,4}
        modpc=5;
        percent = round((iter/total*100)*modpc)/modpc;
        if isempty(h)
            if dispmode == 3
                h = waitbar(percent/100,['Progress: ' num2str(percent) '% done']);
            else
                h = waitbar(percent/100,['Progress: ' num2str(percent) '% done'],'CreateCancelBtn',@closereq);
            end
        elseif (~mod(percent,modpc))
            waitbar(percent/100,h,['Progress: ' num2str(percent) '% done']);
        end
end

function closereq(varargin)
shh = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on');
currFig = get(0,'CurrentFigure');
set(0,'ShowHiddenHandles',shh);
delete(currFig);