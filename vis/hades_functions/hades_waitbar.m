function progf= waitbar(varargin)

% USAGE: 
% waitbar(prog)                  % Update progress bar.
% waitbar(text)                  % New message text.
% waitbar(prog, text)            % Update progress bar and new message text.
% waitbar(text, arg1, arg2, ...) % New message text, works like sprintf(...).
% waitbar()                      % Deletes the waitbar window.
% 
% And for compatibility with the original Matlab waitbar:
% waitbar(prog, handle)
% handle= waitbar(prog, handle, text)
% handle= waitbar(prog, text, <property name/value list>)
% 
% A Matlab waitbar compatible replacement that is faster and more functional. 
% 
% ARGUMENTS:
% prog         (scalar, 0 <= prog <= 1) Progress value.
% text         (string) Sets/changes a message.
% handle       (numeric) Is ignored -- for compatibility with original waitbar.
% callback     (string) Is ignored -- for compatibility with original waitbar.
% <prop. ...>  (list) Sets figure properties to respective value. Special properties:
%              'DelayPeriod'   (default: 5.0)        Time until window is shown, [0, 60]. 
%              'LingerPeriod'  (default: 1.0)        Time until window is hidden, [0, 60]. 
%              'MinUpdateTime' (default: 0.5)        Minimum time between updates, [0, 10].
%              'Name'          (default: 'Progress') Progress bar window name.
%              'BarColor'      (default: 'b')        Progress bar colour. 
%              'CreateCancelBtn'                     For compatibility. Ignored. 
% 
% FUNCTIONALITY, like Matlab's original waitbar, but with some additions:
% 1.  It stays on top of other figures. 
% 2.  Minimal, really minimal, execution time overhead (some 10% of original waitbar's).
% 3.  Simpler to use, simpler calling and no figure handle to pass along.
% 4.  Only one waitbar window, so no old ones left around.
% 5.  Remaining time estimated and presented.
% 6.  Delay of 5 seconds before the waitbar is shown (no flashing at quick tasks).
% 7.  Window position is remembered during a Matlab session.
% 8.  Window is hidden 2 seconds after reaching prog == 1.0. If another task follows, 
%     there is time to continue progress reporting with no window flashing.
% 9.  If there is no window environment, textual feedback is used.
% 10. Information on process start and cpu time used. 
% 
% BUGS: If you find any bugs, please let me know: peder at axensten dot se.
% 
% EXAMPLE 1, simplest possible:
% nEnd= 100;
% for n= 1:nEnd
% pause(0.2); % Do stuff
% waitbar(n/nEnd);
% end
% 
% EXAMPLE 2, you can be more informative:
% nEnd= 100;
% waitbar(0, '(1 of 2) Preparing...', 'DelayPeriod', 0, 'LingerPeriod', 0, ...
%                 'MinUpdateTime', 0, 'Name', 'Other title...', 'BarColor', [0 1 0]);
% for n= 1:nEnd
% pause(0.1); % Prepare stuff
% waitbar(n/nEnd);
% end
% waitbar(0, '(2 of 2) Calculating...');
% for n= 1:nEnd
% pause(0.1); % Calculate stuff
% waitbar(n/nEnd);
% end
% 
% HISTORY:
% Version 1.0, 2006-06-14. 
% Version 1.1, 2006-06-15.
% - Fixed a timing bug and a debug info leak. 
% Version 1.2, 2006-07-16.
% - Added textual feedback, used when no graphical interface is available. 
% - Improved compatibility with original waitbar. 
% - Rewrote help text and comments a bit. 
% Version 1.3, 2006-07-27.
% - Mixed calls to plot-like functions and waitbar now works better.
% Version 1.4, 2006-09-09.
% - Made default values changeable via <prop/value> pairs. See DelayPeriod, LingerPeriod, 
%   MinUpdateTime, Name, and BarColor under "<prop. ...>", above. 
% - Mixing calls to plot-like functions and waitbar should now work. Many thanks to Yair Altman,  
%   who suggested fixes. 
% - Fixed a minor cosmetic bug (placement of bar). 
% Version 1.5, 2006-09-12.
% - Prefixes percentage to window title, e.g. '29% Progress'. 
% - Now the waitbar window is put on top of all figures after every call to it. 
% - Improved the output of the textual interface. 
% Version 1.6, 2006-09-23.
% - Now uses the GNU General Public License. 
% - Better compatibility with Matlab 6.5. 
% Version 1.7, 2006-09-26.
% - Added information on process start and cpu time used.
% - Improved timer handling. One rare warning removed. 
% - The percentage in the window title wasn't always updated. Fixed. Many thanks to Alberto 
%   Schiavone, who found this and other bugs.
% 
% COPYRIGHT:   2005, 2006 Peder Axensten. 
%              This file may be used according to the GNU General Public License.

% KEYWORDS:    wait bar, progress bar, process, time left, ETA
% INSPIRATORS: progressbar (6922), waitwaitbar (10795), Matlab waitbar.
% REQS:	       Matlab 6.5.1 (R13SP1). Might work (but is untested) on older versions.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%%% Persistent variable's name must not conflict with workplace variable names.
	persistent win; 
	
	
	%%% Check if the minimum time between updates is up.
	thistime= clock;
	if(~isempty(win) && (abs(thistime(6)-win.last) <= win.mintime))
		if((nargin == 1) && isnumeric(varargin{1}) && (varargin{1} < 1))
			return;		% Calling waitbar(prog).
		elseif((nargin == 2) && isnumeric(varargin{1}) && (varargin{1} < 1) && isnumeric(varargin{2}))
			return;		% Calling waitbar(prog, handle).
		end
	end
	win.last=		thistime(6);
	if(~isfield(win, 'mintime'))
		win.mintime=	0.5; % [seconds]			% Default minimum update period.
	end
	
	
	%%% A call to waitbar() hides the progress window.	%%%%% No arguments: Hide window!
	if(nargin == 0)
		db=			dbstack;
		progf=		[];
		if(length(db) < 2),	db(2).name= '';	end			% Calling from the command win.
		try 
			if(isempty(win.start) || isempty(findstr(db(2).name, 'timercb')))
				% We want to reset this, in case we actually finished a process and then 
				% call waitbar without a message (previous one would show). If we close 
				% the waitbar window in midprocess, the message text should not be deleted.
				win.msg=			'';
				set(win.fighdl, 'Name', win.figtitle);	% Remove progress from window title.
				
				set(win.fighdl, 'Visible', 'off');
				win.start=			[];
				win.starttime=		'';
				win.vistimer=		handle_timer(win.vistimer, -1, '');
				win.lingertimer=	handle_timer(win.lingertimer, -1, '');
			end
			progf=	win.fighdl;
		catch
		end
		return;
	end
	
	
	%%% Check if we have a waitbar window already, create a new one if there is none.
	try	
		if(~win.text)
			ud=		get(win.fighdl, 'UserData');		% Does win.fighdl exist?
			if(isempty(ud) || (ud ~= 6808795978))
				win=	make_waitbarwin(win); 
			end
		end
	catch	% If we catch an error here, it's because the progress window was closed.
		win=	make_waitbarwin(win); 
	end
    
	
	%%% Are we starting a new process (making the window visible)?
	if(isempty(win.starttime))
		win.starttime=	sprintf('Started at %d:%02d:%02d, cpu time used: ', ...
									thistime(4), thistime(5), floor(thistime(6)));
		win.startcpu=	cputime;
	end
	
	
	%%% We have a waitbar and we want to do something with it. 
	prog= 			varargin{1};						% Progress value. 
	msgstr=			0;									% Default value (no message).
	if(isnumeric(prog) && (numel(prog) == 1))			%%%%% Update waitbar.
		
		
		%% Is there a message string in the call? If so, get it.
		%% Are there any property name/value pairs in the call? If so, get them.
		if((nargin == 3) && isnumeric(varargin{2}) && ischar(varargin{3}))
			msgstr=		varargin{3};		% waitbar(prog, handle, message).
		elseif((nargin == 2) && (isnumeric(varargin{2}) || isempty(varargin{2})))
			% Do nothing, but let through.	% waitbar(prog, handle).
		elseif((nargin >= 2) && ischar(varargin{2}) && (2*floor(nargin/2) == nargin))
			msgstr=		varargin{2};		% waitbar(prog, message, <optional name/value list>).
			if(~win.text && (nargin > 2) && ~strcmp(varargin{3}, 'CreateCancelBtn'))
				win=		set_properties(win, prog, varargin{3:end});
			end
		elseif(nargin >= 2)
			error('Wrong input argument(s), see ''help waitbar'' for correct use. [A]');
		end
		
		%% Create the eta string.
		if((prog < 0) || (prog > 1.001))
			error('The progress value must be in the interval [0,1].');
		elseif(prog >= 1)
			prog=		1;								% Maybe rounding errors.
			etastr=		'100% [0:00:00]';				% No ETA at end.
			
			% Wait to remove waitbar in case another (sub)task follows directly.
			win.lingertimer=	handle_timer(win.lingertimer, win.linger, [mfilename ';']);
			win.start=		[];							% Mark that we're done.
		elseif((prog == 0) || isempty(win.start))
			win.start=	clock;							% (Re)start timing.
			win.head=	true;							% We need to rewrite progress header.
			etastr=		'';								% No ETA at start.
		else
			s=		etime(clock,win.start);				% Time used so far.
			etastr=	sprintf('%3d%% [%s]', floor(100*prog), get_timestr(1 + s*(1/prog - 1)));
		end
		
		%% Update the progress information. 
		if(win.text)
			win=				text_progress(prog, etastr, msgstr, win);
		else
			set(win.etahdl,  'String', sprintf('%s\n', etastr));
			set(win.timehdl, 'String', [win.starttime get_timestr(cputime - win.startcpu)]);
			set(win.fighdl,  'Name',   [num2str(floor(100*prog)) '% ' win.figtitle]);
			set(win.proghdl, 'XData',  [0 prog prog 0]);
			win.winpos=			get(win.fighdl, 'Position');% Save progress window position. 
			if(ischar(msgstr)), set(win.msghdl, 'String', sprintf('%s\n', msgstr));	end
			drawnow;									% Force redraw.
		end
		if(ischar(msgstr)),		win.msg= msgstr;	end
	
	%%% We are calling with message argument(s) only. 
	elseif(ischar(prog))								%%%%% New message.
		progf=		waitbar(0, sprintf(varargin{:}));	% Default progress in this case is 0. 
		return;
	else 
		error('Wrong input argument(s), see ''help waitbar'' for correct use. [B]');
	end
	
	%%% Start delay timer. 
	if(strcmpi(get(win.fighdl, 'Visible'), 'off') && isempty(win.vistimer))
		win.vistimer=	handle_timer(win.vistimer, win.delay, ...
				['try set(' num2str(win.fighdl, 99) ', ''Visible'', ''on''); catch end']);
	end
	
	%%% Put the waitbar window on top of all others.
	children=		allchild(0);
	if(~win.text && (numel(children) > 1) && (children(1) ~= win.fighdl))
		uistack(win.fighdl, 'top');
	end
	
	%%% Return the figure handle (for campatibility reasons only). 
	progf=			win.fighdl;
end


function handle= handle_timer(handle, delay, timerfcn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Handle timers.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	if(~isempty(handle))		% Remove timer. 
		stop(handle);
		delete(handle);
		handle=		[];
	end
	if(delay >= 0)				% Run timer. 
		try		% We might not have the required JAVA engine (textual interafce). 
			handle=	timer('TimerFcn', timerfcn, 'StartDelay', delay, 'Tag', 'waitbarTimer');
			start(handle);
		catch; end
	end
end


function timestr= get_timestr(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Return a time string, given seconds.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	h=			floor(s/3600);						% Hours.
	s=			s - h*3600;
	m=			floor(s/60);						% Minutes.
	s=			s - m*60;							% Seconds.
	timestr=	sprintf('%d:%02d:%02d', h, m, floor(s));
end


function data= set_properties(data, prog, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Read property/value list and set values accordingly.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	if(mod(length(varargin), 2) ~= 0)
		error('Missing a value in property name/value list.');
	end
	for n= 1:2:length(varargin)
		switch(lower(varargin{n}))
			case 'createcancelbtn'	% Ignore...
			case 'delayperiod'
				data.delay=		min(60, max(0, varargin{n+1}));	% In [0, 60]. 
			case 'lingerperiod'
				data.linger=	min(60, max(0, varargin{n+1}));	% In [0, 60]. 
			case 'minupdatetime'
				data.mintime=	min(10, max(0, varargin{n+1}));	% In [0, 10]. 
			case 'name'
				data.figtitle=	varargin{n+1};
			case 'barcolor'
				data.barcolor=	varargin{n+1};
				set(data.proghdl, 'FaceColor', varargin{n+1});
			otherwise
				set(data.fighdl, varargin{n}, varargin{n+1});
		end
	end
end


function data= text_progress(prog, etastr, msgstr, data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Write progress.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	if(ischar(msgstr))
		if(~strcmp(data.msg, msgstr))
			fprintf('%s\r***** %s\n', repmat(' ', 1, data.maxlen+17), msgstr);
		end
		data.header=	true;
	end
	
	if(prog > 0)
		if(data.header)
			width=		get(0,'CommandWindowSize');		% [width height] in chars.
			if    (width(1) >= 110),	divs=	10;
			elseif(width(1) >=  70),	divs=	 5;
			elseif(width(1) >=  60),	divs=	 4;
			elseif(width(1) >=  40),	divs=	 2;
			else						divs=	 1;
			end
			data.maxlen=	divs*10;
		
			data.theprog=		repmat('=', 1, data.maxlen);
			data.theback=		repmat('         :', 1, divs);
			data.theback(end)=	' ';
			data.header=		false;
		end
	
		fprintf('[%s%s] %s\r', data.theprog(1:floor(prog*data.maxlen)), ...
									data.theback(floor(prog*data.maxlen)+1:end), etastr);
		if(prog >= 1), fprintf('%s\r', repmat(' ', 1, data.maxlen+17));	end
	end
end


function data= make_waitbarwin(data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Create default waitbar (gui or text?).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	if(~isfield(data, 'text'))
		data.text=			false;
		data.maxlen=		0;
		data.delay=			5.0; % [seconds]			% Default delay until window is visible.
		data.linger=		1.0; % [seconds]			% Default delay until window is hidden.
		data.figtitle=		'Progress';					% Default bar window title. 
		data.barcolor=		'b';						% Default bar color. 
		data.msg=			'';							% No message yet. 
		data.fighdl=		[];
		data.proghdl=		[];
		data.etahdl=		[];
		data.msghdl=		[];
		data.start=			[];							% Start present bar timing.
		data.starttime=		'';							% Time when progress window was shown.
		data.vistimer=		[];
		data.lingertimer=	[];
	end
	
	warning('off', 'MATLAB:m_warning_end_without_block');	% For Matlab 6.5. 
	
	% Call gui or text
	if(1 == prod(get(0,'ScreenSize')))				%%%%% Textual interface.
		data.text=			true;
		data.header=		true;
	else											%%%%% Graphival interface.
		% Get environment data.
		oldUnits=			get(0, 'Units');
		scrSz=				get(0, 'ScreenSize');
		set(0, 'Units', 'points');
		ppp=				72/get(0, 'ScreenPixelsPerInch');
		set(0, 'Units', oldUnits);
	
		% Calculate position of waitbar window and bar.
		barSze=				[400 16];	% [hor vert]
		barPlc=				[barSze(2) barSze(2)];
		barPos=				[barPlc barSze];
		winSze=				barSze + 2*barPlc + [0 30];
		if(~isfield(data, 'winpos'))
			data.winpos=		[(scrSz(3:4)-winSze)/2 winSze]*ppp;	% [left bottom width height]
		end
	
		% Delete all pre-existing waitbar graphical objects. 
		showhid=			get(0, 'showhid');
		set(0, 'showhid', 'on');
		try delete(findobj('Tag', 'Waitbar_Fig', 'UserData', 6808795978)); catch end
		set(0, 'showhid', showhid);
		
		% Delete timers. 
		data.vistimer=		handle_timer(data.vistimer, -1, '');
		data.lingertimer=	handle_timer(data.lingertimer, -1, '');
		
		% Create the interface objects. 
		data.fighdl=	figure(...
			'Units',				'points', ...
			'Position',				data.winpos, ...
			'NumberTitle',			'off', ...
			'Resize',				'off', ...
			'MenuBar',				'none', ...
			'Visible', 				'off', ...
			'Name', 				data.figtitle, ...
			'CloseRequestFcn',		'waitbar();', ...
			'HandleVisibility',		'off', ...
			'IntegerHandle',		'off', ...
			'Tag',					'Waitbar_Fig', ...
			'UserData', 			6808795978 );
		axhdl=			axes(...
			'Units',				'points', ...
			'Position',				barPos*ppp, ...
			'XLim',					[0 1], ...
			'YLim',					[0 1], ...
			'Box',					'on', ...
			'parent',				data.fighdl, ...
			'ytick',				[], ...
			'xtick',				[] );
		data.proghdl=	patch(...
			'XData',				[0 0 0 0], ...
			'YData',				[0 0 1 1], ...
			'FaceColor',			data.barcolor, ...
			'parent',				axhdl, ...
			'EraseMode',			'none' );
		data.etahdl=	text(1, 1, '', ...
			'parent',				axhdl, ...
			'Interpreter',			'none', ...
			'VerticalAlignment',	'bottom', ...
			'HorizontalAlignment',	'right' );
		data.timehdl=	text(1, -0.2, '', ...
			'parent',				axhdl, ...
			'Interpreter',			'none', ...
			'VerticalAlignment',	'top', ...
			'HorizontalAlignment',	'right', ...
			'FontSize',				9 );
		data.msghdl=	text(0, 1, data.msg, ...
			'parent',				axhdl, ...
			'Interpreter',			'none', ...
			'FontWeight',			'bold', ...
			'VerticalAlignment',	'bottom', ...
			'HorizontalAlignment',	'left' );
	end
end

