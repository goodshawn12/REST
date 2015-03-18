function varargout = nut_timef_stats(varargin)
% NUT_TIMEF_STATS M-file for nut_timef_stats.fig
%

% Last Modified by GUIDE v2.5 09-Dec-2011 18:13:17

% Here and below:
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to nut_timef_stats (see VARARGIN)

if(~strcmp('14',version('-release')))
    warning('off')  % hack so warning for "created by matlab 7.0" is turned off
end
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @nut_timef_stats_OpeningFcn, ...
    'gui_OutputFcn',  @nut_timef_stats_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin & ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
warning('on')  % we want to see all other warnings.


%%---------------------------------------------------
function nut_timef_stats_OpeningFcn(hObject, eventdata, handles, varargin)
% --- Executes just before nut_timef_stats is made visible.

% Choose default command line output for nut_timef_stats
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

global stats

stats=struct('type','ttest1_snpm','docomp',false,'fsel','one','frqstats',false,'corr4multfreq',false,'tsel','one','tsstats',false, ...
    'corr4multtime',true,'covar',[],'whichtail','both','markmode',true,'doavg',true,'dostat',true,'usespatnormed',true);

%%---------------------------------------------------
function nut_timef_stats_CloseRequestFcn(hObject, eventdata, handles, varargin)
% --- Executes just before nut_timef_stats is made visible.

delete(hObject);


%%---------------------------------------------------
function varargout = nut_timef_stats_OutputFcn(hObject, eventdata, handles)
if(isfield(handles,'output'))  
    varargout{1} = handles.output;
end


% --- Executes on button press in nut_check_stats.
function nut_check_stats_Callback(hObject, eventdata, handles)
% hObject    handle to nut_check_stats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global stats

if ( strcmpi(get(handles.nut_menu_isspatnorm,'Checked'),'off') && stats.docomp )
    errordlg('The type of average/statistics you selected is not implemented (yet?) for non-spatially normalized files.')
    return
end

if stats.markmode>0
    
    if stats.markmode==1
        % Load Pointer File
        [voifile,voipath]=uigetfile('*.mat','Load Pointer File...');
        if ~ischar(voifile), return, end
        load([voipath voifile]);
        if ~exist('voi','var'), errordlg('Invalid pointer file.'), return, end
    else
        [voifile,voipath]=uigetfile('*.mat','Load Population File...');        
        if ~ischar(voifile), return, end
        poppath = [voipath voifile];
        pop=load(poppath);
        if isfield(pop,'pop'), pop=pop.pop; end
        if ( ~isfield(pop,'s') || iscell(pop.s) ), errordlg('Not a valid population file.'), return, end
        voi.subjnr = pop.subjnr;  
        voi.condnr = pop.condnr;    clear pop
        if isscalar(voi.condnr), voi.condnr = voi.condnr .* ones(size(voi.subjnr)); end
    end
    
    % Subjects and conditions
    if strcmp(get(handles.nut_edit_subnr,'string'),'all')
        stats.subjsel=unique(voi.subjnr); %voi.subjnr(ismember(voi.subjnr,unique(voi.subjnr)));  % This is like unique but not sorted.
    else
        try
            stats.subjsel=evalin('base',get(handles.nut_edit_subnr,'string'));
        catch
            stats.subjsel=str2num(get(handles.nut_edit_subnr,'string'));
        end
    end
    stats.numsubj=length(stats.subjsel);
    switch stats.type
    case 'ttest2_snpm'
        stats.subjnr = voi.subjnr(ismember(voi.subjnr,stats.subjsel));
        if strcmpi(get(handles.nut_edit_condnr,'string'),'all')
            stats.condnr=voi.condnr;
            set(handles.nut_edit_condnr,'string',int2str(stats.condnr))
        else
            stats.condnr=str2num(get(handles.nut_edit_condnr,'string'));
            if isempty(stats.condnr)
                errordlg('Invalid condition number(s).'), return
            end
            if size(stats.condnr,1)>1, errordlg('No semicolons allowed in condition field for unpaired t-tests.'), return, end
            if ~isequal(length(stats.condnr),length(stats.subjnr))
                sco = ismember(voi.condnr,stats.condnr);
                stats.condnr = voi.condnr(sco);
                stats.subjsel = stats.subjsel(ismember(stats.subjsel,voi.subjnr(sco)));
                stats.subjnr = stats.subjsel;   
                stats.numsubj = length(stats.subjsel);
            end
        end
        if ~isequal(length(stats.condnr),length(stats.subjnr))
            errordlg('You must have only 1 condition for each subject.'), return
        end
        stats.comps=unique(stats.condnr); 
        stats.condsel=stats.comps;
        stats.numtests=1;
    case 'anova_snpm'
        stats.subjnr = voi.subjnr(ismember(voi.subjnr,stats.subjsel));
        if strcmp(get(handles.nut_edit_condnr,'string'),'all')
            stats.condnr=voi.condnr;
            stats.comps=unique(stats.condnr);
            set(handles.nut_edit_condnr,'string',int2str(stats.comps))
        else
            stats.comps = str2num(get(handles.nut_edit_condnr,'string'));   
            if size(stats.comps,1)>1, errordlg('No semicolons allowed in condition field for ANOVA.'), return, end
            stats.condnr = voi.condnr(ismember(voi.condnr,stats.comps));            
        end
        stats.condsel=stats.comps;
        if ~isequal(length(stats.condnr),length(stats.subjnr)), errordlg('Something is wrong with your Pointer File.'), return, end        
        stats.numtests=1;
    otherwise    
        if ~stats.docomp
            if strcmp(get(handles.nut_edit_condnr,'string'),'all')
                stats.condsel=unique(voi.condnr);
                set(handles.nut_edit_condnr,'string',int2str(stats.condsel))
            else
                stats.condsel=str2num(get(handles.nut_edit_condnr,'string'));     
            end
            stats.numtests=length(stats.condsel);
            stats.comps=stats.condsel';
            stats.subjnr=repmat(stats.subjsel,[stats.numtests 1]);
            stats.condnr=repmat(stats.comps,[1 stats.numsubj]);
        else
            stats.comps=str2num(get(handles.nut_edit_condnr,'String'));
            if (isempty(stats.comps) || size(stats.comps,2)~=2)
                errordlg('Invalid condition numbers.'), return
            end
            stats.condsel=unique(stats.comps(:))';
            stats.numtests=size(stats.comps,1);
            stats.subjnr=repmat(stats.subjsel,[stats.numtests 2]);
            stats.condnr=[stats.comps(:,1)*ones(1,stats.numsubj) stats.comps(:,2)*ones(1,stats.numsubj)];
        end
    end
    
    % Prepare path and filenames
    if stats.markmode==1
        stats.files=cell(max(stats.subjsel),max(stats.condsel));
        for c=stats.condsel
            for s=stats.subjsel
                if any(any(stats.subjnr==s & stats.condnr==c))
                v=find(voi.subjnr==s & voi.condnr==c);
                if ~isempty(v)
                    [fpath,ffile,fext]=fileparts(voi.pathnames{v});
                    if strcmpi(get(handles.nut_menu_isspatnorm,'Checked'),'off')    % User preference
                        suffix='.mat';
                    elseif ~isempty(cell2mat(regexp(ffile,{'_spatnorm','_subj','_fcmri'})))     % if files are already spatnorm
                        suffix='.mat';
                    else                            % usually, _spatnorm has to be added
                        suffix='_spatnorm.mat';
                    end
                    stats.files{s,c}=fullfile(fpath,[ffile suffix]);
                    if ~exist(stats.files{s,c},'file'), errordlg(sprintf('Could not find file %s.',stats.files{s,c})), return, end
                else  %if ~strcmp(stats.type,'ttest2_snpm')
                    errordlg(sprintf('Could not find condition %d of subject %d.',c,s)), return
                end
                end
            end
        end
        clear voi f* c s v
        
    else        % stats.markmode==2;
        stats.files = cell(stats.numtests,1);
        for nn=1:stats.numtests
            if nn>1
               [voifile,voipath]=uigetfile('*.mat',sprintf('Load Population File #%d (Cancel if same)...',nn));        
               if ischar(voifile), poppath=[voipath voifile]; end
            end
            stats.files{nn}=poppath;
        end
    end
  
else
    stats.numsubj=str2num(get(handles.nut_edit_subnr,'string'));
    if isempty(stats.numsubj) || length(stats.numsubj)>1
        errordlg('Please enter total number of subjects/files in subjects field.'),return
    end
    stats.subjsel = 1:stats.numsubj;
    stats.subjnr  = stats.subjsel;
    switch stats.type
    case 'ttest2_snpm'
        stats.condnr = str2num(get(handles.nut_edit_condnr,'string'));   
        if isempty(stats.condnr)
            errordlg('Invalid condition numbers. Please specify a condition number in the condition field for each subject.'), return
        end
        if size(stats.condnr,1)>1, errordlg('No semicolons allowed in condition field for unpaired t-tests.'), return, end
        if ~isequal(length(stats.condnr),length(stats.subjnr))
            errordlg('Please specify a condition number in the condition field for each subject.'), return
        end
        stats.comps=unique(stats.condnr); 
        stats.condsel=stats.comps;
    otherwise
        stats.condsel = 1;
        stats.condnr  = ones(1,stats.numsubj);
        stats.comps   = 1;
    end
    stats.numtests= 1;
    stats.files=cell(stats.numsubj,max(stats.condsel));
    for s=1:stats.numsubj
        [beamfile,beampath]=uigetfile('s_beamtf*.mat',sprintf('Load s_beamtf file for subject nr %d...',s));
        if ~ischar(beamfile), return, end
        stats.files{s,stats.condnr(s)}=fullfile(beampath,beamfile);
    end
end

% Covariate of interest for correlation
if strcmp(stats.type,'corr_snpm')
    try
        stats.behavdata = evalin('base',get(handles.nut_edit_corrvar,'string'));
        stats.behavdata = stats.behavdata(:);  % make sure it's a column vector
    catch
        errordlg(sprintf('Could not find variable %s in the Matlab workspace.',get(handles.nut_edit_corrvar,'string'))), return
    end
    if length(stats.behavdata)~=stats.numsubj
        errordlg('Number of elements in correlation covariate does not match the number of subjects.'), return
    end
    if ~all(isfinite(stats.behavdata))
        bad = ~isfinite(stats.behavdata);
        stats.behavdata(bad)= [];
        stats.subjsel(bad)  = [];
        if stats.markmode<2, stats.files(bad,:)  = [];
        else stats.subj2remove=find(bad);
        end
        stats.numsubj       = length(stats.subjsel);
        stats.subjnr        = repmat(stats.subjsel,[stats.numtests 1]);
        clear bad
    end       
end

% Confounding covariates
if any(strcmp(stats.type,{'ttest1_snpm' 'ttest2_snpm'}))
    covar = get(handles.nut_edit_covar,'string');
    if ~isempty(covar) && ~strcmpi(covar,'none')
        try
            stats.covar = evalin('base',covar);
        catch
            errordlg(sprintf('Could not find variable %s in the Matlab workspace.',covar)), return
        end
        if size(stats.covar,1)~=stats.numsubj
            errordlg('Column size of confounding covariate(s) do not match the number of subjects.'), return
        end
        stats.covar = stats.covar - repmat(mean(stats.covar,1),[stats.numsubj 1]);
    end
end
    
% Time and Frequency Parameters
stats.frqband=str2num(get(handles.nut_edit_freqnr,'string'));
stats.frqtxt=int2str(stats.frqband')';
stats.timepts=str2num(get(handles.nut_edit_timepts,'string'));
stats.timetxt=int2str(stats.timepts')';
if strcmp(stats.fsel,'mean') && ~stats.frqstats && isfield(stats,'corr4multfreq')       % clear unnecessary fields
    stats=rmfield(stats,'corr4multfreq');
end
if strcmp(stats.tsel,'sum') && isfield(stats,'corr4multtime')
    stats=rmfield(stats,'corr4multtime');
end

% Number of SnPM permutations
if stats.dostat
    stats.numperm=str2num(get(handles.nut_edit_numperm,'string'));  
    switch stats.type  
        case 'ttest2_snpm'
            tmp=find(diff(sort(stats.condnr)));   % tmp: number of subjects in first condition
            autonumperm = round(exp(gammaln(stats.numsubj+1)-gammaln(tmp+1)-gammaln(stats.numsubj-tmp+1))); clear tmp
        case 'anova_snpm'
            autonumperm = 2^(stats.numsubj-1);
        otherwise  % {'ttest1_snpm' 'corr_snpm'}
            autonumperm = 2^stats.numsubj;
    end
    stats.isaproxperm=false;
    if (isempty(stats.numperm) || length(stats.numperm)>1)
        stats.numperm=autonumperm;
    elseif stats.numperm>autonumperm
        errordlg(sprintf('Maximum number of SnPM permutations is %d.',autonumperm)), return
    elseif stats.numperm<autonumperm
        stats.isaproxperm=true;
    end
    set(handles.nut_edit_numperm,'string',int2str(stats.numperm))  
    if stats.numperm<30, errordlg('Not enough subjects/permutations for SnPM.'), return, end
end

% Update GUI
set(handles.nut_edit_subnr,'string',sprintf('%d ',stats.subjsel)) 

% Enable "Run" if no errors so far
set([handles.nut_run_stats handles.nut_population_button],'Enable','on')

% --- Executes on button press in nut_run_stats.
function nut_run_stats_Callback(hObject, eventdata, handles)
% hObject    handle to nut_run_stats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global stats

save nut_stat_settings stats        % save stats settings, in case not enough memory for GUI window
if stats.docomp
    stats=nut_beamcomp(stats);
end
if stats.tsstats || stats.frqstats
    stats=nut_ttest_freqband(stats);  
end
nut_beamstats(stats);

% --- Executes on button press in nut_population_button.
function nut_population_button_Callback(hObject, eventdata, handles)
% hObject    handle to nut_population_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global stats

save nut_stat_settings stats        % save stats settings, in case not enough memory for GUI window
if stats.docomp
    stats=nut_beamcomp(stats);
end
if stats.tsstats || stats.frqstats
    stats=nut_ttest_freqband(stats);  
end
nut_beampopulation(stats);

% --- Executes on button press in nut_clear_stats.
function nut_clear_stats_Callback(hObject, eventdata, handles)
% hObject    handle to nut_clear_stats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global stats

stats=struct('type','ttest1_snpm','docomp',false,'fsel','one','frqstats',false,'corr4multfreq',false,'tsel','one','tsstats',false, ...
    'corr4multtime',true,'covar',[],'whichtail','both','markmode',true,'doavg',true,'dostat',true,'usespatnormed',true);

set_stattype_string(handles);
set(handles.nut_pop_stattype,'value',2);
set([handles.text_covar handles.nut_edit_corrvar],'Enable','off')
set(handles.nut_edit_numperm,'string','automatic','Enable','on')

set(handles.nut_edit_subnr,'string','all')
set(handles.nut_edit_condnr,'string','all','Enable','on')

set(handles.nut_edit_freqnr,'string','all')
set(handles.nut_radio_freqsel1,'Value',1)
set(handles.nut_check_frqstats,'Value',0)
set(handles.nut_check_corr4multfreq,'Enable','on','Value',0)

set(handles.nut_edit_timepts,'string','all')
set(handles.nut_radio_tone,'Value',1)
set(handles.nut_check_corr4multtime,'Enable','on','Value',1)

set(handles.nut_hypo_posneg,'Value',3,'string',{'Synchronization','Desynchronization','Both'},'Enable','on') 

set_filetype_string(handles);
set(handles.nut_pop_filetype,'Value',2)

set([handles.nut_run_stats handles.nut_population_button],'Enable','off')

set(handles.nut_menu_isspatnorm,'Checked','on')

% --- Executes on selection change in nut_hypo_posneg.
function nut_hypo_posneg_Callback(hObject, eventdata, handles)
% hObject    handle to nut_hypo_posneg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns nut_hypo_posneg contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nut_hypo_posneg

global stats

val=get(hObject,'Value');
switch val
case 1
    stats.whichtail='pos';
case 2
    stats.whichtail='neg';
case 3
    stats.whichtail='both';    
end

% --- Executes during object creation, after setting all properties.
function nut_hypo_posneg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_hypo_posneg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in nut_radio_freqsel1.
function nut_radio_freqsel1_Callback(hObject, eventdata, handles)
% hObject    handle to nut_radio_freqsel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nut_radio_freqsel1
global stats
stats.fsel='one';  
set(handles.nut_check_corr4multfreq,'Enable','on')
stats.corr4multfreq=get(handles.nut_check_corr4multfreq,'Value');

% --- Executes on button press in nut_radio_freqmean.
function nut_radio_freqmean_Callback(hObject, eventdata, handles)
% hObject    handle to nut_radio_freqmean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nut_radio_freqmean
global stats
stats.fsel='mean';   
if ~stats.frqstats, set(handles.nut_check_corr4multfreq,'Enable','off'), end
%stats.corr4multfreq=get(handles.nut_check_corr4multfreq,'Value');

% --- Executes during object creation, after setting all properties.
function nut_pop_stattype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_pop_stattype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in nut_pop_stattype.
function nut_pop_stattype_Callback(hObject, eventdata, handles)
% hObject    handle to nut_pop_stattype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns nut_pop_stattype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nut_pop_stattype
global stats

options=get(hObject,'string');
switch options{get(hObject,'value')}
    case 'Average activation from baseline w/o statistics'
        set([handles.text_covar handles.nut_edit_corrvar],'Enable','off')
        set([handles.text41 handles.nut_edit_numperm],'Enable','off')
        set([handles.nut_radio_tmeanstat handles.text48 handles.nut_radio_tone handles.nut_radio_tsum handles.nut_check_corr4multtime],'Enable','off')
        if ~stats.frqstats && ~stats.tsstats, set(handles.nut_hypo_posneg,'Enable','off'), end
        set(handles.nut_edit_condnr','String','all')     
        set(handles.nut_edit_covar,'Enable','off','string','none')
        set(handles.text57,'Enable','off')
        stats.type='none';
        stats.doavg=true;
        stats.docomp=false;        
    case 'Average activation and SnPM t-tests for one sample'
        set([handles.text_covar handles.nut_edit_corrvar],'Enable','off')        
        set([handles.text41 handles.nut_edit_numperm],'Enable','on') 
        set([handles.nut_radio_tmeanstat handles.text48 handles.nut_radio_tone handles.nut_radio_tsum],'Enable','on')
        if strcmp(stats.tsel,'one'), set(handles.nut_check_corr4multtime,'Enable','on'), end
        set(handles.nut_hypo_posneg,'Enable','on')
        set(handles.nut_edit_condnr','String','all')   
        set([handles.text57 handles.nut_edit_covar],'Enable','on')
        stats.type='ttest1_snpm';  
        stats.doavg=true;
        stats.docomp=false;
    case 'Average difference between condition pairs w/o stats'
        set([handles.text_covar handles.nut_edit_corrvar],'Enable','off')        
        set([handles.text41 handles.nut_edit_numperm],'Enable','off')
        set([handles.nut_radio_tmeanstat handles.text48 handles.nut_radio_tone handles.nut_radio_tsum handles.nut_check_corr4multtime],'Enable','off')        
        if ~stats.frqstats && ~stats.tsstats, set(handles.nut_hypo_posneg,'Enable','off'), end        
        set(handles.nut_edit_condnr,'Enable','on','String','e.g. [1 2;3 4]')   
        set(handles.nut_edit_covar,'Enable','off','string','none')
        set(handles.text57,'Enable','off')
        stats.type='none';
        stats.doavg=true;
        stats.docomp=true;        
    case 'Average difference and SnPM paired t-tests'
        set([handles.text_covar handles.nut_edit_corrvar],'Enable','off')        
        set([handles.text41 handles.nut_edit_numperm],'Enable','on') 
        set([handles.nut_radio_tmeanstat handles.text48 handles.nut_radio_tone handles.nut_radio_tsum],'Enable','on')
        if strcmp(stats.tsel,'one'), set(handles.nut_check_corr4multtime,'Enable','on'), end
        set(handles.nut_hypo_posneg,'Enable','on')     
        set(handles.nut_edit_condnr,'Enable','on','String','e.g. [1 2;3 4]')
        set([handles.text57 handles.nut_edit_covar],'Enable','on')
        stats.type='ttest1_snpm';
        stats.doavg=true;
        stats.docomp=true;              
    case 'SnPM unpaired t-tests'
        set([handles.text_covar handles.nut_edit_corrvar],'Enable','off')        
        set([handles.text41 handles.nut_edit_numperm],'Enable','on')
        set([handles.nut_radio_tmeanstat handles.text48 handles.nut_radio_tone handles.nut_radio_tsum],'Enable','on')
        if strcmp(stats.tsel,'one'), set(handles.nut_check_corr4multtime,'Enable','on'), end
        set(handles.nut_hypo_posneg,'Enable','on')        
        set(handles.nut_edit_condnr,'Enable','on','String','all')
        set([handles.text57 handles.nut_edit_covar],'Enable','on')        
        stats.type='ttest2_snpm';  
        stats.doavg=false;
        stats.docomp=false;
    case 'SnPM within-subject ANOVA'
        set([handles.text_covar handles.nut_edit_corrvar],'Enable','off')        
        set([handles.text41 handles.nut_edit_numperm],'Enable','on')
        set([handles.nut_radio_tmeanstat handles.text48 handles.nut_radio_tone handles.nut_radio_tsum],'Enable','on')
        if strcmp(stats.tsel,'one'), set(handles.nut_check_corr4multtime,'Enable','on'), end
        set(handles.nut_hypo_posneg,'Enable','on')        
        set(handles.nut_edit_condnr,'Enable','on','String','all')
        set(handles.nut_edit_covar,'Enable','off','string','none')
        set(handles.text57,'Enable','off')        
        stats.type='anova_snpm';  
        stats.doavg=false;
        stats.docomp=false;    
    case 'SnPM correlation'
        set([handles.text_covar handles.nut_edit_corrvar],'Enable','on')
        set([handles.text41 handles.nut_edit_numperm],'Enable','on') 
        set([handles.nut_radio_tmeanstat handles.text48 handles.nut_radio_tone handles.nut_radio_tsum],'Enable','on')
        if strcmp(stats.tsel,'one'), set(handles.nut_check_corr4multtime,'Enable','on'), end
        set(handles.nut_hypo_posneg,'Enable','on')
        set(handles.nut_edit_condnr','String','all')    
        set(handles.nut_edit_covar,'Enable','off','string','none')
        set(handles.text57,'Enable','off')        
        stats.type='corr_snpm';  
        stats.doavg=false;
        stats.docomp=false;
    otherwise
        errordlg('Somebody messed up your code!'), return
end
stats.dostat=~strcmp(stats.type,'none');
set_filetype_string(handles);

if ~stats.markmode
    nut_pop_filetype_Callback(handles.nut_pop_filetype,[],handles);
end

% --- Executes on selection change in nut_pop_filetype.
function nut_pop_filetype_Callback(hObject, eventdata, handles)
% hObject    handle to nut_pop_filetype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns nut_pop_filetype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nut_pop_filetype
global stats

options=get(hObject,'string');
switch options{get(hObject,'Value')}
    case 's_beam files'
        stats.markmode=0;
        set(handles.nut_edit_subnr,'String','# of files')
        if ~isempty(strmatch(stats.type, {'none' 'ttest1_snpm' 'snpm' 'corr_snpm'}))
            set(handles.nut_edit_condnr,'String','1','Enable','off')
        end
    case 'Pointer file'
        stats.markmode=1;        
        set(handles.nut_edit_subnr,'string','all')
        set(handles.nut_edit_condnr,'Enable','on')
        nut_pop_stattype_Callback(handles.nut_pop_stattype,[],handles);
    case 'Population file(s)'
        stats.markmode=2;        
        set(handles.nut_edit_subnr,'string','all')
        set(handles.nut_edit_condnr,'Enable','on')
        nut_pop_stattype_Callback(handles.nut_pop_stattype,[],handles);
end
set_stattype_string(handles);


% --- Executes during object creation, after setting all properties.
function nut_pop_filetype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_pop_filetype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in nut_check_corr4multfreq.
function nut_check_corr4multfreq_Callback(hObject, eventdata, handles)
% hObject    handle to nut_check_corr4multfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nut_check_corr4multfreq
global stats
stats.corr4multfreq=get(hObject,'Value');

% --- Executes on button press in nut_check_frqstats.
function nut_check_frqstats_Callback(hObject, eventdata, handles)
% hObject    handle to nut_check_frqstats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nut_check_frqstats
global stats
stats.frqstats=get(hObject,'Value');
if stats.frqstats
    set(handles.nut_hypo_posneg,'Enable','on','Value',1,'string',{'Synchronization','Desynchronization'}) 
    stats.whichtail='pos';
    set(handles.nut_check_corr4multfreq,'Enable','on')
    stats.corr4multfreq=get(handles.nut_check_corr4multfreq,'Value');
else
    if ~stats.tsstats, set(handles.nut_hypo_posneg,'string',{'Synchronization','Desynchronization','Both'}), end
    if ~stats.tsstats && strcmp(stats.type,'none'), set(handles.nut_hypo_posneg,'Enable','off'), end
    if strcmp(stats.fsel,'mean'), set(handles.nut_check_corr4multfreq,'Enable','off'), end    
end

% --- Executes on button press in nut_radio_tone.
function nut_radio_tone_Callback(hObject, eventdata, handles)
% hObject    handle to nut_radio_tone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nut_radio_tone
global stats

stats.tsel='one';
set(handles.nut_check_corr4multtime,'Enable','on','Value',1)
stats.corr4multtime=true;
stats.tsstats=false;
if ~stats.frqstats, set(handles.nut_hypo_posneg,'string',{'Synchronization','Desynchronization','Both'}), end

% --- Executes on button press in nut_radio_tsum.
function nut_radio_tsum_Callback(hObject, eventdata, handles)
% hObject    handle to nut_radio_tsum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nut_radio_tsum
global stats

stats.tsel='sum';
stats.tsstats=false;
if ~stats.frqstats, set(handles.nut_hypo_posneg,'string',{'Synchronization','Desynchronization','Both'}), end
set(handles.nut_check_corr4multtime,'Enable','off')

% --- Executes on button press in nut_check_corr4multtime.
function nut_check_corr4multtime_Callback(hObject, eventdata, handles)
% hObject    handle to nut_check_corr4multtime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nut_check_corr4multtime
global stats
stats.corr4multtime=get(hObject,'Value');


% --- Executes on button press in nut_radio_tmeanstat.
function nut_radio_tmeanstat_Callback(hObject, eventdata, handles)
% hObject    handle to nut_radio_tmeanstat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nut_radio_tmeanstat
global stats

stats.tsel='meanstat';
stats.tsstats=true;
stats.corr4multtime=true;
set(handles.nut_hypo_posneg,'Value',1,'string',{'Synchronization','Desynchronization'}) 
set(handles.nut_check_corr4multtime,'Enable','on','Value',1)

% --------------------------------------------------------------------
function nut_save_settings_Callback(hObject, eventdata, handles)
% hObject    handle to nut_save_settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global stats
save nut_stat_settings stats

msgbox('Settings saved as nut_stat_settings.mat.')

% --------------------------------------------------------------------
function nut_menu_isspatnorm_Callback(hObject, eventdata, handles)
% hObject    handle to nut_menu_isspatnorm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global stats

switch get(hObject,'Checked')
    case 'on'   % User unchecks the item
        answer = questdlg('Warning: Using non spatially normalized data works only if all subjects have the same voxel coordinates. This can be the case in some EEG settings. Are you sure you want to use non-normalized data?', ...
            'Question', 'Yes', 'Cancel', 'Yes');
        if strcmp(answer,'Yes')
            set(hObject,'Checked','off')
            stats.usespatnormed=false;            
        end
    case 'off'  % User checks the item
        set(hObject,'Checked','on')
        stats.usespatnormed=true;
end

%----------------------------------------------------------
function isok=set_stattype_string(handles)
global stats

isok=true;
STS={
'Average activation from baseline w/o statistics';
'Average activation and SnPM t-tests for one sample';
'Average difference between condition pairs w/o stats';
'Average difference and SnPM paired t-tests';
'SnPM unpaired t-tests';
'SnPM within-subject ANOVA';
'SnPM correlation'};

if stats.markmode==0
    STS([3 4 6])=[];
elseif stats.markmode==2
    STS([3 4])=[];
end

currval=get(handles.nut_pop_stattype,'Value');
currstr=get(handles.nut_pop_stattype,'String');
currsel=currstr{currval};

newval=strmatch(currsel,STS);
if isempty(newval)
    msgbox('The selected test does not work with s_beam files. Use pointer files instead.')
    set(handles.nut_pop_filetype,'Value',2)
    stats.markmode=1;
    isok=false;
else
    set(handles.nut_pop_stattype,'String',STS,'Value',newval);
end

%----------------------------------------------------------
function isok=set_filetype_string(handles)
global stats

isok=true;
STS={
's_beam files';
'Pointer file';
'Population file(s)'};

if stats.docomp
    STS = STS(2);
elseif strcmp(stats.type,'anova_snpm')
    STS(1)=[];
end

currval=get(handles.nut_pop_filetype,'Value');
currstr=get(handles.nut_pop_filetype,'String');
currsel=currstr{currval};

newval=strmatch(currsel,STS);
if isempty(newval)
    newval=strmatch('Pointer file',STS);
    stats.markmode=1;
    isok=false;
end
set(handles.nut_pop_filetype,'String',STS,'Value',newval);

function nut_edit_corrvar_Callback(hObject, eventdata, handles)
% hObject    handle to nut_edit_corrvar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_edit_corrvar as text
%        str2double(get(hObject,'String')) returns contents of nut_edit_corrvar as a double


% --- Executes during object creation, after setting all properties.
function nut_edit_corrvar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_edit_corrvar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nut_edit_covar_Callback(hObject, eventdata, handles)
% hObject    handle to nut_edit_covar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_edit_covar as text
%        str2double(get(hObject,'String')) returns contents of nut_edit_covar as a double

% --- Executes during object creation, after setting all properties.
function nut_edit_covar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_edit_covar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nut_edit_subnr_Callback(hObject, eventdata, handles)
% hObject    handle to nut_edit_subnr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_edit_subnr as text
%        str2double(get(hObject,'String')) returns contents of nut_edit_subnr as a double
set(handles.nut_run_stats,'Enable','off')

% --- Executes during object creation, after setting all properties.
function nut_edit_subnr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_edit_subnr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nut_edit_condnr_Callback(hObject, eventdata, handles)
% hObject    handle to nut_edit_condnr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_edit_condnr as text
%        str2double(get(hObject,'String')) returns contents of nut_edit_condnr as a double
set(handles.nut_run_stats,'Enable','off')

% --- Executes during object creation, after setting all properties.
function nut_edit_condnr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_edit_condnr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nut_edit_freqnr_Callback(hObject, eventdata, handles)
% hObject    handle to nut_edit_freqnr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_edit_freqnr as text
%        str2double(get(hObject,'String')) returns contents of nut_edit_freqnr as a double
set(handles.nut_run_stats,'Enable','off')

% --- Executes during object creation, after setting all properties.
function nut_edit_freqnr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_edit_freqnr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function nut_edit_numperm_Callback(hObject, eventdata, handles)
% hObject    handle to nut_edit_numperm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_edit_numperm as text
%        str2double(get(hObject,'String')) returns contents of nut_edit_numperm as a double
set(handles.nut_run_stats,'Enable','off')

% --- Executes during object creation, after setting all properties.
function nut_edit_numperm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_edit_numperm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in nut_edit_timepts.
function nut_edit_timepts_Callback(hObject, eventdata, handles)
% hObject    handle to nut_edit_timepts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.nut_run_stats,'Enable','off')

% --- Executes during object creation, after setting all properties.
function nut_edit_timepts_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_edit_timepts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --------------------------------------------------------------------
function Menu_File_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Menu_Preferences_Callback(hObject, eventdata, handles)
% hObject    handle to Menu_Preferences (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


