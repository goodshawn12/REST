function varargout = nut_results_viewer(varargin)
% NUT_RESULTS_VIEWER M-file for nut_results_viewer.fig
%

% Last Modified by GUIDE v2.5 13-Jun-2012 10:42:59

% Here and below:
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to nut_results_viewer (see VARARGIN)

if(~strcmp('14',version('-release')))
    warning('off')  % hack so warning for "created by matlab 7.0" is turned off
end
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @nut_results_viewer_OpeningFcn, ...
    'gui_OutputFcn',  @nut_results_viewer_OutputFcn, ...
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
warning('off','MATLAB:FINITE:obsoleteFunction'); % but not this one
                                                 % New Matlab versions complain
                                                 % about SPM2

%%---------------------------------------------------
function nut_results_viewer_OpeningFcn(hObject, eventdata, handles, varargin)
% --- Executes just before nut_results_viewer is made visible.

% Choose default command line output for nut_results_viewer
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

global rivets st nuts ndefaults;

% start nutmeg if not already started
if(isempty(nuts))
    nutmeg
end

global beam;

rivets.fig = hObject;
rivets.sliderenable = 'on';
nut_defaults;

if(isempty(varargin)) %occurs if not using matlab 7.0
    global bolts
    if isfield(bolts,'outname')
        default_beamfile=bolts.outname;
    elseif(isfield(nuts,'megfile'));
        default_beamfile = ['s_beam_' nuts.meg.filename '.mat'];
    else
        default_beamfile = '';
    end
    clear global bolts

    [ufile, upath]=uigetfile('*.mat','Select s_beam or pointer file...',default_beamfile);
    if isequal(ufile,0) || isequal(upath,0)
        delete(hObject)
        return;
    end
    upath=fullfile(upath,ufile); clear ufile
else  % if varargin is not empty, it contains path
    upath=varargin{1};
end

barhandle = waitbar(0,'Preparing to display results...');

if exist(upath,'file') || exist([upath '.mat'],'file')
    beam=load(upath);
    if(isfield(beam,'voi'))     % If a pointer file is given as input
        voi=beam.voi;
        beam=[];
    end

    if exist('voi','var')    % if pointer file
        rivets.voi=voi; clear voi
        [beampathname,rivets.beamfilename,dum]=fileparts(rivets.voi.pathnames{1});     % beamfile is always without .mat extension to avoid problems below
        rivets.beampath=fullfile(beampathname,rivets.beamfilename);
        beam=load(rivets.beampath); % load first dataset
        set(handles.nut_getset_patientnr,'string',[cellstr(int2str(unique(rivets.voi.subjnr)'));'New'], ...
            'Value',find(unique(rivets.voi.subjnr)==rivets.voi.subjnr(1)))
        set(handles.nut_getset_conditionnr,'string',[cellstr(int2str(unique(rivets.voi.condnr)'));'New'], ...
            'Value',find(unique(rivets.voi.condnr)==rivets.voi.condnr(1)))
    else                    % if s_beam file
        rivets.beampath=upath;
        [beampathname,rivets.beamfilename,dum] = fileparts(rivets.beampath);           % beamfile is always without .mat extension to avoid problems below
        if isempty(beampathname), beampathname=pwd; rivets.beampath=fullfile(beampathname,rivets.beamfilename);   % make sure the full path is stored
        else, dum=pwd; cd(beampathname), beampathname=pwd; rivets.beampath=fullfile(beampathname,rivets.beamfilename); cd(dum), clear dum
        end
    end
    clear beampathname dum
else
    delete(hObject), delete(barhandle)
    error('NUT_RESULTS_VIEWER: Invalid filename.')
end
set(hObject,'Name',['NUTMEG Results Viewer - ' rivets.beamfilename]);

%%% Legacy compatibility %%%
beam=nut_beam_legacy_compatibility(beam);
if ~isfield(beam,'s')
    delete(hObject), delete(barhandle)
    error('NUT_RESULTS_VIEWER: Not a valid s_beam file or pointer file.')
end
% evoked results did not previously have bands field, now required
if(~isfield(beam,'bands'))
    beam.bands = [0 1000];  % [0 Inf];  Inf will crash the viewer in timef mode!
elseif any(beam.bands(:,2)==0)
    beam.bands(beam.bands(:,2)==0,2)=1000;  % This can happen with highpass filter
end
% This seems to be required now as well...
if ~isfield(beam,'params') || ~isfield(beam.params,'beamformertype')
    beam.params.beamformertype='unknown';
end

% Add some stats corrections if possible
beam = nut_FDR(beam,ndefaults.FDR);
if size(beam.s{1},2)>5, beam = nut_durationcorrect_stats(beam); end

%%% set all experimental flags here... %%%
if( (size(beam.s{1},3) > 1 || length(beam.timepts)<3) && ~strcmp(beam.params.beamformertype,'nut_Saketini') && ~strcmp(beam.params.beamformertype,'nut_Champagne'))

    rivets.timefmode=true;     % display in color coded time-frequency plane
    rivets.normflag=ndefaults.tfbf.normintens;  %option to normalize to 0:1000 or leave s_beam^2 values as is.
    set(handles.nut_norm_box,'Value',ndefaults.tfbf.normintens)
    set(handles.nut_check_logscale,'Enable','on')
    set([handles.nut_view_interval_button handles.nut_view_time_text handles.nut_view_time2_text handles.text13 handles.text14],'Enable','off')
else
    rivets.timefmode=false;     % display line graph
    rivets.normflag=ndefaults.tsbf.normintens;
    set(handles.nut_norm_box,'Value',ndefaults.tsbf.normintens)
    set(handles.nut_check_logscale,'Enable','off')
    set([handles.nut_view_interval_button handles.nut_view_time_text handles.nut_view_time2_text handles.text13 handles.text14],'Enable','on')
end
if ( size(beam.s,2)>2 && (~isfield(beam,'sinfo') || length(beam.sinfo)<3 || strcmpi(beam.sinfo{3},'noise')) ) 
    set(handles.nut_subtractnoise_box,'Enable','on')    
else
    set(handles.nut_subtractnoise_box,'Enable','off')
end
rivets.displaymode = 1; % 1 -> normal blob overlay on orthogonal MRI slices, 2-> overlay on rendered brain
if(ndefaults.vwr.plotweights)
    rivets.plotweights = true;
    rivets.W_fig = figure;
    rivets.W_ax = axes('Parent',rivets.W_fig);
end

% set GUI to default
set(handles.nut_saverender_button,'Value',0)
set(handles.nut_colorscale_menu,'Value',1)
set([handles.nut_maxtime_button handles.nut_maxactivation_button handles.nut_maxvoxel_button],'Enable','on')
if isfield(beam,'sinfo') && any(strcmpi(beam.sinfo{1},{'T' 'F' 'R' 'correlation'}))
    threstring={}; 
elseif isfield(beam,'sinfo') && ~isempty(strfind(beam.sinfo{1},'Coherence'))
    threstring={'Coherence'};
else
    threstring={'Power'}; 
end 
if isfield(beam,'snpm')
    if isfield(beam,'sinfo'), snpmtestlabel = ['SnPM ' beam.sinfo{1}]; else snpmtestlabel='SnPM T'; end
    threstring=[threstring;{snpmtestlabel;'SnPM p uncorr';'SnPM p corr';'SnPM p FDR corr'}]; clear snpmtestlabel
    if isfield(beam.snpm,'p_cluster_corr_pos') || isfield(beam.snpm,'p_cluster_corr_neg'), threstring=[threstring;{'SnPM p cluster corr'}]; end
    if isfield(beam.snpm,'p_duration_corr_pos') || isfield(beam.snpm,'p_duration_corr_neg'), threstring=[threstring;{'SnPM p duration corr'}]; end
end
if isfield(beam,'ttest')
    threstring=[threstring;{'T-test T'}]; 
    fn=fieldnames(beam.ttest); pidx = strmatch('p_',fn); fn=strrep(fn(pidx),'_',' ');
    threstring=[threstring;cellfun(@(x) ['T-test ' x],fn,'UniformOutput',false)]; clear fn pidx
%     if isfield(beam.ttest,'p_cluster_corr'), threstring=[threstring;{'T-test p cluster corr'}]; end    
%     if isfield(beam.ttest,'p_duration_corr'), threstring=[threstring;{'T-test p duration corr'}]; end
end
if isfield(beam,'corr') && isfield(beam.corr,'p_uncorr')
    threstring=[threstring;{'Corr Coeff (r)'}]; 
    fn=fieldnames(beam.corr); pidx = strmatch('p_',fn); fn=strrep(fn(pidx),'_',' ');
    threstring=[threstring;cellfun(@(x) ['Corr ' x],fn,'UniformOutput',false)]; clear fn pidx    
%     if isfield(beam.corr,'p_cluster_corr'), threstring=[threstring;{'Corr p cluster corr'}]; end
%     if isfield(beam.corr,'p_duration_corr'), threstring=[threstring;{'Corr p duration corr'}]; end    
end 
if isfield(beam,'normcdf'), threstring=[threstring;{'p'}]; end
thresval=get(handles.nut_threshold_menu,'Value');
if thresval > length(threstring)    % wrong setting
    thresval=1;
    set(handles.nut_thresholdpos_text,'String','0','Enable','on');
    set(handles.nut_thresholdneg_text,'String','0','Enable','on');
end
set(handles.nut_threshold_menu,'string',threstring,'Value',thresval), clear threstring thresval dum
if isfield(beam,'connectionfile')
    set(handles.nut_fixvoxel,'Enable','on','Value',0)
    if isfield(beam,'rois')
        set([handles.nut_edit_roiradius handles.text43 handles.text44],'Enable','off')
    else
        set([handles.nut_edit_roiradius handles.text43 handles.text44],'Enable','on')
    end
else
    set([handles.nut_fixvoxel handles.nut_edit_roiradius handles.text43 handles.text44],'Enable','off')
end
if( ~rivets.timefmode && ndims(beam.s{1})>3 )
    set([handles.nut_orientated_menu handles.nut_noisetype],'Enable','on')
else
    set([handles.nut_orientated_menu handles.nut_noisetype],'Enable','off')
end
if ( isfield(beam,'labels') && isfield(beam.labels,'contrasts') )
    set(handles.nut_contrastselect_menu,'string',beam.labels.contrasts)
else    % Reset to default in case it has been changed
    set(handles.nut_contrastselect_menu,'string',  {'F-ratio (dB)'; ...
                                                    'F-ratio (raw)'; ...
                                                    '% Change'; ...
                                                    'CTF Pseudo-F'; ...
                                                    'Difference'; ...
                                                    'Difference (noise normalized)'; ...
                                                    'Greater of Active/Control'; ...
                                                    'Active Power'; ...
                                                    'Control Power'; ...
                                                    'Noise Power'; ...
                                                    'Wilcoxon Z'; ...
                                                    'Active/Noise'});
end

waitbar(0.3,barhandle);

spm_image('init',beam.coreg.mripath); % load/reload structural MRI;

% reposition SPM window
set(st.fig,'Units','normalized');
spmwindow_position = get(st.fig,'Position');
spmwindow_position(1:2) = [0.48 0]; % place SPM window at top right
set(st.fig,'Position',spmwindow_position,'Resize','on');
% set(st.fig,'Units','points');  % set units back to pixels so fonts don't go nuts

if(0)  % stub for plotting time series in SPM window
    axialpos = get(st.vols{1}.ax{1}.ax,'Position'); % position of axial view
    handles.nut_ts_axes = axes('position',[(axialpos(1)+axialpos(3)+.05) (axialpos(2)+.02) (axialpos(3)-.05) (axialpos(4)-.05)]);
end

waitbar(0.5,barhandle);

% following jazz is needed because SPM's blob freaks out when given
% negative coordinates, and needs a dilation to know these are big voxels
% furthermore, coords to SPM must be integers.
% also, spm_orthviews will blow up your computer if it encounters singles
rivets.voxelsMRI = double(nut_meg2mri(beam.voxels));
[rivets.voxelsblob,rivets.blob2mri_tfm]=nut_voxels2blob(beam);
% rivets.voxelsMRI = double(nut_coordtfm(beam.voxels,beam.coreg.meg2mri_tfm));

% fixes silly bug if there's only a single time-frequency window
if(size(beam.s{1},1)==1&&length(size(beam.s{1}))==2)%transpose not defined for ND arrays
    beam.s{1}=beam.s{1}';
end

% for some phucked reason, have to initialize plot with a point for
% nut_plotbeamtf to work properly. 'tis phucked, i tells ya'!!!
axes(handles.nut_ts_axes); plot([0 0]);

set(handles.nut_view_time_text,'String',beam.timepts(1));
set(handles.nut_view_time2_text,'String',beam.timepts(end));

global interval_select_handle; % used by ts_axes ButtonDownFcn
interval_select_handle = @interval_button_Callback;

waitbar(0.7,barhandle);
nut_spmfig_setup;
rivets.ts_refresh_handle = @plot_ts;  %% pass handle since technically it's a private function

waitbar(1,barhandle);
close(barhandle);

if(0 && isfield(nuts,'meg'))
    set(handles.nut_megts_axes,'Visible','on');
    axes(handles.nut_megts_axes);
    rivets.meg.data=nuts.meg.data;
    plot(beam.timepts,1e15*rivets.meg.data);
    axis tight; xlabel('Time (ms)'); ylabel('Magnetic Field (fT)');
    set(handles.nut_display_rawmeg_button,'Visible','Off');
else
    set(handles.nut_display_rawmeg_button,'Visible','On');
end

% store SPM's default colormap
rivets.spmmap = get(st.fig,'ColorMap');

tmp = hot(64 + 16);  tmp = tmp([1:64] + 16,:);
rivets.grayhot = [gray(64); tmp]; % SPM's gray-hot colormap
rivets.graygray = [.8*gray(64);gray(64)]; % gray-gray colormap

% gray-gray colormap for bipolar maps
graysymmetric = flipud(gray(32));
graysymmetric = [graysymmetric(1:(end-1),:); gray(32)];
rivets.graygraybipolar = [.8*gray(64); .8*graysymmetric];

% gray-jet colormap
rivets.grayjet = [gray(64);jet(64)];

% hot-symmetric
tmpsquish=tmp(1:2:end,:);
hotsymmetric = [flipud(tmpsquish); tmpsquish];
rivets.grayhotsymmetric=[gray(64);hotsymmetric];

rivets.freqselect = 1;
rivets.timeselect = 1;
% rivets.setscale=0;
rivets.tfwindoworient_all=0;

% set to "Active Power" contrast type for evoked time series with no
% control/baseline period
contrastcontents = get(handles.nut_contrastselect_menu,'String');
contrastvalue =  find(strcmp(contrastcontents,'Active Power'));
if ~isempty(contrastvalue)
    set(handles.nut_contrastselect_menu,'Value',contrastvalue);
else
    set(handles.nut_contrastselect_menu,'Value',1);
end
nut_contrastselect_menu_Callback(hObject,[],handles);

% Deal with pointer file-----------------
if ~isfield(rivets,'voi')   % initialize voi structure, if it doesn't exist
    nut_clear_subjcond_Callback([],[],handles,false);
else                        % clear voi structure, if necessary
    asubj=get(handles.nut_getset_patientnr,'String'); asubj=str2double(asubj{get(handles.nut_getset_patientnr,'Value')});
    acond=get(handles.nut_getset_conditionnr,'String'); acond=str2double(acond{get(handles.nut_getset_conditionnr,'Value')});
    f=find(rivets.voi.subjnr==asubj & rivets.voi.condnr==acond);    % if current s_beam file is not in rivets.voi.pathnames, voi structure is cleared
    if ( isempty(f) || ( ~strcmpi(rivets.beampath,rivets.voi.pathnames{f}) && ~strcmpi([rivets.beampath '.mat'],rivets.voi.pathnames{f}) ) ) 
        nut_clear_subjcond_Callback([],[],handles,false);
    end
    clear f asubj acond
end
nvoi=length(rivets.voi.coords);
if nvoi>0
    avoi=get(handles.nut_getset_voinumber,'Value');
    set(handles.nut_getset_voinumber,'string',[{rivets.voi.coords.label}';'New'])
    if avoi>nvoi, avoi=1; set(handles.nut_getset_voinumber,'Value',avoi); end
    set(handles.nut_voilabel,'String',rivets.voi.coords(avoi).label)
    nut_displayvoi_Callback([],[],handles);     % Display selected voxel
    set([handles.nut_deletevoi handles.nut_displayvoi handles.nut_modifyvoi],'Enable','on')
end
asubj=get(handles.nut_getset_patientnr,'Value');
if (asubj>length(rivets.voi.subjectmaxtf) || rivets.voi.subjectmaxtf(asubj)==0)
    rivets.voi.subjectmaxtf(asubj)=rivets.maxtf;
else
    if rivets.maxtf>rivets.voi.subjectmaxtf(asubj)
        rivets.voi.subjectmaxtf(asubj)=rivets.maxtf;
    elseif rivets.maxtf<rivets.voi.subjectmaxtf(asubj)
        rivets.maxtf=rivets.voi.subjectmaxtf(asubj);
        nut_colorscale_menu_Callback([],[],handles);
    end
end


%%---------------------------------------------------
function varargout = nut_results_viewer_OutputFcn(hObject, eventdata, handles)
if(isfield(handles,'output'))  % prevents crashing if s_beam loading cancelled
    varargout{1} = handles.output;
end

%%---------------------------------------------------
function nut_blobstyle_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
% if ispc
%     set(hObject,'BackgroundColor','white');
% else
%     set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
% end

%%---------------------------------------------------
function nut_blobstyle_Callback(hObject, eventdata, handles)
% --- Executes on selection change in nut_blobstyle.
paint_activation(hObject,eventdata,handles);


%%-----------------------------------
function nut_time_text_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
% if ispc
%     set(hObject,'BackgroundColor','white');
% else
%     set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
% end


%%-----------------------------------------------
function nut_time_text_Callback(hObject, eventdata, handles)

global rivets beam
time = str2double(get(handles.nut_time_text,'String'));
set(handles.nut_time2_text,'String',time);
timept = dsearchn(beam.timepts,time);

rivets.timeselect = timept;
nut_refreshdisplay(handles);

%%-----------------------------------------------
function nut_time2_text_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
% if ispc
%     set(hObject,'BackgroundColor','white');
% else
%     set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
% end

%%--------------------------------------------
function nut_time2_text_Callback(hObject, eventdata, handles)

global rivets beam
time2 = str2double(get(handles.nut_time2_text,'String'));
timept2 = dsearchn(beam.timepts,time2);
%set(handles.nut_time2_text,'String',timept2);

rivets.timeselect = round((rivets.timeselect+timept2)/2);
nut_refreshdisplay(handles);

% ---------------------------------------------
function nut_threshold_menu_Callback(hObject, eventdata, handles)
% hObject    handle to nut_threshold_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns nut_threshold_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nut_threshold_menu
global beam

thmode=get(hObject,'string');
thmode=thmode{get(hObject,'Value')};

if any(strcmp(thmode,{'Power' 'Coherence' 'Corr Coeff (r)'}))
    set(handles.nut_thresholdpos_text,'string','0','Enable','on')
    set(handles.nut_thresholdneg_text,'string','0','Enable','on')
elseif strcmp(thmode,'SnPM T')
    if isfield(beam,'snpm')
        set(handles.nut_thresholdpos_text,'string','0','Enable','on')
        set(handles.nut_thresholdneg_text,'string','0','Enable','on')
    else, errordlg('No SnPM statistics found in current s_beam file.'), return
    end
elseif strncmp(thmode,'SnPM p',6) 
    if isfield(beam,'snpm')
        set([handles.nut_thresholdneg_text handles.nut_thresholdpos_text],'string','1','Enable','off')
        if isfield(beam.snpm,'p_uncorr_pos')
            set(handles.nut_thresholdpos_text,'string','0.05','Enable','on')
        end
        if isfield(beam.snpm,'p_uncorr_neg')
            set(handles.nut_thresholdneg_text,'string','0.05','Enable','on')
        end
    else, errordlg('Could not find SnPM statistics in current s_beam file.'), return
    end   
elseif strcmp(thmode, 'T-test T')
    if isfield(beam,'ttest')
        set(handles.nut_thresholdpos_text,'string','0','Enable','on')
        set(handles.nut_thresholdneg_text,'string','0','Enable','on')
    else, errordlg('No t-tests found in current s_beam file.'), return
    end
elseif strncmp(thmode,'T-test p',8)
    if isfield(beam,'ttest')
        set(handles.nut_thresholdpos_text,'string','','Enable','off')
        set(handles.nut_thresholdneg_text,'string','','Enable','off')
        if ~isempty(strmatch(beam.ttest.tail,{'pos','both'}))
            set(handles.nut_thresholdpos_text,'string','0.05','Enable','on')
            %set(handles.nut_thresholdpos_text,'string',num2str(beam.ttest.cutoff,'%1.5f'),'Enable','on')
        end
        if ~isempty(strmatch(beam.ttest.tail,{'neg','both'}))
            set(handles.nut_thresholdneg_text,'string','0.05','Enable','on')                
            %set(handles.nut_thresholdneg_text,'string',num2str(beam.ttest.cutoff,'%1.5f'),'Enable','on')
        end
    else, errordlg('No t-tests found in current s_beam file.'), return
    end
elseif strcmp(thmode,'p')
    if isfield(beam,'normcdf')
        set(handles.nut_thresholdpos_text,'string','','Enable','off')
        set(handles.nut_thresholdneg_text,'string','','Enable','off')
        if ~isempty(strmatch(beam.normcdf.tail,{'pos','both'}))
            set(handles.nut_thresholdpos_text,'string',num2str(beam.normcdf.cutoff,'%1.5f'),'Enable','on')
        end
        if ~isempty(strmatch(beam.normcdf.tail,{'neg','both'}))
            set(handles.nut_thresholdneg_text,'string',num2str(beam.normcdf.cutoff,'%1.5f'),'Enable','on')
        end
    else, errordlg('No p values found.'), return
    end
elseif strncmp(thmode,'Corr p',6)  
    if isfield(beam,'corr')
        set(handles.nut_thresholdpos_text,'string','','Enable','off')
        set(handles.nut_thresholdneg_text,'string','','Enable','off')
        if ~isempty(strmatch(beam.corr.tail,{'pos','both'}))
            set(handles.nut_thresholdpos_text,'string','0.05','Enable','on')
        end
        if ~isempty(strmatch(beam.corr.tail,{'neg','both'}))
            set(handles.nut_thresholdneg_text,'string','0.05','Enable','on')
        end
    else, errordlg('No correlation stats found.'), return
    end    
end

nut_refreshdisplay(handles);


% ---------------------------------------------
function nut_threshold_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_threshold_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%-----------------------------------------------
function nut_thresholdpos_text_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
% if ispc
%     set(hObject,'BackgroundColor','white');
% else
%     set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
% end


%%-----------------------------------------------
function nut_thresholdpos_text_Callback(hObject, eventdata, handles)
% global timevector timept timept2  ts_axes_obj;
nut_refreshdisplay(handles);


function nut_thresholdneg_text_Callback(hObject, eventdata, handles)
% hObject    handle to nut_thresholdneg_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_thresholdneg_text as text
%        str2double(get(hObject,'String')) returns contents of nut_thresholdneg_text as a double
nut_refreshdisplay(handles);


% --- Executes during object creation, after setting all properties.
function nut_thresholdneg_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_thresholdneg_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%----------------------------------------------------------
function nut_do_threshold(thmode,handles)

global rivets beam

if isempty(thmode) || isnumeric(thmode)
    thmode=get(handles.nut_threshold_menu,'string');
    thmode=thmode{get(handles.nut_threshold_menu,'Value')};
end

ssz=size(rivets.s);
thpos=str2double(get(handles.nut_thresholdpos_text,'String'));
thneg=str2double(get(handles.nut_thresholdneg_text,'String'));
%rivets.threshold = zeros(size(beam.voxels,1),2);
rivets.threshold=false([ssz(1:2) 2]);

switch thmode
    case {'Power' 'Corr Coeff (r)' 'Coherence'}
        rivets.threshold(:,:,1) = (rivets.scalefactor*rivets.s(:,:,rivets.freqselect) > thpos);
        rivets.threshold(:,:,2) = (rivets.scalefactor*rivets.s(:,:,rivets.freqselect) < thneg);

    case {'SnPM T' 'SnPM F' 'SnPM R'}
        mul=ssz./size(beam.snpm.T);  mul=mul(1:2);
        rivets.threshold(:,:,1) = (repmat(beam.snpm.T(:,:,rivets.freqselect),mul) > thpos);
        rivets.threshold(:,:,2) = (repmat(beam.snpm.T(:,:,rivets.freqselect),mul) < thneg);
        %(find(beam.snpm.T(:,min(rivets.timeselect,size(beam.snpm.T,2)),rivets.freqselect) < thpos),1) = Inf;
        %(find(beam.snpm.T(:,min(rivets.timeselect,size(beam.snpm.T,2)),rivets.freqselect) > thneg),2) = -Inf;
        
    case 'T-test T'
        mul=ssz./size(beam.ttest.T);  mul=mul(1:2);
        rivets.threshold(:,:,1) = (repmat(beam.ttest.T(:,:,rivets.freqselect),mul) > thpos);
        rivets.threshold(:,:,2) = (repmat(beam.ttest.T(:,:,rivets.freqselect),mul) < thneg);

    case 'p'
        mul=ssz./size(beam.normcdf.p_uncorr);  mul=mul(1:2);
        if ~isempty(strmatch(beam.normcdf.tail,{'pos','both'}))
            rivets.threshold(:,:,1)=(repmat(beam.normcdf.p_uncorr(:,:,rivets.freqselect),mul) < thpos);
        end
        if ~isempty(strmatch(beam.normcdf.tail,{'neg','both'}))
            rivets.threshold(:,:,2)=(repmat(beam.normcdf.p_uncorr(:,:,rivets.freqselect),mul) < thneg);
        end
    
    otherwise
        
        if ~isempty(strfind(thmode,' p '))  % deal with all the stat possibilities
            [meth,field] = strtok(thmode);
            meth = lower(strrep(meth,'-',''));
            field = strrep(field(2:end),' ','_');
            
            switch meth
                case 'snpm'
                    mul=ssz./size(beam.snpm.T);  mul=mul(1:2);
                    if isfield(beam.snpm,[field '_pos'])
                        %f = find(beam.snpm.p_corr_pos(:,min(rivets.timeselect,size(beam.snpm.T,2)),rivets.freqselect) <= thpos);
                        %rivets.threshold(f,1) = 0;
                        rivets.threshold(:,:,1)=(repmat(beam.snpm.([field '_pos'])(:,:,rivets.freqselect),mul) < thpos);
                    end
                    if isfield(beam.snpm,[field '_neg'])
                        rivets.threshold(:,:,2)=(repmat(beam.snpm.([field '_neg'])(:,:,rivets.freqselect),mul) < thneg);
                    end
                
                case 'ttest'
                    mul=ssz./size(beam.(meth).(field));  mul=mul(1:2);
                    if ~isempty(strmatch(beam.(meth).tail,{'pos','both'}))
                        rivets.threshold(:,:,1)=( (repmat(beam.(meth).(field)(:,:,rivets.freqselect),mul) < thpos) & (repmat(beam.(meth).T(:,:,rivets.freqselect),mul) > 0) );
                    end
                    if ~isempty(strmatch(beam.(meth).tail,{'neg','both'}))
                        rivets.threshold(:,:,2)=( (repmat(beam.(meth).(field)(:,:,rivets.freqselect),mul) < thneg) & (repmat(beam.(meth).T(:,:,rivets.freqselect),mul) < 0) );
                    end
                    
                case 'corr'
                    mul=ssz./size(beam.(meth).(field));  mul=mul(1:2);
                    if ~isempty(strmatch(beam.(meth).tail,{'pos','both'}))
                        rivets.threshold(:,:,1)=( (repmat(beam.(meth).(field)(:,:,rivets.freqselect),mul) < thpos) & (repmat(beam.s{1}(:,:,rivets.freqselect),mul) > 0) );
                    end
                    if ~isempty(strmatch(beam.(meth).tail,{'neg','both'}))
                        rivets.threshold(:,:,2)=( (repmat(beam.(meth).(field)(:,:,rivets.freqselect),mul) < thneg) & (repmat(beam.s{1}(:,:,rivets.freqselect),mul) < 0) );
                    end              
                    
                otherwise
                    errordlg('I do not know what to do with your threshold method %s...',thmode)
            end
            
        else
            errordlg('I do not know what to do with your threshold method %s...',thmode)
        end
end

%%-----------------------------------------------
%function nut_ts_button_Callback(hObject, eventdata, handles)
% --- Executes on button press in nut_ts_button.
%plot_ts;
%return;

%%-----------------------------------------------
function nut_ts_axes_ButtonDownFcn(hObject, eventdata, handles)
% --- Executes on button press in nut_interval_button.

pos = get(handles.nut_ts_axes,'CurrentPoint');
time = pos(1,1);
yselect = pos(1,2); % y = freq for tf spectrograms; y = intensity/amplitude for time series

% FIXME: this should be equivalent to interval_button_Callback
% then we can remove the interval button!!!
switch(get(gcf,'SelectionType'))
    case 'normal' % left click
        nut_tfselect(time,yselect,handles);
    case 'alt' % right click
        nut_tfselect(time,yselect,handles);
    case 'extend' %middle click
        nut_tfselect(time,yselect,handles);
    case 'open' % double click
        nut_tfselect(time,yselect,handles);
end
        

%-------------------------------------------------------------
function nut_tfselect(time,freqselect,handles)
global rivets beam;
timept = dsearchn(beam.timepts,time);

set(handles.nut_time_text,'String',beam.timepts(timept));
set(handles.nut_time2_text,'String',beam.timepts(timept));
rivets.timeselect = timept;

if rivets.timefmode
    freqpt = find(all([freqselect > beam.bands(:,1) freqselect < beam.bands(:,2)],2));
    if isempty(freqpt), return, end
    %rivets.s_beam = rivets.s(:,:,freqpt);
    rivets.freqselect = freqpt;
end

nut_refreshdisplay(handles);
return;

%%-----------------------------------------------
function nut_interval_button_Callback(hObject, eventdata, handles)
% --- Executes on button press in nut_interval_button.
global rivets beam;

% initialize times
% set(handles.nut_time_text,'String',beam.timepts(1));
% set(handles.nut_time2_text,'String',beam.timepts(end));
% plot_ts; % erase any previous interval lines

if rivets.timefmode
    [time,ypos,button]=ginput(1);
    timept = dsearchn(beam.timepts,time);

    set(handles.nut_time_text,'String',time);
    set(handles.nut_time2_text,'String',time);
    plot_ts;

    freqselect = abs(ypos);
    freqpt = find( beam.bands(:,1)<=freqselect & beam.bands(:,2)>=freqselect );
    if isempty(freqpt), return, end

    rivets.freqselect = freqpt;
    rivets.timeselect = timept;
    %rivets.s_beam = rivets.s(:,:,freqpt);
else
    [time,crap,button]=ginput(1);
    set(handles.nut_time_text,'String',time);
    set(handles.nut_time2_text,'String',time);
    plot_ts;
    timept = dsearchn(beam.timepts,time);

    % now for the second selection
    button = 0;
    while(button~=1 && button~=3)
        [time2,ypos,button]=ginput(1);
        if (button==3)   % just use first time point if right button clicked
            timept2 = timept;
            %elseif (button==2)  % if middle button clicked, set threshold
            %    threshold = abs(ypos);
            %    set(handles.nut_threshold_text,'String',threshold);
            %    rivets.threshold = [repmat(threshold,size(beam.voxels,1),1) repmat(-threshold,size(beam.voxels,1),1)];
            %    plot_ts;
        else  % select second time point for other buttons
            set(handles.nut_time2_text,'String',time2);
            plot_ts;
            timept2 = dsearchn(beam.timepts,time2);
        end
    end
    rivets.timeselect = round((timept+timept2)/2);
end

nut_refreshdisplay(handles);

%%-----------------------------------------------
function nut_view_interval_button_Callback(hObject, eventdata, handles)
% --- Executes on button press in nut_interval_button.
global rivets beam;

% initialize times
set(handles.nut_view_time_text,'String',beam.timepts(1));
set(handles.nut_view_time2_text,'String',beam.timepts(end));
plot_ts; % erase any previous interval lines

[time,crap,button]=ginput(1);
timept = dsearchn(beam.timepts,time);
set(handles.nut_view_time_text,'String',timept);
% set(handles.nut_view_time2_text,'String',beam.timepts(end));

% now for the second selection
button = 0;
while(button~=1 && button~=3)
    [time2,ypos,button]=ginput(1);
    if (button==3)   % just use first time point if right button clicked
        timept2 = timept;
        %     elseif (button==2)  % if middle button clicked, set threshold
        %         threshold = ypos;
        %         set(handles.nut_thresholdpos_text,'String',threshold);
    else  % select second time point for other buttons
        timept2 = dsearchn(beam.timepts,time2);
    end
end

set(handles.nut_view_time2_text,'String',timept2);
plot_ts;

%%---------------------------------------------------
function nut_animate_button_Callback(hObject, eventdata, handles)
% --- Executes on clicking of nut_animate_button button.
% retrieve selected coordinates, and inverse transform to MEG mm
global rivets beam st;
persistent stopflag;

slices = false;

switch(get(hObject,'String'))
    case 'Animate'
        orthview = 3; % 1 selects axial, 2 coronal, 3 sagittal

        stepsize = 1;
        if(get(handles.nut_savemovie_box,'Value'))
            kubrick = true;
        else
            kubrick = false;
        end

        if(kubrick) % save movie mode
            if(ispc)
                errordlg('Sorry, you can''t do this with Windows yet. You might consider Linux. :)');
                return
            end

            status = unix('which convert');
            if(status) % if status is nonzero, "convert" failed and probably is not installed
                errordlg('Sorry, this feature requires the convert function from ImageMagick (http://www.imagemagick.org/)');
                return
            end
            if ~exist('movie','dir'), mkdir movie, end
        end
        stopflag = false;  % set stopflag to false so we can slide
        set(hObject,'String','Stop Animation');

        if(0)  % moves time series to SPM window... unfortunately, SPM's blob colorbar craps all over it
            set(handles.nut_ts_axes,'Units','Normalized');
            set(handles.nut_ts_axes,'Position',[0.55 0.45 0.4 0.25],'Parent',st.fig); % move time series to SPM window
        elseif(isempty(st.vols{1})) % for brain rendering...
            figure(st.fig);
            %texthandle = text(1,1,'0.0 ms','FontSize',14,'HorizontalAlignment','right');
        else % so for now, let's just add a box to indicate time
            imagewidth = get(st.vols{1}.ax{1}.d,'XData');
            %texthandle = xlabel(st.vols{1}.ax{orthview}.ax,'0.0 ms','FontSize',14,'HorizontalAlignment','right','Position',[imagewidth(2) -15 1]);
        end

        framenum = 0;
        % make sure figures export correctly:
        set(st.fig,'PaperPositionMode','auto');
        set(st.fig,'PaperOrientation','portrait');

        MEGpos = (spm_orthviews('pos'))';
        % figure out which MEG voxel is closest, then plot the time series for it
        if(isfield(st.vols{1},'premul'))
            MEGvoxelindex = dsearchn(nut_coordtfm(rivets.voxelsMRI,st.vols{1}.premul),MEGpos);
        else % st.vols{1}.premul may not exist when displaying on rendered brain
            MEGvoxelindex = dsearchn(rivets.voxelsMRI,MEGpos);
        end

        if(kubrick & slices)
            tempfig=figure;
            set(st.vols{1}.ax{orthview}.ax,'Units','pixels');

            for f=1:4
                nut_tfselect(beam.timepts(1),mean(beam.bands(f,1:2)),handles);
                copyh(f) = copyobj(st.vols{1}.ax{orthview}.ax,tempfig);
                pos = get(copyh(f),'Position');
                set(copyh(f),'Position',[(f-1)*(pos(3)+25) 50 pos(3) pos(4)], 'ButtonDownFcn',[]);
                copykids = get(copyh(f),'Children');
                set(copykids,'DeleteFcn',[]); % kill kids without killing twins
                title(copyh(f),[num2str(beam.bands(f,1)) '-' num2str(beam.bands(f,2)) ' Hz'],'FontSize',24,'Color',[1 1 0])
            end
            %             cbar = copyobj(rivets.cbar,tempfig);
            set(tempfig,'InvertHardCopy','off','PaperPositionMode','auto','Color',[0 0 0],'Colormap',get(st.fig,'Colormap'),'Position',[1 50 4*(pos(3)+25) pos(4)+100]);

        end

        %         timezero = dsearchn(beam.timepts,0);

        %         for i = 1:stepsize:length(beam.timepts)
        for i = 1:stepsize:length(beam.timepts)
            % someday computers will be fast enough that you'd want to pause
            % between frames... but not as of Apr. 2003!
            % pause(0.1);

            set(handles.nut_time_text,'String',beam.timepts(i));
            set(handles.nut_time2_text,'String',beam.timepts(i));
            framenum = framenum + 1;
            %             plot_ts;
            %             paint_activation(hObject,eventdata,handles);
            f = rivets.freqselect;
            nut_tfselect(beam.timepts(i),mean(beam.bands(f,1:2)),handles);
            texthandle = annotation('textbox',[0 0 1 0.5],'string',[num2str(beam.timepts(i),'%0.1f') ' ms'],'Color',[1 1 0],'FontSize',24,'horizontalalignment','center','verticalalignment','middle');

            if(kubrick)
                if(slices)
                    % save as /tmp/frame0.###
                    delete(copyh);
                    for f=1:4
                        nut_tfselect(beam.timepts(i),mean(beam.bands(f,1:2)),handles);
                        copyh(f) = copyobj(st.vols{1}.ax{orthview}.ax,tempfig);
                        pos = get(copyh(f),'Position');
                        set(copyh(f),'Position',[(f-1)*(pos(3)+25) 50 pos(3) pos(4)], 'ButtonDownFcn',[]);
                        copykids = get(copyh(f),'Children');
                        set(copykids,'DeleteFcn',[]); % kill kids without killing twins
                        title(copyh(f),[num2str(beam.bands(f,1)) '-' num2str(beam.bands(f,2)) ' Hz'],'FontSize',24,'Color',[1 1 0])
                    end
                    %                 copyh = copyobj(st.vols{1}.ax{3}.ax,tempfig);
                    %                 pos = get(copyh,'Position')
                    %                 set(copyh,'Position',[1 1 pos(3) pos(4)], 'ButtonDownFcn',[]);
                    %                 copykids = get(copyh,'Children');
                    %                 set(copykids,'DeleteFcn',[]); % kill kids without killing twins
                    % saveas(tempfig,['movie/MEGframe' num2str(framenum,'%0.3d')],'tiff')
                    print(tempfig,'-dtiff','-painters','-r72',['movie/MEGframe' num2str(framenum,'%0.3d') '.tiff']);
                else
                    drawnow
                    print(st.fig,'-dtiff','-painters','-r72',['movie/MEGframe' num2str(framenum,'%0.3d') '.tiff']);
                end
            end
            delete(texthandle);
            if(stopflag) % break loop when stop button is pressed
                return
            end
        end

        if kubrick
            [moviepath,moviename] = fileparts(rivets.beampath);
            moviepath = fullfile(moviepath,[moviename '.mpg']);
            watchon;
            disp('Please wait... piecing movie together...');
            if rivets.displaymode==1
                unix(sprintf('convert -crop -0-300 -delay %d movie/MEGframe*.tiff %s',round(100/beam.srate),moviepath));
            else
                unix(sprintf('convert -crop -0+422 -delay %d movie/MEGframe*.tiff %s',30,moviepath));     %round(100/beam.srate)
            end
            disp(['Movie saved as ' moviepath]);
            status = unix('which mplayer');
            if(status) % if status is nonzero, "gmplayer" failed and probably is not installed
                warning('MPEG not displayed since mplayer not found. (http://www.mplayerhq.hu)');
            else
                unix(['mplayer -noaspect -loop 0 ' moviepath '&']);
            end
            %             delete(fullfile(tempdir,'MEGframe*.tiff'));
            watchoff;
        end
    case 'Stop Animation'
        stopflag = true;
        set(hObject,'String','Animate');
        return
end

set(hObject,'String','Animate');


%%-----------------------------------------------
function nut_results_viewer_CloseRequestFcn(hObject, eventdata, handles)
% --- Executes during object deletion, before destroying properties.
global rivets st beam corrs;

spm_image('init',beam.coreg.mripath); % reload structural MRI;
spm_figure('ColorMap','gray'); % reset original colormap
if(isfield(corrs,'fig') && ishandle(corrs.fig)), delete(corrs.fig), end
clear global beam rivets corrs
delete(hObject);


%%------------------------------------------
function nut_ERP_Overlay_button_Callback(hObject, eventdata, handles)
% --- Executes on button press in nut_ERP_Overlay_button.
global rivets beam;
uiload;

prompt   = 'Which electrode?';
title    = 'Select ERP';
lines = 1;
def     = {'1'};
answer   = inputdlg(prompt,title,lines,def);

axes(handles.nut_ts_axes); hold on;
time = -225+.4992511:(0.4992511):450;
plot(time,erpdat3(str2num(cell2mat(answer)),:).^2/3600,'r');
hold off;
return;

%%------------------------------------------
function nut_saverender_button_Callback(hObject, eventdata, handles)
% --- Executes on button press in nut_saverender_button.
% threshold = str2double(get(handles.nut_thresholdpos_text,'String'));
% nut_create_spmimage;

global rivets beam
if(get(hObject,'Value'))
    if(~isfield(beam.coreg,'brainrender_path'))
        if strcmp(spm('ver'),'SPM2')
            beam.coreg.brainrender_path = spm_get(1,'render*.mat','Render file',fullfile(spm('Dir'),'rend'));
        elseif ( strcmp(spm('ver'),'SPM8') || strcmp(spm('ver'),'SPM8b') )
            beam.coreg.brainrender_path = spm_select(1,'mesh','Render file',[],fullfile(spm('Dir'),'rend'));
        end
    end
    rivets.displaymode=2;
    set([handles.nut_maxtime_button handles.nut_maxactivation_button handles.nut_maxvoxel_button],'Enable','off')
    nut_render_meg(beam.coreg.brainrender_path);
else
    rivets.displaymode=1;
    set([handles.nut_maxtime_button handles.nut_maxactivation_button handles.nut_maxvoxel_button],'Enable','on')
    spm_image('init',beam.coreg.mripath); % load/reload structural MRI;
    nut_spmfig_setup
    plot_ts; % erase any previous interval lines
    paint_activation(hObject,eventdata,handles)
end


%%------------------------------------------
function activation_bpf_Callback(hObject, eventdata, handles)
% --- Executes on button press in activation_bpf.
plot_ts;
return;


%%--------------------------------------------
function tmp = nut_calc_orientations(handles,contrastselect)

global rivets beam

% if size(beam.s{1},4)>1
% if( ~rivets.timefmode && ndims(beam.s{1})>3 )   
%     if rivets.timefmode==1
%         if rivets.tfwindoworient_all % orientation based on all windows
%             disp('please code me')
%         else  % orientation based each window separately (with average of control?)
%             disp('please code me too')
%         end
%     elseif rivets.timefmode==0
if size(beam.s{1},3)>1
    error('how come rivets.timefmode=0 if you''ve got multiple frequency bands?')
end
orientcontents=get(handles.nut_orientated_menu,'String');
orienttype=orientcontents{get(handles.nut_orientated_menu,'Value')};
cond1=str2double(get(handles.nut_condcontrast1_text,'String'));
cond2=str2double(get(handles.nut_condcontrast2_text,'String'));

switch orienttype
    case 'Power All Orientations'
        switch contrastselect
            case 'Active Power'

                stmp=zeros(size(beam.s{cond1},1),size(beam.s{cond1},2),size(beam.s{cond1},3),'single');
                for ii=1:size(beam.s{cond1},4)
                    stmp=stmp+beam.s{cond1}(:,:,:,ii).^2;
                end
                tmp.sa=stmp;
                clear stmp;
                if size(beam.s,2)>1
                    stmp=zeros(size(beam.s{cond2},1),size(beam.s{cond2},2),size(beam.s{cond2},3),'single');
                    for ii=1:size(beam.s{cond2},4)
                        stmp=stmp+beam.s{cond2}(:,:,:,ii).^2;
                    end
                    tmp.sc=stmp;
                    clear stmp;
                else
                    tmp.sc=0;
                end
                if size(beam.s,2)>3 
                    warning('using beam.s{3} as noise estimate')
                    stmp=zeros(size(beam.s{3},1),size(beam.s{3},2),size(beam.s{3},3),'single');
                    for ii=1:size(beam.s{3},4)
                        stmp=stmp+beam.s{3}(:,:,:,ii).^2;
                    end
                    clear noise
                    noise=stmp;
                    clear stmp;
                end
                %                 rivets.minblob = min(rivets.s(:));
                %                 rivets.maxblob = max(rivets.s(:));
                %                 rivets.scalefactor = 1000/max(abs(rivets.maxblob),abs(rivets.minblob));
                rivets.display_s_II_flag = false;
            case 'Control Power'
                stmp=zeros(size(beam.s{cond2},1),size(beam.s{cond2},2),size(beam.s{cond2},3),'single');
                for ii=1:size(beam.s{cond2},4)
                    stmp=stmp+beam.s{cond2}(:,:,:,ii).^2;
                end
                tmp.sa=stmp;
                clear stmp;
            otherwise
                switch get(handles.nut_noisetype,'Value')
                    case 1
                        [u,s]=svd(beam.params.TotC);
                        noisecov=s(end,end);
                    case 2
                        noisecov=beam.params.Czz;
                end
                if length(beam.W)==3 % beam.W{1} is active, beam.W{2} is control, and beam.W{3} is ave of both
                    activewin=dsearchn(beam.timepts,beam.params.preprocessing.bf_timeinterval');
                    controlwin=dsearchn(beam.timepts,beam.params.preprocessing.bf_ptimeinterval');
                    rivets.timeselect=1;
                    %                     if size(beam.W{3},3)==1 % 1 orientation
                    %                         %             pt = (diag(beam.W{3}'*beam.params.Rzz1*beam.W{3})-diag(beam.W{3}'*beam.params.Czz*beam.W{3}))./(2*diag(beam.W{3}'*beam.params.Czz*beam.W{3}));
                    %                         %             pt = (diag(beam.W{3}'*beam.params.Rzz1*beam.W{3})-diag(beam.W{3}'*beam.params.Czz*beam.W{3}))./(2*diag(beam.W{3}'*sig*beam.W{3}));
                    %                         tmp.sa=squeeze(mean(beam.s{1}(:,activewin).^2,2));
                    %                         tmp.sc=squeeze(mean(beam.s{2}(:,controlwin).^2,2));
                    % %                         pt = (diag(beam.W{3}'*beam.params.Rzz1*beam.W{3})-diag(beam.W{3}'*beam.params.Czz*beam.W{3}))./(2*diag(beam.W{3}'*noise*beam.W{3}));
                    %                     else % 2 or 3 orientations

                    rivets=rmfield(rivets,'s'); % clearing for memory reasons
                    disp('you might need more than 2GB RAM for this part');
                    for ii=1:size(beam.W{3},2) % loop over voxels
                        apower(:,ii)=diag(squeeze(beam.W{3}(:,ii,:))'*beam.params.Rzz1*squeeze(beam.W{3}(:,ii,:)));
                        cpower(:,ii)=diag(squeeze(beam.W{3}(:,ii,:))'*beam.params.Czz*squeeze(beam.W{3}(:,ii,:)));
%                                 tmp.sa(ii,:)=diag(squeeze(beam.W{3}(:,ii,:))'*beam.params.Rzz1*squeeze(beam.W{3}(:,ii,:)));
%                                 tmp.sc(ii,:)=diag(squeeze(beam.W{3}(:,ii,:))'*beam.params.Czz*squeeze(beam.W{3}(:,ii,:)));
                        %                 anoise(:,ii)=diag(squeeze(beam.W{3}(:,ii,:))'*beam.params.Rzz1*squeeze(beam.W{3}(:,ii,:))); %wrong!
                        %                 cnoise(:,ii)=diag(squeeze(beam.W{3}(:,ii,:))'*beam.params.Czz*squeeze(beam.W{3}(:,ii,:)));
                        %                             npower(:,ii)=diag(squeeze(beam.W{3}(:,ii,:))'*noise*squeeze(beam.W{3}(:,ii,:)));
                        noise(:,ii)=diag(squeeze(beam.W{3}(:,ii,:))'*noisecov*squeeze(beam.W{3}(:,ii,:)));
                        %         for ii=1:size(beam.W{3},2)

                    end
                    if size(apower,2)>3
                        apower=apower';
                        cpower=cpower';
                        noise=noise';
                    end
                    tmp.sa=sum(apower,2);
                    tmp.sc=sum(cpower,2);
                    tmp.n=sum(noise,2);
                    if size(tmp.sa,2)>1
                        tmp.sa=tmp.sa';
                        tmp.sc=tmp.sc';
                        tmp.n=tmp.n';
                    end
                    %             pt = ((sum(apower')-sum(cpower'))./(sum(anoise')+sum(cnoise')))';
                    %                         pt = ((sum(apower')-sum(cpower'))./(2*sum(npower')))';
                    %                     end
                end
        end
case 'Amplitude: Constant Orientation'
        rivets.display_s_II_flag = true;
        %         if(~isfield(beam,'s_II'))  % only do this the first time cuz it's slower than the line at the DMV
        %         commented out line above, since s_II might exist from a different contrast above; -jz
        disp('Please wait... computing orientations...');
        disp('todo: orientation based on active only, or active + control?')
        %             if(all(beam.s_ph(:)==0))
        %                 rivets.s_II = beam.s_th;
        %                 rivets.s_perp = beam.s_ph;
        %             else
        %                 switch(5)
        if size(beam.s{1},4)==3
            rivets=rmfield(rivets,'s');
            disp('this is a slow for-loop over voxels')
            for jj=1:size(beam.s{1},1) % help! do you know how to make this not a for-loop?
                S=cov(squeeze(beam.s{1}(jj,:,1,:)));
                % [uu,dd]=svd(inv(S)); %sarang trusts eig over svd for consistency
                [vv,dd]=eig(inv(S));
                [jnk,ord]=sort(diag(dd),'descend');
                rotinv=vv(:,ord)';
                tmp.sa(jj,:,1,:)=flipdim((rotinv*squeeze(beam.s{1}(jj,:,1,:))')',2);
%                         tmp.sc(jj,:,1,:)=(rotinv*squeeze(beam.s{2}(jj,:,1,:))')';
%                         tmp.n(jj,:,1,:)=(rotinv*squeeze(beam.s{3}(jj,:,1,:))')';
            end
        elseif size(beam.s{1},4)==2
            %                     case 1  % compute rho over Rzz interval -- original method
            t_int = 1:size(beam.s{1},2);  % use all available time points
            %    t_int=210:270
            rhoselectneg = find(mean(beam.s{1}(:,t_int,1,1),2)./mean(beam.s{1}(:,t_int,1,2),2)<0);
            % rhoselectneg = find(mean(beam.s_th(:,t_int)./beam.s_ph(:,t_int),2)<0);
            rho = atan(sqrt(mean(beam.s{1}(:,t_int,1,1).^2,2)./mean(beam.s{1}(:,t_int,1,2).^2,2)));
            rho(rhoselectneg) = -rho(rhoselectneg);
            rho(rivets.currMEGvoxelindex),

            cosrho = repmat(cos(rho),1,size(beam.s{1},2));
            sinrho = repmat(sin(rho),1,size(beam.s{1},2));
            tmp.sa(:,:,1,1) = beam.s{1}(:,:,1,2).*cosrho + beam.s{1}(:,:,1,1).*sinrho;
            tmp.sa(:,:,1,2) = beam.s{1}(:,:,1,1).*cosrho - beam.s{1}(:,:,1,2).*sinrho;
%                     tmp.sc(:,:,1,1) = beam.s{2}(:,:,1,2).*cosrho + beam.s{2}(:,:,1,1).*sinrho;
%                     tmp.sc(:,:,1,2) = beam.s{2}(:,:,1,1).*cosrho - beam.s{2}(:,:,1,2).*sinrho;
%                     tmp.n(:,:,1,1) = beam.s{3}(:,:,1,2).*cosrho + beam.s{3}(:,:,1,1).*sinrho;
%                     tmp.n(:,:,1,2) = beam.s{3}(:,:,1,1).*cosrho - beam.s{3}(:,:,1,2).*sinrho;
            clear cosrho sinrho s_ph s_th;  % done with these large matrices
        elseif size(beam.s{1},4)==1
            tmp.sa=beam.s{1};
        else
            disp('orientations on crack')
        end
    case 'Amplitude: Instant Orientation'
        % case 2  % compute rho for each time point -- this is crap: need to do something to make polarities consistent across time
        if size(beam.s{1},4)==3
            tmp.sa=zeros(size(beam.s{1}));
            disp('this is a very slow nested for-loop over voxels and time')
            disp('and crashes a bit too often, not sure if possible to do? fixme?')
            for jj=1:size(beam.s{1},1) % help! do you know how to make this not a for-loop?
                for kk=1:(size(beam.s{1},2)-2)
                S=cov(squeeze(beam.s{1}(jj,kk:(kk+2),1,:)));
                [vv,dd]=eig(inv(S));
                [jnk,ord]=sort(diag(dd),'descend');
                rotinv=vv(:,ord)';
%                         [uu,dd]=svd(inv(S));                        
%                         [jnk,ord]=sort(diag(dd),'descend');
%                         rotinv=vv(:,ord)'; %help! this was empirically determined to do this reorder (rather than just inv(v)), but is it true always??
                tmp.sa(jj,kk,1,:)=(rotinv*squeeze(beam.s{1}(jj,kk,1,:)))';
%                         tmp.sc(jj,:,1,:)=(rotinv*squeeze(beam.s{2}(jj,kk,1,:))')';
%                         tmp.n(jj,:,1,:)=(rotinv*squeeze(beam.s{3}(jj,kk,1,:))')';
                end
                tmp.sa(jj,(end-1):end,1,:)=(rotinv*squeeze(beam.s{1}(jj,(kk+1):(kk+2),1,:))')';
            end
        elseif size(beam.s{1},4)==2                
            tmp.sa = zeros(size(beam.s{1}));
            rivets.rho = (atan(beam.s{1}(:,:,1,1) ./ beam.s{1}(:,:,1,2)));
            tmp.sa(:,:,1,1) = beam.s{1}(:,:,1,2).*cos(rivets.rho) + beam.s{1}(:,:,1,1).*sin(rivets.rho);
            tmp.sa(:,:,1,2) = beam.s{1}(:,:,1,1).*cos(rivets.rho) - beam.s{1}(:,:,1,2).*sin(rivets.rho);
            clear s_ph s_th;  % done with these large matrices
        elseif size(beam.s{1},4)==1
            tmp.sa=beam.s{1};
        else
            disp('orientations on crack')
        end
    case 'Amplitude: Mob Orientation'
        %                     case 3  % sarang's fancy method (that doesn't make sense, but could be useful for other purposes)
        %%%%%%%% TODO: We can improve on this by computing orientation
        %%%%%%%% using moving average of rho (i.e., using 20ms sliding
        %%%%%%%% windows)
        tmp.sa=zeros(size(beam.s{1}));
        disp('this does 2 different things, depending if you have 2 or 3 orientations')
        disp('for 3 orient, not doing the mob polarity as it should for 2 orient..ya wanna code it?')
        if size(beam.s{1},4)==3
            disp('this is a very slow nested for-loop over voxels and time')
            disp('and crashes a bit too often, not sure if possible to do? fixme?')
            for jj=1:size(beam.s{1},1) % help! do you know how to make this not a for-loop?
                for kk=1:(size(beam.s{1},2)-20)
                S=cov(squeeze(beam.s{1}(jj,kk:(kk+20),1,:)));
                [vv,dd]=eig(inv(S));
                [jnk,ord]=sort(diag(dd),'descend');
                rotinv=vv(:,ord)';
%                         [uu,dd]=svd(inv(S));                        
%                         [jnk,ord]=sort(diag(dd),'descend');
%                         rotinv=vv(:,ord)'; %help! this was empirically determined to do this reorder (rather than just inv(v)), but is it true always??
                tmp.sa(jj,kk,1,:)=(rotinv*squeeze(beam.s{1}(jj,kk,1,:)))';
%                         tmp.sc(jj,:,1,:)=(rotinv*squeeze(beam.s{2}(jj,kk,1,:))')';
%                         tmp.n(jj,:,1,:)=(rotinv*squeeze(beam.s{3}(jj,kk,1,:))')';
                end
                tmp.sa(jj,(end-19):end,1,:)=(rotinv*squeeze(beam.s{1}(jj,(kk+1):(kk+20),1,:))')';
            end                    
        elseif size(beam.s{1},4)==2
            tanrho = beam.s{1}(:,:,1,1) ./ beam.s{1}(:,:,1,2);
            % polarities much noisier than magnitudes... so let's
            % decide polarity by mob justice -- 20 neighbors on either side...
            signrho = sign(filtfilt(ones(20,1),1,sign(tanrho')));
            rho2 = atan(signrho .* abs(tanrho'));
            % moving window average... more mob justice
            windowsize=20;
            rho = (filtfilt(ones(1,windowsize),1,rho2)/(windowsize*windowsize))';
            tmp.sa(:,:,1,1) = beam.s{1}(:,:,1,2).*cos(rho) + beam.s{1}(:,:,1,1).*sin(rho);
            tmp.sa(:,:,1,2) = beam.s{1}(:,:,1,1).*cos(rho) - beam.s{1}(:,:,1,2).*sin(rho);
            % clear s_ph s_th;  % done with these large matrices
        elseif size(beam.s{1},4)==1
            tmp.sa=beam.s{1};
        else
            disp('orientations on crack')
        end
    case 'Amplitude: Original Orientation'
        %                     case 4
        tmp.sa=beam.s{1};
        %                 end
        %             end
        %         end

        %%%% normalize s_II but keep activation_baseline
        %                 min_s = min(rivets.s_II(:));
        %                 max_s = max(rivets.s_II(:));
        %                 if(rivets.normflag)
        %                     rivets.maxblob = max(abs(max_s),abs(min_s));
        %                     rivets.scalefactor = 1000/rivets.maxblob;
        %                     rivets.minblob = -rivets.maxblob;
        %                 else
        %                     rivets.scalefactor = 1;
        %                 end

        % rivets.minblob=min_s;
        % rivets.maxblob=max_s;

        %%%%% commented out, probably
        %%%%% not needed anymore, but could be useful at some point???  %%%%%
        %
        %                 if(get(handles.nut_statsenabled_box,'Value'))  % EXPERIMENTAL: statistical thresholding
        %                     s_beam = rivets.s_II;
        %                     timepre=str2num(get(handles.Tstat_prestim_time,'String'));
        %                     prestim = find(beam.timewindow<=timepre);
        %                     if length(prestim)==0
        %                         disp('you need prestimulus data to compute a T value')
        %                     else
        %
        %                         sigma = sqrt(mean(s_beam(:,prestim).^2,2) - mean(abs(s_beam(:,prestim)),2).^2);
        %                         T = (abs(s_beam(:,prestim)) - repmat(mean(abs(s_beam(:,prestim)),2),[1 size(prestim,1)])) ./ repmat(sigma,[1 size(prestim,1)]);
        %                         Tmax = sort(max(T,[],2));
        %                         p = floor(0.99*size(Tmax,1));
        %                         Tthresh = Tmax(p);
        %                         rivets.threshold(:,1) = rivets.scalefactor*Tthresh*sigma + mean(abs(s_beam(:,prestim)),2);
        %                         rivets.threshold(:,2) = -rivets.threshold(:,1);
        %
        %                         %             threshold = round(threshold/max(abs(rivets.maxblob),abs(rivets.minblob))*1000);
        %                     end
        %                 end

    case 'Multi Time Series'
        % just hack it so it's like another orientation
        tmp.sa = beam.s{1};
        for ii=2:size(beam.s,1)
            tmp.sa(:,:,:,ii) = beam.s{ii,1};
        end
end
%     end
% else
%     tmp.sa=beam.s{1};
%     if size(beam.s,2)>1
%         tmp.sc=beam.s{2};
%     else
%         tmp.sc=zeros(size(tmp.sa),'single');
%     end
%end

% --- Executes on selection change in nut_contrastselect_menu.
function nut_contrastselect_menu_Callback(hObject, eventdata, handles)

global beam rivets

contrastcontents = get(handles.nut_contrastselect_menu,'String');
if ~iscell(contrastcontents), contrastcontents = {contrastcontents}; end    % This can happen if only 1 entry in the menu.
contrastselect = contrastcontents{get(handles.nut_contrastselect_menu,'Value')};

if(get(handles.nut_subtractnoise_box,'Value')==1)
    if size(beam.s,2)<3
        noise=0;
    else
        noise=beam.s{3};
    end
else
    noise=0;
end

% if (0 && isfield(beam,'W'))
%     switch get(handles.nut_noisetype,'Value')
%         case 1
%             [u,s]=svd(beam.params.TotC);
%             noisecov=s(end,end);
%         case 2
%             noisecov=beam.params.Czz;
%     end
%     if size(beam.W)==3
%     if size(beam.W{3},3)==1
%         %             pt = (diag(beam.W{3}'*beam.params.Rzz1*beam.W{3})-diag(beam.W{3}'*beam.params.Czz*beam.W{3}))./(2*diag(beam.W{3}'*beam.params.Czz*beam.W{3}));
%         %             pt = (diag(beam.W{3}'*beam.params.Rzz1*beam.W{3})-diag(beam.W{3}'*beam.params.Czz*beam.W{3}))./(2*diag(beam.W{3}'*sig*beam.W{3}));
%         pt = (diag(beam.W{3}'*beam.params.Rzz1*beam.W{3})-diag(beam.W{3}'*beam.params.Czz*beam.W{3}))./(2*diag(beam.W{3}'*noise*beam.W{3}));
%     else
%         for ii=1:size(beam.W{3},2)
%             apower(:,ii)=diag(squeeze(beam.W{3}(:,ii,:))'*beam.params.Rzz1*squeeze(beam.W{3}(:,ii,:)));
%             cpower(:,ii)=diag(squeeze(beam.W{3}(:,ii,:))'*beam.params.Czz*squeeze(beam.W{3}(:,ii,:)));
%             %                 anoise(:,ii)=diag(squeeze(beam.W{3}(:,ii,:))'*beam.params.Rzz1*squeeze(beam.W{3}(:,ii,:))); %wrong!
%             %                 cnoise(:,ii)=diag(squeeze(beam.W{3}(:,ii,:))'*beam.params.Czz*squeeze(beam.W{3}(:,ii,:)));
%             npower(:,ii)=diag(squeeze(beam.W{3}(:,ii,:))'*noise*squeeze(beam.W{3}(:,ii,:)));
%         end
%         %             pt = ((sum(apower')-sum(cpower'))./(sum(anoise')+sum(cnoise')))';
%         pt = ((sum(apower')-sum(cpower'))./(2*sum(npower')))';
%     end
%     end
% end

cond1=str2double(get(handles.nut_condcontrast1_text,'String'));
cond2=str2double(get(handles.nut_condcontrast2_text,'String'));
if cond1>size(beam.s,2)
    cond1=1;
    set(handles.nut_condcontrast1_text,'String','1');
end
if cond2>size(beam.s,2)
    cond2=2;
    set(handles.nut_condcontrast2_text,'String','2');
end

noorientations = ( rivets.timefmode || ndims(beam.s{1})<4 );
if ~noorientations
    tmp = nut_calc_orientations(handles,contrastselect);
end

switch contrastselect
    case 'F-ratio (dB)'
        if noorientations
            rivets.s = 10*log10(abs((beam.s{cond1}-noise)./(beam.s{cond2}-noise)));
        else
            rivets.s = 10*log10(abs((tmp.sa-noise)./(tmp.sc-noise)));
        end
    case 'F-ratio (raw)'
        if noorientations
            rivets.s = (beam.s{cond1}-noise)./(beam.s{cond2}-noise);
        else
            rivets.s = (tmp.sa-noise)./(tmp.sc-noise);
        end
    case '% Change'
        if noorientations
            rivets.s = 100*(beam.s{cond1} - beam.s{cond2})./(beam.s{cond2}-noise);
        else
            rivets.s = 100*(tmp.sa - tmp.sc)./(tmp.sc-noise);
        end
    case 'CTF Pseudo-F'
        if noorientations
            F = (beam.s{cond1}-noise)./(beam.s{cond2}-noise);
        else
            F = (tmp.sa-noise)./(tmp.sc-noise);
        end
        selectpos = find(F>=1);
        selectneg = find(F<1 & F>0);
        % selectzero = find(rivets.s==0);
        rivets.s(selectpos) = F(selectpos) - 1;
        rivets.s(selectneg) = 1-(1./F(selectneg));
    case 'Difference'
        if noorientations
            rivets.s = beam.s{cond1} - beam.s{cond2};
        else
            rivets.s = tmp.sa - tmp.sc;
        end
    case 'Difference (noise normalized)'
        if noorientations
            rivets.s = (beam.s{cond1} - beam.s{cond2})./noise;
        else
            rivets.s = (tmp.sa - tmp.sc)./noise;
        end
    case 'Greater of Active/Control'
        if noorientations
            rivets.s = beam.s{cond1};
            select = find(beam.s{cond2}>beam.s{cond1});
            rivets.s(select)=-beam.s{cond2}(select);
        else
            rivets.s = tmp.sa;
            select = find(tmp.sc>tmp.sa);
            rivets.s(select)=-tmp.sc(select);
        end
    case 'Active Power'
        if noorientations
            rivets.s = beam.s{cond1}-noise;
        else
            rivets.s = tmp.sa-noise;
        end
    case 'Control Power'
        if size(beam.s,2)<2, errordlg('This contrast is not available for your dataset.'), return, end
        rivets.s = beam.s{cond2}-noise;
        %         rivets.s = 10*log10(beam.s{2}./repmat(mean(beam.s{2},2),[1 41 1])) -noise;
    case 'Noise Power'
        rivets.s = noise;
    case 'Wilcoxon Z'
        if ~isfield(beam,'z'), errordlg('This contrast is not available for your dataset.'), return, end
        rivets.s = beam.z;
    case 'Active/Noise'
        if ~any(noise(:)), errordlg('This contrast is not available for your dataset.'), return, end
        if noorientations
            rivets.s = beam.s{1}./noise;
        else
            rivets.s = tmp.sa./noise;
        end
    case 'Imaginary Coherence'
        rivets.s = beam.s{1};
    case 'Magnitude Squared Coherence'
        rivets.s = beam.s{2};
    case 'Real Coherence'
        rivets.s = beam.s{3};
    case 'Abs Imaginary Coherence'
        rivets.s = abs(beam.s{1});
    case 't-test std-dev normalised'
        if isfield(beam,'sasd')
            rivets.s = (beam.s{cond1}-beam.s{cond2})./((beam.ssd{cond1}+beam.ssd{cond2})/2);
        else
            disp('you need to create beam.sasd and beam.scsd first')
        end
    case 'z of active'
        if isfield(beam,'sasd')
            rivets.s = (beam.s{cond1})./(beam.ssd{cond1});
        else
            disp('you need to create beam.sasd and beam.scsd first')
        end
    case 'z of control'
        if isfield(beam,'scsd')
            rivets.s = (beam.s{cond2})./(beam.ssd{cond2});
        else
            disp('you need to create beam.sasd and beam.scsd first')
        end
    otherwise
        error('Something funky this way comes...');
end
clear tmp

% rmfield rivets sa, sc, n
rivets.maxtf = max(abs(rivets.s(:)));

rivets.minblob=-rivets.maxtf;
rivets.maxblob=rivets.maxtf;

if(rivets.normflag)
    rivets.scalefactor = 1000/max(abs(rivets.maxblob),abs(rivets.minblob));
else
    rivets.scalefactor = 1;
end

if ~isempty(hObject)    % if user has clicked on menu we need to refresh, otherwise this function is called by another callback
    nut_tfselect(beam.timepts(rivets.timeselect),mean(beam.bands(rivets.freqselect,:)),handles);
end

%%------------------------------------------
function nut_lpf_text_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
% if ispc
%     set(hObject,'BackgroundColor','white');
% else
%     set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
% end


%%------------------------------------------
function nut_lpf_text_Callback(hObject, eventdata, handles)
plot_ts;
return;


%%------------------------------------------
function nut_hpf_text_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
% if ispc
%     set(hObject,'BackgroundColor','white');
% else
%     set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
% end


%%------------------------------------------
function nut_hpf_text_Callback(hObject, eventdata, handles)
plot_ts;

%%------------------------------------------
function nut_maxactivation_button_Callback(hObject, eventdata, handles)
% --- Executes on button press in nut_maxactivation_button.
global rivets st beam;

cursorpos = spm_orthviews('pos')';
%[MEGvoxelindex,distances] = dsearchn(rivets.voxelsMRI,nut_coordtfm(cursorpos,inv(st.vols{1}.premul)));
distances = nut_rownorm(nut_coord_diff(rivets.voxelsMRI,cursorpos));
exclude_radius = str2double(get(handles.nut_exclude_radius_text,'String'));
within_radius = str2double(get(handles.nut_within_radius_text,'String'));
searchindices = find(distances>=exclude_radius & distances<=within_radius);
if isempty(searchindices), return, end

s_beam_interval = rivets.s(searchindices,:,rivets.freqselect,1);

% only search voxels which are above threshold
s_beam_interval( ~any( rivets.threshold(searchindices,:,:),3 ) ) = NaN;
if all(isnan(s_beam_interval(:))), return, end
searchindices2 = find(any(isfinite(rivets.threshold(searchindices,:)),2));
s_beam_interval=s_beam_interval(searchindices2,:);
searchindices=searchindices(searchindices2);

actcontents = get(handles.nut_activation_menu,'String');
actselect = actcontents{get(handles.nut_activation_menu,'Value')};

switch(actselect)
    case 'abs max'
        [a,maxvoxels]=max(abs(s_beam_interval));
        [b,maxtime]=max(a);
    case 'max'
        [a,maxvoxels]=max(s_beam_interval);
        [b,maxtime]=max(a);
    case 'min'
        [a,maxvoxels]=min(s_beam_interval);
        [b,maxtime]=min(a);
    otherwise
        error('well, this is confusing.')
end

MEGvoxelindex=searchindices(maxvoxels(maxtime));
%     spm_orthviews('Reposition',nut_coordtfm(rivets.voxelsMRI(MEGvoxelindex,:),st.vols{1}.premul))
nut_reposition(nut_coordtfm(rivets.voxelsMRI(MEGvoxelindex,:),st.vols{1}.premul));

timept = maxtime;
timept2 = maxtime;
time = beam.timepts(timept);
time2 = beam.timepts(timept2);
set(handles.nut_time_text,'String',time);
set(handles.nut_time2_text,'String',time2);

rivets.timeselect = timept;
% nut_do_threshold([],handles);
paint_activation(hObject,eventdata,handles);

%%------------------------------------------
function nut_maxvoxel_button_Callback(hObject, eventdata, handles)
% --- Executes on button press in nut_maxvoxel_button.
global rivets st beam;
time = str2double(get(handles.nut_time_text,'String'));
time2 = str2double(get(handles.nut_time2_text,'String'));
timept = dsearchn(beam.timepts,time);
timept2 = dsearchn(beam.timepts,time2);

cursorpos = spm_orthviews('pos')';
%[MEGvoxelindex,distances] = dsearchn(rivets.voxelsMRI,nut_coordtfm(cursorpos,inv(st.vols{1}.premul)));
distances = nut_rownorm(nut_coord_diff(rivets.voxelsMRI,cursorpos));
exclude_radius = str2double(get(handles.nut_exclude_radius_text,'String'));
within_radius = str2double(get(handles.nut_within_radius_text,'String'));
searchindices = find(distances>=exclude_radius & distances<=within_radius);
if isempty(searchindices), return, end

s_beam_interval = mean(rivets.s(searchindices,timept:timept2,rivets.freqselect,1),2);

% only search voxels which are above threshold
s_beam_interval( ~any( rivets.threshold(searchindices,rivets.timeselect,:),3 ) ) = NaN;
if all(isnan(s_beam_interval)), return, end
searchindices2 = find( ~isinf(rivets.threshold(searchindices,1)) );
s_beam_interval=s_beam_interval(searchindices2);
searchindices=searchindices(searchindices2);


actcontents = get(handles.nut_activation_menu,'String');
actselect = actcontents{get(handles.nut_activation_menu,'Value')};

switch(actselect)
    case 'abs max'
        [a,maxvoxel]=max(abs(s_beam_interval));
    case 'max'
        [a,maxvoxel]=max(s_beam_interval);
    case 'min'
        [a,maxvoxel]=min(s_beam_interval);
    otherwise
        error('well, this is confusing.')
end

MEGvoxelindex=searchindices(maxvoxel);
nut_reposition(nut_coordtfm(rivets.voxelsMRI(MEGvoxelindex,:),st.vols{1}.premul))

%%------------------------------------------
function nut_maxtime_button_Callback(hObject, eventdata, handles)
% --- Executes on button press in nut_maxtime_button.
global rivets beam;
cursorpos = spm_orthviews('pos')';
MEGvoxelindex = dsearchn(rivets.voxelsMRI,cursorpos);

s_beam_interval = rivets.s(MEGvoxelindex,:,rivets.freqselect,1);
if get(handles.nut_threshold_menu,'Value')>1        % only search time points which are statistically significant
    s_beam_interval( find( ~any( rivets.threshold(MEGvoxelindex,:,:),3 ) ) ) = NaN;
    if all(isnan(s_beam_interval)), return, end
end

actcontents = get(handles.nut_activation_menu,'String');
actselect = actcontents{get(handles.nut_activation_menu,'Value')};

switch(actselect)
    case 'abs max'
        [a,maxtime]=max(abs(s_beam_interval),[],2);
    case 'max'
        [a,maxtime]=max(s_beam_interval,[],2);
    case 'min'
        [a,maxtime]=min(s_beam_interval,[],2);
    otherwise
        error('well, this is confusing.')
end

timept = maxtime;
timept2 = maxtime;
time = beam.timepts(timept);
time2 = beam.timepts(timept2);
set(handles.nut_time_text,'String',time);
set(handles.nut_time2_text,'String',time2);

rivets.timeselect = timept;
nut_do_threshold([],handles);
paint_activation(hObject,eventdata,handles);


%%------------------------------------------
function nut_applyfilter_button_Callback(hObject, eventdata, handles)
% --- Executes on button press in nut_applyfilter_button.
nut_refreshdisplay(handles);

%%-------------------------------------------
function plot_ts(hObject)
% refresh time series plot (input argument used only if "apply to volume" button is pressed)

global rivets st beam ndefaults;

handles = guidata(rivets.fig);

apply_all_flg = 0;

%%%% set properties for selected display style
stylecontents = get(handles.nut_style_menu,'String');
styleselect = stylecontents{get(handles.nut_style_menu,'Value')};

switch styleselect
    case 'Normal Style'
%         bf_ts_color = [0 0 1];
        bf_ts_color = [0 0 1; 0 1 0; 1 0 0];
        bf_ts2_color = [0 1 0];
        bf_ts_width = 0.5;
        bf_ts2_width = bf_ts_width;
        bf_ts_xhair_color = [1 0 0];
        bf_ts_xhair_width = 0.5;
        % meg_ts_color = [0 0 1];
        meg_ts_width = 0.5;
        meg_ts_xhair_color = bf_ts_xhair_color;
        meg_ts_xhair_width = bf_ts_xhair_width;
        spm_xhair_color = [0 1 0];
        spm_xhair_width = 0.5;
        if(rivets.minblob < 0)
            set(st.fig,'ColorMap',rivets.grayjet);
            if rivets.timefmode, set(rivets.fig,'ColorMap',jet(64)); end
        else
            set(st.fig,'ColorMap',rivets.grayhot);
            if rivets.timefmode, set(rivets.fig,'ColorMap',hot(64)); end
        end
    case 'Poster Style (Color)'
%         bf_ts_color = [1 0 0];
        bf_ts_color = [1 0 0; 0 .8 0; 0 0 .8];
        bf_ts2_color = [0 .8 0];
        bf_ts_width = 2;
        bf_ts2_width = bf_ts_width;
        bf_ts_xhair_color = [0 .7 1];
        bf_ts_xhair_width = 2;
        % meg_ts_color = [0 0 1];
        meg_ts_width = 2;
        meg_ts_xhair_color = bf_ts_xhair_color;
        meg_ts_xhair_width = bf_ts_xhair_width;
        spm_xhair_color = [0 1 0];
        spm_xhair_width = 1;
        if(rivets.minblob < 0)
            set(st.fig,'ColorMap',rivets.grayjet);
            if rivets.timefmode, set(rivets.fig,'ColorMap',jet(64)); end
        else
            set(st.fig,'ColorMap',rivets.grayhot);
            if rivets.timefmode, set(rivets.fig,'ColorMap',hot(64)); end
        end
    case 'Poster Style (Color2)'
        set([st.slider{:}],'Visible','off');
%         bf_ts_color = [1 0 0];
        bf_ts_color = [1 0 0; 0 .8 0; 0 0 .8];
        bf_ts2_color = [0 .8 0];
        bf_ts_width = 2;
        bf_ts2_width = bf_ts_width;
        bf_ts_xhair_color = [0 .7 1];
        bf_ts_xhair_width = 2;
        % meg_ts_color = [0 0 1];
        meg_ts_width = 2;
        meg_ts_xhair_color = bf_ts_xhair_color;
        meg_ts_xhair_width = bf_ts_xhair_width;
        spm_xhair_color = [0 1 0];
        spm_xhair_width = 0.5;
        if(rivets.minblob < 0)
            set(st.fig,'ColorMap',rivets.grayhotsymmetric);

        else
            set(st.fig,'ColorMap',rivets.grayhot);
        end
    case 'Poster Style (Monochrome)'
%         bf_ts_color = [0 0 0];
        bf_ts_color = [0 0 0; .3 .3 .3; .6 .6 .6];
        bf_ts2_color = [0 0 0];
        bf_ts_width = 2;
        bf_ts2_width = bf_ts_width;
        bf_ts_xhair_color = [0 0 0];
        bf_ts_xhair_width = 2;
        meg_ts_color = [0 0 0];
        meg_ts_width = 2;
        meg_ts_xhair_color = bf_ts_xhair_color;
        meg_ts_xhair_width = bf_ts_xhair_width;
        spm_xhair_color = [.9 .9 .9];
        spm_xhair_width = 0.5;
        if(rivets.minblob < 0)
            set(st.fig,'ColorMap',rivets.graygraybipolar);
        else
            set(st.fig,'ColorMap',rivets.graygray);
        end
    case 'Presentation Style'
%         bf_ts_color = [1 .5 0];
        bf_ts_color = [1 .5 0; .5 1 0; 0 .5 1];
        bf_ts2_color = [1 .5 0];
        bf_ts_width = 2;
        bf_ts2_width = bf_ts_width;
        bf_ts_xhair_color = [1 1 0];
        bf_ts_xhair_width = 2;
        % meg_ts_color = [0 0 1];
        meg_ts_width = 2;
        meg_ts_xhair_color = bf_ts_xhair_color;
        meg_ts_xhair_width = bf_ts_xhair_width;
        spm_xhair_color = [1 1 0];
        spm_xhair_width = 2;
        if(rivets.minblob < 0)
            set(st.fig,'ColorMap',rivets.grayjet);
            if rivets.timefmode, set(rivets.fig,'ColorMap',jet(64)); end
        else
            set(st.fig,'ColorMap',rivets.grayhot);
            if rivets.timefmode, set(rivets.fig,'ColorMap',hot(64)); end
        end
    case {'Paper Style (Color)'}
%         bf_ts_color = [1 0 0];
        bf_ts_color = [1 0 0; 0 1 0; 0 0 1];
        bf_ts2_color = [1 0 0];
        bf_ts_width = 1;
        bf_ts2_width = bf_ts_width;
        bf_ts_xhair_color = [0 0 1];
        bf_ts_xhair_width = 1;
        % meg_ts_color = [0 0 1];
        meg_ts_width = 1;
        meg_ts_xhair_color = bf_ts_xhair_color;
        meg_ts_xhair_width = bf_ts_xhair_width;
        spm_xhair_color = [0 1 0];
        spm_xhair_width = 0.5;
        if(rivets.minblob < 0)
            set(st.fig,'ColorMap',rivets.grayjet);
            if rivets.timefmode, set(rivets.fig,'ColorMap',jet(64)); end
        else
            set(st.fig,'ColorMap',rivets.grayhot);
            if rivets.timefmode, set(rivets.fig,'ColorMap',hot(64)); end
        end
    case 'Paper Style (Monochrome)'
%         bf_ts_color = [0 0 0];
        bf_ts_color = [0 0 0; .3 .3 .3; .6 .6 .6];
        bf_ts2_color = [0 0 0];
        bf_ts_width = 1;
        bf_ts2_width = bf_ts_width;
        bf_ts_xhair_color = [0 0 0];
        bf_ts_xhair_width = 1;
        meg_ts_color = [0 0 0];
        meg_ts_width = 1;
        meg_ts_xhair_color = bf_ts_xhair_color;
        meg_ts_xhair_width = bf_ts_xhair_width;
        spm_xhair_color = [0.9 0.9 0.9];
        spm_xhair_width = 0.5;
        if(rivets.minblob < 0)
            set(st.fig,'ColorMap',rivets.graygraybipolar);
            if rivets.timefmode, set(rivets.fig,'ColorMap',rivets.graygraybipolar(65:end,:)); end
        else
            set(st.fig,'ColorMap',rivets.graygray);
            if rivets.timefmode, set(rivets.fig,'ColorMap',0.8*gray(64)); end
        end
    case 'SAM Style'
%         bf_ts_color = [0 0 1];
        bf_ts_color = [0 0 1; 0 1 0; 1 0 0];
        bf_ts2_color = [0 1 0];
        bf_ts_width = 0.5;
        bf_ts2_width = bf_ts_width;
        bf_ts_xhair_color = [1 0 0];
        bf_ts_xhair_width = 0.5;
        % meg_ts_color = [0 0 1];
        meg_ts_width = 0.5;
        meg_ts_xhair_color = bf_ts_xhair_color;
        meg_ts_xhair_width = bf_ts_xhair_width;
        spm_xhair_color = [0 1 0];
        spm_xhair_width = 0.5;
end

%nut_set_sliders;  % jfh
cursorpos = spm_orthviews('pos')';

% st.vols{1} is empty for brain rendering mode
if(isempty(st.vols{1}))
    MEGvoxelindex = dsearchn(rivets.voxelsMRI,cursorpos);
else
    for i=1:3
        set([st.vols{1}.ax{i}.lx st.vols{1}.ax{i}.ly],'Color',spm_xhair_color,'LineWidth',spm_xhair_width);
    end
    MEGvoxelindex = dsearchn(rivets.voxelsMRI,nut_coordtfm(cursorpos,inv(st.vols{1}.premul)));
end
rivets.currMEGvoxelindex=MEGvoxelindex;

if(ndefaults.vwr.plotweights)
    nut_plotweights(squeeze(beam.W(:,2,MEGvoxelindex)));
end

% Scaling
%--------
view_scale = str2num(get(handles.nut_view_scale_text,'String'));    % str2double does not work here, because values may be empty
view_scale2 = str2num(get(handles.nut_view_scale2_text,'String'));

if isempty(view_scale), view_scale = rivets.minblob; end
if isempty(view_scale2), view_scale2 = rivets.maxblob; end
%if view_scale > beam.timepts(end), view_scale = beam.timepts(end); end
%if view_scale2 < beam.timepts(1), view_scale2 = beam.timepts(1); end
%if view_scale2 > beam.timepts(end), view_scale2 = beam.timepts(end); end

if(view_scale2 < view_scale)
    view_scaletmp = view_scale2;
    view_scale2 = view_scale;
    view_scale = view_scaletmp;
end
    %set(handles.nut_view_scale_text,'String',view_scale);
    %set(handles.nut_view_scale2_text,'String',view_scale2);
    %     view_scalept = dsearchn(beam.timepts,view_scale);
    %     view_scalept2 = dsearchn(beam.timepts,view_scale2);

%     if rivets.setscale
rivets.scalemin=view_scale;
rivets.scalemax=view_scale2;
%     else
        % rivets.scalemin=rivets.minblob*rivets.scalefactor;
        % rivets.scalemax=rivets.maxblob*rivets.scalefactor;
%     end

% Timef-plot
%-----------
if rivets.timefmode
    nut_plotbeamtf(MEGvoxelindex);
end

% Time Axis
%----------
time = str2double(get(handles.nut_time_text,'String'));
time2 = str2double(get(handles.nut_time2_text,'String'));

if time < beam.timepts(1), time = beam.timepts(1); end
if time > beam.timepts(end), time = beam.timepts(end); end
if time2 < beam.timepts(1), time2 = beam.timepts(1); end
if time2 > beam.timepts(end), time2 = beam.timepts(end); end

if(time2 < time)
    timetmp = time2;
    time2 = time;
    time = timetmp;
end

set(handles.nut_time_text,'String',time);
set(handles.nut_time2_text,'String',time2);

timept = dsearchn(beam.timepts,time);
timept2 = dsearchn(beam.timepts,time2);
if size(rivets.s,2)==1 % occurs when timefmode=0, and contrastselect is not Active Power
    timept=1;timept2=1;
end

% Some other stuff
% ----------------
if(isfield(st,'beamin')) % may not exist for brain renderings...
    beam_val = sum(squeeze(mean(rivets.s(MEGvoxelindex,timept:timept2,rivets.freqselect,:),2)));
    set(st.beamin,'String',sprintf('%g',rivets.scalefactor*beam_val)); clear beam_val
end

if(0 & isfield(beam,'snpm'))  % print p-values on command line...
    if(isfield(beam.snpm,'p_corr_pos'))
        beam.snpm.p_corr_pos(MEGvoxelindex,rivets.timeselect,rivets.freqselect)
    end
    if(isfield(beam.snpm,'p_corr_neg'))
        beam.snpm.p_corr_neg(MEGvoxelindex,rivets.timeselect,rivets.freqselect)
    end
end

if(0)
    % plot on MEG time series
    if(strcmp(get(handles.nut_megts_axes,'Visible'),'on'))
        axisvalues = axis(handles.nut_megts_axes);
        line_max=axisvalues(4);
        line_min=axisvalues(3);
        axes(handles.nut_megts_axes); hold on
        plot([time time],[line_min line_max],'Color',meg_ts_xhair_color,'LineWidth',meg_ts_xhair_width);
        plot([time2 time2],[line_min line_max],'Color',meg_ts_xhair_color,'LineWidth',meg_ts_xhair_width);
        hold off
        axes(handles.nut_ts_axes);
    end

    % plot on beamformer time series
    axisvalues = axis;
    line_max=axisvalues(4);
    line_min=axisvalues(3);

    threshold = rivets.threshold(MEGvoxelindex,rivets.timeselect,1);
    thresholdneg = rivets.threshold(MEGvoxelindex,rivets.timeselect,2);
    hold on
    time1_handle = plot([time time],[line_min line_max],'Color',bf_ts_xhair_color,'LineWidth',bf_ts_xhair_width);
    time2_handle = plot([time2 time2],[line_min line_max],'Color',bf_ts_xhair_color,'LineWidth',bf_ts_xhair_width);
    thresh_handle(1) = plot([beam.timepts(1) beam.timepts(end)],[threshold threshold],'Color',bf_ts_xhair_color,'LineWidth',bf_ts_xhair_width);
    thresh_handle(2) = plot([beam.timepts(1) beam.timepts(end)],[thresholdneg thresholdneg],'Color',bf_ts_xhair_color,'LineWidth',bf_ts_xhair_width);
    hold off
end

% Line-plot
% ---------
if ~rivets.timefmode
    view_time = str2double(get(handles.nut_view_time_text,'String'));
    view_time2 = str2double(get(handles.nut_view_time2_text,'String'));

    if view_time < beam.timepts(1), view_time = beam.timepts(1); end
    if view_time > beam.timepts(end), view_time = beam.timepts(end); end
    if view_time2 < beam.timepts(1), view_time2 = beam.timepts(1); end
    if view_time2 > beam.timepts(end), view_time2 = beam.timepts(end); end

    if(view_time2 < view_time)
        view_timetmp = view_time2;
        view_time2 = view_time;
        view_time = view_timetmp;
    end

    set(handles.nut_view_time_text,'String',view_time);
    set(handles.nut_view_time2_text,'String',view_time2);

    view_timept = dsearchn(beam.timepts,view_time);
    view_timept2 = dsearchn(beam.timepts,view_time2);
    if size(rivets.s,2)==1 % occurs when timefmode=0, and contrastselect is not Active Power
        view_timept=1;view_timept2=1;
    end

    fs = beam.srate;
    % these values pertain to the sbeam data, not the meg data (which uses beam.params.preprocessing.*)
    %notch = get(handles.nut_contrastselect_menu,'Value');
    %bpf = get(handles.activation_bpf,'Value');
    %tmp = round(str2num(get(handles.nut_lpf_text,'String'))*10)/10; if tmp < 0, tmp = 0; end, if tmp >= fs/2, tmp = fix(fs/2)-0.1; end, set(handles.nut_lpf_text,'String',num2str(tmp));
    %bpf_low_cutoff = tmp;
    %tmp = round(str2num(get(handles.nut_hpf_text,'String'))*10)/10; if tmp < 0, tmp = 0; end, if tmp >= fs/2, tmp = fix(fs/2)-0.1; end, set(handles.nut_hpf_text,'String',num2str(tmp));
    %bpf_high_cutoff = tmp;

    if apply_all_flg
        ndx = 1:size(rivets.s,1);
    else
        ndx = MEGvoxelindex;
    end

    data = squeeze(rivets.s(ndx,:,rivets.freqselect,:));
    if size(rivets.s,4)==1
        data=data';
    end
    
    if isfield(rivets,'meg')
        if isfield(rivets.meg,'yc_full')
            megdata=rivets.meg.yc_full';
        elseif isfield(rivets.meg,'data')
            megdata=nut_filter(rivets.meg.data,beam.timewindow,fs,beam.params.preprocessing.baseline,beam.params.preprocessing.notch,beam.params.preprocessing.bpf,beam.params.preprocessing.bpf_low_cutoff,beam.params.preprocessing.bpf_high_cutoff);
        end
        axes(handles.nut_megts_axes);
        meg_ts_handle = plot(beam.timewindow(view_timept:view_timept2),1e15*megdata(view_timept:view_timept2,:),'LineWidth',meg_ts_width);
        if(exist('meg_ts_color','var'))
            set(meg_ts_handle,'Color',meg_ts_color);
        end
        axis tight; xlabel('Time (ms)'); ylabel('Magnetic Field (fT)');
    end

    %data = nut_filter(data,beam.timewindow,fs,1,notch,bpf,bpf_low_cutoff,bpf_high_cutoff);
       
    axes(handles.nut_ts_axes);
    %plot(beam.timewindow,data(:,ndx));
    % jz wonders: is this data2 idea defunct?
    if(exist('data2','var'))
        %         plot(beam.timepts(view_timept:view_timept2),rivets.scalefactor*data2(view_timept:view_timept2,ndx,:),'Color',bf_ts2_color,'LineWidth',bf_ts_width);
        %         hold on;
        %         plot(beam.timepts(view_timept:view_timept2),rivets.scalefactor*data(view_timept:view_timept2,ndx,:),'Color',bf_ts_color,'LineWidth',bf_ts_width);
        plot(beam.timepts(view_timept:view_timept2),rivets.scalefactor*data2(view_timept:view_timept2,:),'Color',bf_ts2_color,'LineWidth',bf_ts_width);
        hold on;
        plot(beam.timepts(view_timept:view_timept2),rivets.scalefactor*data(view_timept:view_timept2,:),'Color',bf_ts_color,'LineWidth',bf_ts_width);
    else
        %         plot(beam.timepts(view_timept:view_timept2),rivets.scalefactor*data(view_timept:view_timept2,ndx,:),'Color',bf_ts_color,'LineWidth',bf_ts_width);
        for ii=1:size(data,2)
            plot(beam.timepts(view_timept:view_timept2),rivets.scalefactor*data(view_timept:view_timept2,ii),'Color',bf_ts_color(ii,:),'LineWidth',bf_ts_width);
            hold on;
        end
        hold off;
    end
    % axis([beam.timewindow(view_timept) beam.timewindow(view_timept2) 0 1000]);
    axis tight;
    if isfield(beam,'labels')
        xlabel(beam.labels.xaxis)
        ylabel(beam.labels.yaxis)
        temp = axis(handles.nut_ts_axes);
        axis(handles.nut_ts_axes,[temp(1) temp(2) rivets.scalemin-1e-3 rivets.scalemax+1e-3]); % rescale axes to min:max        
    else
        xlabel('Time (ms)');
        if rivets.normflag
            ylabel('Normalized Intensity (0-1000)');
            temp = axis(handles.nut_ts_axes);
            axis(handles.nut_ts_axes,[temp(1) temp(2) rivets.scalemin rivets.scalemax]); % rescale axes to min:max
        else
            ylabel('Intensity');
            temp = axis(handles.nut_ts_axes);
            %         axis(handles.nut_ts_axes,[temp(1) temp(2) rivets.minblob rivets.maxblob]); % rescale axes to min:max
            axis(handles.nut_ts_axes,[temp(1) temp(2) rivets.scalemin rivets.scalemax]); % rescale axes to min:max
        end
    end

    % plot on MEG time series
    if(strcmp(get(handles.nut_megts_axes,'Visible'),'on'))
        axisvalues = axis(handles.nut_megts_axes);
        line_max=axisvalues(4);
        line_min=axisvalues(3);
        axes(handles.nut_megts_axes); hold on
        plot([time time],[line_min line_max],'Color',meg_ts_xhair_color,'LineWidth',meg_ts_xhair_width);
        plot([time2 time2],[line_min line_max],'Color',meg_ts_xhair_color,'LineWidth',meg_ts_xhair_width);
        hold off
        axes(handles.nut_ts_axes);
    end

    % plot on beamformer time series
    axisvalues = axis;
    line_max=axisvalues(4);
    line_min=axisvalues(3);

    hold on
    time1_handle = plot([time time],[line_min line_max],'Color',bf_ts_xhair_color,'LineWidth',bf_ts_xhair_width);
    time2_handle = plot([time2 time2],[line_min line_max],'Color',bf_ts_xhair_color,'LineWidth',bf_ts_xhair_width);
    if get(handles.nut_threshold_menu,'Value')==1
        threshold = rivets.threshold(MEGvoxelindex,rivets.timeselect,1);
        thresholdneg = rivets.threshold(MEGvoxelindex,rivets.timeselect,2);
        thresh_handle(1) = plot([beam.timepts(1) beam.timepts(end)],[threshold threshold],'Color',bf_ts_xhair_color,'LineWidth',bf_ts_xhair_width);
        thresh_handle(2) = plot([beam.timepts(1) beam.timepts(end)],[thresholdneg thresholdneg],'Color',bf_ts_xhair_color,'LineWidth',bf_ts_xhair_width);
    end
    hold off
end

switch styleselect
    case 'Normal Style'
        set(handles.nut_ts_axes,'FontUnits','points','LineWidth',2,'FontName','Helvetica','FontWeight','bold','FontSize',12,'Color',[1 1 1],'XColor',[0 0 0],'YColor',[0 0 0]);  %changes time series axes numbers font
        set(handles.nut_megts_axes,'FontUnits','points','LineWidth',2,'FontName','Helvetica','FontWeight','bold','FontSize',12,'Color',[1 1 1],'XColor',[0 0 0],'YColor',[0 0 0]);  %changes time series axes numbers font
        set(handles.nut_results_viewer,'Color',[1 1 1]);
        %         set(handles.nut_results_viewer,'Color',[.93 .93 .92]);  %color of background of whole figure
        set(cell2mat(get(handles.nut_ts_axes,{'Xlabel','Ylabel'})),'FontUnits','points','FontName','Helvetica','FontWeight','bold','FontSize',12,'Color',[0 0 0]);  %changes time series axes label
        set(cell2mat(get(handles.nut_megts_axes,{'Xlabel','Ylabel'})),'FontUnits','points','FontName','Helvetica','FontWeight','bold','FontSize',12,'Color',[0 0 0]);  %changes time series axes label
        % set(handles.nut_blobstyle,'Value',1);
    case 'Poster Style (Color)'
        set(handles.nut_ts_axes,'FontUnits','points','LineWidth',2,'FontName','Helvetica','FontWeight','bold','FontSize',18,'Color',[1 1 1],'XColor',[0 0 0],'YColor',[0 0 0]);%changes time series axes numbers font
        set(handles.nut_megts_axes,'FontUnits','points','LineWidth',2,'FontName','Helvetica','FontWeight','bold','FontSize',18,'Color',[1 1 1],'XColor',[0 0 0],'YColor',[0 0 0]);%changes time series axes numbers font
        set(handles.nut_results_viewer,'Color',[1 1 1]);
        set(cell2mat(get(handles.nut_ts_axes,{'Xlabel','Ylabel'})),'FontUnits','points','FontName','Helvetica','FontWeight','bold','FontSize',18);%changes time series axes label
        set(cell2mat(get(handles.nut_megts_axes,{'Xlabel','Ylabel'})),'FontUnits','points','FontName','Helvetica','FontWeight','bold','FontSize',18);%changes time series axes label
        % set(handles.nut_blobstyle,'Value',2);
    case 'Poster Style (Color2)'
        set(handles.nut_ts_axes,'FontUnits','points','LineWidth',2,'FontName','Helvetica','FontWeight','bold','FontSize',18,'Color',[1 1 1],'XColor',[0 0 0],'YColor',[0 0 0]);%changes time series axes numbers font
        set(handles.nut_megts_axes,'FontUnits','points','LineWidth',2,'FontName','Helvetica','FontWeight','bold','FontSize',18,'Color',[1 1 1],'XColor',[0 0 0],'YColor',[0 0 0]);%changes time series axes numbers font
        set(handles.nut_results_viewer,'Color',[1 1 1]);
        set(cell2mat(get(handles.nut_ts_axes,{'Xlabel','Ylabel'})),'FontUnits','points','FontName','Helvetica','FontWeight','bold','FontSize',18);%changes time series axes label
        set(cell2mat(get(handles.nut_megts_axes,{'Xlabel','Ylabel'})),'FontUnits','points','FontName','Helvetica','FontWeight','bold','FontSize',18);%changes time series axes label
        set(handles.nut_blobstyle,'Value',1);
    case 'Poster Style (Monochrome)'
        set(handles.nut_ts_axes,'FontUnits','points','LineWidth',2,'FontName','Helvetica','FontWeight','bold','FontSize',18,'Color',[1 1 1],'XColor',[0 0 0],'YColor',[0 0 0]);%changes time series axes numbers font
        set(handles.nut_megts_axes,'FontUnits','points','LineWidth',2,'FontName','Helvetica','FontWeight','bold','FontSize',18,'Color',[1 1 1],'XColor',[0 0 0],'YColor',[0 0 0]);%changes time series axes numbers font
        set(handles.nut_results_viewer,'Color',[1 1 1]);
        set(cell2mat(get(handles.nut_ts_axes,{'Xlabel','Ylabel'})),'FontUnits','points','FontName','Helvetica','FontWeight','bold','FontSize',18);%changes time series axes label
        set(cell2mat(get(handles.nut_megts_axes,{'Xlabel','Ylabel'})),'FontUnits','points','FontName','Helvetica','FontWeight','bold','FontSize',18);%changes time series axes label
        set(handles.nut_blobstyle,'Value',1);
    case 'Presentation Style'
        set(handles.nut_ts_axes,'FontUnits','points','LineWidth',2,'FontName','Helvetica','FontWeight','bold','FontSize',15,'Color',[0 0 0],'XColor',[1 1 0],'YColor',[1 1 0]);  %changes time series axes numbers font
        set(handles.nut_megts_axes,'FontUnits','points','LineWidth',2,'FontName','Helvetica','FontWeight','bold','FontSize',15,'Color',[0 0 0],'XColor',[1 1 0],'YColor',[1 1 0]);  %changes time series axes numbers font
        set(handles.nut_results_viewer,'Color',[0 0 0]);
        set(cell2mat(get(handles.nut_ts_axes,{'Xlabel','Ylabel'})),'FontUnits','points','FontName','Helvetica','FontWeight','bold','FontSize',15,'Color',[1 1 0]);  %changes time series axes label
        set(cell2mat(get(handles.nut_megts_axes,{'Xlabel','Ylabel'})),'FontUnits','points','FontName','Helvetica','FontWeight','bold','FontSize',15,'Color',[1 1 0]);  %changes time series axes label
        %         set(handles.nut_blobstyle,'Value',2);
    case {'Paper Style (Color)'}
        set(handles.nut_ts_axes,'FontUnits','points','LineWidth',1,'FontName','Times','FontWeight','bold','FontSize',15,'Color',[1 1 1],'XColor',[0 0 0],'YColor',[0 0 0]);  %changes time series axes numbers font
        set(handles.nut_megts_axes,'FontUnits','points','LineWidth',1,'FontName','Times','FontWeight','bold','FontSize',15,'Color',[1 1 1],'XColor',[0 0 0],'YColor',[0 0 0]);  %changes time series axes numbers font
        set(handles.nut_results_viewer,'Color',[1 1 1]);
        set(cell2mat(get(handles.nut_ts_axes,{'Xlabel','Ylabel'})),'FontUnits','points','FontName','Times','FontWeight','bold','FontSize',15);  %changes time series axes label
        set(cell2mat(get(handles.nut_megts_axes,{'Xlabel','Ylabel'})),'FontUnits','points','FontName','Times','FontWeight','bold','FontSize',15);  %changes time series axes label
        set(handles.nut_blobstyle,'Value',2);
    case 'Paper Style (Monochrome)'
        set(handles.nut_ts_axes,'FontUnits','points','LineWidth',1,'FontName','Times','FontWeight','bold','FontSize',15,'Color',[1 1 1],'XColor',[0 0 0],'YColor',[0 0 0]);  %changes time series axes numbers font
        set(handles.nut_megts_axes,'FontUnits','points','LineWidth',1,'FontName','Times','FontWeight','bold','FontSize',15,'Color',[1 1 1],'XColor',[0 0 0],'YColor',[0 0 0]);  %changes time series axes numbers font
        set(handles.nut_results_viewer,'Color',[1 1 1]);
        set(cell2mat(get(handles.nut_ts_axes,{'Xlabel','Ylabel'})),'FontUnits','points','FontName','Times','FontWeight','bold','FontSize',15);  %changes time series axes label
        set(cell2mat(get(handles.nut_megts_axes,{'Xlabel','Ylabel'})),'FontUnits','points','FontName','Times','FontWeight','bold','FontSize',15);  %changes time series axes label
        set(handles.nut_blobstyle,'Value',1);
    case 'SAM Style'   % jfh
        set(handles.nut_ts_axes,'FontUnits','points','LineWidth',2,'FontName','Helvetica','FontWeight','bold','FontSize',12,'Color',[1 1 1],'XColor',[0 0 0],'YColor',[0 0 0]);  %changes time series axes numbers font
        set(handles.nut_megts_axes,'FontUnits','points','LineWidth',2,'FontName','Helvetica','FontWeight','bold','FontSize',12,'Color',[1 1 1],'XColor',[0 0 0],'YColor',[0 0 0]);  %changes time series axes numbers font
        set(handles.nut_results_viewer,'Color',[.93 .93 .92]);  %color of background of whole figure
        set(cell2mat(get(handles.nut_ts_axes,{'Xlabel','Ylabel'})),'FontUnits','points','FontName','Helvetica','FontWeight','bold','FontSize',12,'Color',[0 0 0]);  %changes time series axes label
        set(cell2mat(get(handles.nut_megts_axes,{'Xlabel','Ylabel'})),'FontUnits','points','FontName','Helvetica','FontWeight','bold','FontSize',12,'Color',[0 0 0]);  %changes time series axes label

        mytargintens = 0.7;
        myncolsteps = 64;

        mycolmax = rivets.maxblob;
        mycolmin = rivets.minblob;
        mycolvals = mycolmin:((mycolmax-mycolmin)/(myncolsteps-1)):mycolmax;
        [duh,zvalidx] = min(abs(mycolvals));
        mynnegcolsteps = zvalidx;
        mynposcolsteps = myncolsteps - mynnegcolsteps;

        one2targ = 1:(-(1-mytargintens)/(mynnegcolsteps-1)):mytargintens;
        one2zero = 1:(-1/(mynnegcolsteps-1)):0;
        mycolmapBlue = [one2zero' one2zero' one2targ'];

        targ2one = mytargintens:((1-mytargintens)/(mynposcolsteps-1)):1;
        zero2one =  0:(1/(mynposcolsteps-1)):1;
        mycolmapRed  = [targ2one' zero2one' zero2one'];

        mycolmap = [mycolmapBlue;mycolmapRed];

        set(st.fig,'ColorMap',[gray(64);mycolmap]);
        set(rivets.fig,'ColorMap',mycolmap);
end

set(handles.nut_ts_axes,'ButtonDownFcn','nut_results_viewer(''nut_ts_axes_ButtonDownFcn'',gcbf,[],guidata(gcbf))');
% we have to put this here, because MATLAB kicks it out each time we plot something
% in the axes...

return;

%%--------------------------------------
function paint_activation(hObject,eventdata,handles)
% get information needed to paint blobs and punt to nut_view_beamforming_activations

global rivets st beam;

time = str2double(get(handles.nut_time_text,'String'));
time2 = str2double(get(handles.nut_time2_text,'String'));
timept = dsearchn(beam.timepts,time);   % find nearest index in timevector to given time
timept2 = dsearchn(beam.timepts,time2);
if size(rivets.s,2)==1 % occurs when timefmode=0, and contrastselect is not Active Power
    timept=1;timept2=1;
end
% threshold = str2double(get(handles.nut_thresholdpos_text,'String'));
% nut_view_beamforming_activation(timept:timept2,get(findobj('Tag','nut_blobstyle'),'Value'), threshold/1000*beam.maxblob);
nut_view_beamforming_activation(timept:timept2,get(handles.nut_blobstyle,'Value'), rivets.threshold);
nut_image('shopos');
%set_spm_readouts; % jfh


%%------------------------------------
function nut_colorscale_menu_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
% if ispc
%     set(hObject,'BackgroundColor','white');
% else
%     set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
% end


%%-----------------------------------------------
function nut_colorscale_menu_Callback(hObject, eventdata, handles)
% --- Executes on selection change in nut_colorscale_menu.
global rivets

switch get(handles.nut_colorscale_menu,'Value')
    case 1  % Abs Max
        rivets.minblob=-rivets.maxtf;
        rivets.maxblob=rivets.maxtf;
        %rivets.scalefactor = 1000/max(abs(rivets.maxblob),abs(rivets.minblob));
    case 2  % Min-Max
        rivets.minblob=min(rivets.s(:));
        rivets.maxblob=max(rivets.s(:));
    case 3  % 0-Max
        rivets.minblob=0;
        rivets.maxblob=rivets.maxtf;
    case 4  % Min-0
        rivets.minblob=-rivets.maxtf;
        rivets.maxblob=0;
    case 5  % Abs Max for current band
        maxinband = max(max(abs(rivets.s(:,:,rivets.freqselect))));
        rivets.minblob=-maxinband;
        rivets.maxblob=maxinband;
end

if ~isempty(hObject)    % if user has clicked on menu we need to refresh display, otherwise this function is called by nut_refreshdisplay
    nut_refreshdisplay(handles);
end


% --- Executes on selection change in nut_style_menu.
function nut_style_menu_Callback(hObject, eventdata, handles)
% hObject    handle to nut_style_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns nut_style_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nut_style_menu
plot_ts
paint_activation(hObject,eventdata,handles)

% --- Executes during object creation, after setting all properties.
function nut_style_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_style_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
% if ispc
%     set(hObject,'BackgroundColor','white');
% else
%     set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
% end


% --------------------------------------------------------------------
function nut_exportEPS_menu_Callback(hObject, eventdata, handles)
% hObject    handle to nut_exportEPS_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[EPSfilename, EPSpathname]=uiputfile('*.eps','Save yer EPS...');
% [EPSfilename, EPSpathname]=uiputfile('*.tif','Save yer tif...');
if isequal(EPSfilename,0)|isequal(EPSpathname,0)
    return;
else
    set(gcbf,'PaperPositionMode','auto');
    set(gcbf,'PaperOrientation','portrait');
    set(gcbf,'InvertHardcopy','off');       % keep background color of figure
    set(handles.nut_megts_axes,'Visible','off');  % MEG time series are a pain in the ass when manipulating EPS files...
    set(get(handles.nut_megts_axes,'Children'),'Visible','off');  % MEG time series are a pain in the ass when manipulating EPS files...
    print(gcbf,'-depsc2','-adobecset','-r600','-painters','-noui',fullfile(EPSpathname,EPSfilename));
end

% --------------------------------------------------------------------
function nut_modifycoreg_menu_Callback(hObject, eventdata, handles)
% hObject    handle to nut_modifycoreg_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global nuts;
%
% clear global nuts;
global rivets beam st
% nuts.coreg.mripath = beam.coreg.mripath;
% nuts.coreg.meg2mri_tfm = beam.coreg.meg2mri_tfm;
beam.coreg = nut_CoregistrationTool(beam.coreg);
% uiwait(coregfig);

nut_spmfig_setup
rivets.ts_refresh_handle = @plot_ts;

% for i=1:3
%     SPM_axes_obj(i)=st.vols{1}.ax{i}.ax;
% end
%
% % setup time series refresh to occur when clicking on new MRI
% set(SPM_axes_obj,'ButtonDownFcn',rivets.SPM_axes_ButtonDownFcn);
% following jazz is needed because SPM's blob freaks out when given
% negative coordinates, and needs a dilation to know these are big voxels
% furthermore, coords to SPM must be integers.
rivets.voxelsMRI = nut_meg2mri(beam.voxels);
[rivets.voxelsblob,rivets.blob2mri_tfm]=nut_voxels2blob(beam);

paint_activation(hObject,eventdata,handles);

%nut_slider_set_inc; % jfh


% --------------------------------------------------------------------
function nut_activation_special_menu_Callback(hObject, eventdata, handles)
% hObject    handle to nut_activation_special_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in nut_snbeam_button.
function nut_snbeam_button_Callback(hObject, eventdata, handles)
% hObject    handle to nut_snbeam_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global rivets

settings=nut_normsettings_gui;
if isempty(settings), return, end   % User pressed cancel
if settings.all                     % Normalize all datasets
    iserr=nut_autonorm(settings.newvoxsize);
else                                % Normalize current dataset only
    set(handles.nut_thresholdpos_text,'String','0'); %we want to renorm all voxels.
    set(handles.nut_thresholdneg_text,'String','0');
    %rivets.threshold=zeros(size(rivets.voxelsMRI,1),2);
    nut_do_threshold(1,handles);
    plot_ts;
    paint_activation(hObject,eventdata,handles);
    iserr=nut_normalize_beam(rivets.beampath,[],settings.newvoxsize);
end
if ~iserr
    msgbox('Dataset(s) have been spatially normalized. Click on "Display spat norm" to view the result.')
end

% --- Executes on button press in nut_display_snbeam.
function nut_display_snbeam_Callback(hObject, eventdata, handles)
% hObject    handle to nut_display_snbeam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nut_display_snbeam
global rivets

switch get(hObject,'Value')
    case 1
        [fpath,ffile,fext]=fileparts(rivets.beampath);
        norm_sbeamname=fullfile(fpath,[ffile '_spatnorm.mat']);
        if ~exist(norm_sbeamname,'file')
            errordlg('You have to create a spatial normalization first (click on "Spatially normalize".')
            set(hObject,'Value',0)
        else
            nut_results_viewer(norm_sbeamname);
        end
    case 0
        sbeamfile=rivets.beampath;
        f=strfind(sbeamfile,'_spatnorm');
        if ~isempty(f)
            sbeamfile(f:f+8)=[];
            nut_results_viewer(sbeamfile);
        end
end
%clear global rivets beam
%nut_results_viewer(norm_sbeamname)
%uiwait

% global rivets beam
% rivets=rivets2; beam=beam2; clear rivets2 beam2
% spm_image('init',beam.coreg.mripath); % load/reload structural MRI
% nut_spmfig_setup
% %plot_ts; % erase any previous interval lines
% paint_activation(hObject,eventdata,handles)
% axes(handles.nut_ts_axes), colormap(jet)    % for some reason, it returns with a gray colormap sometimes...

%-------------------------------------------------------
function nut_reshaping_blobs(xyz,t,mat);
global rivets st
rcp      = round(xyz);
dim      = max(rcp,[],2)';
off      = rcp(1,:) + dim(1)*(rcp(2,:)-1 + dim(2)*(rcp(3,:)-1));
vol      = zeros(dim)+NaN;
vol(off) = t;
vol      = reshape(vol,dim);
st.vols{1}.blobs=cell(1,1);
% if st.mode == 0,
%     axpos = get(st.vols{1}.ax{2}.ax,'Position');
% else,
%     axpos = get(st.vols{1}.ax{1}.ax,'Position');
% end;
% ax = axes('Parent',st.fig,'Position',[(axpos(1)+axpos(3)+0.1) (axpos(2)+0.005) 0.05 (axpos(4)-0.01)],'Box','on');
% mx = max([eps max(t)]);
% mn = min([0 min(t)]);
%image([0 1],[mn mx],[1:64]' + 64,'Parent',ax);
%set(ax,'YDir','normal','XTickLabel',[]);
%st.vols{1}.blobs{1} = struct('vol',vol,'mat',mat,'cbar',ax,'max',mx, 'min',mn);
st.vols{1}.blobs{1} = struct('vol',vol,'mat',mat);




% --- Executes during object creation, after setting all properties.
function nut_view_time_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_view_time_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
% if ispc
%     set(hObject,'BackgroundColor','white');
% else
%     set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
% end



function nut_view_time_text_Callback(hObject, eventdata, handles)
% hObject    handle to nut_view_time_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_view_time_text as text
%        str2double(get(hObject,'String')) returns contents of nut_view_time_text as a double

nut_refreshdisplay(handles);


% --- Executes during object creation, after setting all properties.
function nut_view_time2_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_view_time2_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
% if ispc
%     set(hObject,'BackgroundColor','white');
% else
%     set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
% end



function nut_view_time2_text_Callback(hObject, eventdata, handles)
% hObject    handle to nut_view_time2_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_view_time2_text as text
%        str2double(get(hObject,'String')) returns contents of nut_view_time2_text as a double

nut_refreshdisplay(handles);


% --- Executes on button press in nut_norm_box.
function nut_norm_box_Callback(hObject, eventdata, handles)
% hObject    handle to nut_norm_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nut_norm_box
global rivets beam

if(get(hObject,'Value'))
    rivets.normflag = true;
    rivets.scalefactor = 1000/max(abs(rivets.maxblob),abs(rivets.minblob));
else
    rivets.normflag = false;
    rivets.scalefactor = 1;
end
if(get(handles.nut_threshold_menu,'Value')==1)
    thpos=str2double(get(handles.nut_thresholdpos_text,'String'))*rivets.scalefactor;
    thneg=str2double(get(handles.nut_thresholdneg_text,'String'))*rivets.scalefactor;
    set(handles.nut_thresholdpos_text,'String',num2str(thpos));
    set(handles.nut_thresholdneg_text,'String',num2str(thneg));
end
if ~isempty(get(handles.nut_view_scale_text,'String')) %rivets.setscale
    scneg=str2double(get(handles.nut_view_scale_text,'String'))*rivets.scalefactor;      
    set(handles.nut_view_scale_text,'String',num2str(scneg));
end
if ~isempty(get(handles.nut_view_scale2_text,'String')) 
    scpos=str2double(get(handles.nut_view_scale2_text,'String'))*rivets.scalefactor;  
    set(handles.nut_view_scale2_text,'String',num2str(scpos));
end
nut_refreshdisplay(handles);


% --- Executes on button press in nut_display_rawmeg_button.
function nut_display_rawmeg_button_Callback(hObject, eventdata, handles)
% hObject    handle to nut_display_rawmeg_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global rivets nuts
nut_importmeg;
set(handles.nut_display_rawmeg_button,'Visible','Off');
rivets.meg.data=nuts.meg.data;
plot_ts;


% --------------------------------------------------------------------
function MRIcroROI_menu_Callback(hObject, eventdata, handles)
% hObject    handle to MRIcroROI_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global st;

warning('repositioning to MRIcro ROI is totally experimental and probably inaccurate!');

[ROIfile, ROIpath] = uigetfile('*.roi','Select MRIcro ROI...');
if isequal(ROIfile,0)|isequal(ROIpath,0)
    return;
end
ROIfid = fopen(strcat(ROIpath,ROIfile),'r');
ROIshite = fread(ROIfid,'char');  % sometimes first 6 entries is FUBAR, let's hope last few are ok
pos=[ROIshite(end-3);ROIshite(end-2);ROIshite(1)];
pos = [258-pos(1);pos(2);pos(3)]; %%%%% 258 seems to agree better????
tmp = st.vols{1}.premul*st.vols{1}.mat;
tmp(1:3,:)*[pos ; 1];
pos = nut_coordtfm(pos,st.vols{1}.premul*st.vols{1}.mat);
% spm_orthviews('Reposition',pos); % place cursor on coords
nut_reposition(pos);
plot_ts; % refresh time series plot


% --- Executes on selection change in nut_activation_menu.
function nut_activation_menu_Callback(hObject, eventdata, handles)
% hObject    handle to nut_activation_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns nut_activation_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nut_activation_menu


% --- Executes during object creation, after setting all properties.
function nut_activation_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_activation_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --------------------------------------------------------------------
function nut_savebeam_menu_Callback(hObject, eventdata, handles)
% hObject    handle to nut_savebeam_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global beam rivets
[filename, pathname] = uiputfile('*.mat', 'Save Beamformer Volume...',rivets.beamfilename);  
if isequal(filename,0) | isequal(pathname,0)
    return;
else
    filename = fullfile(pathname,filename);
end

% if length(beam.s)==2
%     beam.s(2)=[];
% end
save(filename,'beam');
% if length(beam.s)==1
%     beam.s{2}=ones(size(beam.s{1}));
% end

% --- Executes on button press in nut_subtractnoise_box.
function nut_subtractnoise_box_Callback(hObject, eventdata, handles)
% hObject    handle to nut_subtractnoise_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nut_subtractnoise_box
nut_contrastselect_menu_Callback([],[],handles); % pick contrast method and create initial T-F plot

function nut_exclude_radius_text_Callback(hObject, eventdata, handles)
% hObject    handle to nut_exclude_radius_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_exclude_radius_text as text
%        str2double(get(hObject,'String')) returns contents of nut_exclude_radius_text as a double


% --- Executes during object creation, after setting all properties.
function nut_exclude_radius_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_exclude_radius_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end


function nut_within_radius_text_Callback(hObject, eventdata, handles)
% hObject    handle to nut_within_radius_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_within_radius_text as text
%        str2double(get(hObject,'String')) returns contents of nut_within_radius_text as a double


% --- Executes during object creation, after setting all properties.
function nut_within_radius_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_within_radius_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in nut_check_logscale.
function nut_check_logscale_Callback(hObject, eventdata, handles)
% hObject    handle to nut_check_logscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nut_check_logscale
plot_ts;

%======================
% Voxel Marker Section
%---------------------

% --- Executes on selection change in nut_getset_voinumber.
function nut_getset_voinumber_Callback(hObject, eventdata, handles)
% hObject    handle to nut_getset_voinumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns nut_getset_voinumber contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nut_getset_voinumber

global rivets

avoi=get(hObject,'Value');
nvoi=length(rivets.voi.coords);

if avoi>nvoi
    set([handles.nut_deletevoi handles.nut_displayvoi handles.nut_modifyvoi],'Enable','off')
    set(handles.nut_addvoi,'String','Add to all')
    set(handles.nut_voilabel,'String','')
else
    set(handles.nut_addvoi,'String','Insert to all')
    set([handles.nut_deletevoi handles.nut_displayvoi handles.nut_modifyvoi],'Enable','on')
    set(handles.nut_voilabel,'String',rivets.voi.coords(avoi).label)
end

% --- Executes during object creation, after setting all properties.
function nut_getset_voinumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_getset_voinumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in nut_displayvoi.
function nut_displayvoi_Callback(hObject, eventdata, handles)

global rivets

avoi=get(handles.nut_getset_voinumber,'Value');
ssubj=get(handles.nut_getset_patientnr,'String');
asubj=str2num(ssubj{get(handles.nut_getset_patientnr,'Value')});
scond=get(handles.nut_getset_conditionnr,'String');
acond=str2num(scond{get(handles.nut_getset_conditionnr,'Value')});
idx=find(rivets.voi.subjnr==asubj & rivets.voi.condnr==acond);

nut_reposition(rivets.voi.coords(avoi).MRI(idx,:))


% --- Executes on button press in nut_addvoi.
function nut_addvoi_Callback(hObject, eventdata, handles)
% hObject    handle to nut_addvoi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global beam rivets

avoi=get(handles.nut_getset_voinumber,'Value');
voi2replace=find([1:length(rivets.voi.coords)]>=avoi);

if ~isempty(voi2replace)
    rivets.voi.coords(voi2replace+1)=rivets.voi.coords(voi2replace);
end

nsubjcond=length(rivets.voi.subjnr);

MRI=spm_orthviews('pos')';
MEG=nut_mri2meg(MRI)./10;
rivets.voi.coords(avoi).MRI   = repmat(MRI,[nsubjcond 1]);
rivets.voi.coords(avoi).MEG   = repmat(MEG,[nsubjcond 1]);
index=rivets.currMEGvoxelindex;
rivets.voi.coords(avoi).index = repmat(index,[1 nsubjcond]);
rivets.voi.coords(avoi).label = get(handles.nut_voilabel,'String');

nvoi=length(rivets.voi.coords);
set(handles.nut_getset_voinumber,'string',[{rivets.voi.coords.label}';'New'], ...
    'Value',nvoi+1)
set(handles.nut_addvoi,'String','Add to all')
set([handles.nut_deletevoi handles.nut_displayvoi handles.nut_modifyvoi],'Enable','off')
set(handles.nut_voilabel,'String','')


% --- Executes on button press in nut_deletevoi.
function nut_deletevoi_Callback(hObject, eventdata, handles)
% hObject    handle to nut_deletevoi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global rivets

avoi=get(handles.nut_getset_voinumber,'Value');
rivets.voi.coords(avoi)=[];
nvoi=length(rivets.voi.coords);

if nvoi>0
    set(handles.nut_getset_voinumber,'string',[{rivets.voi.coords.label}';'New'])
    set(handles.nut_getset_voinumber,'Value',nvoi+1)
else
    set(handles.nut_getset_voinumber,'string','New')
    set(handles.nut_getset_voinumber,'Value',1)
end
set([handles.nut_deletevoi handles.nut_displayvoi handles.nut_modifyvoi],'Enable','off')
set(handles.nut_addvoi,'String','Add to all')
set(handles.nut_voilabel,'String','')

function voilabel_Callback(hObject, eventdata, handles)
% hObject    handle to voilabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of voilabel as text
%        str2double(get(hObject,'String')) returns contents of voilabel as a double


% --- Executes during object creation, after setting all properties.
function voilabel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to voilabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in nut_savevoi.
function nut_savevoi_Callback(hObject, eventdata, handles)
% hObject    handle to nut_savevoi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global rivets
voi=rivets.voi;
[filename,pathname]=uiputfile('*.mat','Save Pointer File as...');
if ischar(filename)
    save([pathname filename],'voi')
    msgbox(['Pointer file saved as ' filename])
end

% --- Executes on button press in nut_loadvoi.
function nut_loadvoi_Callback(hObject, eventdata, handles)
% hObject    handle to nut_loadvoi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global rivets

[filename,pathname]=uigetfile('*.mat','Load Pointer File...');
if ~ischar(filename), return, end
load([pathname filename]);
rivets.voi=voi; clear voi
nvoi=length(rivets.voi.coords);
if nvoi>0
    set(handles.nut_getset_voinumber,'string',[{rivets.voi.coords.label}';'New'])
    set(handles.nut_getset_voinumber,'Value',1)
    set(handles.nut_voilabel,'String',rivets.voi.coords(1).label)
    set([handles.nut_deletevoi handles.nut_displayvoi handles.nut_modifyvoi],'Enable','on')
end

idx=strmatch(rivets.beampath,rivets.voi.pathnames);
if ~isempty(idx)
    set(handles.nut_getset_patientnr,'string',[cellstr(int2str(unique(rivets.voi.subjnr)'));'New'])
    set(handles.nut_getset_conditionnr,'string',[cellstr(int2str(unique(rivets.voi.condnr)'));'New'])
    set(handles.nut_getset_patientnr,'Value',find(unique(rivets.voi.subjnr)==rivets.voi.subjnr(idx)))
    set(handles.nut_getset_conditionnr,'Value',find(unique(rivets.voi.condnr)==rivets.voi.condnr(idx)))
    if rivets.maxtf < rivets.voi.subjectmaxtf(rivets.voi.subjnr(idx))
        rivets.maxtf=rivets.voi.subjectmaxtf(rivets.voi.subjnr(idx));
        rivets.minblob=-rivets.maxtf; rivets.maxblob=rivets.maxtf;
        if rivets.timefmode
            caxis(handles.nut_ts_axes,[rivets.minblob rivets.maxblob]); caxis(rivets.cbar,[rivets.minblob rivets.maxblob]);
        end
        paint_activation(hObject,eventdata,handles);
    elseif rivets.maxtf > rivets.voi.subjectmaxtf(rivets.voi.subjnr(idx))
        rivets.voi.subjectmaxtf(rivets.voi.subjnr(idx))=rivets.maxtf;
    end
else
    nsubj=max(rivets.voi.subjnr)+1;
    msgbox(sprintf('The loaded pointer file does not know the current beamforming data. I therefore assume that it corresponds to the data of a new subject nr %d.',nsubj))
    rivets.voi.filenames={rivets.voi.filenames{:} rivets.beamfilename};
    rivets.voi.pathnames={rivets.voi.pathnames{:} rivets.beampath};
    rivets.voi.subjnr=[rivets.voi.subjnr nsubj];
    rivets.voi.condnr=[rivets.voi.condnr 1];
    rivets.voi.subjectmaxtf=[rivets.voi.subjectmaxtf rivets.maxtf];
    idx=find(rivets.voi.subjnr==nsubj-1 & rivets.voi.condnr==1);
    for k=1:nvoi;
        rivets.voi.coords(k).MRI(end+1,:)=rivets.voi.coords(k).MRI(idx,:);
        rivets.voi.coords(k).MEG(end+1,:)=rivets.voi.coords(k).MEG(idx,:);
        rivets.voi.coords(k).index(:,end+1)=rivets.voi.coords(k).index(:,idx);
    end
    set(handles.nut_getset_patientnr,'string',[cellstr(int2str(unique(rivets.voi.subjnr)'));'New'])
    set(handles.nut_getset_conditionnr,'string',[cellstr(int2str(unique(rivets.voi.condnr)'));'New'])
    set(handles.nut_getset_patientnr,'Value',find(unique(rivets.voi.subjnr)==nsubj))
    set(handles.nut_getset_conditionnr,'Value',1)
end

% --- Executes on button press in nut_voi_exportctf.
function nut_voi_exportctf_Callback(hObject, eventdata, handles)
% hObject    handle to nut_voi_exportctf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~exist('ctf_write_sam_targets')
    errordlg('This function requires the CTF toolbox.')
    return
end

global rivets

answer=questdlg('Which voxel markers to you want to export?','Export CTF SAM targets','Current subject and condition only','All subjects and conditions','Current subject and condition only');
switch answer
    case 'Current subject and condition only'
        asubj=unique(rivets.voi.subjnr); asubj=asubj(get(handles.nut_getset_patientnr,'Value'));
        acond=unique(rivets.voi.condnr); acond=acond(get(handles.nut_getset_conditionnr,'Value'));
        %avoi =get(handles.nut_getset_voinumber,'Value');
        idx=find(rivets.voi.subjnr==asubj & rivets.voi.condnr==acond);
        if length(idx)>1
            errordlg('Corrupted subject and condition numbers.'), return
        end
        vertices=cat(3,rivets.voi.coords.MEG);
        vertices=squeeze(vertices(idx,:,:))';
        if size(vertices,2)==1, vertices=vertices'; end     % prevents bug when only 1 VOI exported
        dum=pwd; cd(fileparts(rivets.voi.pathnames{idx}))
        [filename,pathname]=uiputfile('*','Save CTF SAM targets as...');
        if ischar(filename)
            ctf_write_sam_targets(vertices,[pathname filename]);
            msgbox(['SAM targets of current subject and condition saved as ' filename])
        end
        cd(dum)
    case 'All subjects and conditions'
        filename=inputdlg({'Save CTF SAM targets as: '},'',1,{'VOI'});
        if isempty(filename), return, end
        vertices=cat(3,rivets.voi.coords.MEG);
        for k=1:length(rivets.voi.pathnames)
            currvert=squeeze(vertices(k,:,:))';
            if size(currvert,2)==1, currvert=currvert'; end     % prevents bug when only 1 VOI exported
            ctf_write_sam_targets(currvert,sprintf('%s/%s',fileparts(rivets.voi.pathnames{k}),filename{:}));
        end
        msgbox(['SAM targets saved as ' filename{:} ' in SAM folder of each subject and condition.'])
end


% --- Executes on selection change in nut_getset_patientnr.
function nut_getset_patientnr_Callback(hObject, eventdata, handles)
% hObject    handle to nut_getset_patientnr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns nut_getset_patientnr contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nut_getset_patientnr
txt  = get(hObject,'String');
asubj= get(hObject,'Value');
if strcmpi(txt{asubj},'new')
    txt   = txt(1:end-1);
    asubj = int2str(max(str2num(char(txt)))+1);
    set(hObject,'string',[txt;asubj;'New'])
end


% --- Executes during object creation, after setting all properties.
function nut_getset_patientnr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_getset_patientnr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in nut_getset_conditionnr.
function nut_getset_conditionnr_Callback(hObject, eventdata, handles)
% hObject    handle to nut_getset_conditionnr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns nut_getset_conditionnr contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nut_getset_conditionnr
txt  = get(hObject,'String');
acond= get(hObject,'Value');
if strcmpi(txt{acond},'new')
    txt   = txt(1:end-1);
    acond = int2str(max(str2num(char(txt)))+1);
    set(hObject,'string',[txt;acond;'New'])
end

% --- Executes during object creation, after setting all properties.
function nut_getset_conditionnr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_getset_conditionnr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in nut_modifyvoi.
function nut_modifyvoi_Callback(hObject, eventdata, handles)
% hObject    handle to nut_modifyvoi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global rivets

avoi=get(handles.nut_getset_voinumber,'Value');
asubj=unique(rivets.voi.subjnr); asubj=asubj(get(handles.nut_getset_patientnr,'Value'));
acond=unique(rivets.voi.condnr); acond=acond(get(handles.nut_getset_conditionnr,'Value'));
idx=find(rivets.voi.subjnr==asubj & rivets.voi.condnr==acond);
if ~isempty(idx)
    MRI=spm_orthviews('pos')';
    MEG=nut_mri2meg(MRI);
    rivets.voi.coords(avoi).MRI(idx,:)   = MRI;
    rivets.voi.coords(avoi).MEG(idx,:)   = MEG;
    rivets.voi.coords(avoi).index(:,idx) = rivets.currMEGvoxelindex;
    rivets.voi.coords(avoi).label        = get(handles.nut_voilabel,'String');
end
set(handles.nut_getset_voinumber,'string',[{rivets.voi.coords.label}';'New'])

% --- Executes on button press in nut_load_subjcond.
function nut_load_subjcond_Callback(hObject, eventdata, handles)
% hObject    handle to nut_load_subjcond (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global rivets

exsubj=unique(rivets.voi.subjnr);
valsubj=get(handles.nut_getset_patientnr,'Value');
dsubj=valsubj-length(exsubj);

excond=unique(rivets.voi.condnr);
valcond=get(handles.nut_getset_conditionnr,'Value');
dcond=valcond-length(excond);

if dsubj>0, asubj=max(exsubj)+dsubj;
else, asubj=exsubj(valsubj);
end
if dcond>0, acond=max(excond)+dcond;
else acond=excond(valcond);
end
idx=find(rivets.voi.subjnr==asubj & rivets.voi.condnr==acond);
if ~isempty(idx)
    if strcmp(rivets.voi.pathnames{idx},rivets.beampath)
        errordlg('File is already open.')
    else
        nut_results_viewer(rivets.voi.pathnames{idx})
    end
else
    [filename,pathname]=uigetfile('s_beam*.mat','Load beamforming activations...');
    filename=filename(1:end-4);     % Remove .mat extension
    if ~ischar(filename), return, end
    pathname=fullfile(pathname,filename);
    rivets.voi.condnr=[rivets.voi.condnr acond];
    rivets.voi.subjnr=[rivets.voi.subjnr asubj];
    rivets.voi.filenames={rivets.voi.filenames{:} filename};
    rivets.voi.pathnames={rivets.voi.pathnames{:} pathname};
    for k=1:length(rivets.voi.coords);
        rivets.voi.coords(k).MRI(end+1,:)=rivets.voi.coords(k).MRI(end,:);
        rivets.voi.coords(k).MEG(end+1,:)=rivets.voi.coords(k).MEG(end,:);
        rivets.voi.coords(k).index(:,end+1)=rivets.voi.coords(k).index(:,end);
    end
    nut_results_viewer(pathname)
end


% --- Executes on button press in nut_clear_subjcond.
function nut_clear_subjcond_Callback(hObject, eventdata, handles, doquest)
% hObject    handle to nut_clear_subjcond (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global rivets beam
if nargin<4 || doquest
    isyes=questdlg('Are you sure you want to remove the markers of all subjects, conditions, and voxels from workspace?','Clear','Yes, go ahead','Oh no','Oh no');
else
    isyes='Yes, go ahead';
end

if strcmp(isyes,'Yes, go ahead')
    %    rivets.voi=struct('subjnr',1,'condnr',1,'filenames',{{rivets.beamfilename}}, ...
    %    'pathnames',{{rivets.beampath}},'coords',[],'subjectmaxtf',max(abs(beam.s{1}(:))));
    rivets.voi=struct('subjnr',1,'condnr',1,'filenames',{{rivets.beamfilename}}, ...
        'pathnames',{{rivets.beampath}},'coords',[],'subjectmaxtf',rivets.maxtf);
    set(handles.nut_getset_patientnr,'string',{'1';'New'},'Value',1)
    set(handles.nut_getset_conditionnr,'string',{'1';'New'},'Value',1)
    set(handles.nut_getset_voinumber,'string','New','Value',1)
    set([handles.nut_deletevoi handles.nut_displayvoi handles.nut_modifyvoi],'Enable','off')
    set(handles.nut_addvoi,'String','Add to all')
    set(handles.nut_voilabel,'String','Label')
end



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in nut_velec_button.
function nut_velec_button_Callback(hObject, eventdata, handles)
% hObject    handle to nut_velec_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% load raw meg? or filtered?
% load GP1motorbeta
load GP1motor
[MEGvoxelindex,distances] = dsearchn(rivets.voxelsMRI,nut_coordtfm(cursorpos,inv(st.vols{1}.premul)));
Lp = nuts.Lp(:,:,MEGvoxelindex);
beam.params.cn=1;
Rall = nut_tfcov(nuts.meg,[-1000 1000]);
[W,eta] = nut_SAM(Lp,Rall,beam.params.cn);
S = zeros(size(nuts.meg.data,1),size(nuts.meg.data,3),'single');
for ii=1:size(nuts.meg.data,3)
    S(:,ii) = nuts.meg.data(:,:,ii)*W;
end
save motor_virtualelec.mat S
ss = reshape(S,1,size(S,1)*size(S,2));
[ersp,itc,powbase,times,freqs,erspboot] = timef(ss,3000,[-1000 1500],1200,0,'baseline',-500,'maxfreq',150);


% [u,s,v]=svd(Rall);
% sig=s(end,end);
% noise = sum(W.*(sig*W))';



function [W,eta] = nut_SAM(Lp, Rall, cn)

InvRall = inv(Rall);

flags.LCMVcn = cn;
flags.wn = false;
flags.dualstate = true;
flags

[W, eta] = nut_LCMV_Beamformer(Lp,InvRall, flags);


% --------------------------------------------------------------------
function nut_export_ris_Callback(hObject, eventdata, handles)
% hObject    handle to nut_export_ris (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function nut_export_activations_Callback(hObject, eventdata, handles)
% hObject    handle to nut_export_activations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function nut_view_scale_text_Callback(hObject, eventdata, handles)
% hObject    handle to nut_view_scale_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_view_scale_text as text
%        str2double(get(hObject,'String')) returns contents of nut_view_scale_text as a double
global rivets
% rivets.setscale = ~isempty(get(hObject,'string'));  % if user deletes value, scale automatically
plot_ts;
paint_activation(hObject,eventdata,handles);
return;


% --- Executes during object creation, after setting all properties.
function nut_view_scale_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_view_scale_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nut_view_scale2_text_Callback(hObject, eventdata, handles)
% hObject    handle to nut_view_scale2_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_view_scale2_text as text
%        str2double(get(hObject,'String')) returns contents of nut_view_scale2_text as a double
global rivets
% rivets.setscale = ~isempty(get(hObject,'string'));  % if user deletes value, scale automatically
plot_ts;
paint_activation(hObject,eventdata,handles);
return;


% --- Executes during object creation, after setting all properties.
function nut_view_scale2_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_view_scale2_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in nut_view_scale_button.
function nut_view_scale_button_Callback(hObject, eventdata, handles)
% hObject    handle to nut_view_scale_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global rivets beam;

% initialize times
set(handles.nut_view_scale_text,'String',beam.timepts(1));
set(handles.nut_view_scale2_text,'String',beam.timepts(end));
plot_ts; % erase any previous interval lines

[time,crap,button]=ginput(1);
timept = dsearchn(beam.timepts,time);
set(handles.nut_view_time_text,'String',timept);
% set(handles.nut_view_time2_text,'String',beam.timepts(end));

% now for the second selection
button = 0;
while(button~=1 && button~=3)
    [time2,ypos,button]=ginput(1);
    if (button==3)   % just use first time point if right button clicked
        timept2 = timept;
        %     elseif (button==2)  % if middle button clicked, set threshold
        %         threshold = ypos;
        %         set(handles.nut_thresholdpos_text,'String',threshold);
    else  % select second time point for other buttons
        timept2 = dsearchn(beam.timepts,time2);
    end
end

set(handles.nut_view_time2_text,'String',timept2);
plot_ts;
return;


% --- Executes on button press in activation_bpf.
function checkbox20_Callback(hObject, eventdata, handles)
% hObject    handle to activation_bpf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of activation_bpf



function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to nut_lpf_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_lpf_text as text
%        str2double(get(hObject,'String')) returns contents of nut_lpf_text as a double


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_lpf_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to nut_hpf_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_hpf_text as text
%        str2double(get(hObject,'String')) returns contents of nut_hpf_text as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_hpf_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% % --- Executes on button press in nut_applyfilter_button.
% function pushbutton50_Callback(hObject, eventdata, handles)
% % hObject    handle to nut_applyfilter_button (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
%

% --- Executes on selection change in activation_notch.
function activation_notch_Callback(hObject, eventdata, handles)
% hObject    handle to activation_notch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns activation_notch contents as cell array
%        contents{get(hObject,'Value')} returns selected item from activation_notch
plot_ts;
return;


% --- Executes during object creation, after setting all properties.
function activation_notch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to activation_notch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in nut_orientated_menu.
function nut_orientated_menu_Callback(hObject, eventdata, handles)
% hObject    handle to nut_orientated_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns nut_orientated_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nut_orientated_menu
nut_contrastselect_menu_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function nut_orientated_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_orientated_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in nut_save_ts_roi_button.
function nut_save_ts_roi_button_Callback(hObject, eventdata, handles)
% hObject    handle to nut_save_ts_roi_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%% note: this is from activation_viewer, probably obsolete.  %%%

% --- Executes on button press in nut_save_ts_roi_button.
make_ps=false;

dir_listing = dir('*.roi');
roifiles=cell(size(dir_listing));
[roifiles{:}]=deal(dir_listing.name);

global rivets st beam;

for i=1:length(roifiles)
    ROIfid = fopen(roifiles{i},'r');
    ROIshite = fread(ROIfid,'char');  % sometimes first 6 entries is FUBAR, let's hope last few are ok
    pos=[ROIshite(end-3);ROIshite(end-2);ROIshite(1)];
    pos = [258-pos(1);pos(2);pos(3)]; %%%%% TEMP. HACK TO GO FROM MARYAM UNITS TO SARANG UNITS
    tmp = st.vols{1}.premul*st.vols{1}.mat;
    pos = tmp(1:3,:)*[pos ; 1];
    %     spm_orthviews('Reposition',pos); % place cursor on coords
    nut_reposition(pos);
    plot_ts; % refresh time series plot


    cursorpos = spm_orthviews('pos')';
    MEGvoxelindex = dsearchn(beam.voxelsblob,cursorpos);
    if(beam.display_s_II_flag)
        MEGts = squeeze(rivets.s_II(MEGvoxelindex,:));
    else
        MEGts = squeeze(rivets.s_beam(MEGvoxelindex,:));
    end
    if(roifiles{i}(8) == '_')
        if(make_ps)
            spm_print;
            unix(['mv spm2.ps meg_' roifiles{i}(6:7) '.ps']);
        else
            save(['meg_' roifiles{i}(6:7) '.mat'],'timevector','MEGts');
        end
    else
        if(make_ps)
            spm_print;
            unix(['mv spm2.ps meg_' roifiles{i}(6:8) '.ps']);
        else
            save(['meg_' roifiles{i}(6:8) '.mat'],'timevector','MEGts');
        end
    end
end
return;


% --- Executes on selection change in nut_noisetype.
function nut_noisetype_Callback(hObject, eventdata, handles)
% hObject    handle to nut_noisetype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns nut_noisetype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nut_noisetype

nut_contrastselect_menu_Callback(hObject,eventdata,handles);


% --- Executes during object creation, after setting all properties.
function nut_noisetype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_noisetype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in nut_export2DV3D.
function nut_export2DV3D_Callback(hObject, eventdata, handles)
% hObject    handle to nut_export2DV3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global beam rivets ndefaults

nut_mat2img(rivets.timeselect,'nut2dv3d-tmp.img');
[status,result]=system(['cd ' ndefaults.dv3dpath '; python -m RunDV3D.py -base ' beam.coreg.mripath ' -overlay ' pwd '/nut2dv3d-tmp.img']);
if(status)
    warning('python and/or RunDV3D could not be successfully run... are your system and python paths properly set?')
    warning(result)
end


% --- Executes on button press in nut_fixvoxel.
function nut_fixvoxel_Callback(hObject, eventdata, handles)
% hObject    handle to nut_fixvoxel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global rivets beam

if get(hObject,'Value')
    hw = waitbar(0,'Seeding selected voxel/ROI...');
    load(beam.connectionfile);
    
    % Legacy compatibility
    if exist('IC','var')
        CC = IC; clear IC        
        CC.method='icohere';
    elseif exist('CC','var') 
        if ~isfield(CC,'method'), CC.method='ccohere'; end
    else
        error('NUT_RESULTS_VIEWER: Invalid connection file.')
    end
    if size(CC.coh,1)==size(CC.frq,1)     
        CC.coh=permute(CC.coh,[2 3 1]); 
    end    
    [nc,nt,nf] = size(CC.coh);
    if nt==1 && nf>1 && size(beam.s{1},2)==nf  % if spectrogram
        CC.coh=permute(CC.coh,[1 3 2]);
        nt=nf; nf=1;
    end
    nv=size(beam.voxels,1);
    
    useroi=isfield(beam,'rois');
    if useroi
        waitbar(.3,hw);
        nr = length(beam.R.goodroi);
        sv = rivets.currMEGvoxelindex;
        sr = find(beam.R.roi2voxel_tfm(sv,:),1);
        ss = find(beam.R.roi2voxel_tfm(:,sr));
        nss=length(ss);
        g=(any(CC.comps==sr,2));
        [II{1:length(beam.s)}] = deal(zeros(nr,nt,nf));
        switch CC.method
            case 'ccohere'
                II{1}(setdiff(1:nr,sr),:,:) = abs(imag(CC.coh(g,:,:)));     II{1}(isnan(II{1}))=0;            
                II{2}(setdiff(1:nr,sr),:,:) = abs(CC.coh(g,:,:)).^2;        II{2}(isnan(II{2}))=0;
                II{3}(setdiff(1:nr,sr),:,:) = abs(real(CC.coh(g,:,:)));     II{3}(isnan(II{3}))=0;               
            otherwise
                II{1}(setdiff(1:nr,sr),:,:) = abs(CC.coh(g,:,:));           II{1}(isnan(II{1}))=0;
        end
        
    else
    
        RAD=str2num(get(handles.nut_edit_roiradius,'string')); %#ok<ST2NM>
        if( ~isempty(RAD) && RAD>0 )
            waitbar(.4,hw);

            CC.coh = fcm_fisherZtf(CC.coh,CC.method,1);         % Fisher Z-transform
            waitbar(.75,hw);

            distances = nut_rownorm(nut_coord_diff(beam.voxels,beam.voxels(rivets.currMEGvoxelindex,:)));
            ss = find(distances<RAD); % seed voxels
            nss= length(ss);
            T = zeros(nv,nt,nf);
            for k=1:nss
                g = (any(CC.comps==ss(k),2));
                temp = nan(nv,nt,nf);
                temp(setdiff(1:nv,ss(k)),:,:) = CC.coh(g,:,:);
                T = T + temp;
            end
            T(ss,:,:)=[];

            % Average
            T = T / nss;
            waitbar(.9,hw);

            % Inverse Fisher Z-transform and extraction of selected component
            T = fcm_fisherZtf(T,CC.method,-1);
            switch CC.method
                case {'ccohere' 'nccohere'}
                    II{1} = imag(T);
                    II{2} = abs(T).^2;
                    II{3} = real(T);
                otherwise
                    II = {T};
            end
            clear T

        else
            waitbar(.5,hw);
            ss = rivets.currMEGvoxelindex;
            g=(any(CC.comps==ss,2));
            nss=1;
            switch CC.method
                case {'ccohere' 'nccohere'}
                    II{1} = abs(imag(CC.coh(g,:,:)));
                    II{2} = abs(CC.coh(g,:,:)).^2;
                    II{3} = abs(real(CC.coh(g,:,:)));                    
                otherwise
                    II = {abs(CC.coh(g,:,:))};
            end
        end
    end
    
    for k=1:length(II)
        if useroi
            for k3=1:nf
                beam.s{k}(:,:,k3) = beam.R.roi2voxel_tfm * II{k}(:,:,k3);
            end
        else
            beam.s{k}(setdiff(1:nv,ss),:,:) = II{k};
        end
        maxconn = max(beam.s{k},[],1);
        beam.s{k}(ss,:,:)=repmat(maxconn,[nss 1 1]);
    end
    set(handles.nut_snbeam_button,'Enable','off')    %Make sure spatnorm is not created with wrong data 
    waitbar(1,hw);
    close(hw);    
    
else
    beam=load(rivets.beampath);  % reset
    beam=nut_beam_legacy_compatibility(beam);
    set(handles.nut_snbeam_button,'Enable','on')   
end

nut_contrastselect_menu_Callback([],[],handles);
nut_refreshdisplay(handles);

    
%-------------------------------
function nut_refreshdisplay(handles)

global rivets beam

nut_do_threshold([],handles);
nut_colorscale_menu_Callback([],[],handles);
switch(rivets.displaymode)
    case 1
        plot_ts; % erase any previous interval lines
        paint_activation([],[],handles)
    case 2
        spm_image('init',beam.coreg.mripath); % load/reload structural MRI;
        nut_spmfig_setup
        plot_ts; % erase any previous interval lines
        paint_activation([],[],handles)
        nut_render_meg(beam.coreg.brainrender_path);
end

function nut_edit_roiradius_Callback(hObject, eventdata, handles)
% hObject    handle to nut_edit_roiradius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_edit_roiradius as text
%        str2double(get(hObject,'String')) returns contents of nut_edit_roiradius as a double


% --- Executes during object creation, after setting all properties.
function nut_edit_roiradius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_edit_roiradius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nut_condcontrast1_text_Callback(hObject, eventdata, handles)
% hObject    handle to nut_condcontrast1_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_condcontrast1_text as text
%        str2double(get(hObject,'String')) returns contents of nut_condcontrast1_text as a double


% --- Executes during object creation, after setting all properties.
function nut_condcontrast1_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_condcontrast1_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nut_condcontrast2_text_Callback(hObject, eventdata, handles)
% hObject    handle to nut_condcontrast2_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_condcontrast2_text as text
%        str2double(get(hObject,'String')) returns contents of nut_condcontrast2_text as a double


% --- Executes during object creation, after setting all properties.
function nut_condcontrast2_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_condcontrast2_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
