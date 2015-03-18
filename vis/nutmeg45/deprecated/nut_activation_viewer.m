function varargout = nut_activation_viewer(varargin)
% NUT_ACTIVATION_VIEWER M-file for nut_activation_viewer.fig
%

% Last Modified by GUIDE v2.5 04-Aug-2008 15:04:32

% Here and below:
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to nut_activation_viewer (see VARARGIN)

if(~strcmp('14',version('-release')))
    warning('off')  % hack so warning for "created by matlab 7.0" is turned off
end
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @nut_activation_viewer_OpeningFcn, ...
    'gui_OutputFcn',  @nut_activation_viewer_OutputFcn, ...
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
function nut_activation_viewer_OpeningFcn(hObject, eventdata, handles, varargin)
% --- Executes just before nut_activation_viewer is made visible.

% Choose default command line output for nut_activation_viewer
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

global rivets st nuts beam ndefaults;

rivets.fig = hObject;

if(isempty(varargin)) %occurs if not using matlab 7.0
    global bolts
    if isfield(bolts,'outname')
        default_beamfile=bolts.outname; 
    elseif(isfield(nuts,'megfile'));
        default_beamfile = ['s_beam_' nuts.meg.filename '.mat'];
    else
        default_beamfile = '';
    end
    clear global bolts;pack;

    [beamfilename, beampathname]=uigetfile('s_beam*.mat','Select MEG reconstruction file...',default_beamfile);
    if isequal(beamfilename,0)|isequal(beampathname,0)
        delete(hObject);
        return;
    else
        beampath = fullfile(beampathname,beamfilename);
    end
else  % if varargin is not empty, it contains beampath
    beampath=varargin{1};
    [beampathname,beamfilename] = fileparts(beampath);
end

rivets.beamfilename = beamfilename;

set(hObject,'Name',['NUTMEG Activation Viewer - ' rivets.beamfilename]);

barhandle = waitbar(0,['Preparing to display results...']);

beam=load(beampath); % should load "beam" structure
if isfield(beam,'beam')
    beam=beam.beam;
end
rivets.beampath = beampath; % remember, load(beampath) will clear anything previously in the beam structure!
clear beampath beampathname beamfilename;
% set all experimental flags here...
rivets.normflag=true;  %option to normalize to 0:1000 or leave s_beam^2 values as is.  
rivets.oldnorm=false; 
rivets.display_s_II_flag=false;  % assume we just want s_beam to start
rivets.display_s_perp_flag=false; % does nothing unless s_II is enabled
rivets.displaymode = 1; % 1 -> normal blob overlay on orthogonal MRI slices, 2-> overlay on rendered brain
if(ndefaults.vwr.plotweights)
    rivets.plotweights = true;
    rivets.W_fig = figure;
    rivets.W_ax = axes('Parent',rivets.W_fig);
end

waitbar(0.3,barhandle);

spm_image('init',beam.coreg.mripath); % load/reload structural MRI; clears VOI emphasis from nut_view_scan_region.m

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

% fit to desired MEG voxel grid, shift such that origin is (0,0,0)mm
% rivets.voxelsMRI = nut_coordtfm(beam.voxels,beam.coreg.meg2mri_tfm);
rivets.voxelsMRI = nut_meg2mri(double(beam.voxels));

% following jazz is needed because SPM's blob freaks out when given
% negative coordinates, and needs a dilation to know these are big voxels
% furthermore, coords to SPM must be integers.
if(beam.voxelsize(1,1)<0)
    beam.voxelsize = abs(beam.voxelsize);   % to account for negative voxelsizes (radiological???)
    warndlg('negative voxel size??? whassup with that?!');
end

if(round(rivets.voxelsMRI)-rivets.voxelsMRI<1e-4)
    rivets.blob2mri_tfm = [ beam.voxelsize(1)                0             0 beam.voxelsize(1)*floor(st.bb(1,1)/beam.voxelsize(1))-30
                                       0 beam.voxelsize(2)                 0 beam.voxelsize(2)*floor(st.bb(1,2)/beam.voxelsize(2))-30
                                       0                 0 beam.voxelsize(3) beam.voxelsize(3)*floor(st.bb(1,3)/beam.voxelsize(3))-30
                                       0                 0                 0                                                     1 ];
    rivets.voxelsblob = nut_coordtfm(rivets.voxelsMRI,inv(rivets.blob2mri_tfm));
elseif(round(beam.voxels)-beam.voxels<1e-12)
    rivets.blob2meg_tfm = [ beam.voxelsize(1)           0                  0   min(beam.voxels(:))-beam.voxelsize(1)-30
                                       0 beam.voxelsize(2)                 0      min(beam.voxels(:))-beam.voxelsize(2)-30
                                       0                 0 beam.voxelsize(3)      min(beam.voxels(:))-beam.voxelsize(3)-30
                                       0                 0                 0                                                     1 ];
    rivets.voxelsblob = double(nut_coordtfm(beam.voxels,inv(rivets.blob2meg_tfm)));
    rivets.blob2mri_tfm=beam.coreg.meg2mri_tfm*rivets.blob2meg_tfm;
else
    warning('all your voxels are not in either MRI or MEG space??')
    return
end

rivets.blob2mri_tfm=double(rivets.blob2mri_tfm);

clear voxels crap;

% compatibility and legacy support...
if(isfield(beam,'s'))
    beam.s_th = beam.s{1};
    if length(beam.s)==2
        beam.s_ph=beam.s{2};
    else
        beam.s_ph=0;
    end
%     beam.timewindow = beam.timepts; % we should make beam.timepts the new way, but later...
elseif(~isfield(beam,'s_th'))
    % only have weights, so need to load MEG data set and compute s_th and s_ph
    nut_importmeg(beam.meg.filename);
    %data = nuts.meg.data(:,nuts.meg.goodchannels,:);
    meg = nut_filter(nuts.meg.data(:,beam.meg.goodchannels,:),nuts.meg.latency,nuts.meg.srate,beam.params.preprocessing.baseline,beam.params.preprocessing.notch,beam.params.preprocessing.bpf,beam.params.preprocessing.bpf_low_cutoff,beam.params.preprocessing.bpf_high_cutoff);

    beam.s_th = transpose(squeeze(beam.W(:,1,:)))*meg';
    beam.s_ph = transpose(squeeze(beam.W(:,2,:)))*meg';
end
% We have s_th and s_ph, so we must compute s_beam (power) here:
if(isfield(beam,'s_z'))
    rivets.s_beam = double(beam.s_th.^2 + beam.s_ph.^2 + beam.s_z.^2);
elseif(isfield(beam,'s_th'))
    rivets.s_beam = double(beam.s_th.^2 + beam.s_ph.^2);
end

rivets.minblob=min(rivets.s_beam(:));
rivets.maxblob=max(rivets.s_beam(:));
if(rivets.normflag)
    rivets.scalefactor = 1000/max(abs(rivets.maxblob),abs(rivets.minblob));
else
    rivets.scalefactor = 1;
end

set(handles.nut_view_time_text,'String',beam.timewindow(1));
set(handles.nut_view_time2_text,'String',beam.timewindow(end));

global interval_select_handle; % used by ts_axes ButtonDownFcn
interval_select_handle = @interval_button_Callback;

waitbar(0.7,barhandle);
nut_spmfig_setup;
rivets.ts_refresh_handle = @plot_ts;  %% pass handle since technically it's a private function

waitbar(1,barhandle);
close(barhandle);

if(isfield(nuts,'meg'))
    set(handles.nut_megts_axes,'Visible','on');
    axes(handles.nut_megts_axes);
    if ~isfield(nuts,'preprocessing')
        nut_preprocessing_defaults
    end
    rivets.meg.data=nut_filter(nuts.meg.data,nuts.meg.latency,nuts.meg.srate,nuts.preprocessing.baseline,nuts.preprocessing.notch,nuts.preprocessing.bpf,nuts.preprocessing.bpf_low_cutoff,nuts.preprocessing.bpf_high_cutoff);
%     rivets.meg.data=mean(nuts.meg.data,3);
%     plot(beam.timewindow,1e15*rivets.meg.data);
    if isfield(nuts.meg,'eegflag') & nuts.meg.eegflag    % NUTEEG mod
        plot(nuts.meg.latency,1e6*rivets.meg.data);
        ylabel('Electric Potential (uV)');
    else
        plot(nuts.meg.latency,1e15*rivets.meg.data);
        ylabel('Magnetic Field (fT)');
    end
    axis tight; xlabel('Time (ms)');
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
rivets.graygraybipolar = [.8*gray(64); graysymmetric];

rivets.grayjet = [gray(64);jet(64)]; % gray-jet colormap

tmpsquish=tmp(1:2:end,:);
hotsymmetric = [flipud(tmpsquish); tmpsquish];
rivets.grayhotsymmetric=[gray(64);hotsymmetric];

% set initial display to max activation
rivets.threshold = [repmat(500,size(beam.voxels,1),1) repmat(-500,size(beam.voxels,1),1)];
nut_maxactivation_button_Callback(hObject,eventdata,handles);
% nut_style_menu_Callback(hObject,eventdata,handles)



%%---------------------------------------------------
function varargout = nut_activation_viewer_OutputFcn(hObject, eventdata, handles)
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


function nut_spmfig_setup
global st rivets beam
%%%%%%%%%%%%% new code to add boxes for MNI coords, MEG coords, and
%%%%%%%%%%%%% beamformer intensity values
fg = spm_figure('GetWin','Graphics');
WS = spm('winscale');

% reposition MRI intensity box
inlabel = findobj('String','Intensity:','HorizontalAlignment','center');

% get rid of some stuff, sweep other stuff under the rug
delete(findobj('String','Crosshair Position'));
delete(findobj('ToolTipString','move crosshairs to origin'));
delete(inlabel);
set(st.in,'Visible','off');

beaminlabel = uicontrol(fg,'Style','Text','Position',[60 350 120 020].*WS,'String','Activation Intensity:');
st.beamin = uicontrol(fg,'Style','edit', 'Position',[175 350  85 020].*WS,'String','');

uicontrol(fg,'Style','Text', 'Position',[75 255 35 020].*WS,'String','MNI:');
st.mnip = uicontrol(fg,'Style','edit', 'Position',[110 255 135 020].*WS,'String','','Callback','nut_image(''setposmni'')','ToolTipString','move crosshairs to MNI mm coordinates');
st.mnilabel=uicontrol('Style','text','BackgroundColor',[1 1 1],'Units','normalized','Position',[.65 .55 .3 .1],'FontSize',12);

uicontrol(fg,'Style','Text', 'Position',[75 315 35 020].*WS,'String','MEG:');
st.megp = uicontrol(fg,'Style','edit', 'Position',[110 315 135 020].*WS,'String',sprintf('%.1f %.1f %.1f',nut_mri2meg(spm_orthviews('pos')')),'Callback','nut_image(''setposmeg'')','ToolTipString','move crosshairs to MEG mm coordinates');

set(st.mp,'Callback','spm_image(''setposmm''); nut_image(''shopos'');');
set(st.vp,'Callback','spm_image(''setposvx''); nut_image(''shopos'');');

%add sliders for scrolling through MRI
ax_pos = get(st.vols{1}.ax{1}.ax,'Position');
cor_pos = get(st.vols{1}.ax{2}.ax,'Position');
sag_pos = get(st.vols{1}.ax{3}.ax,'Position');
slidercallback1 = 'global st; pos=str2num(get(st.vp,''String'')); set(st.vp,''String'',num2str([st.vols{1}.dim(1)*get(gcbo,''Value'') pos(2) pos(3)])); spm_image(''setposvx''); nut_image(''shopos'')';
slidercallback2 = 'global st; pos=str2num(get(st.vp,''String'')); set(st.vp,''String'',num2str([pos(1) st.vols{1}.dim(2)*get(gcbo,''Value'') pos(3)])); spm_image(''setposvx''); nut_image(''shopos'')';
slidercallback3 = 'global st; pos=str2num(get(st.vp,''String'')); set(st.vp,''String'',num2str([pos(1) pos(2) st.vols{1}.dim(3)*get(gcbo,''Value'')])); spm_image(''setposvx''); nut_image(''shopos'')';
dim_mm = st.vols{1}.dim .* diag(st.vols{1}.mat)';
st.slider{1}=uicontrol(fg,'Style','Slider','Callback',slidercallback1,'SliderStep',[beam.voxelsize(1)/abs(dim_mm(1)) .1],'Units','normalized','Position',[ax_pos(1) ax_pos(2)+ax_pos(4) ax_pos(3) .015]);
st.slider{2}=uicontrol(fg,'Style','Slider','Callback',slidercallback2,'SliderStep',[beam.voxelsize(2)/dim_mm(2) .1],'Units','normalized','Position',[ax_pos(1)+ax_pos(3) ax_pos(2) .02 ax_pos(4)]);
st.slider{3}=uicontrol(fg,'Style','Slider','Callback',slidercallback3,'SliderStep',[beam.voxelsize(3)/dim_mm(3) .1],'Units','normalized','Position',[cor_pos(1)+cor_pos(3) cor_pos(2) .02 cor_pos(4)]);


%                 set(st.slider{i},'Value',max((posmrimm(i)-st.bb(1,i))/dim(i),0));

% tell time series to refresh after clicking around in SPM volume
rivets.ts_refresh_handle = @plot_ts;  %% pass handle since technically it's a private function
for i=1:3  % step through three orthogonal views
    SPM_axes_obj(i) = st.vols{1}.ax{i}.ax;
end


SPM_axes_ButtonDownFcn = [get(SPM_axes_obj(1),'ButtonDownFcn') 'nut_image(''shopos'');'];
set(SPM_axes_obj,'ButtonDownFcn',SPM_axes_ButtonDownFcn);

if(isfield(beam.coreg,'norm_mripath') || ~isempty(regexp(beam.coreg.mripath,[filesep 'w[^' filesep ']+\.img'])) || ~isempty(findstr(beam.coreg.mripath,['templates' filesep 'T1.img'])))
    load('ihb_DataBase.cdb','-mat');
    rivets.MNIdb = MNIdb;
    [meshx,meshy,meshz]=ndgrid(MNIdb.minX:MNIdb.voxX:MNIdb.maxX,MNIdb.minY:MNIdb.voxY:MNIdb.maxY,MNIdb.minZ:MNIdb.voxZ:MNIdb.maxZ);
    rivets.MNIdb.coords = [meshx(:) meshy(:) meshz(:)];
end


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
% time2 = str2double(get(handles.nut_time2_text,'String'));
set(handles.nut_time2_text,'String',num2str(time));

switch rivets.displaymode
    case 1
        plot_ts;
        paint_activation(hObject,eventdata,handles);
    case 2
        spm_image('init',beam.coreg.mripath); % load/reload structural MRI;
        nut_spmfig_setup
        plot_ts; % erase any previous interval lines
        paint_activation([],[],handles)
        nut_render_meg(beam.coreg.brainrender_path);
end
return


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
switch(rivets.displaymode)
    case 1
        plot_ts;
        paint_activation(hObject,eventdata,handles);
    case 2
        spm_image('init',beam.coreg.mripath); % load/reload structural MRI;
        nut_spmfig_setup
        plot_ts; % erase any previous interval lines
        paint_activation([],[],handles)
        nut_render_meg(beam.coreg.brainrender_path);
end



%%-----------------------------------------------
function nut_threshold_text_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
% if ispc
%     set(hObject,'BackgroundColor','white');
% else
%     set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
% end


%%-----------------------------------------------
function nut_threshold_text_Callback(hObject, eventdata, handles)
% global timevector timept timept2  ts_axes_obj; 
global rivets beam
rivets.threshold(:,1) = repmat(str2double(get(hObject,'String')),size(beam.voxels,1),1);
switch(rivets.displaymode)
    case 1
        plot_ts;
        paint_activation(hObject,eventdata,handles);
    case 2
        spm_image('init',beam.coreg.mripath); % load/reload structural MRI;
        nut_spmfig_setup
        plot_ts; % erase any previous interval lines
        paint_activation([],[],handles)
        nut_render_meg(beam.coreg.brainrender_path);
end
return;


%%-----------------------------------------------
function nut_ts_button_Callback(hObject, eventdata, handles)
% --- Executes on button press in nut_ts_button.
plot_ts;
return;


%%-----------------------------------------------
function nut_interval_button_Callback(hObject, eventdata, handles)
% --- Executes on button press in nut_interval_button.
global rivets beam;

% initialize times
set(handles.nut_time_text,'String',beam.timewindow(1));
set(handles.nut_time2_text,'String',beam.timewindow(end));
plot_ts; % erase any previous interval lines

[time,crap,button]=ginput(1);

set(handles.nut_time_text,'String',time);
set(handles.nut_time2_text,'String',time);
plot_ts;

button = 0;
% threshold = 0;
% now for the second selection
while(button~=1 && button~=3)
    [time2,ypos,button]=ginput(1);
    timept = dsearchn(beam.timewindow,time);
    if (button==3)   % just use first time point if right button clicked
        timept2 = timept;
    elseif (button==2)  % if middle button clicked, set threshold
        threshold = abs(ypos);
        set(handles.nut_threshold_text,'String',threshold);
        rivets.threshold = [repmat(threshold,size(beam.voxels,1),1) repmat(-threshold,size(beam.voxels,1),1)];
		plot_ts;
    else  % select second time point for other buttons
        set(handles.nut_time2_text,'String',time2);
        plot_ts;
    end
end

switch(rivets.displaymode)
    case 1
        % initialize times
        % plot_ts; % erase any previous interval lines

        paint_activation([],[],handles)
        
    case 2
        spm_image('init',beam.coreg.mripath); % load/reload structural MRI; 
        nut_spmfig_setup
        plot_ts; % erase any previous interval lines
        paint_activation([],[],handles)
        nut_render_meg(beam.coreg.brainrender_path);
end

% paint_activation(hObject,eventdata,handles)

return;

%%-----------------------------------------------
function nut_view_interval_button_Callback(hObject, eventdata, handles)
% --- Executes on button press in nut_interval_button.
global rivets beam;

% initialize times
set(handles.nut_view_time_text,'String',beam.timewindow(1));
set(handles.nut_view_time2_text,'String',beam.timewindow(end));
plot_ts; % erase any previous interval lines

[time,crap,button]=ginput(1);

set(handles.nut_view_time_text,'String',time);
% set(handles.nut_view_time2_text,'String',beam.timewindow(end));

button = 0;
% now for the second selection
while(button~=1 && button~=3)
    [time2,ypos,button]=ginput(1);
    timept = dsearchn(beam.timewindow,time);
    if (button==3)   % just use first time point if right button clicked
        timept2 = timept;
%     elseif (button==2)  % if middle button clicked, set threshold
%         threshold = ypos;
%         set(handles.nut_threshold_text,'String',threshold);
    else  % select second time point for other buttons
        timept2 = dsearchn(beam.timewindow,time2);
    end
end

set(handles.nut_view_time2_text,'String',time2);
plot_ts;

return;


%%---------------------------------------------------
function nut_animate_button_Callback(hObject, eventdata, handles)
% --- Executes on clicking of nut_animate_button button.
% retrieve selected coordinates, and inverse transform to MEG mm
global rivets beam st;
persistent stopflag;

switch(get(hObject,'String'))
    case 'Animate'
        stepsize = 5;
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
        end
        stopflag = false;  % set stopflag to false so we can slide
        set(hObject,'String','Stop Animation');
        
        if(0)  % moves time series to SPM window... unfortunately, SPM's blob colorbar craps all over it
            set(handles.nut_ts_axes,'Units','Normalized');
            set(handles.nut_ts_axes,'Position',[0.55 0.45 0.4 0.25],'Parent',st.fig); % move time series to SPM window
        elseif(rivets.displaymode==2) % for brain rendering...
            figure(st.fig);
            texthandle = text(1,1,'0.0 ms','FontSize',14,'HorizontalAlignment','right');
        else % so for now, let's just add a box to indicate time
            imagewidth = get(st.vols{1}.ax{1}.d,'XData');
            texthandle = xlabel(st.vols{1}.ax{3}.ax,'0.0 ms','FontSize',14,'HorizontalAlignment','right','Position',[imagewidth(2) -15 1]);
        end

        framenum = 0;
        % make sure figures export correctly:
        set(st.fig,'PaperPositionMode','auto');
        set(st.fig,'PaperOrientation','portrait');
        
        MEGpos = (spm_orthviews('pos'))';
        % figure out which MEG voxel is closest, then plot the time series for it
        if rivets.displaymode==1
            MEGvoxelindex = dsearchn(nut_coordtfm(rivets.voxelsMRI,st.vols{1}.premul),MEGpos);
        else
            MEGvoxelindex = dsearchn(rivets.voxelsMRI,MEGpos);
        end


        timezero = dsearchn(beam.timewindow,0);

        for i = timezero:stepsize:length(beam.timewindow)
            % someday computers will be fast enough that you'd want to pause
            % between frames... but not as of Apr. 2003!
            % pause(0.1);

            set(handles.nut_time_text,'String',beam.timewindow(i));
            set(handles.nut_time2_text,'String',beam.timewindow(i));
%             set(texthandle,'String',[num2str(beam.timewindow(i),'%0.1f') ' ms']);

            framenum = framenum + 1;
            switch(rivets.displaymode)
                case 1
                    plot_ts;
                    paint_activation([],[],handles)
                case 2
                    spm_image('init',beam.coreg.mripath); % load/reload structural MRI;
                    nut_spmfig_setup
                    plot_ts; % erase any previous interval lines
                    paint_activation([],[],handles)
                    nut_render_meg(beam.coreg.brainrender_path);
            end
            %             paint_activation(hObject,eventdata,handles);
            texthandle = text(-200,-100,[num2str(beam.timewindow(i),'%0.1f') ' ms'],'Color',[1 1 0],'FontSize',24);

            if(kubrick)
                % save as /tmp/frame0.###
                print(st.fig,'-dtiff','-r72','-noui',['movie/MEGframe' num2str(framenum,'%0.3d') '.tiff']);
            end
            
            if(stopflag) % break loop when stop button is pressed
                delete(texthandle);
                return
            end
        end
        delete(texthandle);
        if(kubrick)
            [moviepath,moviename] = fileparts(rivets.beampath);
            moviepath = fullfile(moviepath,[moviename '.mpg']);
            watchon;
            disp('Please wait... piecing movie together...');
            if rivets.displaymode==1
            unix(['convert -crop -0-300 ' 'movie/MEGframe*.tiff ' moviepath]);
            else
            unix(['convert -crop -0+400 ' 'movie/MEGframe*.tiff ' moviepath]);
            end
            disp(['Movie saved as ' moviepath]);
            status = unix('which gmplayer');
            if(status) % if status is nonzero, "gmplayer" failed and probably is not installed
                warning('MPEG not displayed since gmplayer not found. (http://www.mplayerhq.hu)');
            else
                unix(['gmplayer -fps 10 -noaspect ' moviepath '&']);
            end
            delete(fullfile(tempdir,'MEGframe*.tiff'));
            watchoff;
        end
    case 'Stop Animation'
        stopflag = true;
        set(hObject,'String','Animate');
        return
end

set(hObject,'String','Animate');


%%-----------------------------------------------
function nut_activation_viewer_CloseRequestFcn(hObject, eventdata, handles)
% --- Executes during object deletion, before destroying properties.
global rivets st beam;
if(exist('st','var'))
    spm_image('init',beam.coreg.mripath); % reload structural MRI; clears VOI emphasis from nut_view_scan_region.m
    spm_figure('ColorMap','gray'); % reset original colormap
end
clear global beam rivets
delete(hObject);

%%------------------------------------------------
function nut_save_nut_ts_button_Callback(hObject, eventdata, handles)
% --- Executes on button press in nut_save_nut_ts_button.
error('FIXME -- this is so old there''s no way it still works!');
global rivets beam;
% global timevector s_beam voxels display_s_II_flag s_II; % beamx beamy;
cursorpos = spm_orthviews('pos')';
MEGvoxelindex = dsearchn(voxelsMRI,cursorpos);
if(rivets.display_s_II_flag)
    MEGts = squeeze(rivets.s_II(MEGvoxelindex,:));
else    
    MEGts = squeeze(s_beam(MEGvoxelindex,:));
end
uisave({'beam.timewindow', 'MEGts'});
return;


%%------------------------------------------------
function nut_MRIcroROI_button_Callback(hObject, eventdata, handles)
% --- Executes on button press in MRIcroROI_button.
global rivets st;

warning('repositioning to MRIcroROI is totally experimental and probably inaccurate!');

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
% threshold = str2double(get(handles.nut_threshold_text,'String'));
% nut_create_spmimage;
global rivets beam
if(get(hObject,'Value'))
    if(~isfield(beam.coreg,'brainrender_path'))
        beam.coreg.brainrender_path = spm_get(1,'render*.mat','Render file',fullfile(spm('Dir'),'rend'));
    end
    rivets.displaymode=2;
    nut_render_meg(beam.coreg.brainrender_path);
else
    rivets.displaymode=1;
    spm_image('init',beam.coreg.mripath); % load/reload structural MRI;
    nut_spmfig_setup
    plot_ts; % erase any previous interval lines
    paint_activation(hObject,eventdata,handles)
end


%%------------------------------------------
function nut_activation_notch_box_Callback(hObject, eventdata, handles)
% --- Executes on button press in nut_activation_notch_box.
plot_ts;
return;

%%------------------------------------------
function activation_bpf_Callback(hObject, eventdata, handles)
% --- Executes on button press in activation_bpf.
plot_ts;
return;


% --- Executes on selection change in activation_notch.
function activation_notch_Callback(hObject, eventdata, handles)
plot_ts;
return;


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

actcontents = get(handles.nut_activation_menu,'String');
actselect = actcontents{get(handles.nut_activation_menu,'Value')};

cursorpos = spm_orthviews('pos')';
[MEGvoxelindex,distances] = dsearchn(rivets.voxelsMRI,nut_coordtfm(cursorpos,inv(st.vols{1}.premul)));
distances = nut_rownorm(nut_coord_diff(rivets.voxelsMRI,cursorpos));
exclude_radius = str2double(get(handles.nut_exclude_radius_text,'String'));
within_radius = str2double(get(handles.nut_within_radius_text,'String'));
searchindices1 = find(distances >= exclude_radius);
searchindices2 = find(distances <= within_radius);
searchindices = intersect(searchindices1,searchindices2);

switch(get(handles.nut_orientations_menu,'Value'))
    case 1
        s_beam_interval = rivets.s_beam(searchindices,:);
    case 2
        s_beam_interval = rivets.s_II(searchindices,:);
    case 3
        s_beam_interval = beam.s{rivets.tsflag}(searchindices,:);
    case 4
        s_beam_interval = rivets.s_beam(searchindices,:);
    case 5
        s_beam_interval = rivets.s_beam(searchindices,:);
    case 6
        s_beam_interval = rivets.s_beam(searchindices,:);
    case 7
        s_beam_interval = rivets.s_beam(searchindices,:);
    otherwise
        errordlg('Operation not yet supported.');
end

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
% spm_orthviews('Reposition',nut_coordtfm(rivets.voxelsMRI(MEGvoxelindex,:),st.vols{1}.premul))
nut_reposition(nut_coordtfm(rivets.voxelsMRI(MEGvoxelindex,:),st.vols{1}.premul))

timept = maxtime;
timept2 = maxtime;
time = beam.timewindow(timept);
time2 = beam.timewindow(timept2);
set(handles.nut_time_text,'String',time);
set(handles.nut_time2_text,'String',time2);

plot_ts;
paint_activation(hObject,eventdata,handles);


%%------------------------------------------
function nut_maxvoxel_button_Callback(hObject, eventdata, handles)
% --- Executes on button press in nut_maxvoxel_button.
global rivets st beam;
time = str2double(get(handles.nut_time_text,'String'));
time2 = str2double(get(handles.nut_time2_text,'String'));
timept = dsearchn(beam.timewindow,time);
timept2 = dsearchn(beam.timewindow,time2);

switch(get(handles.nut_orientations_menu,'Value'))
    case 1
        s_beam_interval = mean(rivets.s_beam(:,timept:timept2),2);
    case 2
        s_beam_interval = mean(rivets.s_II(:,timept:timept2),2);
    case 3
        s_beam_interval = mean(beam.s{rivets.tsflag}(:,timept:timept2),2);
    case 4
        s_beam_interval = rivets.s_beam(:,1);
    case 5
        s_beam_interval = rivets.s_beam(:,1);
    case 6
        s_beam_interval = rivets.s_beam(:,1);
    case 7
        s_beam_interval = rivets.s_beam(:,1);
    otherwise
        errordlg('Operation not yet supported.');
end



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

MEGvoxelindex=maxvoxel;
nut_reposition(nut_coordtfm(rivets.voxelsMRI(MEGvoxelindex,:),st.vols{1}.premul))

plot_ts;
paint_activation(hObject,eventdata,handles);
return;


%%------------------------------------------
function nut_maxtime_button_Callback(hObject, eventdata, handles)
% --- Executes on button press in nut_maxtime_button.
global rivets beam;
cursorpos = spm_orthviews('pos')';
MEGvoxelindex = dsearchn(rivets.voxelsMRI,cursorpos);

actcontents = get(handles.nut_activation_menu,'String');
actselect = actcontents{get(handles.nut_activation_menu,'Value')};

switch(get(handles.nut_orientations_menu,'Value'))
    case 1
        s_beam_interval = rivets.s_beam(MEGvoxelindex,:);
    case 2
        s_beam_interval = rivets.s_II(MEGvoxelindex,:);
    case 3
        s_beam_interval = beam.s{rivets.tsflag}(MEGvoxelindex,:);
    case 4
        s_beam_interval = rivets.s_beam(MEGvoxelindex,:);
    case 5
        s_beam_interval = rivets.s_beam(MEGvoxelindex,:);
    case 6
        s_beam_interval = rivets.s_beam(MEGvoxelindex,:);
    case 7
        s_beam_interval = rivets.s_beam(MEGvoxelindex,:);
    otherwise
        errordlg('Operation not yet supported.');
end


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
time = beam.timewindow(timept);
time2 = beam.timewindow(timept2);
set(handles.nut_time_text,'String',time);
set(handles.nut_time2_text,'String',time2);

plot_ts;
paint_activation(hObject,eventdata,handles);


%%------------------------------------------
function nut_applyfilter_button_Callback(hObject, eventdata, handles)
% --- Executes on button press in nut_applyfilter_button.
plot_ts;
% repaint blobs and time series with new filtered version
paint_activation(hObject,eventdata,handles);

return;

%%-------------------------------------------
function plot_ts(hObject)
% refresh time series plot (input argument used only if "apply to volume" button is pressed)

global rivets st beam ndefaults;

handles = guidata(rivets.fig);


apply_all_flg = 0;
if exist('hObject','var')
   if strcmp(lower(get(hObject,'Tag')),'nut_applyfilter_button')
      apply_all_flg = 1;
   end
end

rivets.tsflag=1;
% set(handles.nut_ts_axes,'ButtonDownFcn','rivets.tsflag=handles.nut_ts_axes')


%%%% set properties for selected display style
stylecontents = get(handles.nut_style_menu,'String');
styleselect = stylecontents{get(handles.nut_style_menu,'Value')};

switch styleselect
    case 'Normal Style'
        switch(rivets.displaymode)
                case 1
                    set([st.slider{:}],'Visible','on');
            case 2
        end
        bf_ts_color = [0 0 1];
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
        else
            set(st.fig,'ColorMap',rivets.grayhot);
        end
    case 'Poster Style (Color)'
        set([st.slider{:}],'Visible','off');
        bf_ts_color = [1 0 0];
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
             set(st.fig,'ColorMap',rivets.grayjet);
        else
            set(st.fig,'ColorMap',rivets.grayhot);
        end
    case 'Poster Style (Color2)'
        set([st.slider{:}],'Visible','off');
        bf_ts_color = [1 0 0];
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
        bf_ts_color = [0 0 0];
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
        set([st.slider{:}],'Visible','off');
        bf_ts_color = [1 .5 0];
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
        else
            set(st.fig,'ColorMap',rivets.grayhot);
        end
    case {'Paper Style (Color)'}
        set([st.slider{:}],'Visible','off');
        bf_ts_color = [1 0 0];
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
        else
            set(st.fig,'ColorMap',rivets.grayhot);
        end
    case 'Paper Style (Monochrome)'
        set([st.slider{:}],'Visible','off');
        bf_ts_color = [0 0 0];
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
        else
            set(st.fig,'ColorMap',rivets.graygray);
        end
    case 'SAM Style'
        set([st.slider{:}],'Visible','off');
        bf_ts_color = [0 0 1];
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
switch(rivets.displaymode)
    case 1
        for i=1:3
            set([st.vols{1}.ax{i}.lx st.vols{1}.ax{i}.ly],'Color',spm_xhair_color,'LineWidth',spm_xhair_width);
        end
        MEGvoxelindex = dsearchn(rivets.voxelsMRI,nut_coordtfm(cursorpos,inv(st.vols{1}.premul)));
    case 2
        MEGvoxelindex = dsearchn(rivets.voxelsMRI,cursorpos);
end

if(ndefaults.vwr.plotweights)
    nut_plotweights(squeeze(beam.W(:,2,MEGvoxelindex)));
end


% mydata = rivets.s_beam(MEGvoxelindex,:)';  % jfh
% mydatval = mydata(1);  % jfh
% set(findobj('Tag','nut_crosshair_val'),'String',sprintf('%0.3f',mydatval));  % jfh
% set_spm_readouts;  % jfh

time = str2double(get(handles.nut_time_text,'String'));
time2 = str2double(get(handles.nut_time2_text,'String'));

if time < beam.timewindow(1), time = beam.timewindow(1); end
if time > beam.timewindow(end), time = beam.timewindow(end); end
if time2 < beam.timewindow(1), time2 = beam.timewindow(1); end
if time2 > beam.timewindow(end), time2 = beam.timewindow(end); end

if(time2 < time)
    timetmp = time2;
    time2 = time;
    time = timetmp;
end

set(handles.nut_time_text,'String',time);
set(handles.nut_time2_text,'String',time2);

timept = dsearchn(beam.timewindow,time);
timept2 = dsearchn(beam.timewindow,time2);

if rivets.display_s_II_flag
    beam_val = mean(rivets.s_II(MEGvoxelindex,timept:timept2));
else
    beam_val = mean(rivets.s_beam(MEGvoxelindex,timept:timept2));
end

switch(rivets.displaymode)
    case 1
set(st.beamin,'String',sprintf('%g',rivets.scalefactor*beam_val));
    case 2
end



view_time = str2double(get(handles.nut_view_time_text,'String'));
view_time2 = str2double(get(handles.nut_view_time2_text,'String'));

if view_time < beam.timewindow(1), view_time = beam.timewindow(1); end
if view_time > beam.timewindow(end), view_time = beam.timewindow(end); end
if view_time2 < beam.timewindow(1), view_time2 = beam.timewindow(1); end
if view_time2 > beam.timewindow(end), view_time2 = beam.timewindow(end); end

if(view_time2 < view_time)
    view_timetmp = view_time2;
    view_time2 = view_time;
    view_time = view_timetmp;
end

set(handles.nut_view_time_text,'String',view_time);
set(handles.nut_view_time2_text,'String',view_time2);

view_timept = dsearchn(beam.timewindow,view_time);
view_timept2 = dsearchn(beam.timewindow,view_time2);

fs = beam.srate;
% these values pertain to the sbeam data, not the meg data (which uses beam.params.preprocessing.*)
notch = get(handles.activation_notch,'Value');
bpf = get(handles.activation_bpf,'Value');
tmp = round(str2num(get(handles.nut_lpf_text,'String'))*10)/10; if tmp < 0, tmp = 0; end, if tmp >= fs/2, tmp = fix(fs/2)-0.1; end, set(handles.nut_lpf_text,'String',num2str(tmp));
bpf_low_cutoff = tmp;
tmp = round(str2num(get(handles.nut_hpf_text,'String'))*10)/10; if tmp < 0, tmp = 0; end, if tmp >= fs/2, tmp = fix(fs/2)-0.1; end, set(handles.nut_hpf_text,'String',num2str(tmp));
bpf_high_cutoff = tmp;

if apply_all_flg
   ndx = 1:size(rivets.s_beam,1);
else
   ndx = MEGvoxelindex;
end

% start with unmodified good channels
if rivets.display_s_II_flag
   data = rivets.s_II(ndx,:)';
   if(rivets.display_s_perp_flag)
       data2 = rivets.s_perp(ndx,:)';
       data2 = nut_filter(data2,beam.timewindow,fs,1,notch,bpf,bpf_low_cutoff,bpf_high_cutoff);
   end
else
   data = rivets.s_beam(ndx,:)';
end
if isfield(rivets,'meg')
    if isfield(rivets.meg,'yc_full')
       megdata=rivets.meg.yc_full';
    elseif isfield(rivets.meg,'data')
       megdata=nut_filter(rivets.meg.data,beam.timewindow,fs,beam.params.preprocessing.baseline,beam.params.preprocessing.notch,beam.params.preprocessing.bpf,beam.params.preprocessing.bpf_low_cutoff,beam.params.preprocessing.bpf_high_cutoff);
%        megdata=rivets.meg.data;
    end
    axes(handles.nut_megts_axes);
    
    global nuts % NUTEEG mod
    if isfield(nuts.meg,'eegflag') & nuts.meg.eegflag
        meg_ts_handle = plot(beam.timewindow(view_timept:view_timept2),1e6*megdata(view_timept:view_timept2,:),'LineWidth',meg_ts_width);
        ylabel('Electric Potential (uV)');
    else
        meg_ts_handle = plot(beam.timewindow(view_timept:view_timept2),1e15*megdata(view_timept:view_timept2,:),'LineWidth',meg_ts_width);
        ylabel('Magnetic Field (fT)');
    end
    if(exist('meg_ts_color','var'))
        set(meg_ts_handle,'Color',meg_ts_color);
    end
    axis tight; xlabel('Time (ms)');
end

data = nut_filter(data,beam.timewindow,fs,1,notch,bpf,bpf_low_cutoff,bpf_high_cutoff);

% apply filters to whole volume flag
if apply_all_flg
   ndx = MEGvoxelindex;
   rivets.s_beam = data';
   
   %%%% need to rescale s_beam
   warning('I can''t imagine this is still correct.')
   %%%% TODO: recode this to match new normalization scheme
   min_s = min(rivets.s_beam(:));
   max_s = max(rivets.s_beam(:));
   if rivets.normflag
       rivets.minblob = 0; % for compatibility with nut_view_beamforming_activation
       rivets.maxblob = 1000; % ditto
   else
       rivets.minblob=min_s;
       rivets.maxblob=max_s;
   end
else
   ndx = 1;
end

axes(handles.nut_ts_axes);
%plot(beam.timewindow,data(:,ndx));
if(exist('data2','var'))
    plot(beam.timewindow(view_timept:view_timept2),rivets.scalefactor*data2(view_timept:view_timept2,ndx),'Color',bf_ts2_color,'LineWidth',bf_ts_width);
    hold on;
    plot(beam.timewindow(view_timept:view_timept2),rivets.scalefactor*data(view_timept:view_timept2,ndx),'Color',bf_ts_color,'LineWidth',bf_ts_width);
else
    plot(beam.timewindow(view_timept:view_timept2),rivets.scalefactor*data(view_timept:view_timept2,ndx),'Color',bf_ts_color,'LineWidth',bf_ts_width);
end
% axis([beam.timewindow(view_timept) beam.timewindow(view_timept2) 0 1000]);
axis tight; 
xlabel('Time (ms)');

if rivets.normflag
    if rivets.display_s_II_flag
        ylabel('Normalized Intensity (-1000 to 1000)');
        temp = axis(handles.nut_ts_axes);
        axis(handles.nut_ts_axes,[temp(1) temp(2) -1000 1000]); % rescale axes to -1000:1000
    else
        ylabel('Normalized Intensity (0-1000)');
    end
else
    ylabel('Intensity');
    temp = axis(handles.nut_ts_axes);
    axis(handles.nut_ts_axes,[temp(1) temp(2) rivets.minblob rivets.maxblob]); % rescale axes to min:max
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
% if(get(handles.nut_statsenabled_box,'Value'))  % EXPERIMENTAL: statistical thresholding
%     threshold = rivets.thresh(MEGvoxelindex);
% else
%     threshold = str2double(get(handles.nut_threshold_text,'String'));
% end

threshold = rivets.threshold(MEGvoxelindex,1);
thresholdneg = rivets.threshold(MEGvoxelindex,2);
hold on
time1_handle = plot([time time],[line_min line_max],'Color',bf_ts_xhair_color,'LineWidth',bf_ts_xhair_width);
time2_handle = plot([time2 time2],[line_min line_max],'Color',bf_ts_xhair_color,'LineWidth',bf_ts_xhair_width);
thresh_handle(1) = plot([beam.timewindow(1) beam.timewindow(end)],[threshold threshold],'Color',bf_ts_xhair_color,'LineWidth',bf_ts_xhair_width);
thresh_handle(2) = plot([beam.timewindow(1) beam.timewindow(end)],[thresholdneg thresholdneg],'Color',bf_ts_xhair_color,'LineWidth',bf_ts_xhair_width);
hold off

switch styleselect
    case 'Normal Style'
        set(handles.nut_ts_axes,'FontUnits','points','LineWidth',2,'FontName','Helvetica','FontWeight','bold','FontSize',12,'Color',[1 1 1],'XColor',[0 0 0],'YColor',[0 0 0]);  %changes time series axes numbers font
        set(handles.nut_megts_axes,'FontUnits','points','LineWidth',2,'FontName','Helvetica','FontWeight','bold','FontSize',12,'Color',[1 1 1],'XColor',[0 0 0],'YColor',[0 0 0]);  %changes time series axes numbers font
        set(handles.nut_activation_viewer,'Color',[.93 .93 .92]);  %color of background of whole figure
        set(cell2mat(get(handles.nut_ts_axes,{'Xlabel','Ylabel'})),'FontUnits','points','FontName','Helvetica','FontWeight','bold','FontSize',12,'Color',[0 0 0]);  %changes time series axes label
        set(cell2mat(get(handles.nut_megts_axes,{'Xlabel','Ylabel'})),'FontUnits','points','FontName','Helvetica','FontWeight','bold','FontSize',12,'Color',[0 0 0]);  %changes time series axes label
        % set(handles.nut_blobstyle,'Value',1);
        
        

    case 'Poster Style (Color)'
        set(handles.nut_ts_axes,'FontUnits','points','LineWidth',2,'FontName','Helvetica','FontWeight','bold','FontSize',18,'Color',[1 1 1],'XColor',[0 0 0],'YColor',[0 0 0]);%changes time series axes numbers font
        set(handles.nut_megts_axes,'FontUnits','points','LineWidth',2,'FontName','Helvetica','FontWeight','bold','FontSize',18,'Color',[1 1 1],'XColor',[0 0 0],'YColor',[0 0 0]);%changes time series axes numbers font
        set(handles.nut_activation_viewer,'Color',[1 1 1]);
        set(cell2mat(get(handles.nut_ts_axes,{'Xlabel','Ylabel'})),'FontUnits','points','FontName','Helvetica','FontWeight','bold','FontSize',18);%changes time series axes label
        set(cell2mat(get(handles.nut_megts_axes,{'Xlabel','Ylabel'})),'FontUnits','points','FontName','Helvetica','FontWeight','bold','FontSize',18);%changes time series axes label
        set(handles.nut_blobstyle,'Value',1);
    case 'Poster Style (Color2)'
        set(handles.nut_ts_axes,'FontUnits','points','LineWidth',2,'FontName','Helvetica','FontWeight','bold','FontSize',18,'Color',[1 1 1],'XColor',[0 0 0],'YColor',[0 0 0]);%changes time series axes numbers font
        set(handles.nut_megts_axes,'FontUnits','points','LineWidth',2,'FontName','Helvetica','FontWeight','bold','FontSize',18,'Color',[1 1 1],'XColor',[0 0 0],'YColor',[0 0 0]);%changes time series axes numbers font
        set(handles.nut_activation_viewer,'Color',[1 1 1]);
        set(cell2mat(get(handles.nut_ts_axes,{'Xlabel','Ylabel'})),'FontUnits','points','FontName','Helvetica','FontWeight','bold','FontSize',18);%changes time series axes label
        set(cell2mat(get(handles.nut_megts_axes,{'Xlabel','Ylabel'})),'FontUnits','points','FontName','Helvetica','FontWeight','bold','FontSize',18);%changes time series axes label
        set(handles.nut_blobstyle,'Value',1);
    case 'Poster Style (Monochrome)'
        set(handles.nut_ts_axes,'FontUnits','points','LineWidth',2,'FontName','Helvetica','FontWeight','bold','FontSize',18,'Color',[1 1 1],'XColor',[0 0 0],'YColor',[0 0 0]);%changes time series axes numbers font
        set(handles.nut_megts_axes,'FontUnits','points','LineWidth',2,'FontName','Helvetica','FontWeight','bold','FontSize',18,'Color',[1 1 1],'XColor',[0 0 0],'YColor',[0 0 0]);%changes time series axes numbers font
        set(handles.nut_activation_viewer,'Color',[1 1 1]);
        set(cell2mat(get(handles.nut_ts_axes,{'Xlabel','Ylabel'})),'FontUnits','points','FontName','Helvetica','FontWeight','bold','FontSize',18);%changes time series axes label
        set(cell2mat(get(handles.nut_megts_axes,{'Xlabel','Ylabel'})),'FontUnits','points','FontName','Helvetica','FontWeight','bold','FontSize',18);%changes time series axes label
        set(handles.nut_blobstyle,'Value',1);
    case 'Presentation Style'
        set(handles.nut_ts_axes,'FontUnits','points','LineWidth',2,'FontName','Helvetica','FontWeight','bold','FontSize',15,'Color',[0 0 0],'XColor',[1 1 0],'YColor',[1 1 0]);  %changes time series axes numbers font
        set(handles.nut_megts_axes,'FontUnits','points','LineWidth',2,'FontName','Helvetica','FontWeight','bold','FontSize',15,'Color',[0 0 0],'XColor',[1 1 0],'YColor',[1 1 0]);  %changes time series axes numbers font
        set(handles.nut_activation_viewer,'Color',[0 0 0]);
        set(cell2mat(get(handles.nut_ts_axes,{'Xlabel','Ylabel'})),'FontUnits','points','FontName','Helvetica','FontWeight','bold','FontSize',15,'Color',[1 1 0]);  %changes time series axes label
        set(cell2mat(get(handles.nut_megts_axes,{'Xlabel','Ylabel'})),'FontUnits','points','FontName','Helvetica','FontWeight','bold','FontSize',15,'Color',[1 1 0]);  %changes time series axes label
        set(handles.nut_blobstyle,'Value',2);
    case {'Paper Style (Color)'}
        set(handles.nut_ts_axes,'FontUnits','points','LineWidth',1,'FontName','Times','FontWeight','normal','FontSize',15,'Color',[1 1 1],'XColor',[0 0 0],'YColor',[0 0 0]);  %changes time series axes numbers font
        set(handles.nut_megts_axes,'FontUnits','points','LineWidth',1,'FontName','Times','FontWeight','normal','FontSize',15,'Color',[1 1 1],'XColor',[0 0 0],'YColor',[0 0 0]);  %changes time series axes numbers font
        set(handles.nut_activation_viewer,'Color',[1 1 1]);
        set(cell2mat(get(handles.nut_ts_axes,{'Xlabel','Ylabel'})),'FontUnits','points','FontName','Times','FontWeight','normal','FontSize',15);  %changes time series axes label
        set(cell2mat(get(handles.nut_megts_axes,{'Xlabel','Ylabel'})),'FontUnits','points','FontName','Times','FontWeight','normal','FontSize',15);  %changes time series axes label
        set(handles.nut_blobstyle,'Value',2);
    case 'Paper Style (Monochrome)'
        set(handles.nut_ts_axes,'FontUnits','points','LineWidth',1,'FontName','Times','FontWeight','normal','FontSize',15,'Color',[1 1 1],'XColor',[0 0 0],'YColor',[0 0 0]);  %changes time series axes numbers font
        set(handles.nut_megts_axes,'FontUnits','points','LineWidth',1,'FontName','Times','FontWeight','normal','FontSize',15,'Color',[1 1 1],'XColor',[0 0 0],'YColor',[0 0 0]);  %changes time series axes numbers font
        set(handles.nut_activation_viewer,'Color',[1 1 1]);
        set(cell2mat(get(handles.nut_ts_axes,{'Xlabel','Ylabel'})),'FontUnits','points','FontName','Times','FontWeight','normal','FontSize',15);  %changes time series axes label
        set(cell2mat(get(handles.nut_megts_axes,{'Xlabel','Ylabel'})),'FontUnits','points','FontName','Times','FontWeight','normal','FontSize',15);  %changes time series axes label
        set(handles.nut_blobstyle,'Value',1);
    case 'SAM Style'   % jfh
        set(handles.nut_ts_axes,'FontUnits','points','LineWidth',2,'FontName','Helvetica','FontWeight','bold','FontSize',12,'Color',[1 1 1],'XColor',[0 0 0],'YColor',[0 0 0]);  %changes time series axes numbers font
        set(handles.nut_megts_axes,'FontUnits','points','LineWidth',2,'FontName','Helvetica','FontWeight','bold','FontSize',12,'Color',[1 1 1],'XColor',[0 0 0],'YColor',[0 0 0]);  %changes time series axes numbers font
        set(handles.nut_activation_viewer,'Color',[.93 .93 .92]);  %color of background of whole figure
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
end


return;


%%--------------------------------------
function paint_activation(hObject,eventdata,handles)
% get information needed to paint blobs and punt to nut_view_beamforming_activations

global rivets st beam;
time = str2double(get(handles.nut_time_text,'String'));
time2 = str2double(get(handles.nut_time2_text,'String'));
timept = dsearchn(beam.timewindow,time);   % find nearest index in timevector to given time
timept2 = dsearchn(beam.timewindow,time2);
% threshold = str2double(get(handles.nut_threshold_text,'String'));
% nut_view_beamforming_activation(timept:timept2,get(findobj('Tag','nut_blobstyle'),'Value'), threshold/1000*beam.maxblob);
nut_view_beamforming_activation(timept:timept2,get(handles.nut_blobstyle,'Value'), rivets.threshold);
% nut_image('shopos');
%set_spm_readouts; % jfh

%%--------------------------------------
function nut_save_ts_roi_button_Callback(hObject, eventdata, handles)
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


%%------------------------------------
function nut_orientations_menu_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
% if ispc
%     set(hObject,'BackgroundColor','white');
% else
%     set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
% end


%%-----------------------------------------------
function nut_orientations_menu_Callback(hObject, eventdata, handles)
% --- Executes on selection change in nut_orientations_menu.
global rivets st beam;

switch get(handles.nut_orientations_menu,'Value')
    case 1
        if(isfield(beam,'s_z'))
            rivets.s_beam = double(beam.s_th.^2 + beam.s_ph.^2 + beam.s_z.^2);
        elseif(isfield(beam,'s_th'))
            rivets.s_beam = double(beam.s_th.^2 + beam.s_ph.^2);
        end
        rivets.minblob = min(rivets.s_beam(:));
        rivets.maxblob = max(rivets.s_beam(:));
        rivets.scalefactor = 1000/max(abs(rivets.maxblob),abs(rivets.minblob));
        rivets.display_s_II_flag = false;
    case 2
        rivets.display_s_II_flag = true;
        if(~isfield(beam,'s_II'))  % only do this the first time cuz it's slower than the line at the DMV
            disp('Please wait... computing s_II and s_perp...');
            if(all(beam.s_ph(:)==0))
                rivets.s_II = beam.s_th;
                rivets.s_perp = beam.s_ph;
            else
                %%%%%%%% TODO: We can improve on this by computing orientation
                %%%%%%%% using moving average of rho (i.e., using 20ms sliding
                %%%%%%%% windows)
                switch(1)
                    case 1  % compute rho over Rzz interval -- original method
                        t_int = 1:size(beam.s_th,2);  % use all available time points
                        %    t_int=210:270
                        rhoselectneg = find(mean(beam.s_th(:,t_int),2)./mean(beam.s_ph(:,t_int),2)<0);
                        % rhoselectneg = find(mean(beam.s_th(:,t_int)./beam.s_ph(:,t_int),2)<0);
                        rho = atan(sqrt(mean(beam.s_th(:,t_int).^2,2)./mean(beam.s_ph(:,t_int).^2,2)));
                        rho(rhoselectneg) = -rho(rhoselectneg);

                        cosrho = repmat(cos(rho),1,size(beam.s_ph,2));
                        sinrho = repmat(sin(rho),1,size(beam.s_ph,2));
                        rivets.s_II = beam.s_ph.*cosrho + beam.s_th.*sinrho;
                        rivets.s_perp = beam.s_th.*cosrho - beam.s_ph.*sinrho;

                        clear cosrho sinrho s_ph s_th;  % done with these large matrices

                    case 2  % compute rho for each time point -- this is crap... need to do something to make polarities consistent across time
                        % this is not crap if you actually expect a rotating dipole!  -jz
                        rivets.s_II = zeros(size(beam.s_th));
                        rivets.s_perp = zeros(size(beam.s_th));
                        rivets.rho = (atan(beam.s_th ./ beam.s_ph));

                        rivets.s_II = beam.s_ph.*cos(rivets.rho) + beam.s_th.*sin(rivets.rho);
                        rivets.s_perp = beam.s_th.*cos(rivets.rho) - beam.s_ph.*sin(rivets.rho);
                        clear s_ph s_th;  % done with these large matrices
                    case 3  % sarang's fancy method (that doesn't make sense, but could be useful for other purposes)
                        rivets.s_II = zeros(size(beam.s_th));
                        rivets.s_perp = zeros(size(beam.s_th));
                        tanrho = beam.s_th ./ beam.s_ph;

                        % polarities much noisier than magnitudes... so let's
                        % decide polarity by mob justice -- 20 neighbors on either side...
                        signrho = sign(filtfilt(ones(20,1),1,sign(tanrho')));

                        rho2 = atan(signrho .* abs(tanrho'));
                        % moving window average... more mob justice
                        windowsize=20;
                        rho = (filtfilt(ones(1,windowsize),1,rho2)/(windowsize*windowsize))';


                        rivets.s_II = beam.s_ph.*cos(rho) + beam.s_th.*sin(rho);
                        rivets.s_perp = beam.s_th.*cos(rho) - beam.s_ph.*sin(rho);
                        % clear s_ph s_th;  % done with these large matrices
                    case 4
                        rivets.s_II = beam.s_th;
                        rivets.s_perp = beam.s_ph;
                end
            end
        end

        %%%% normalize s_II but keep activation_baseline
        min_s = min(rivets.s_II(:));
        max_s = max(rivets.s_II(:));
        if(rivets.normflag)
            rivets.maxblob = max(abs(max_s),abs(min_s));
            rivets.scalefactor = 1000/rivets.maxblob;
            rivets.minblob = -rivets.maxblob;
        else
            rivets.scalefactor = 1;
        end

        % rivets.minblob=min_s;
        % rivets.maxblob=max_s;

        if(get(handles.nut_statsenabled_box,'Value'))  % EXPERIMENTAL: statistical thresholding
            s_beam = rivets.s_II;
            timepre=str2num(get(handles.Tstat_prestim_time,'String'));
            prestim = find(beam.timewindow<=timepre);
            if length(prestim)==0
                disp('you need prestimulus data to compute a T value')
            else

                sigma = sqrt(mean(s_beam(:,prestim).^2,2) - mean(abs(s_beam(:,prestim)),2).^2);
                T = (abs(s_beam(:,prestim)) - repmat(mean(abs(s_beam(:,prestim)),2),[1 size(prestim,1)])) ./ repmat(sigma,[1 size(prestim,1)]);
                Tmax = sort(max(T,[],2));
                p = floor(0.99*size(Tmax,1));
                Tthresh = Tmax(p);
                rivets.threshold(:,1) = rivets.scalefactor*Tthresh*sigma + mean(abs(s_beam(:,prestim)),2);
                rivets.threshold(:,2) = -rivets.threshold(:,1);

                %             threshold = round(threshold/max(abs(rivets.maxblob),abs(rivets.minblob))*1000);
            end
        end
    case 3
        disp('sorry this may not be working yet');
    case 4  % psuedo-T
        switch get(handles.nut_noisetype,'Value')
            case 1
                [u,s]=svd(beam.params.TotC);
                noise=s(end,end);
            case 2
                noise=beam.params.Czz;
        end
        rivets.display_s_II_flag = false;
        if size(beam.W{3},3)==1
%             pt = (diag(beam.W{3}'*beam.params.Rzz1*beam.W{3})-diag(beam.W{3}'*beam.params.Czz*beam.W{3}))./(2*diag(beam.W{3}'*beam.params.Czz*beam.W{3}));
%             pt = (diag(beam.W{3}'*beam.params.Rzz1*beam.W{3})-diag(beam.W{3}'*beam.params.Czz*beam.W{3}))./(2*diag(beam.W{3}'*sig*beam.W{3}));
            pt = (diag(beam.W{3}'*beam.params.Rzz1*beam.W{3})-diag(beam.W{3}'*beam.params.Czz*beam.W{3}))./(2*diag(beam.W{3}'*noise*beam.W{3}));
        else
            for ii=1:size(beam.W{3},2)
                apower(:,ii)=diag(squeeze(beam.W{3}(:,ii,:))'*beam.params.Rzz1*squeeze(beam.W{3}(:,ii,:)));
                cpower(:,ii)=diag(squeeze(beam.W{3}(:,ii,:))'*beam.params.Czz*squeeze(beam.W{3}(:,ii,:)));
%                 anoise(:,ii)=diag(squeeze(beam.W{3}(:,ii,:))'*beam.params.Rzz1*squeeze(beam.W{3}(:,ii,:))); %wrong!
%                 cnoise(:,ii)=diag(squeeze(beam.W{3}(:,ii,:))'*beam.params.Czz*squeeze(beam.W{3}(:,ii,:)));
                npower(:,ii)=diag(squeeze(beam.W{3}(:,ii,:))'*noise*squeeze(beam.W{3}(:,ii,:)));
            end
%             pt = ((sum(apower')-sum(cpower'))./(sum(anoise')+sum(cnoise')))';
            pt = ((sum(apower')-sum(cpower'))./(2*sum(npower')))';
        end
            
        rivets.s_beam=repmat(pt,1,length(beam.timewindow));
        min_s = min(rivets.s_beam(:));
        max_s = max(rivets.s_beam(:));
        rivets.maxblob = max(abs(max_s),abs(min_s));
        rivets.minblob = -rivets.maxblob;
        if(rivets.normflag)
            rivets.scalefactor = 1000/rivets.maxblob;
        else
            rivets.scalefactor = 1;
        end
    case 5 % powerdiff (might be psuedo-T if lead field or weight is normalized)
        rivets.display_s_II_flag = false;
        if size(beam.W{3},3)==1
        pt = diag(beam.W{3}'*beam.params.Rzz1*beam.W{3})-diag(beam.W{3}'*beam.params.Czz*beam.W{3});
        else
            for ii=1:size(beam.W{3},2)
                apower(:,ii)=diag(squeeze(beam.W{3}(:,ii,:))'*beam.params.Rzz1*squeeze(beam.W{3}(:,ii,:)));
                cpower(:,ii)=diag(squeeze(beam.W{3}(:,ii,:))'*beam.params.Czz*squeeze(beam.W{3}(:,ii,:)));
            end
            pt = (sum(apower')-sum(cpower'))';
        end
        rivets.s_beam=repmat(pt,1,length(beam.timewindow));
        min_s = min(rivets.s_beam(:));
        max_s = max(rivets.s_beam(:));
        rivets.maxblob = max(abs(max_s),abs(min_s));
        rivets.minblob = -rivets.maxblob;
        if(rivets.normflag)
            rivets.scalefactor = 1000/rivets.maxblob;
        else
            rivets.scalefactor = 1;
        end
    case 6  % psuedo-T, with separate active and control weights
        % fixme: is this right?  why not active cov matrix in denominator?
        % fixme2: is sum of power ok?  
        rivets.display_s_II_flag = false;
        switch get(handles.nut_noisetype,'Value')
            case 1
                [u,s]=svd(beam.params.TotC);
                noise=s(end,end);
            case 2
                noise=beam.params.Czz;
        end
        if size(beam.W{3},3)==1  
%             pt = (diag(beam.W{1}'*beam.params.Rzz1*beam.W{1})-diag(beam.W{2}'*beam.params.Czz*beam.W{2}))./(diag(beam.W{1}'*beam.params.Czz*beam.W{1})+diag(beam.W{2}'*beam.params.Czz*beam.W{2}));
            pt = (diag(beam.W{1}'*beam.params.Rzz1*beam.W{1})-diag(beam.W{2}'*beam.params.Czz*beam.W{2}))./(diag(beam.W{1}'*noise*beam.W{1})+diag(beam.W{2}'*noise*beam.W{2}));
        else
            for ii=1:size(beam.W{3},2)
                apower(:,ii)=diag(squeeze(beam.W{1}(:,ii,:))'*beam.params.Rzz1*squeeze(beam.W{1}(:,ii,:)));
                cpower(:,ii)=diag(squeeze(beam.W{2}(:,ii,:))'*beam.params.Czz*squeeze(beam.W{2}(:,ii,:)));
%                 anoise(:,ii)=diag(squeeze(beam.W{1}(:,ii,:))'*beam.params.Rzz1*squeeze(beam.W{1}(:,ii,:)));
%                 cnoise(:,ii)=diag(squeeze(beam.W{2}(:,ii,:))'*beam.params.Czz*squeeze(beam.W{2}(:,ii,:)));
                anpower(:,ii)=diag(squeeze(beam.W{1}(:,ii,:))'*noise*squeeze(beam.W{1}(:,ii,:)));
                cnpower(:,ii)=diag(squeeze(beam.W{2}(:,ii,:))'*noise*squeeze(beam.W{2}(:,ii,:)));
            end
            pt = ((sum(apower')-sum(cpower'))./(sum(anpower')+sum(cnpower')))';
        end
            
        rivets.s_beam=repmat(pt,1,length(beam.timewindow));
        min_s = min(rivets.s_beam(:));
        max_s = max(rivets.s_beam(:));
        rivets.maxblob = max(abs(max_s),abs(min_s));
        rivets.minblob = -rivets.maxblob;
        if(rivets.normflag)
            rivets.scalefactor = 1000/rivets.maxblob;
        else
            rivets.scalefactor = 1;
        end
    case 7 % powerdiff (might be psuedo-T if lead field or weight is normalized), separate active and control weights
        rivets.display_s_II_flag = false;
        if size(beam.W{3},3)==1
            pt = diag(beam.W{1}'*beam.params.Rzz1*beam.W{1})-diag(beam.W{2}'*beam.params.Czz*beam.W{2});
        else
            for ii=1:size(beam.W{3},2)
                apower(:,ii)=diag(squeeze(beam.W{1}(:,ii,:))'*beam.params.Rzz1*squeeze(beam.W{1}(:,ii,:)));
                cpower(:,ii)=diag(squeeze(beam.W{2}(:,ii,:))'*beam.params.Czz*squeeze(beam.W{2}(:,ii,:)));
            end
            pt = (sum(apower')-sum(cpower'))';
        end
        rivets.s_beam=repmat(pt,1,length(beam.timewindow));
        min_s = min(rivets.s_beam(:));
        max_s = max(rivets.s_beam(:));
        rivets.maxblob = max(abs(max_s),abs(min_s));
        rivets.minblob = -rivets.maxblob;
        if(rivets.normflag)
            rivets.scalefactor = 1000/rivets.maxblob;
        else
            rivets.scalefactor = 1;
        end
end

plot_ts;
paint_activation(hObject,eventdata,handles);


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

% for i=1:3
%     SPM_axes_obj(i)=st.vols{1}.ax{i}.ax;
% end
% 
% % setup time series refresh to occur when clicking on new MRI
% set(SPM_axes_obj,'ButtonDownFcn',rivets.SPM_axes_ButtonDownFcn);

rivets.voxelsMRI = nut_meg2mri(beam.voxels);

rivets.blobtransform = [ beam.voxelsize(1)                0                  0 beam.voxelsize(1)*floor(st.bb(1,1)/beam.voxelsize(1))
                                       0 beam.voxelsize(2)                 0 beam.voxelsize(2)*floor(st.bb(1,2)/beam.voxelsize(2))
                                       0                 0 beam.voxelsize(3) beam.voxelsize(3)*floor(st.bb(1,3)/beam.voxelsize(3))
                                       0                 0                 0                                                     1 ];

rivets.voxelsblob = nut_coordtfm(rivets.voxelsMRI,inv(rivets.blobtransform));
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
global rivets beam
set(handles.nut_threshold_text,'String','0'); %we want to renorm all voxels.
rivets.threshold = zeros(size(rivets.voxelsMRI,1),2);
plot_ts;
paint_activation(hObject,eventdata,handles);
if isfield(beam,'coreg') && isfield(beam.coreg,'norm_mripath')
    norm_sbeamname=nut_normalize_beam(rivets.beampath,beam.coreg.norm_mripath);
else
    norm_sbeamname=nut_normalize_beam(rivets.beampath);
end
clear global rivets
nut_activation_viewer_OpeningFcn(handles.nut_activation_viewer, eventdata, handles, norm_sbeamname)

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
plot_ts;
paint_activation(hObject,eventdata,handles);
return;


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
global rivets beam;

plot_ts;
paint_activation(hObject,eventdata,handles);
return


% --- Executes on button press in nut_norm_box.
function nut_norm_box_Callback(hObject, eventdata, handles)
% hObject    handle to nut_norm_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nut_norm_box
global rivets beam

threshold = rivets.threshold;

if(get(hObject,'Value'))
    rivets.normflag = true;
    rivets.scalefactor = 1000/max(abs(rivets.maxblob),abs(rivets.minblob));
    threshold = round(threshold*rivets.scalefactor);
    %max(abs(rivets.maxblob),abs(rivets.minblob))*1000);
else
    rivets.normflag = false;
    threshold = threshold/rivets.scalefactor;
    %1000*max(abs(rivets.maxblob),abs(rivets.minblob));
    rivets.scalefactor = 1;
end
set(handles.nut_threshold_text,'String',num2str(threshold(1)));
rivets.threshold = threshold;
plot_ts;
paint_activation(hObject,eventdata,handles);


% --- Executes on button press in nut_display_rawmeg_button.
function nut_display_rawmeg_button_Callback(hObject, eventdata, handles)
% hObject    handle to nut_display_rawmeg_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global rivets nuts
nut_importmeg;
set(handles.nut_display_rawmeg_button,'Visible','Off');
rivets.meg.data=nut_filter(nuts.meg.data,nuts.meg.latency,nuts.meg.srate,nuts.preprocessing.baseline,nuts.preprocessing.notch,nuts.preprocessing.bpf,nuts.preprocessing.bpf_low_cutoff,nuts.preprocessing.bpf_high_cutoff);
plot_ts;


% --------------------------------------------------------------------
function MRIcroROI_menu_Callback(hObject, eventdata, handles)
% hObject    handle to MRIcroROI_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global rivets st;

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


% --- Executes on button press in nut_statsenabled_box.
function nut_statsenabled_box_Callback(hObject, eventdata, handles)
% hObject    handle to nut_statsenabled_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nut_statsenabled_box

global rivets beam

if(get(hObject,'Value'))
    set(handles.nut_threshold_text,'Enable','off');
    set(handles.nut_orientations_menu,'Value',2)
    nut_orientations_menu_Callback(handles.nut_orientations_menu,eventdata,handles);
%     if(~isfield(rivets,'threshold'))
%         s_beam = rivets.s_beam;
%         prestim = find(beam.timewindow<=0);
%         %     for r = 1:size(s_beam,1)
%         %         sigma(r) = sqrt(mean(s_beam(r,prestim).^2) - mean(s_beam(r,prestim)).^2);
%         %         T(r,:) = (abs(s_beam(r,prestim)) - mean(abs(s_beam(r,prestim)))) / sigma(r);
%         %     end
% 
%         sigma = sqrt(mean(s_beam(:,prestim).^2,2) - mean(s_beam(:,prestim),2).^2);
%         T = (abs(s_beam(:,prestim)) - repmat(mean(abs(s_beam(:,prestim)),2),[1 size(prestim,1)])) ./ repmat(sigma,[1 size(prestim,1)]);
%         Tmax = sort(max(T,[],2));
%         p = floor(0.99*size(Tmax,1));
%         Tthresh = Tmax(p);
%         rivets.threshold = ( Tthresh*sigma + mean(abs(s_beam(:,prestim)),2) ) ;% .^2;
%     end
else
    set(handles.nut_threshold_text,'Enable','on');
end    


% --------------------------------------------------------------------
function nut_savebeam_menu_Callback(hObject, eventdata, handles)
% hObject    handle to nut_savebeam_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global beam rivets
[filename, pathname] = uiputfile('*.mat', 'Save Beamformer Volume...',rivets.beamfilename);  %optional as long as saving VOIvoxels to nuts
if isequal(filename,0) | isequal(pathname,0)
    return;
else
    filename = fullfile(pathname,filename);
end

save(filename,'beam');



function nut_exclude_radius_text_Callback(hObject, eventdata, handles)
% hObject    handle to nut_exclude_radius_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_exclude_radius_text as text
%        str2double(get(hObject,'String')) returns contents of nut_exclude_radius_text as a double


% --- Executes on button press in checkbox14.
function checkbox14_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox14


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


% --------------------------------------------------------------------
function nut_sbeam2ana_menu_Callback(hObject, eventdata, handles)
% hObject    handle to nut_sbeam2ana_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global beam rivets
%global voxelsize voxels s_beam timept timept2 pathname filename megfile
[megimageoutname,megimageoutpath]=uiputfile('*.img','create a name for your img of activations');  %%% FIXME! megfile is a very misleading variable name...
beam.megimageoutname=[megimageoutpath megimageoutname];
megvoxsize=beam.voxelsize;
% structpath=pathname;
% structname=filename;

% Vs=spm_vol([structpath,structname]);
Vs=spm_vol(beam.coreg.mripath);
timept=str2num(get(handles.nut_time_text,'String'));
timept2=str2num(get(handles.nut_time2_text,'String'));
timeofint=dsearchn(beam.timewindow,timept):dsearchn(beam.timewindow,timept2);
spixdim=Vs.private.hdr.dime.pixdim;

Vm.fname=beam.megimageoutname;
xdim=floor(Vs.dim(1)*abs(spixdim(2))/megvoxsize(1));
ydim=floor(Vs.dim(2)*abs(spixdim(3))/megvoxsize(2));
zdim=floor(Vs.dim(3)*abs(spixdim(4))/megvoxsize(3));
Vm.dim=[xdim ydim zdim Vs.dim(4)];

Vs_origmat=[spixdim(2) 0 0 -(Vs.dim(1)+1)/2*spixdim(2); 0 spixdim(3) 0 -(Vs.dim(2)+1)/2*spixdim(3); 0 0 spixdim(4) -(Vs.dim(3)+1)/2*spixdim(4); 0 0 0 1];
Vs_rot=Vs.mat*inv(Vs_origmat);

Vm.mat=[megvoxsize(1) 0 0 -(xdim+1)/2*megvoxsize(1);0 megvoxsize(2) 0 -(ydim+1)/2*megvoxsize(2);0 0 megvoxsize(3) -(zdim+1)/2*megvoxsize(3);0 0 0 1];
Vm.mat=Vs_rot*Vm.mat;
Vm.pinfo=Vs.pinfo;
Vm.descrip=[];

%create header file *.hdr
Vout=spm_create_vol(Vm);
iVsmat=inv(Vs.mat);
iVmmat=inv(Vm.mat);
megcoord=round(nut_coordtfm(round(rivets.voxelsMRI),iVmmat)');
if length(timeofint)>1
    active=mean(rivets.s_beam(:,timeofint),2);
else
    active=rivets.s_beam(:,timeofint);
end

active=active-min(active);
active=active/(max(active))*1000;

tmp   = zeros(xdim,ydim,zdim);
for ii=1:length(active),
    tmp(megcoord(1,ii),megcoord(2,ii),megcoord(3,ii))=active(ii);
end

for jj=1:zdim,
    Vspm=spm_write_plane(Vm,squeeze(tmp(:,:,jj)),jj);
end
Vspm = spm_close_vol(Vspm);



function nut_thresholdneg_text_Callback(hObject, eventdata, handles)
% hObject    handle to nut_thresholdneg_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nut_thresholdneg_text as text
%        str2double(get(hObject,'String')) returns contents of nut_thresholdneg_text as a double
global rivets beam
rivets.threshold(:,2) = repmat(str2double(get(hObject,'String')),size(beam.voxels,1),1);
plot_ts;
paint_activation(hObject,eventdata,handles);
return;


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



function Tstat_prestim_time_Callback(hObject, eventdata, handles)
% hObject    handle to Tstat_prestim_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tstat_prestim_time as text
%        str2double(get(hObject,'String')) returns contents of Tstat_prestim_time as a double


% --- Executes during object creation, after setting all properties.
function Tstat_prestim_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tstat_prestim_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in nut_noisetype.
function nut_noisetype_Callback(hObject, eventdata, handles)
% hObject    handle to nut_noisetype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns nut_noisetype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nut_noisetype

nut_orientations_menu_Callback(hObject,eventdata,handles);


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


