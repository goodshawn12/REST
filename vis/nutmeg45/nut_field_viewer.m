function varargout = nut_field_viewer(varargin)
% nut_field_viewer M-file for nut_field_viewer.fig
%

% Last Modified by GUIDE v2.5 30-May-2009 19:37:03

% Here and below:
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to nut_field_viewer (see VARARGIN)

if(~strcmp('14',version('-release')))
    warning('off')  % hack so warning for "created by matlab 7.0" is turned off
end
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @nut_field_viewer_OpeningFcn, ...
    'gui_OutputFcn',  @nut_field_viewer_OutputFcn, ...
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
function nut_field_viewer_OpeningFcn(hObject, eventdata, handles, varargin)
% --- Executes just before nut_field_viewer is made visible.

% Choose default command line output for nut_field_viewer
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

global rivets st nuts ndefaults;

rivets.fig = hObject;

set(hObject,'Name',['NUTMEG Field Viewer - ' nuts.meg.filename]);

spm_image('init',nuts.coreg.mripath); % load/reload structural MRI; clears VOI emphasis from nut_view_scan_region.m

% reposition SPM window
set(st.fig,'Units','normalized');
spmwindow_position = get(st.fig,'Position');
spmwindow_position(1:2) = [0.48 0]; % place SPM window at top right
set(st.fig,'Position',spmwindow_position,'Resize','on');
% set(st.fig,'Units','points');  % set units back to pixels so fonts don't go nuts

% fit to desired MEG voxel grid, shift such that origin is (0,0,0)mm
% rivets.voxelsMRI = nut_coordtfm(beam.voxels,beam.coreg.meg2mri_tfm);
rivets.voxelsMRI = nut_meg2mri(double(nuts.voxels));
% rivets.voxelsMRI = double(nuts.voxels);

% following jazz is needed because SPM's blob freaks out when given
% negative coordinates, and needs a dilation to know these are big voxels
% furthermore, coords to SPM must be integers.
if(nuts.voxelsize(1,1)<0)
    nuts.voxelsize = abs(nuts.voxelsize);   % to account for negative voxelsizes (radiological???)
    warndlg('negative voxel size??? whassup with that?!');
end

    rivets.blob2mri_tfm = [ nuts.voxelsize(1)                0             0 nuts.voxelsize(1)*floor(st.bb(1,1)/nuts.voxelsize(1))-30
                                       0 nuts.voxelsize(2)                 0 nuts.voxelsize(2)*floor(st.bb(1,2)/nuts.voxelsize(2))-30
                                       0                 0 nuts.voxelsize(3) nuts.voxelsize(3)*floor(st.bb(1,3)/nuts.voxelsize(3))-30
                                       0                 0                 0                                                     1 ];
    rivets.voxelsblob = nut_coordtfm(rivets.voxelsMRI,inv(rivets.blob2mri_tfm));

rivets.blob2mri_tfm=double(rivets.blob2mri_tfm);

clear voxels crap;

rivets.ts_refresh_handle = @plot_ts;  %% pass handle since technically it's a private function
rivets.sensorCallback_handle = @paint_activation;

nut_spmfig_setup;

rivets.displaymode = 1;
rivets.threshold = [0 0];
rivets.currentSensor = 1;

set(handles.nut_sensorlist_menu,'String',nuts.meg.sensor_labels(nuts.meg.goodchannels));

plot_ts;
paint_activation(hObject,eventdata,handles);


% nut_maxactivation_button_Callback(hObject,eventdata,handles);



%%---------------------------------------------------
function varargout = nut_field_viewer_OutputFcn(hObject, eventdata, handles)
if(isfield(handles,'output'))  % prevents crashing if s_beam loading cancelled
    varargout{1} = handles.output;
end


%%------------------------------------------
function nut_maxvoxel_button_Callback(hObject, eventdata, handles)
% --- Executes on button press in nut_maxvoxel_button.
global rivets st nuts

[a,MEGvoxelindex]=max(abs(rivets.Lmax(rivets.currentSensor,:)));
nut_reposition(nut_coordtfm(rivets.voxelsMRI(MEGvoxelindex,:),st.vols{1}.premul))

plot_ts;
paint_activation(hObject,eventdata,handles);
return;


%%------------------------------------------
function nut_maxsensor_button_Callback(hObject, eventdata, handles)
% --- Executes on button press in nut_maxsensor_button.
global rivets nuts
cursorpos = spm_orthviews('pos')';
MEGvoxelindex = dsearchn(rivets.voxelsMRI,cursorpos);

[a,rivets.currentSensor]=max(abs(rivets.Lmax(:,MEGvoxelindex)));

plot_ts;
paint_activation(hObject,eventdata,handles);




%%-------------------------------------------
function plot_ts(hObject)
% refresh time series plot (input argument used only if "apply to volume" button is pressed)

global nuts rivets st beam ndefaults;

handles = guidata(rivets.fig);

% disp(['plot_ts: ' nuts.meg.sensor_labels{rivets.currentSensor}]);

set(handles.nut_sensorlist_menu,'Value',rivets.currentSensor);


lfmode=get(handles.nut_selectLcomp_menu,'String');
lfmode=lfmode{get(handles.nut_selectLcomp_menu,'Value')};

switch(lfmode)
    case 'Projected L'
        [rivets.Lmax,rivets.Lmin] = findLscalar(nuts.Lp(nuts.meg.goodchannels,:,:));
    case 'Lx'
        rivets.Lmax = squeeze(nuts.Lp(nuts.meg.goodchannels,1,:));
    case 'Ly'
        rivets.Lmax = squeeze(nuts.Lp(nuts.meg.goodchannels,2,:));
    case 'Lz'
        rivets.Lmax = squeeze(nuts.Lp(nuts.meg.goodchannels,3,:));
    otherwise
        error('wonkiness.');
end        


apply_all_flg = 0;
if exist('hObject','var')
   if strcmp(lower(get(hObject,'Tag')),'nut_applyfilter_button')
      apply_all_flg = 1;
   end
end

rivets.tsflag=1;
% set(handles.nut_ts_axes,'ButtonDownFcn','rivets.tsflag=handles.nut_ts_axes')


%%%% set properties for selected display style
%    case 'Normal Style'
%         switch(rivets.displaymode)
%                 case 1
%                     set([st.slider{:}],'Visible','on');
%         end
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
        
%         if(1)
%             warning('colormap test!!!');
%         end

        rivets.grayjet = [gray(64);jet(64)]; % gray-jet colormap
        set(st.fig,'ColorMap',rivets.grayjet);


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

global nuts
axes(handles.nut_ts_axes);

% hack to correspond coil positions to final lead field rows
if isfield(nuts.meg,'chanmixMtx')
for ii=1:size(nuts.meg.chanmixMtx{1},1)
    [jnk,tmp]=find(nuts.meg.chanmixMtx{1}(ii,:)==1);
    if(~isempty(tmp))
        coilselect(ii)=tmp(1);
    end
end
coilselect = coilselect(nuts.meg.goodchannels);
else
    coilselect = nuts.meg.goodchannels;
end

rivets.sensormap = nut_plot_bfield(nuts.meg.sensorCoord(coilselect,:),rivets.Lmax(:,MEGvoxelindex)); axis equal; axis off
plot_sensors
set(st.beamin,'String',sprintf('%g',rivets.Lmax(rivets.currentSensor,MEGvoxelindex)));

ratio = abs(rivets.Lmax(rivets.currentSensor,MEGvoxelindex))/max(abs(rivets.Lmax(rivets.currentSensor,:)))


function plot_sensors
global nuts rivets

% disp(['plot_sensors: ' nuts.meg.sensor_labels{rivets.currentSensor}]);


handles = guidata(rivets.fig);

axes(handles.nut_ts_axes);
hold on
for ii=1:size(rivets.sensormap,1)
    rivets.ploth(ii) = plot(-rivets.sensormap(ii,2),-rivets.sensormap(ii,1),'ow','ButtonDownFcn',['global rivets; rivets.currentSensor = ' num2str(ii) '; feval(rivets.ts_refresh_handle); feval(rivets.sensorCallback_handle);']);
%    set(rivets.ploth(ii),'Color','white','LineWidth',0.5,'LineStyle','o','ButtonDownFcn',['global rivets; rivets.currentSensor = ' num2str(ii) '; feval(rivets.sensorCallback_handle);']);
end
set(rivets.ploth(rivets.currentSensor),'Color','red','LineWidth',1.5,'LineStyle','x','ButtonDownFcn','');
hold off

return;


%%--------------------------------------
function paint_activation(hObject,eventdata,handles)
% get information needed to paint blobs and punt to nut_view_beamforming_activations

global nuts rivets st

% disp(['paint_activation: ' nuts.meg.sensor_labels{rivets.currentSensor}]);



%threshold = [0 0];

sensor = rivets.currentSensor;
select = (rivets.Lmax(sensor,:) > rivets.threshold(:,1)) | (rivets.Lmax(sensor,:) < rivets.threshold(:,2));


spm_orthviews('rmblobs',1);
spm_orthviews('addblobs',1,rivets.voxelsblob(select,:)',rivets.Lmax(sensor,select),rivets.blob2mri_tfm);
%st.vols{1}.blobs{1}.max = rivets.maxblob*rivets.scalefactor;
%st.vols{1}.blobs{1}.min = rivets.minblob*rivets.scalefactor;
maxblob = max(abs(st.vols{1}.blobs{1}.max),abs(st.vols{1}.blobs{1}.min));
st.vols{1}.blobs{1}.max = maxblob;
st.vols{1}.blobs{1}.min = -maxblob;


% we shouldn't have to do this -- even after above hack, we need to force
% the colorbar to reflect the full range as well
image([0 1],[st.vols{1}.blobs{1}.min st.vols{1}.blobs{1}.max],[1:64]' + 64,'Parent',st.vols{1}.blobs{1}.cbar);
set(st.vols{1}.blobs{1}.cbar,'YDir','normal','XTickLabel',[]);
spm_orthviews('Redraw')

nut_image('shopos');



function nut_thresholdneg_text_Callback(hObject, eventdata, handles)
global rivets nuts
rivets.threshold(2) = str2double(get(hObject,'String'));
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



function nut_threshold_text_Callback(hObject, eventdata, handles)
global rivets nuts
rivets.threshold(1) = str2double(get(hObject,'String'));
plot_ts;
paint_activation(hObject,eventdata,handles);
return;


% --- Executes during object creation, after setting all properties.
function nut_threshold_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_threshold_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [Lmax,Lmin]=findLscalar(Lp)

Lmax = zeros(size(Lp,1),size(Lp,3));
Lmin = zeros(size(Lp,1),size(Lp,3));
for ii=1:size(Lp,3)
    [u,s,v]=svd(Lp(:,:,ii)'*Lp(:,:,ii));
    Lmax(:,ii)=Lp(:,:,ii)*v(:,1);
    Lmin(:,ii)=Lp(:,:,ii)*v(:,3);
end


% --- Executes on selection change in nut_selectLcomp_menu.
function nut_selectLcomp_menu_Callback(hObject, eventdata, handles)
% hObject    handle to nut_selectLcomp_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns nut_selectLcomp_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nut_selectLcomp_menu


% --- Executes during object creation, after setting all properties.
function nut_selectLcomp_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_selectLcomp_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in nut_sensorlist_menu.
function nut_sensorlist_menu_Callback(hObject, eventdata, handles)
% hObject    handle to nut_sensorlist_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns nut_sensorlist_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nut_sensorlist_menu

global nuts rivets

% disp(['sensorlist_menu_Callback1: ' nuts.meg.sensor_labels{rivets.currentSensor}]);

rivets.currentSensor = get(hObject,'Value');

% disp(['sensorlist_menu_Callback2: ' nuts.meg.sensor_labels{rivets.currentSensor}]);

plot_ts;

% disp(['sensorlist_menu_Callback3: ' nuts.meg.sensor_labels{rivets.currentSensor}]);

paint_activation;

% disp(['sensorlist_menu_Callback4: ' nuts.meg.sensor_labels{rivets.currentSensor}]);

% --- Executes during object creation, after setting all properties.
function nut_sensorlist_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_sensorlist_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


