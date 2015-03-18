function varargout = nut_CoregistrationTool(varargin)
% NUT_COREGISTRATIONTOOL M-file for nut_CoregistrationTool.fig
%      NUT_COREGISTRATIONTOOL, by itself, creates a new NUT_COREGISTRATIONTOOL or raises the existing
%      singleton*.
%
%      H = NUT_COREGISTRATIONTOOL returns the handle to a new NUT_COREGISTRATIONTOOL or the handle to
%      the existing singleton*.
%
%      NUT_COREGISTRATIONTOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NUT_COREGISTRATIONTOOL.M with the given input arguments.
%
%      NUT_COREGISTRATIONTOOL('Property','Value',...) creates a new NUT_COREGISTRATIONTOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before nut_CoregistrationTool_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to nut_CoregistrationTool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help nut_CoregistrationTool

% Last Modified by GUIDE v2.5 16-Jun-2009 20:58:46

warning('off')  % hack so warning for "created by matlab 7.0" is turned off
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @nut_CoregistrationTool_OpeningFcn, ...
                   'gui_OutputFcn',  @nut_CoregistrationTool_OutputFcn, ...
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


% --- Executes just before nut_CoregistrationTool is made visible.
function nut_CoregistrationTool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to nut_CoregistrationTool (see VARARGIN)

% Choose default command line output for nut_CoregistrationTool
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes nut_CoregistrationTool wait for user response (see UIRESUME)
% uiwait(handles.nut_CoregistrationTool);

global coreg;
coreg = varargin{1};

% fill in previously selected information, if available
if isfield(coreg, 'mripath')
    if length(coreg.mripath)>35
        imagestring=['...' coreg.mripath(end-30:end)];
        set(handles.nut_ImageNameText,'String',imagestring);
    else
        set(handles.nut_ImageNameText,'String',coreg.mripath);
    end
end
if isfield(coreg,'norm_mripath')
    if length(coreg.norm_mripath)>40
        imagestring=['...' coreg.norm_mripath(end-35:end)];
        set(handles.nut_normimage_nametext,'String',imagestring);
    else
        set(handles.nut_normimage_nametext,'String',coreg.norm_mripath);
    end
end
if isfield(coreg,'brainrender_path')
    if length(coreg.brainrender_path)>40
        imagestring=['...' coreg.brainrender_path(end-35:end)];
        set(handles.nut_renderedbrain_nametext,'String',imagestring);
    else
        set(handles.nut_renderedbrain_nametext,'String',coreg.brainrender_path);
    end
end

% if(isfield(coreg,'fiducials_mri_mm'))
%     set(handles.nut_left_text,'String',['MRI: ' num2str(coreg.fiducials_mri_mm(1,:),'%.1f ')]);
%     set(handles.nut_right_text,'String',['MRI: ' num2str(coreg.fiducials_mri_mm(2,:),'%.1f ')]);
%     set(handles.nut_nose_text,'String',['MRI: ' num2str(coreg.fiducials_mri_mm(3,:),'%.1f ')]);
% end



global screws
screws.coreg.handles = handles;
nut_coreg_enabler



% --- Outputs from this function are returned to the command line.
function varargout = nut_CoregistrationTool_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;

global coreg screws
uiwait(handles.output);
clear global screws
varargout{1} = coreg;

% --- Executes during object creation, after setting all properties.
function nut_orientation_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nut_orientation_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on selection change in nut_orientation_menu.
function nut_orientation_menu_Callback(hObject, eventdata, handles)
% hObject    handle to nut_orientation_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns nut_orientation_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nut_orientation_menu
global coreg
coreg.orientation=get(hObject,'Value');
nut_coreg_enabler;


% --- Executes on button press in nut_coreg_firstpass.
function nut_coreg_firstpass_Callback(hObject, eventdata, handles)
% hObject    handle to nut_coreg_firstpass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in nut_lsc.
function nut_lsc_Callback(hObject, eventdata, handles)
% hObject    handle to nut_lsc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in nut_showlsc.
function nut_showlsc_Callback(hObject, eventdata, handles)
% hObject    handle to nut_showlsc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in nut_nose.
function nut_nose_Callback(hObject, eventdata, handles)
% hObject    handle to nut_nose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in nut_normMRI.
function nut_normMRI_Callback(hObject, eventdata, handles)
% hObject    handle to nut_normMRI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global coreg screws rivets
if get(screws.coreg.handles.nut_orientation_menu,'Value')==2
    errordlg('you must use a neurological MRI for the whole spatial normalization thing to work here');
end
[mrifile, norm_mripath]=uigetfile('w*.img;w*.nii','Please Select Normed MRI...');
if isequal(mrifile,0)|isequal(norm_mripath,0)
    return;
end

h=waitbar(0.2,'Loading normalized image volume...');

coreg.norm_mripath = fullfile(norm_mripath,mrifile);
if length(coreg.norm_mripath)>40
    imagestring=['...' coreg.norm_mripath(end-35:end)];
    set(findobj('Tag','nut_normimage_nametext'),'String',imagestring);
else
    set(findobj('Tag','nut_normimage_nametext'),'String',coreg.norm_mripath);
end
snmat=[norm_mripath mrifile(2:end-4) '_sn.mat'];
yimg=[norm_mripath 'y_' mrifile(2:end)];
if ~exist([yimg(1:end-4) '.nii'],'file') && ~exist([yimg(1:end-4) '.img'],'file')
    if strncmp(spm('ver'),'SPM8',4)
        cwd=pwd;
        cd(norm_mripath); %ugly hack because SPM8 wants to write files to 'pwd'
        spm_jobman('initcfg');

        % generate deformation field (MNI -> individual MRI)
        matlabbatch{1}.spm.util.defs.comp{1}.sn2def.matname = {snmat};
        matlabbatch{1}.spm.util.defs.comp{1}.sn2def.vox = NaN(1,3);
        matlabbatch{1}.spm.util.defs.comp{1}.sn2def.bb = NaN(2,3);
        matlabbatch{1}.spm.util.defs.ofname = mrifile(2:(end-4));
        matlabbatch{1}.spm.util.defs.fnames = '';
        matlabbatch{1}.spm.util.defs.interp = 1;
        
        % generate inverse deformation field (individual MRI -> MNI)
        matlabbatch{2}.spm.util.defs.comp{1}.inv.comp{1}.sn2def.matname = {snmat};
        matlabbatch{2}.spm.util.defs.comp{1}.inv.comp{1}.sn2def.vox = NaN(1,3);
        matlabbatch{2}.spm.util.defs.comp{1}.inv.comp{1}.sn2def.bb = NaN(2,3);
        matlabbatch{2}.spm.util.defs.comp{1}.inv.space = {coreg.norm_mripath};
        matlabbatch{2}.spm.util.defs.ofname = ['i' mrifile(2:(end-4))];
        matlabbatch{2}.spm.util.defs.fnames = '';
        matlabbatch{2}.spm.util.defs.interp = 1;

        spm_jobman('run',matlabbatch);
        cd(cwd);
    elseif(strcmp(spm('ver'),'SPM2'))
        nut_spm_sn2def('def',snmat);  % you'll need this later for beam coord on-the-fly transform
        iyimg=[norm_mripath 'iy_' mrifile(2:end)];

        if exist(iyimg,'file')~=2
            inputimg=[norm_mripath mrifile(2:end)];
            nut_spm_invdef(yimg,inputimg);
        end
    else
        errordlg(['NUTMEG is not compatible with ' spm('ver') '. Please place SPM2 or SPM8 in your path, which may be downloaded from http://www.fil.ion.ucl.ac.uk/spm'])
        return
    end
    waitbar(0.6,h);
end

if(~isfield(rivets,'TalDB'))
%     load('ihb_DataBase.cdb','-mat'); % load MNI labels
%     rivets.MNIdb = MNIdb;
%     [meshx,meshy,meshz]=ndgrid(MNIdb.minX:MNIdb.voxX:MNIdb.maxX,MNIdb.minY:MNIdb.voxY:MNIdb.maxY,MNIdb.minZ:MNIdb.voxZ:MNIdb.maxZ);
%     rivets.MNIdb.coords = [meshx(:) meshy(:) meshz(:)];
%     
%     rivets.MNIdm = nut_loadMNI;
    rivets.TalDB = load('talairachDB');
end

waitbar(1,h);
close(h);


% % --- Executes when user attempts to close nut_CoregistrationTool.
% function nut_CoregistrationTool_CloseRequestFcn(hObject, eventdata, handles)
% % hObject    handle to nut_CoregistrationTool (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hint: delete(hObject) closes the figure
% clear global coreg screws
% delete(hObject);
% 

% --- Executes on button press in nut_load_lsc.
function nut_load_lsc_Callback(hObject, eventdata, handles)
% hObject    handle to nut_load_lsc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in nut_renderedbrain.
function nut_renderedbrain_Callback(hObject, eventdata, handles)
% hObject    handle to nut_renderedbrain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global coreg screws
if get(screws.coreg.handles.nut_orientation_menu,'Value')==2
    errordlg('you must use a neurological MRI for the whole rendered brain bit to work here');
end
[brainrender_file, brainrender_path]=uigetfile('render_*.mat','Please Select SPM-Rendered Brain...');
if isequal(brainrender_file,0)|isequal(brainrender_path,0)
    return;
end

coreg.brainrender_path = fullfile(brainrender_path,brainrender_file);
if length(coreg.brainrender_path)>40
    imagestring=['...' coreg.brainrender_path(end-35:end)];
    set(findobj('Tag','nut_renderedbrain_nametext'),'String',imagestring);
else
    set(findobj('Tag','nut_renderedbrain_nametext'),'String',coreg.brainrender_mripath);
end
