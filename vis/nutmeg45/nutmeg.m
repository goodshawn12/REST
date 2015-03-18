function varargout = nutmeg(varargin)
% NUTMEG M-file for nutmeg.fig
%   A standard .fig-associated m-file. Launches the GUI, returns a handle,
%   handles callbacks. Main entrance point to the NUTMEG toolbox.
%

% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)
% varargout  cell array for returning output args (see VARARGOUT);

warning('off')  % hack so warning for "created by matlab 7.0" is turned off
% --Here initialization code starts.
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @nutmeg_OpeningFcn, ...
                   'gui_OutputFcn',  @nutmeg_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin & ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1}); 
end
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% Here initialization ends. 'Do not edit' they ask.
warning('on')  % we want to see all other warnings.


%%---------------------------------------------------
function nutmeg_OpeningFcn(hObject, eventdata, handles, varargin)
% --- Executes just before nutmeg is made visible.
handles.output = hObject;
guidata(hObject, handles);
set(handles.nut_logo_axes, 'Units', 'pixels');
banner = imread(['nutmeg.png']); % Read the image file banner.jpg
info = imfinfo(['nutmeg.png']); % Determine the size of the image file
position = get(handles.nut_logo_axes, 'Position');
axes(handles.nut_logo_axes);
bannerhandle=image(banner);
set(handles.nut_logo_axes, ...
    'Visible', 'off', ...
    'Units', 'pixels');
set(bannerhandle,'ButtonDownFcn','nut_about');

watchon;

% compatibility checks -----------------------
if(~exist('spm','file'))
    errordlg('NUTMEG requires SPM2 or SPM8 in your MATLAB path. They may be downloaded from http://www.fil.ion.ucl.ac.uk/spm')
    return
end

if(strcmp(spm('ver'),'SPM8b') || strcmp(spm('ver'),'SPM8'))
    disp('Be patient, don''t click anything yet, SPM8 takes a moment to open its MRI viewer.')
elseif(strcmp(spm('ver'),'SPM5'))
    errordlg(['NUTMEG is not compatible with SPM5. SPM2 or SPM8 may be downloaded from http://www.fil.ion.ucl.ac.uk/spm'])
    return
elseif(strcmp(spm('ver'),'SPM2'))
    % anything special to say?
else
    errordlg(['NUTMEG is not compatible with ' spm('ver') '. Please place SPM2 or SPM8 in your path, which may be downloaded from http://www.fil.ion.ucl.ac.uk/spm'])
    return
end

whichcell2mat = which('cell2mat');
if(strfind(whichcell2mat,'eeglab'))
    errordlg(['There is a conflict with EEGLAB''s version of cell2mat. Please remove ' whichcell2mat ' from your path.']);
    return
end

% hack to get mcc to recognize certain important functions exist
if(0)
    nut_select_VOI;
    nut_Default_Beamformer;
end
% end compatibility checks -------------------

% set up paths
nut_setpath

% nutmegpath = fileparts(which('nutmeg'));
% beamformerpath = fullfile(nutmegpath,'beamformers'); % path to nutmeg/beamformers
% externalpath = fullfile(nutmegpath,'external'); % path to nutmeg/external
% denoiserpath = fullfile(nutmegpath,'denoisers'); % path to nutmeg/denoisers
% dataimporterpath = fullfile(nutmegpath,'data_importers');
% lfpath = fullfile(nutmegpath,'leadfield_obtainers');
% templatespath = fullfile(nutmegpath,'templates');
% nuteegpath = fullfile(nutmegpath,'nuteeg');
% nuteegutilpath = fullfile(nuteegpath,'util');
% nuteegutilsphpath = fullfile(nuteegutilpath,'Sphere_Tessellation');
% addpath(beamformerpath,externalpath,denoiserpath,dataimporterpath,...
%     lfpath,templatespath,nuteegpath,nuteegutilpath,nuteegutilsphpath);

nut_defaults;

if ( ~isempty(varargin) && isstruct(varargin{1}) )  % if nuts structure given as input
    global nuts screws st defaults
    nuts = varargin{1};
    spm('Defaults','fMRI');
    set(handles.nut_megfile,'String','(none loaded)');
    nut_enabler(handles);
    nut_refresh_image;
else
    nut_reset;   % flush everything and set up proper globals, etc. 
    global nuts
end

nuts.fig = hObject;

nut_spmfig_setup;

monpos = get(0,'MonitorPositions');
figpos = get(nuts.fig,'Position');
set(nuts.fig,'Position',[figpos(1) figpos(1)+monpos(1,2) figpos(3) figpos(4)]);

watchoff;




%%---------------------------------------------------
function varargout = nutmeg_OutputFcn(hObject, eventdata, handles)
% --- Outputs from this function are returned to the command line.
varargout{1} = handles.output;

