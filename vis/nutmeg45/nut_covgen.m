function varargout = nut_covgen(varargin)
% NUT_COVGEN M-file for nut_covgen.fig
%   A standard .fig-associated m-file. Launches the GUI, returns a handle,
%   handles callbacks. Main entrance point to the NUT_COVGEN toolbox.
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
                   'gui_OpeningFcn', @nut_covgen_OpeningFcn, ...
                   'gui_OutputFcn',  @nut_covgen_OutputFcn, ...
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
function nut_covgen_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to crap (see VARARGIN)

% Choose default command line output for crap
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes crap wait for user response (see UIRESUME)
% uiwait(handles.figure1);

global nuts

% load with existing dataset info, if any
if(isfield(nuts,'megfile'))
    set([handles.nut_activeds, handles.nut_controlds],'String',nuts.meg.filename);
end
if(isfield(nuts,'preprocessing'))
    % set active time window
    set(handles.nut_activestart_text,'String',nuts.preprocessing.bf_timeinterval(1));
    set(handles.nut_activeend_text,'String',nuts.preprocessing.bf_timeinterval(2));
    % set control time window
    set(handles.nut_controlstart_text,'String',nuts.preprocessing.bf_timeinterval(1));
    set(handles.nut_controlend_text,'String',nuts.preprocessing.bf_timeinterval(2));
end
    



%%---------------------------------------------------
function varargout = nut_covgen_OutputFcn(hObject, eventdata, handles)
% --- Outputs from this function are returned to the command line.
varargout{1} = handles.output;


% --------------------------------------------------------------------
function nut_load_Rcon_menu_Callback(hObject, eventdata, handles)
% hObject    handle to nut_load_Rcon_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function nut_special_menu_Callback(hObject, eventdata, handles)
% hObject    handle to nut_special_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function nut_help_menu_Callback(hObject, eventdata, handles)
% hObject    handle to nut_help_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in nut_selectactive_button.
function nut_selectactive_button_Callback(hObject, eventdata, handles)
% hObject    handle to nut_selectactive_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global nuts ndefaults;

if exist('megdata_fullpath','var')
    [megdata_path,megdata_filename,megdata_ext]=fileparts(megdata_fullpath);
    megdata_path = [megdata_path '/'];  % to maintain consistency with uigetfile output
    megdata_filename = [megdata_filename megdata_ext];
else
    [megdata_filename, megdata_path]=uigetfile('*.bin;*.meg4;*.txt','Select MEG data file (BTi,CTF,KIT)...');
    if isequal(megdata_filename,0)|isequal(megdata_path,0)
        return;
    end
end

if(strcmp(megdata_path((end-3):(end-1)),'.ds'))
    megsys = 'CTF'
    nuts.meg.filename = megdata_path(1:(end-1));
elseif(strcmp(megdata_filename((end-2):end),'.ds'))
    megsys = 'CTF'
    nuts.meg.filename = megdata_fullpath;
else  % otherwise it's something else
    [crap,morecrap,ext] = fileparts(megdata_filename);
    switch(ext)
        case '.bin'
            megsys = 'BTi'
        case '.txt'
            megsys = 'KIT'
        case '.mat'
            megsys = 'matfile'
        otherwise
            errordlg('Sorry, unknown data type.');
    end
    nuts.meg.filename = fullfile(megdata_path,megdata_filename);
end

[crap,megfile,ext] = fileparts(nuts.meg.filename);
clear crap morecrap;
megfile = [megfile ext];

switch(hObject)
    case handles.nut_selectactive_button
        set([handles.nut_activeds handles.nut_controlds],'String',megfile);
    case handles.nut_selectcontrol_button
        set(handles.nut_controlds,'String',megfile);
    otherwise
        warning('hmm. something''s wonky here.');
end


% --- Executes on button press in nut_selectcontrol_button.
function nut_selectcontrol_button_Callback(hObject, eventdata, handles)
% hObject    handle to nut_selectcontrol_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nut_selectactive_button_Callback(hObject,eventdata,handles);


% --- Executes on button press in nut_gencov_button.
function nut_gencov_button_Callback(hObject, eventdata, handles)
% hObject    handle to nut_gencov_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global nuts ndefaults

activeds = get(handles.nut_activeds,'String');
[crap,morecrap,ext] = fileparts(activeds);
switch(ext)
    case '.ds'
        megsys = 'CTF'
    case '.bin'
        megsys = 'BTi'
    case '.txt'
        megsys = 'KIT'
    case '.mat'
        megsys = 'matfile'
    otherwise
        errordlg('Sorry, unknown data type.');
end

switch(megsys)
    case 'CTF'
        currentdir = pwd;  % remember current directory so we can get back to it later

        [crap, filename, ext] = fileparts(activeds);
        clear crap;

        ctf = ctf_read_res4(activeds);
        set(handles.nut_activestart_text,'String','0');
        set(handles.nut_activeend_text,'String',num2str(ctf.setup.end_msec));
        set(handles.nut_controlstart_text,'String',num2str(ctf.setup.start_msec));
        set(handles.nut_controlend_text,'String','0');
        if ctf.setup.number_trials > 3
            trials = 1:ctf.setup.number_trials;
            
            dstime = 'all';  % for RAM conservation, may want to adjust this

            R = zeros(length(ctf.sensor.index.all_sens));
            datamean = zeros(ctf.setup.number_samples, length(ctf.sensor.index.all_sens));
            conserveRAM = false
            if(conserveRAM)
                tic
                for trialidx=trials
                    ctf = ctf_read(activeds,'meg',dstime,trialidx,[],0);
                    %%%%%%%%%%%%%%% filtering would go here %%%%%%%%%%%%%%%%
                    datamean = datamean + ctf.data;
                    R = R + ctf.data' * ctf.data;
                end
                toc
            else
                % ramneeded: *4 for single precision, *8 for double precision
                ramneeded = length(ctf.sensor.index.all_sens) * ctf.setup.number_samples * length(trials) * 4;
                chunks = ceil(ramneeded/(ndefaults.RAM*2^30)); % convert GB to bytes by multiplying 2^30
                chunklength = ceil(length(trials)/chunks);
                for ii=1:chunks
                    ctf = ctf_read(activeds,'meg',dstime,(1:chunklength) + (ii-1)*chunklength,[],0);
                    %%%%%%%%%%%%%%% filtering would go here %%%%%%%%%%%%%%%%
                    datamean = datamean + sum(ctf.data,3);
                    for jj=1:chunklength
                        R = R + ctf.data(:,:,jj)' * ctf.data(:,:,jj);
                    end
                end
            end
            
        else
            % if number of trials is 1, 2, or 3, dataset is most likely
            % averaged, and only first "trial" contains relevant info
            trials = 1:1;
        end
        nuts.meg{1}.R = R;
        nuts.meg{1}.datamean = datamean;
        
        

    case 'BTi'
    case 'KIT'
    case 'Neuromag'
    otherwise  % user cancelled
        return
end


%         nuts.preprocessing.denoisertype = bolts.denoise{get(handles.nut_whichdenoiser,'Value')};
%         denoise_str = strrep(['nut_' nuts.preprocessing.denoisertype],' ','_');
%         feval(denoise_str,handles,nuts.preprocessing.signalspace,timept1,timept2); % update bolts.meg, bolts.params.Rzz1 (bolts.meg should be updated from nuts.meg.data every time Beamformer Tool is opened otherwise the denoising will be applied multiple times)
% 
%         % call code that computes inverse of autocorrelation/covariance matrix
%         nuts.preprocessing.invtype = bolts.invtype{get(handles.nut_which_invtype,'Value')};
%         bolts.params.InvRzz = nut_inv(bolts.params.Rzz, nuts.preprocessing.invtype); % invert Rzz matrix and store in bolts


% --- Executes on button press in pushbutton73.
function pushbutton73_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton73 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu20.
function popupmenu20_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu20 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu20


% --- Executes during object creation, after setting all properties.
function popupmenu20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


