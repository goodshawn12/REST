
%   HADES_MAIN starts the HADES GUI and manages the Main window and the
%   callback to the other functions.
% 
%   Copyright (C) 2010, Cristina Campi, Annalisa Pascarella, Michele Piana, Alberto Sorrentino



function[]= hades_main(action,varargin)
if nargin<1,
   clear global;
   action='gui';
end;

feval(action,varargin{:});

return;

function[]=gui()
clear global
%%% definition of the global variables
global pf_gui;
global pf;
global fwd;
global b;
global data;

%%% setting the fields of the struct pf
pf.parameters=[];       %%% it can value 1 or 0 if the parameters are saved or not, respectively
pf.subjects_dir=[];     %%% path of the subjects directory 
pf.results_dir=[];      %%% path of the result directory     
pf.subject=[];          %%% name of the experiment  
pf.g_matrix =[];        %%% lead field matrix
pf.g_matrix_dir=[];     %%% lead field matrix directory
pf.vertices = [];       %%% vertices matrix
pf.vertices_dir = [];   %%% vertices matrix directory
pf.orient = [];         %%% orientations matrix (no mandatory)
pf.orient_dir = [];     %%% orientations matrix directory (no mandatory)  
pf.cm_evol=[];          %%% radius [cm] used for computing the evolution matrix. It is not necessary if the 
                        %%% evolution matrix is already avalaible 
pf.evol = [];           %%% evolution matrix
pf.evol_dir=[];         %%% evolution matrix directory
pf.ssp = [];            %%% ssp matrix (not mandatory)
pf.ssp_dir = [];        %%% ssp matrix directory (not mandatory)
pf.cov_matrix=[];       %%% noise covariance matrix (not mandatory)
pf.cov_matrix_dir=[];   %%% noise covariance matrix directory (not mandatory)
pf.sigma=[];            %%% mean variance of noise during the pre-stimulus
pf.sigma_par=1;         %%% discrepancy parameter
pf.np=[];               %%% number of particles
pf.zero_time=[];        %%% first time point [s] of the data
pf.final_time=[];       %%% final time point [s] of the data
pf.freq=[];             %%% frequency
pf.time_interval=[];    %%% vector containig the time interval
pf.dataset_number=[];   %%% number of the dataset that will be analyzed
pf.t1=[];               %%% the first time point of the analysis
pf.t2=[];               %%% the last time point of the analysis
pf.data = [];           %%% data matrix 
pf.data_dir =[];        %%% data matrix directory
pf.p1=[];               %%% first pre-stimulus time point
pf.p2=[];               %%% last pre-stimulus time point
pf.autosave=0;          %%% parameter for auto-save the particle fitler results
pf.matrix_pw=[];            %%% auxiliary matrix for the noise covariance  


%%% the fields of the struct fwd contain information about the forward problem
fwd.check=0;
fwd.fig=[];
fwd.delta='';
fwd.edge_cm='';
fwd.str_fwd_fif=[];
fwd.str_lead_field=[];
fwd.str_vert=[];
fwd.str_orient=[];
fwd.name_fwd_fif=[];
fwd.location_fwd_fif=[];
fwd.name_lead_field=[];
fwd.location_lead_field=[];
fwd.name_vertices=[];
fwd.location_vertices=[];
fwd.name_orient=[];
fwd.location_orient=[];
fwd.name_versors=[];
fwd.location_versors=[];
fwd.name_sensors=[];
fwd.location_sensors=[];
fwd.name_orient=[];
fwd.location_orient=[];

%%% the fields of the struct b contain information about data
b.fig=[];
b.str_b=[];
b.name=[];
b.dataset_number=[];
b.cmd_dataset=[];
b.sigma=[];
b.str_sigma=[];
b.p1=[];
b.p2=[];
b.name_b=[];
b.location_b=[];
b.t1=[];
b.t2=[];
b.str_t1=[];
b.str_t2=[];
b.zero_time=[];
b.final_time=[];
b.freq=[];
b.label_ptimes_fif=[];
b.label_times_dat=[];
b.data_show=[];
b.fif_tot=[];
b.dataset_tot=[];


data=[];

scrsz = get(0,'ScreenSize');

pf_gui.fig=figure('HandleVisibility','Callback','Menubar','none',...
    'Name','HADES main', 'NumberTitle','off', ...
        'Visible','on', 'BackingStore','off');


pos=get(pf_gui.fig,'pos');
pos(1,3)=950;
pos(1,4)=500;
pos(1,1)=(scrsz(1,3)-pos(1,3))/2;
pos(1,2)=(scrsz(1,4)-pos(1,4))/2;

set(pf_gui.fig,'pos',pos);

         
    
pf_gui.panel1=uipanel('parent',pf_gui.fig,'Title','Input files','units','pixel',...
    'Position',[10 210 930 290],'FontWeight','bold');
%%% Command button for setting the subject directory
pf_gui.cmd_subject_dir=uicontrol(pf_gui.panel1,'style', 'pushbutton',...
        'units', 'pixel',...
        'position',[10 240 150 30],...
        'enable','on',...
        'callback', 'hades_main(''load_subject_dir'');' ,...
        'string', 'subject directory');
pf_gui.str_subjects_dir=uicontrol(pf_gui.panel1,'style', 'edit',...
        'units', 'pixel',...
        'position',[170 240  750 30],...
        'visible','off',...
        'enable','off',...
        'string', '');
%%% Command button for loading the forward problem
pf_gui.cmd_fwd_yes = uicontrol(pf_gui.panel1,'style', 'pushbutton',...
        'units', 'pixel',...
         'position', [10 200 150 30],...
        'string', 'load lead field',...
        'enable','on',...
        'callback', 'hades_main(''fwd_yes'');');
pf_gui.str_fwd_yes = uicontrol(pf_gui.panel1,'style', 'edit',...
        'units', 'pixel',...
         'position', [170 200 750 30],...
        'string', '',...
        'enable','off',...
        'visible','off');
% pf_gui.cmd_fwd_no = uicontrol(pf_gui.panel1,'style', 'pushbutton',...
%           'string', 'compute lead field',...
%           'units', 'pixel',...
%          'position', [170 200 150 30],...      
%         'enable','on',...
%         'callback', 'hades_main(''fwd_no'');'); 
%%% Command button and string for loading the evolution matrix  
pf_gui.cmd_evol=uicontrol(pf_gui.panel1,'style', 'pushbutton',...
        'units', 'pixel',...
        'position', [10 160 150 30],...
        'string', 'load neighbours',...
        'enable','on',...
        'callback', 'hades_main(''load_evol'');');
pf_gui.str_evol=uicontrol(pf_gui.panel1,'style', 'edit',...
        'units', 'pixel',...
        'string', '',...
        'position', [170 160 750 30],...
        'enable','off',...
        'visible','off');  
%%% string for computing the evolution matrix
pf_gui.label_compute = uicontrol(pf_gui.panel1,'style', 'text',...
        'units', 'pixel',...
        'string', 'or compute neighbours',...
        'position', [10 137 150 15]);
pf_gui.flag_evol = uicontrol(pf_gui.panel1,'style', 'checkbox',...
        'units', 'pixel',...
        'string', ' ',...
        'enable','on',...
        'position', [170 130 30 30]); 
pf_gui.label_within = uicontrol(pf_gui.panel1,'style', 'text',...
        'units', 'pixel',...
        'string', 'within',...
        'position', [190 137 40 15]);
pf_gui.str_cm_evol = uicontrol(pf_gui.panel1,'style', 'edit',...
        'units', 'pixel',...
        'string', ' ',...
        'position', [240 130 30 30]);  
 pf_gui.label_cm = uicontrol(pf_gui.panel1,'style', 'text',...
        'units', 'pixel',...
        'string', 'cm',...
        'position', [270 137 30 15]);
 
%%% Command button for loading the ssp matrix
pf_gui.cmd_ssp =uicontrol(pf_gui.panel1,'style', 'pushbutton',...
        'units', 'pixel',...
        'position', [10 90 150 30],...
        'enable','on',...
        'callback', 'hades_main(''load_ssp'');' ,...
        'string', 'load ssp');
pf_gui.str_ssp = uicontrol(pf_gui.panel1,'style', 'edit',...
        'units', 'pixel',...
        'string', '',...
        'enable','off',...
        'visible','off',...
        'position', [170 90 590 30]);
pf_gui.no_ssp =uicontrol(pf_gui.panel1,'style', 'pushbutton',...
        'units', 'pixel',...
        'position', [770 90 150 30],...
        'enable','on',...
        'callback', 'hades_main(''remove_ssp'');' ,...
        'string', 'remove ssp');   
%%% Command button for loading the covariance noise matrix
pf_gui.cmd_cov_matrix =uicontrol(pf_gui.panel1,'style', 'pushbutton',...
        'units', 'pixel',...
        'position',  [10 50 150 30],...
        'enable','on',...
        'callback', 'hades_main(''load_cov_matrix'');' ,...
        'string', 'load noise covariance');
pf_gui.str_cov_matrix = uicontrol(pf_gui.panel1,'style', 'edit',...
        'units', 'pixel',...
        'string', '',...
        'enable','off',...    
        'visible','off',...    
        'position', [170 50 590 30]);
pf_gui.no_cov =uicontrol(pf_gui.panel1,'style', 'pushbutton',...
        'units', 'pixel',...
        'position', [770 50 150 30],...
        'enable','on',...
        'callback', 'hades_main(''remove_cov'');' ,...
        'string', 'remove noise covariance');   
   
%%% Command button for loading the data 
pf_gui.cmd_boot_no = uicontrol(pf_gui.panel1,'style', 'pushbutton',...
        'units', 'pixel',...
        'string', 'load data',...
        'position', [10 10 150 30],...
        'enable','on',...
        'callback', 'hades_main(''boot_no'');'); 
pf_gui.str_boot_no = uicontrol(pf_gui.panel1,'style', 'edit',...
        'units', 'pixel',...
        'string', '',...
        'position', [170 10 750 30],...
        'visible','off',...
        'enable','off');     
pf_gui.panel2=uipanel('parent',pf_gui.fig,'Title','PF parameters','units','pixel',...
    'Position',[10 150 930 50],'FontWeight','bold');    
%%% String for setting the number of particles  
pf_gui.label_np = uicontrol(pf_gui.panel2,'style', 'text',...
        'units', 'pixel',...
        'string', 'number of particles',...
        'position', [10 10 150 15]);
pf_gui.str_np = uicontrol(pf_gui.panel2,'style', 'edit',...
        'units', 'pixel',...
        'string', '20000',...
        'position', [170 5 100 30]);
pf_gui.label_sigma_par = uicontrol(pf_gui.panel2,'style', 'text',...
        'units', 'pixel',...
        'string', 'discrepancy',...
        'position', [405 10 80 15]);
pf_gui.str_sigma_par = uicontrol(pf_gui.panel2,'style', 'edit',...
        'units', 'pixel',...
        'string', '1',...
        'position', [495 5 100 30]);
 
pf_gui.panel3=uipanel('parent',pf_gui.fig,'Title','Output file','units','pixel',...
    'Position',[10 40 930 100],'FontWeight','bold');       
%%% Command button and string for setting the result directory 
pf_gui.cmd_results_dir=uicontrol(pf_gui.panel3,'style', 'pushbutton',...
        'units', 'pixel',...
        'position',[10 50 150 30],...
        'enable','on',...
        'callback', 'hades_main(''load_results_dir'');' ,...
        'string', 'results directory');
pf_gui.str_results_dir=uicontrol(pf_gui.panel3,'style', 'edit',...
        'units', 'pixel',...
        'position',[170 50 750 30],...
        'visible','off',...
        'enable','off',...
        'string', '');
%%% Command button and string for loading the experiment
pf_gui.label_subject = uicontrol(pf_gui.panel3,'style', 'text',...
        'string', 'subject/experiment',...
        'units', 'pixel',...
        'position', [10 20 150 15]);
pf_gui.str_subject = uicontrol(pf_gui.panel3,'style', 'edit',...
        'string', ' ',...
        'visible','on',...
        'enable','on',...
        'units', 'pixel',...
        'position', [170 10 635 30]);
pf_gui.label_autosave = uicontrol(pf_gui.panel3,'style', 'text',...
        'string', 'auto-save',...
        'units', 'pixel',...
        'position', [815 20 75 15]);    
pf_gui.flag_autosave = uicontrol(pf_gui.panel3,'style', 'checkbox',...
        'units', 'pixel',...
        'string', ' ',...
        'enable','on',...
        'position', [890 10 30 40]); 
 

%%% Command button for saving the current parameters contained in the
%%% struct pf
pf_gui.cmd_save=uicontrol(pf_gui.fig,'style', 'pushbutton',...
        'units', 'pixel',...
        'position', [30 5 150 30],...
        'enable','on',...
        'callback', 'hades_main(''save_parameters'');',...
        'string', 'save settings');  
pf_gui.cmd_parameters =uicontrol(pf_gui.fig,'style', 'pushbutton',...
        'units', 'pixel',...
        'position',[220 5 150 30],...
        'enable','on',...
        'callback', 'hades_load_parameters;' ,...
        'string', 'load settings');
%%% Command button for runnig the particle filter
pf_gui.cmd_GO=uicontrol(pf_gui.fig,'style', 'pushbutton',...
        'units', 'pixel',...
        'position', [410  5 150 30],...
        'enable','on',...
        'callback', 'hades_main(''GO'');',...
        'string', 'run PF ');
%%% Command button for the post processing
pf_gui.cmd_post_processing=uicontrol(pf_gui.fig,'style', 'pushbutton',...
        'units', 'pixel',...
        'position', [600 5 150 30],...
        'enable','on',...
        'callback', 'hades_main(''post_processing'');',...
        'string', 'load results');
%%% exit    
pf_gui.cmd_exit=uicontrol(pf_gui.fig,'style', 'pushbutton',...
        'units', 'pixel',...
        'position', [790  5 150 30],...
        'enable','on',...
        'callback', 'hades_main(''exit'');',...
        'string', 'exit');
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% functions for loading the subject directory
function[]=load_subject_dir()
global pf;
global pf_gui;
start_path=pwd;
pf.subjects_dir=uigetdir(start_path, 'Choose the subject directory');
if pf.subjects_dir == 0
else
    set(pf_gui.str_subjects_dir,'visible','on','string',strtrim(pf.subjects_dir));
    if strcmp(get(pf_gui.str_results_dir,'string'),'')==1
        set(pf_gui.str_results_dir,'visible','on','string',strtrim(pf.subjects_dir));
        pf.results_dir=strtrim(pf.subjects_dir);
    end
    cd(pf.subjects_dir);
end
return

%%% functions for loading the results directory
function[]=load_results_dir()
global pf;
global pf_gui;
start_path=pwd;
pf.results_dir=uigetdir(start_path, 'Choose the folder where the results will be saved');
if pf.results_dir == 0
else
set(pf_gui.str_results_dir,'visible','on','string',strtrim(pf.results_dir));
end
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% function for loading the lead field
function[]=fwd_yes()
global pf;
global fwd;
fwd.check=0;
if isempty(findobj('Name','Loading the lead field','-and','-depth',6))==0
    %%% show a window for loading the lead field if it is already
    figure(fwd.fig)
else         
    scrsz = get(0,'ScreenSize');
    fwd.fig=figure('HandleVisibility','Callback','Menubar','none',...
        'Name','Loading the lead field', 'NumberTitle','off', ...
            'Visible','on', 'BackingStore','off');
    pos(1,3)=500;
    pos(1,4)=300;
    pos(1,1)=(scrsz(1,3)-pos(1,3))/2;
    pos(1,2)=(scrsz(1,4)-pos(1,4))/2;
    set(fwd.fig,'pos',pos);
    %%% Command button for loading a .fif file for the lead field 
    fwd.cmd_fwd_fif=uicontrol(fwd.fig,'style', 'pushbutton',...
            'units', 'pixel',...            
            'position', [10 240 100 30],...
            'enable','on',...
            'callback', 'hades_main(''load_fwd_fif'');' ,...
            'string', 'load fif');
    fwd.str_fwd_fif=uicontrol(fwd.fig,'style', 'edit','units', 'pixel',...
            'position', [110 240 400 30],...
            'enable','off');
    if isempty(fwd.location_fwd_fif)==0
        set(fwd.str_fwd_fif,'string', fullfile(fwd.location_fwd_fif,fwd.name_fwd_fif));
    elseif isempty(fwd.location_fwd_fif)==1 && isempty(pf.g_matrix_dir)==0
        if strcmp(pf.g_matrix_dir(end-3:end),'.fif')==1
            set(fwd.str_fwd_fif,'string', pf.g_matrix_dir);
        end
    end
    fwd.label_field=uicontrol(fwd.fig,'style', 'text',...
            'units', 'pixel',...
            'string', 'or',...
            'position', [10 210 40 15]);
    %%% Command button for loading the lead field matrix  
    fwd.cmd_lead_field=uicontrol(fwd.fig,'style', 'pushbutton',...
            'units', 'pixel',...
            'position', [10 160 100 30],...
            'enable','on',...
            'callback', 'hades_main(''load_lead_field'');' ,...
            'string', 'lead field');
    fwd.str_lead_field=uicontrol(fwd.fig,'style', 'edit','units', 'pixel',...
            'position', [110 160 400 30],...
            'enable','off',...
            'visible','on');   
    if isempty(fwd.location_lead_field)==0
        set(fwd.str_lead_field,'string', fullfile(fwd.location_lead_field,fwd.name_lead_field));
    elseif isempty(fwd.location_lead_field)==1 && isempty(pf.g_matrix_dir)==0
        if strcmp(pf.g_matrix_dir(end-3:end),'.dat')==1
            set(fwd.str_lead_field,'string', pf.g_matrix_dir);
        else
            set(fwd.str_lead_field,'string', '');
        end
    else
        set(fwd.str_lead_field,'string', '');
    end        
    %%% Command button for loading the vertices matrix 
    fwd.cmd_lead_vert=uicontrol(fwd.fig,'style', 'pushbutton','units', 'pixel',...
            'position', [10 130 100 30],...
            'enable','on',...
            'callback', 'hades_main(''load_vertices'');' ,...
            'string', 'vertices'); 
    fwd.str_vert=uicontrol(fwd.fig,'style', 'edit','units', 'pixel',...
            'position', [110 130 400 30],...
            'visible','on',...
            'enable','off',...
            'string', fullfile(fwd.location_vertices,fwd.name_vertices)); 
    if isempty(fwd.location_vertices)==0
        set(fwd.str_vert,'string', fullfile(fwd.location_vertices,fwd.name_vertices));
    elseif isempty(fwd.location_vertices)==1 && isempty(pf.vertices_dir)==0
        set(fwd.str_vert,'string', pf.vertices_dir);
    else
        set(fwd.str_vert,'string', '');
    end         
    %%% Command button for loading the orientations matrix      
    fwd.cmd_lead_orient=uicontrol(fwd.fig,'style', 'pushbutton','units', 'pixel',...
            'position', [10 100 100 30],...
            'enable','on',...
            'callback', 'hades_main(''load_orient'');' ,...
            'string', 'orientations');  
    fwd.str_orient=uicontrol(fwd.fig,'style', 'edit','units', 'pixel',...
            'position', [110 100 400 30],...
            'visible','on',...
            'enable','off',...
            'string', fullfile(fwd.location_orient,fwd.name_orient)); 
    if isempty(fwd.location_orient)==0
        set(fwd.str_orient,'string', fullfile(fwd.location_orient,fwd.name_orient));
    elseif isempty(fwd.location_orient)==1 && isempty(pf.orient_dir)==0
        set(fwd.str_orient,'string', pf.orient_dir);
    else
        set(fwd.str_orient,'string', '');
    end     
    fwd.remove_orient=uicontrol(fwd.fig,'style', 'pushbutton','units', 'pixel',...
            'position', [10 70 150 30],...
            'enable','on',...
            'callback', 'hades_main(''remove_orient'');' ,...
            'string', 'remove orientations');  
    %%% Command button for canceling this window    
    fwd.cmd_exit=uicontrol(fwd.fig,'style', 'pushbutton',...
            'units', 'pixel',...
            'position', [350 20 150 30],...
            'enable','on',...
            'callback', 'hades_main(''cancel_fwd'');',...
            'string', 'cancel');
    %%% Command button for storing in the struct pf the loaded file(s)     
    fwd.cmd_return=uicontrol(fwd.fig,'style', 'pushbutton',...
            'units', 'pixel',...
            'position', [150 20 150 30],...
            'enable','on',...
            'callback', 'hades_main(''done_fwd_yes'');',...
            'string', 'ok');
end
return;

%%% function for loading the .fif file containing the lead field
function[]=load_fwd_fif()
global fwd;
global pf;
if isempty(pf.subjects_dir)==0
    if pf.subjects_dir==0
    else
        cd(pf.subjects_dir);
    end
end
[fwd.name_fwd_fif,fwd.location_fwd_fif]=uigetfile('*.fif', 'Choose the <subject>-fwd.fif file');
if  fwd.name_fwd_fif==0
else
    set(fwd.str_fwd_fif,'visible','on','string',fullfile(fwd.location_fwd_fif,fwd.name_fwd_fif));
end
return;

%%% function for loading the lead field matrix
function[]=load_lead_field()
global fwd;
global pf;
if isempty(pf.subjects_dir)==0
    if pf.subjects_dir==0
    else
        cd(pf.subjects_dir);
    end
end
[fwd.name_lead_field,fwd.location_lead_field]=uigetfile('*.dat;*.mat', 'Choose the lead field file');
if  fwd.name_lead_field==0
else
    str_lead_field=fullfile(fwd.location_lead_field,fwd.name_lead_field);
    set(fwd.str_lead_field,'visible','on','string',str_lead_field);
end
return;
%%% function for loading the vertices matrix
function[]=load_vertices()
global fwd;
global pf;
if isempty(pf.subjects_dir)==0
    cd(pf.subjects_dir);
end
[fwd.name_vertices,fwd.location_vertices]=uigetfile('*.dat;*.mat', 'Choose the vertices file');
if  fwd.name_vertices==0
else
    str_vert=fullfile(fwd.location_vertices,fwd.name_vertices);
    set(fwd.str_vert,'visible','on','string',str_vert);
end
return;
%%% function for loading the orientations matrix
function[]=load_orient()
global fwd;
global pf;
if isempty(pf.subjects_dir)==0
    cd(pf.subjects_dir);
end
[fwd.name_orient,fwd.location_orient]=uigetfile('*.dat;*.mat', 'Choose the orientations file');
if  fwd.name_orient==0
    set(fwd.str_orient,'visible','on','string','');
else
    str_orient=fullfile(fwd.location_orient,fwd.name_orient);
    set(fwd.str_orient,'visible','on','string',str_orient);
end
return;
function[]=remove_orient()
global fwd;
global pf;
    set(fwd.str_orient,'visible','off','string','');
    fwd.location_orient=[];
    pf.orient_dir=[];
return;

%%% function for storing in the struct pf the loaded file(s)
function[]=done_fwd_yes()
global fwd;
global pf_gui;
name_fwd_fif=fwd.name_fwd_fif;
name_lead_field=fwd.name_lead_field;
str_fwd_fif=get(fwd.str_fwd_fif,'string');
str_lead_field=get(fwd.str_lead_field,'string');
if isempty(str_fwd_fif)==0 && isempty(str_lead_field)==1
    if name_fwd_fif==0
     msgbox('Please load the lead field','Error','error'); 
    else    
        hades_load_fwd_fif(str_fwd_fif);
        set(pf_gui.str_fwd_yes,'string',str_fwd_fif,'visible','on');
    end
elseif isempty(str_fwd_fif)==1 && isempty(str_lead_field)==0 
   if name_lead_field~=0
        hades_load_fwd_dat(str_lead_field);
        set(pf_gui.str_fwd_yes,'string',str_lead_field,'visible','on');
   else
        msgbox('Please load the lead field','Error','error');
   end
elseif (isempty(str_lead_field)==1 && isempty(str_fwd_fif)==1) 
    msgbox('Please load the lead field','Error','error'); 
elseif (strcmp(name_lead_field,'0')==1 && strcmp(name_fwd_fif,'0')==1)
    msgbox('Please load the lead field','Error','error'); 
else
    choice_fwd = questdlg('Please choose the lead field matrix you want to load', ...
            'Too many lead fields!', ...
            '.fif','.dat/.mat','.fif');
            % Handle response
        switch choice_fwd
        case '.fif'
            hades_load_fwd_fif(get(fwd.str_fwd_fif,'string'));	
            set(pf_gui.str_fwd_yes,'string',str_fwd_fif,'visible','on');
        case '.dat/.mat'
            hades_load_fwd_dat(get(fwd.str_lead_field,'string'));
            set(pf_gui.str_fwd_yes,'string',str_lead_field,'visible','on');
        end    
end
return

function[]=cancel_fwd()
global fwd;
global pf_gui;
global pf;
if isempty(pf.g_matrix_dir)==1
    set(pf_gui.str_fwd_yes,'string','','visible','off');
else
    set(pf_gui.str_fwd_yes,'string',pf.g_matrix_dir);
    ext=pf.g_matrix_dir(end-3:end);
    if strcmp(ext,'.fif')
        set(fwd.str_fwd_fif,'string',pf.g_matrix_dir);
    else
        set(fwd.str_lead_field,'string',pf.g_matrix_dir);
    end
end
set(fwd.str_vert,'string',pf.vertices_dir);
set(fwd.str_orient,'string',pf.orient_dir);
close

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% function for loading the evolution matrix
function[]=load_evol()
global pf_gui;
global pf;

if isempty(pf.subjects_dir)==0
    cd(pf.subjects_dir);
end
[evol_name,evol_location]=uigetfile('*.dat;*.mat', 'Choose the neighbours matrix');
if  evol_name==0
else
    str_evol=fullfile(evol_location,evol_name);
    set(pf_gui.str_evol,'visible','on','string',str_evol);
    pf.evol_dir=str_evol;
    ext = str_evol(end-3:end);
    if strcmp(ext,'.dat')==1
        pf.evol=load(pf.evol_dir);
    elseif strcmp(ext,'.mat')==1
        evol_mat=load(str_evol);        
        pf.evol=getfield(evol_mat,cell2mat(fieldnames(evol_mat)));
    end   
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% functions for loading the ssp and noise covariance matrices
function[]=load_ssp()
global pf_gui;
global pf;
if isempty(pf.subjects_dir)==0
    cd(pf.subjects_dir);
end
[ssp_name,ssp_location]=uigetfile('*.dat;*.mat', 'Choose the ssp matrix');
set(pf_gui.str_ssp,'visible','on');
if  ssp_name==0
else
    str_ssp=fullfile(ssp_location,ssp_name);
    set(pf_gui.str_ssp,'string',str_ssp);
    pf.ssp_dir=str_ssp;
    ext = str_ssp(end-3:end);
    if strcmp(ext,'.dat')==1
        pf.ssp=fullfile(str_ssp);
    elseif strcmp(ext,'.mat')==1
        ssp_mat=load(str_ssp);        
        pf.ssp=getfield(ssp_mat,cell2mat(fieldnames(ssp_mat)));
    end 
end
return

function[]=remove_ssp()
global pf_gui;
global pf;

set(pf_gui.str_ssp,'visible','off');
pf.ssp=[];
pf.ssp_dir=[];

return


function[]=load_cov_matrix()
global pf_gui;
global pf;
if strcmp(get(pf_gui.str_subjects_dir,'string'),'')==0
    cd(get(pf_gui.str_subjects_dir,'string'));
end
[cov_matrix_name,cov_matrix_location]=uigetfile('*.dat;*.mat', 'Choose the noise covariance matrix');
set(pf_gui.str_cov_matrix,'visible','on');
if  cov_matrix_name==0
else
    str_cov_matrix=fullfile(cov_matrix_location,cov_matrix_name);
    set(pf_gui.str_cov_matrix,'string', str_cov_matrix);
    pf.cov_matrix_dir= str_cov_matrix;
    ext = str_cov_matrix(end-3:end);
    if strcmp(ext,'.dat')==1
        pf.cov_matrix=load(str_cov_matrix);
    elseif strcmp(ext,'.mat')==1
        cov_matrix_mat=load(str_cov_matrix);        
        pf.cov_matrix=getfield(cov_matrix_mat,cell2mat(fieldnames(cov_matrix_mat)));
    end 
    n_ch = size(pf.cov_matrix,1);
    par_reg = 0.1;
    variances = diag(pf.cov_matrix);
    % MATRIX COV DIAGONALE 
    pf.cov_matrix_diag = diag(diag(pf.cov_matrix)); % matrix diagonale 
    [v_cov d_cov] = svd(pf.cov_matrix_diag);
    d_cov_reg = d_cov + par_reg*mean(variances)*eye(size(d_cov));
    d_i = 1./(diag(d_cov_reg)); %D^-1
    pf.inv_cov_matrix = v_cov*diag(d_i)*v_cov'; % SIGMA^-1 = V*D^-1*V'
    pf.matrix_pw = diag(sqrt(d_i))*v_cov';  % D^-1/2*V'   
end
return

function[]=remove_cov()
global pf_gui;
global pf;

set(pf_gui.str_cov_matrix,'visible','off');
pf.cov_matrix=[];
pf.cov_matrix_dir=[];
pf.matrix_pw=[];
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% functions for loading the magnetic field
function[]=boot_no()
global pf;
global b;
if isempty(findobj('Name','Loading the magnetic field','-and','-depth',6))==0
   figure(b.fig);
else
    scrsz = get(0,'ScreenSize');
    b.fig=figure('HandleVisibility','Callback','Menubar','none',...
        'Name','Loading the magnetic field', 'NumberTitle','off','Visible','on', 'BackingStore','off');
    pos=get(b.fig,'pos');
    pos(1,3)=500;
    pos(1,4)=300;
    pos(1,1)=(scrsz(1,3)-pos(1,3))/2;
    pos(1,2)=(scrsz(1,4)-pos(1,4))/2;set(b.fig,'pos',pos);
    %%% command button for loading the magnetic field
    b.cmd_b=uicontrol(b.fig,'style', 'pushbutton','units', 'pixel',...
            'position', [10 260 150 30],...
            'enable','on',...
            'callback', 'hades_main(''load_magnetic_field'');' ,...
            'string', 'load the magnetic field');
    b.str_b = uicontrol(b.fig,'style', 'edit',...
            'units', 'pixel',...
            'enable','off',...
            'string', '',...
            'position', [170 260 260 30]); 
    if isempty(pf.data_dir)==0
        set(b.str_b,'string',num2str(pf.data_dir));
    end
    %%% command button for showing the magnetic field     
    b.cmd_show_b=uicontrol(b.fig,'style', 'pushbutton',...
            'units', 'pixel',...
            'position', [450 260 50 30],...
            'enable','on',...
            'callback', 'hades_main(''show_b'');',...
            'string', 'show');
    %%% strings for the initial time point, the final time poin t and the frequency    
    b.label_times_dat = uicontrol(b.fig,'style', 'text',...
            'visible','off',...
            'units', 'pixel',...
            'string', 'initial time [ms]',...
             'position', [10 220 100 30]);
    b.zero_time = uicontrol(b.fig,'style', 'edit',...
            'visible','off',...
            'units', 'pixel',...
            'string', pf.zero_time,...
            'position', [120 220 50 30]);
    b.label_freq_dat = uicontrol(b.fig,'style', 'text',...
            'visible','off',...
            'units', 'pixel',...
            'string', 'frequency [Hz]',...
             'position', [180 220 100 30]);
    b.freq= uicontrol(b.fig,'style', 'edit',...
        'visible','off',...
            'units', 'pixel',...
            'string',pf.freq,...
            'position', [290 220 50 30]);
   %%% strings for the dataset number     
   b.label_ptimes_fif = uicontrol(b.fig,'style', 'text',...
            'units', 'pixel',...
              'visible','off',...
            'string', ' you loaded the dataset number',...
             'position', [10 220 170 30]);
    b.dataset_number = uicontrol(b.fig,'style', 'edit',...
            'units', 'pixel',...
            'visible','off',...
            'string', pf.dataset_number,...
            'enable','off',...
            'position', [190 220 50 30]);     
    b.cmd_dataset = uicontrol(b.fig,'style', 'pushbutton',...
            'units', 'pixel',...
            'visible','off',...
            'string', 'change dataset',...
            'enable','off',...
            'callback','hades_main(''load_data_fif'');',...
            'position', [250 220 150 30]);    
    aux_data=get(b.str_b,'string');    
    if strcmp(aux_data,'')==0
        ext = aux_data(end-3:end);
        if strcmp(ext,'.dat')==1 || strcmp(ext,'.mat')==1
            set(b.label_times_dat,'visible','on'); 
            set(b.label_freq_dat,'visible','on'); 
            set(b.zero_time,'visible','on'); 
            set(b.final_time,'visible','on');
            set(b.freq,'visible','on');
            set(b.dataset_number,'visible','off');  
            set(b.label_ptimes_fif,'visible','off');
            set(b.cmd_dataset,'visible','off'); 
        elseif strcmp(ext,'.fif')==1
            if isempty(pf.parameters)==1
                set(b.cmd_dataset,'visible','on','enable','on');
            else
                set(b.cmd_dataset,'visible','on','enable','off');
            end
            set(b.label_ptimes_fif,'visible','on');
            set(b.dataset_number,'visible','on')
            set(b.zero_time,'visible','off'); 
            set(b.final_time,'visible','off');
            set(b.freq,'visible','off');  
            set(b.label_freq_dat,'visible','off');  
            set(b.label_times_dat,'visible','off');  
        end
    end
    %%% strings for the pre stimulus interval
    b.label_prestimulus = uicontrol(b.fig,'style', 'text',...
            'units', 'pixel',...
            'string', 'pre-stimulus interval from [ms] ',...
             'position', [10 180 170 30]);
    b.p1 = uicontrol(b.fig,'style', 'edit',...
            'units', 'pixel',...
            'position', [190 180 50 30]);
    if isempty(pf.data)==0
        set(b.p1,'string',pf.p1);
    else
        set(b.p1,'string','');
    end    
    b.label_prestimulus_to = uicontrol(b.fig,'style', 'text',...
            'units', 'pixel',...
            'string', 'to',...
             'position', [250 180 50 30]);
    b.p2 = uicontrol(b.fig,'style', 'edit',...
            'units', 'pixel',...
            'position', [310 180 50 30]);
    if isempty(pf.data)==0
        set(b.p2,'string',pf.p2);
    else
        set(b.p2,'string','');
    end    
    %%% strings for select the time interval you want to analyze
    b.label_times = uicontrol(b.fig,'style', 'text',...
            'units', 'pixel',...
            'string', 'analyze from [ms]',...
            'position', [10 140 170 30]);
    b.str_t1 = uicontrol(b.fig,'style', 'edit',...
            'units', 'pixel',...
            'position', [190 140 50 30]);
    if isempty(pf.data)==0
        set(b.str_t1,'string',1000*pf.time_interval(pf.t1));
    else
        set(b.str_t1,'string','');
    end    
    b.label_times_to = uicontrol(b.fig,'style', 'text',...
            'units', 'pixel',...
            'string', 'to',...
             'position', [250 140 50 30]);    
    b.str_t2 = uicontrol(b.fig,'style', 'edit',...
            'units', 'pixel',...
            'position', [310 140 50 30]);  
    if isempty(pf.data)==0
        set(b.str_t2,'string',1000*pf.time_interval(pf.t2));
    else
        set(b.str_t2,'string','');
    end    
    %%% command button for showing the magnetic field and the mean noise variance during the prestimulus    
    b.cmd_show_sigma=uicontrol(b.fig,'style', 'pushbutton',...
            'units', 'pixel',...
            'position', [10 100 245 30],...
            'enable','on',...
            'callback', 'hades_main(''show_sigma'');',...
            'string', 'show the pre-stimulus standard deviation');
    %%% command button for computing the covariance noise from the prestimulus    
    b.cmd_compute_pw=uicontrol(b.fig,'style', 'pushbutton',...
            'units', 'pixel',...
            'position', [265 100 225 30],...
            'enable','on',...
            'callback', 'hades_main(''compute_pw'');',...
            'string', 'compute noise covariance matrix');    
    %%% string for showing and changing the value of the mean noise variance    
    b.label_sigma=uicontrol(b.fig,'style', 'text',...
                    'units', 'pixel',...
                    'string', 'the pre-stimulus standard deviation is',...
                    'position', [10 60 170 30]);
    b.str_sigma=uicontrol(b.fig,'style', 'edit',...
                    'units', 'pixel',...
                    'string', '',...
                    'enable','off',...
                    'position', [190 60 150 30]); 
    if isempty(pf.data)==0
        set(b.str_sigma,'string',num2str(pf.sigma));
    end
    b.cmd_update_sigma=uicontrol(b.fig,'style', 'pushbutton',...
                    'units', 'pixel',...
                    'string', 'update',...
                    'enable','on',...
                    'callback','hades_main(''load_prestimulus_interval'');',...
                    'position', [350 60 100 30]);           
    %%% commund button for storing in the struct pf the loaded quantities            
    b.cmd_return=uicontrol(b.fig,'style', 'pushbutton',...
            'units', 'pixel',...
            'position', [150 20 150 30],...
            'enable','on',...
            'callback', 'hades_main(''done_boot_no'');',...
            'string', 'ok');
    %%% command button for cancel this window    
    b.cmd_exit=uicontrol(b.fig,'style', 'pushbutton',...
            'units', 'pixel',...
            'position', [350 20 150 30],...
            'enable','on',...
            'callback', 'hades_main(''cancel_data'');',...
            'string', 'cancel');
end
return;
%%% function for loading the magnetic field
function[]=load_magnetic_field()
global b;
global pf;
global data;
if isempty(pf.subjects_dir)==0
    if exist(pf.subjects_dir)==0
        %msgbox(strcat(['The directoty ',pf.subjects_dir,' does not exist']),'Error','error');
    else
        cd(pf.subjects_dir)
    end
else
end
[b.name_b,b.location_b]=uigetfile('*.fif;*.dat;*.mat', 'Load the magnetic field');
if  b.name_b==0
    b_name=get(b.str_b,'string');
else
    b_name=fullfile(b.location_b,b.name_b);
    set(b.str_b,'visible','on','string',fullfile(b.location_b,b.name_b));
end
ext =  b_name(end-3:end);
pf.matrix_pw=[];
if strcmp(ext,'.dat')==1 || strcmp(ext,'.mat')==1
    set(b.label_times_dat,'visible','on'); 
    set(b.label_freq_dat,'visible','on'); 
    set(b.zero_time,'visible','on','string',''); 
    set(b.final_time,'visible','on','string','');
    set(b.freq,'visible','on','string','');
    set(b.dataset_number,'visible','off');    
    set(b.cmd_dataset,'visible','off'); 
    set(b.label_ptimes_fif,'visible','off');
    set(b.p1,'string','');
    set(b.p2,'string','');
    set(b.str_t1,'string','');
    set(b.str_t2,'string','');
    set(b.str_sigma,'string','');
    if strcmp(ext,'.dat')==1
        data=load(b_name);
    elseif strcmp(ext,'.mat')==1
        data_mat=load(b_name);      
        data=getfield(data_mat,cell2mat(fieldnames(data_mat)));
        pf.parameters=[];
    end       
elseif strcmp(ext,'.fif')==1
    set(b.label_ptimes_fif,'visible','on');
    set(b.dataset_number,'visible','on','string','');
    set(b.cmd_dataset,'visible','on','enable','on');
    set(b.zero_time,'visible','off'); 
    set(b.final_time,'visible','off');
    set(b.freq,'visible','off');  
    set(b.label_times_dat,'visible','off'); 
    set(b.label_freq_dat,'visible','off'); 
    set(b.p1,'string','');
    set(b.p2,'string','');
    set(b.str_t1,'string','');
    set(b.str_t2,'string','');
    set(b.str_sigma,'string','');
    b.fif_tot=fiff_read_evoked_all(b_name);
    b.dataset_tot=size(b.fif_tot.evoked,2);
    pf.parameters=[];
   load_data_fif
end

return;

function[]=load_data_fif
global b_fwd;
global b;
global pf;
scrsz = get(0,'ScreenSize');
b_fwd.fig=figure('HandleVisibility','Callback','Menubar','none',...
    'Name','Loading the dataset ', 'NumberTitle','off','Visible','on', 'BackingStore','off','WindowStyle','modal');
    pos=get(b_fwd.fig,'pos');
    pos(1,3)=300;
    pos(1,4)=80;
    pos(1,1)=(scrsz(1,3)-pos(1,3))/2;
    pos(1,2)=(scrsz(1,4)-pos(1,4))/2;
    set(b_fwd.fig,'pos',pos); 
    clear comment
    comment=[];
    for t=1:b.dataset_tot
        if t==1
            comment=[comment,b.fif_tot.evoked(1,t).comment];
        else
            comment=[comment,'|',b.fif_tot.evoked(1,t).comment];
        end
    end
    b_fwd.label= uicontrol(b_fwd.fig,'style', 'text',...
        'units', 'pixel',...
        'string', 'Select the dataset',...
         'position', [10 40 100 20]);
     if isempty(pf.dataset_number)==1
         val=1;
     else
         val=pf.dataset_number;
     end
    b_fwd.popup=uicontrol(b_fwd.fig,'style', 'popup',...
        'units', 'pixel',...
        'String', comment,...
        'value',val,...
        'enable','on',...
         'position', [110 30 190 30],...
        'callback','hades_main(''set_dataset_fif'');');
return
function[]=set_dataset_fif()
global b_fwd;
global pf;
global b;
global data;
val = get(b_fwd.popup,'Value');
pf.dataset_number=val;
set(b.dataset_number,'string',val);
aux_b_fif=fiff_read_evoked(fullfile(b.location_b,b.name_b),pf.dataset_number);
str_ch=aux_b_fif.info.ch_names;
for t=1:size(str_ch,2)-1
    aux1=cell2mat(str_ch(1,t));
    a1=aux1(1:3);
    aux2=cell2mat(str_ch(1,t+1));
    a2=aux2(1:3);
    if a1~=a2
        break
    end
end
data=aux_b_fif.evoked.epochs(1:t,:)';
pf.zero_time=double(aux_b_fif.evoked.first);
pf.freq = aux_b_fif.info.sfreq;
for i=1:size(data,1)
    pf.time_interval(i)=pf.zero_time/1000+(i-1)*(1/pf.freq);
end
pf.final_time=1000*pf.time_interval(end);
close

return    

%%% function for showing the magnetic field
function[]=show_b()
global b;
global pf;
global data;
if isempty(findobj('Name','Magnetic field','-and','-depth',6))==0
    figure(b.plot)
else
    if isempty(data)==1 
      msgbox('You have to load the magnetic field','Error', 'error');
    else
        b_name=get(b.str_b,'string');
        ext = b_name(end-3:end);
        if strcmp(ext,'.dat')==1 || strcmp(ext,'.mat')==1 
            b.plot=figure;
            set(b.plot,'Name','Magnetic field');
            plot(data);
            xlabel('samples');
            ylabel('magnetic field [T]');
         elseif strcmp(ext,'.fif')==1       
            b.plot=figure;
            set(b.plot,'Name','Magnetic field');
            plot(1000*pf.time_interval(1:size(data,1)),data);
            xlabel('time [ms]')
            ylabel('magnetic field [T]');
        end
    end
end
                                     
                    
return

function[]=load_prestimulus_interval()
global pf;
global b;
global data;
if isempty(data)==1
    msgbox( 'You have to load a magnetic field!','Error','error');
        return
else
    b_name=get(b.str_b,'string');
    ext =  b_name(end-3:end);
    if strcmp(ext,'.dat')==1 || strcmp(ext,'.mat')==1 
        pf.zero_time=str2num(get(b.zero_time,'string'));
        pf.freq = str2num(get(b.freq,'string'));
        if isempty(pf.zero_time)==1 || isempty(pf.freq)==1
            msgbox( 'You have to load the a correct initial time point and frequency!','Error','error');
            return
        else
            for t=1:size(data,1)
                pf.time_interval(t)=pf.zero_time/1000+(t-1)*(1/pf.freq);
            end
            pf.final_time=1000*pf.time_interval(end);
        end
    end
    pf.p1 = str2num(get(b.p1,'string'));
    pf.p2 = str2num(get(b.p2,'string'));
    if isempty(pf.p1)==1 || isempty(pf.p2)==1
        msgbox( 'You have to load the a correct time point!','Error','error')
        return
    else 
        if pf.zero_time>pf.p1
            msgbox('The first prestimulus time point has to be bigger than first time point',...
              'Error','error')
            return
        elseif pf.p1>pf.p2
            msgbox('The first prestimulus time point has to be smaller than the second one',...
              'Error','error')
            return
        elseif pf.p2>pf.final_time
              msgbox( 'The second prestimulus time point is bigger than the last time point',...
                     'Error','error');
                 return
        elseif pf.p1<=pf.p2 && pf.p2<=pf.final_time && pf.zero_time<=pf.p1
            aux1=repmat(pf.p1,1,size(pf.time_interval,2));
            [v1,prestim1]=min(abs(1000*pf.time_interval-aux1));
            aux2=repmat(pf.p2,1,size(pf.time_interval,2));
            [v2,prestim2]=min(abs(1000*pf.time_interval-aux2));
            b.sigma = sqrt(mean(var(data(prestim1:prestim2,:))));
            set(b.str_sigma,'string',num2str(b.sigma));    
        end
    end
end
return
        

function[]=show_sigma()
global b;
global pf;
global data;
b_name=get(b.str_b,'string');
if (isempty(data)==1 && isempty(b_name)==1) || strcmp(get(b.p1,'string'),'')==1 || strcmp(get(b.p2,'string'),'')==1
        msgbox('You have to load the magnetic field and the prestimulus interval','Error', 'error');
        return
elseif (isempty(data)==0 && isempty(b_name)==0) || strcmp(get(b.p1,'string'),'')==0 || strcmp(get(b.p2,'string'),'')==0
      load_prestimulus_interval
     
    if isempty(findobj('Name','Magnetic field','-and','-depth',6))==1
        b.plot=figure;
        set(b.plot,'Name','Magnetic field');
        plot(1000*pf.time_interval(1:size(data,1)),data);
        hold on
        plot(1000*pf.time_interval(1:size(data,1)),3*b.sigma,'ko','MarkerFaceColor','k'); 
        plot(1000*pf.time_interval(1:size(data,1)),-3*b.sigma,'ko','MarkerFaceColor','k');
        xlabel('time [ms]');
        ylabel('magnetic field [T]');
    else
        figure(b.plot);
        plot(1000*pf.time_interval(1:size(data,1)),data);
        hold on
        plot(1000*pf.time_interval(1:size(data,1)),3*b.sigma,'ko','MarkerFaceColor','k'); 
        plot(1000*pf.time_interval(1:size(data,1)),-3*b.sigma,'ko','MarkerFaceColor','k');
        xlabel('time [ms]');
    end
        set(b.str_sigma,'string',num2str(b.sigma));      
else
    msgbox( 'You forget something','Error', 'error');
    return
end
return

function[]=compute_pw()
global b;
global pf;
global data;
b_name=get(b.str_b,'string');
if (isempty(data)==1 && isempty(b_name)==1) || strcmp(get(b.p1,'string'),'')==1 || strcmp(get(b.p2,'string'),'')==1
        msgbox('You have to load the magnetic field and the prestimulus interval','Error', 'error');
        return
elseif (isempty(data)==0 && isempty(b_name)==0) || strcmp(get(b.p1,'string'),'')==0 || strcmp(get(b.p2,'string'),'')==0
      load_prestimulus_interval
      
      aux1=repmat(pf.p1,1,size(pf.time_interval,2));
      [aux_prestim1,prestim1]=min(abs(1000*pf.time_interval-aux1));
      aux2=repmat(pf.p2,1,size(pf.time_interval,2));
      [aux_prestim2,prestim2]=min(abs(1000*pf.time_interval-aux2));
      mean_prestimulus = repmat(mean(data(prestim1:prestim2,:)),size(data(prestim1:prestim2,:),1),1);
      pf.cov_matrix=(data(prestim1:prestim2,:) - mean_prestimulus)'*(data(prestim1:prestim2,:)-mean_prestimulus)/(prestim2-prestim1);
      n_ch = size(pf.cov_matrix,1);      
      par_reg = 0.1;
      variances = diag(pf.cov_matrix);      
      % MATRIX COV DIAGONALE
      pf.cov_matrix_diag = diag(diag(pf.cov_matrix)); % matrix diagonale     
      [v_cov d_cov] = svd(pf.cov_matrix_diag);      
      d_cov_reg = d_cov + par_reg*mean(variances)*eye(size(d_cov));      
      d_i = 1./(diag(d_cov_reg)); %D^-1      
      pf.inv_cov_matrix = v_cov*diag(d_i)*v_cov'; % SIGMA^-1 = V*D^-1*V'
      pf.matrix_pw= diag(sqrt(d_i))*v_cov';  % D^-1/2*V'
end
return
%%% function for storing the set values in the struct pf
function[]=done_boot_no()
global pf;
global pf_gui;
global b;
global data;
if isempty(data)==1 
  msgbox( 'You have to load the magnetic field','Error', 'error');
  return
else
    load_prestimulus_interval
    pf.sigma=b.sigma;
    %%%%check on the time interval
    b.t1 = str2num(get(b.str_t1,'string'));
    b.t2 = str2num(get(b.str_t2,'string'));
    if isempty(b.t1)==0 && isempty(b.t2)==0
        aux1=repmat(b.t1,1,size(pf.time_interval,2));
        [v1,b_step1]=min(abs(1000*pf.time_interval-aux1));
        aux2=repmat(b.t2,1,size(pf.time_interval,2));
        [v2,b_step2]=min(abs(1000*pf.time_interval-aux2)); 
        pf.t1 = b_step1;
        pf.t2 = b_step2;
        if b.t1>b.t2
            msgbox('Attention, t1>t2','Error','error');
            return
        end
        if b.t2>pf.final_time
            msgbox( 'Attention, the final time is bigger than the time duration of the data','Error','error');
            return
        end
        if pf.zero_time>b.t1
             msgbox( 'The first time has to be bigger than the first time point',...
              'Error','error');
            return
        end
        if  pf.t1<=pf.t2 && b.t2<=pf.final_time && b.t1>=pf.zero_time
            close(b.fig)
            msgbox('The magnetic field and time interval were loaded')
            if isempty(pf.data_dir)==1
                pf.data_dir= fullfile(b.location_b,b.name_b);
            end
            pf.data=data;
            set(pf_gui.str_boot_no,'visible','on','string',pf.data_dir);
            return
        end
    else
        msgbox( 'Please insert the time interval','Error','error');
        return
    end
end
return

function[]=cancel_data()
global b;
global pf_gui;
global pf;
global data;
data=pf.data;
if isempty(pf.data_dir)==1
    set(pf_gui.str_boot_no,'string','','visible','off');
    set(b.str_b,'string','');
else
    set(pf_gui.str_boot_no,'string',num2str(pf.data_dir));
    set(b.str_b,'string',num2str(pf.data_dir));
end
close
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% function for saving the data present in the GUI

function[]=save_parameters()
global pf_gui;
global pf;
if isempty(pf.results_dir)==0
    cd(pf.results_dir);
end
choice_par = questdlg('Please choose the input parameters you want to save', ...
'Saving parameters!', ...
'only input directories','input directories and data','only input directories');
% Handle response
switch choice_par
case 'only input directories'
    [save_name,save_path] = uiputfile('*.mat','Save these set of directories');
    if save_path == 0 
        return
    end
    pf_aux.parameters=1;
    pf_aux.subjects_dir=pf.subjects_dir;
    pf_aux.results_dir=pf.results_dir;      %%% path of the result directory     
    pf_aux.subject=pf.subject;              %%% name of the experiment  
    pf_aux.g_matrix_dir=pf.g_matrix_dir;    %%% lead field matrix directory
    pf_aux.vertices_dir=pf.vertices_dir;    %%% vertices matrix directory
    pf_aux.orient_dir=pf.orient_dir;        %%% orientations matrix directory (no mandatory)  
    pf_aux.cm_evol=str2num(get(pf_gui.str_cm_evol,'string'));              %%% radius [cm] used for computing the evolution matrix. It is not necessary if the
    pf_aux.evol_dir=pf.evol_dir;            %%% evolution matrix directory
    pf_aux.ssp_dir=pf.ssp_dir;              %%% ssp matrix directory (not mandatory)
    pf_aux.cov_matrix_dir=pf.cov_matrix_dir;%%% noise covariance matrix directory (not mandatory)
    pf_aux.sigma=pf.sigma;              %%% mean variance of noise during the pre-stimulus
    if strcmp(get(pf_gui.str_sigma_par,'string'),'')==1
        pf_aux.sigma_par=1;
    else
        pf_aux.sigma_par=(str2num(get(pf_gui.str_sigma_par,'string')));
    end
    pf_aux.np=floor(str2num(get(pf_gui.str_np,'string')));      %%% number of particles
    pf_aux.zero_time=pf.zero_time;                              %%% first time point [s] of the data
    pf_aux.final_time=pf.final_time;                            %%% final time point [s] of the data
    pf_aux.freq=pf.freq;                                        %%% frequency
    pf_aux.time_interval=pf.time_interval;                      %%% vector containig the time interval
    pf_aux.dataset_number=pf.dataset_number;                    %%% number of the dataset that will be analyzed
    pf_aux.t1=pf.t1;                                            %%% the first time point of the analysis
    pf_aux.t2=pf.t2;                                            %%% the last time point of the analysis
    pf_aux.data_dir =pf.data_dir;                               %%% data matrix directory
    pf_aux.p1=pf.p1;                                            %%% first pre-stimulus time point
    pf_aux.p2=pf.p2;  
    pf_aux.autosave=pf.autosave;
    if (get(pf_gui.flag_autosave,'Value')) == (get(pf_gui.flag_autosave,'Max'))
        pf_aux.autosave =1;
    else
        pf_aux.autosave =0;
    end
    %saves the gui data

    save(fullfile(save_path,save_name),'pf_aux','-mat');
case 'input directories and data'
    [save_name,save_path] = uiputfile('*.mat','Save these set of paramters');
    if save_path == 0 
        return
    end
    pf.parameters=1;
    pf.subject = get(pf_gui.str_subject,'string');
    pf.cm_evol = str2num(get(pf_gui.str_cm_evol,'string'));
    pf.np=floor(str2num(get(pf_gui.str_np,'string')));
    if strcmp(get(pf_gui.str_sigma_par,'string'),'')==1
        pf.sigma_par=1;
    else
        pf.sigma_par=(str2num(get(pf_gui.str_sigma_par,'string')));
    end
    if (get(pf_gui.flag_autosave,'Value')) == (get(pf_gui.flag_autosave,'Max'))
        pf.autosave =1;
    else
        pf.autosave =0;
    end
    %saves the gui data
    save(fullfile(save_path,save_name),'pf','-mat');
end
return;



%%% function for running the particle filter 
function[]=GO()
global pf_gui;
global pf;
global b;
global string_ris;
(get(pf_gui.flag_autosave,'Value'));
(get(pf_gui.flag_autosave,'Max'));
if (get(pf_gui.flag_autosave,'Value')) == (get(pf_gui.flag_autosave,'Max'))
    pf.autosave =1;
else
    pf.autosave =0;
end
n_check=0;
[n_check,err_msg]=hades_check;
if n_check==9   
    string_ris = fullfile(pf.results_dir,pf.subject);
    if strcmp(string_ris(end-3:end),'.mat')==0
       string_ris = strcat(string_ris,'.mat');
    end
    if pf.autosave==0
        hades_pf
    else
        if exist(string_ris)
            choice_save_ris = questdlg(strcat(['The results are going to be saved on the already existing file ',string_ris,'. Do you want to overwrite it?']), ...
            'Saving results...', 'yes','no','yes');
            % Handle response
            switch choice_save_ris
            case 'yes'
                hades_pf
            case 'no'
                [save_name_ris_new,save_path_ris_new] = uiputfile('*.mat','Choose the file where save the results');
                if save_path_ris_new == 0 
                    return
                end
                string_ris=fullfile(save_path_ris_new,save_name_ris_new);
                hades_pf
            end
        else
            hades_pf
        end
    end          
else
    err_msg;
end
return;
%%% function for running the post preocessing
function[]=post_processing()
global pf;
global ris;
global fwd;
global gui_vis;
global pf_gui;
if isempty(pf.results_dir)==0
    if pf.results_dir~=0
        cd(pf.results_dir)
    end
end
if isfield(ris,'subject')==1
else
    ris.subject='';
end
if isempty(findobj('Name',strcat(['Results ',ris.subject]),'-and','-depth',6))==0  
    figure(gui_vis.fig);
    else 
    [subject_name,subject_location]=uigetfile('*.mat', 'Choose the subject/experiment file');
    if  subject_name==0
    else
        loaded_ris=load(fullfile(subject_location,subject_name));
       
        if isfield(loaded_ris,'ris')==0 
            msgbox('The loaded .m file contains unexpected data','Error','error');
            return
        else
            ris=loaded_ris.ris;
            set(pf_gui.str_subject,'string',ris.subject,'visible','on');
            set(pf_gui.str_boot_no,'string',ris.data,'visible','on');
            set(pf_gui.str_ssp,'string',ris.ssp,'visible','on');
            set(pf_gui.str_cov_matrix,'string',ris.cov_matrix,'visible','on');
            set(pf_gui.str_np,'string',ris.particles);
            set(pf_gui.str_sigma_par,'string',ris.sigma_par);
            if isempty(pf.vertices)==1 
                fwd_yes
                while (fwd.check==0) && isempty(findobj('Name','Loading the lead field','-and','-depth',6))==0
                    fwd_yes
                end
                if fwd.check==1
                    close
                    hades_visualization
                end                       
            else
                hades_visualization
            end
        end
    end
end
return
%%% close all the open windows of the GUI
function[]=exit()
global ris;
       
obj_b=[];
obj_b=findobj('Name','Loading the magnetic field','-and','-depth',6);
if isempty(obj_b)==1
else
    close (obj_b);
end
obj_fwd=[];
obj_fwd=findobj('Name','Loading the lead field','-and','-depth',6);
if isempty(obj_fwd)==1
else
    close (obj_fwd);
end
obj_ris=[];
if isfield(ris,'subject')==1
    obj_fwd=findobj('Name',strcat(['Results ',ris.subject]),'-and','-depth',6);
    if isempty(obj_fwd)==1
    else
        close (obj_fwd);
    end
end
obj_err=[];
obj_err=findobj('Name','Error','-and','-depth',6);
if isempty(obj_err)==1
else
    close (obj_err);
end
close
return

