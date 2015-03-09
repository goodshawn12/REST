function[]=hades_load_parameters()
%%% HADES_LOAD_PARAMETERS loads a .mat file containing settings previously
%%% saved.
%   Copyright (C) 2010, Cristina Campi, Annalisa Pascarella, Michele Piana, Alberto
%   Sorrentino
global pf;
global pf_aux;
global fwd;
global pf_gui;
global b;
global data;
[par_name,par_location]=uigetfile('*.mat', 'Choose the .mat file containing the parameters');
if par_name==0;
    return
else
    loaded_par=load(fullfile(par_location,par_name)); %%% load the .mat file containing the parameters 
    %%% associated with an experiment
    if isfield(loaded_par,'pf')==0 && isfield(loaded_par,'pf_aux')==0
        msgbox('The loaded .m file contains unexpected data','Error','error');
        return
    else
        msgbox('Loading the parameters');
        if isfield(loaded_par,'pf')==1
            pf=loaded_par.pf;
        elseif isfield(loaded_par,'pf_aux')==1
            pf=loaded_par.pf_aux;
            pf.parameters=1;
            pf.g_matrix =[];        %%% lead field matrix
            pf.vertices = [];       %%% vertices matrix
            pf.orient = [];         %%% orientations matrix (no mandatory)   
            pf.evol = [];           %%% evolution matrix
            pf.ssp = [];            %%% ssp matrix (not mandatory)
            pf.cov_matrix=[];       %%% noise covariance matrix (not mandatory)
            pf.data = [];           %%% data matrix 
        end
        if pf.autosave ==1 
            set(pf_gui.flag_autosave,'Value',1);
        end
        if  isempty(pf.subjects_dir)==0
            set(pf_gui.str_subjects_dir,'visible','on','string',pf.subjects_dir);
        end
        if isempty(pf.results_dir)==0
            set(pf_gui.str_results_dir,'visible','on','string',pf.results_dir);
        end
        if isempty(pf.subject)==0
            set(pf_gui.str_subject,'visible','on','string',pf.subject);
        end
        %%% set the lead field
        if isempty(pf.g_matrix)==1
            if isempty(pf.g_matrix_dir)==0
                if exist(pf.g_matrix_dir)==0
                    msgbox(strcat(['The file ',pf.g_matrix_dir,' does not exist']),'Error','error');
                else
                    ext=pf.g_matrix_dir(end-3:end);
                    if strcmp(ext,'.fif')==1
                      fwd.str_fwd_fif = [];
                        set(fwd.str_fwd_fif,'visible','on','string',pf.g_matrix_dir);
                        set(pf_gui.str_fwd_yes,'visible','on','string',pf.g_matrix_dir);
                        if exist('mne_read_forward_solution.m')==0
                            msgbox('The file mne_read_forward_solution.m is not in the path','Error','error');
                        else
                            fwd_aux = mne_read_forward_solution(pf.g_matrix_dir);
                            pf.g_matrix=fwd_aux.sol.data;
                            pf.vertices=fwd_aux.source_rr;
                            if size(fwd_aux.source_rr,1)==size(fwd_aux.source_nn,1)
                                pf.orient=fwd_aux.source_nn;
                            end
                        end
                    elseif strcmp(ext,'.dat')==1
                        set(fwd.str_lead_field,'visible','on','string',pf.g_matrix_dir);
                        set(pf_gui.str_fwd_yes,'visible','on','string',pf.g_matrix_dir);
                        pf.g_matrix=load(pf.g_matrix_dir);
                    elseif strcmp(ext,'.mat')==1
                        set(fwd.str_lead_field,'visible','on','string',pf.g_matrix_dir);
                        set(pf_gui.str_fwd_yes,'visible','on','string',pf.g_matrix_dir);
                        field_mat=load(pf.g_matrix_dir);   
                        pf.g_matrix=getfield(field_mat,cell2mat(fieldnames(field_mat)));
                    end  
                end
            end
        else
          fwd.str_fwd_fif = [];
            set(fwd.str_fwd_fif,'visible','on','string',pf.g_matrix_dir);
            set(pf_gui.str_fwd_yes,'visible','on','string',pf.g_matrix_dir);
        end
        if isempty(pf.vertices)==1
            if isempty(pf.vertices_dir)==0
                set(fwd.str_vert,'visible','on','string',pf.vertices_dir);
                if exist(pf.vertices_dir)==0
                    msgbox(strcat(['The file ',pf.vertices_dir,' does not exist']),'Error','error');
                else
                    ext =  pf.vertices_dir(end-3:end);
                    if strcmp(ext,'.dat')==1
                        pf.vertices=load(pf.vertices_dir);
                    elseif strcmp(ext,'.mat')==1
                        vert_mat=load(pf.vertices_dir);        
                        pf.vertices=getfield(vert_mat,cell2mat(fieldnames(vert_mat)));
                    end  
                end
            end
        end
        if isempty(pf.orient)==1
            if isempty(pf.orient_dir)==0
                set(fwd.str_orient,'visible','on','string',pf.orient_dir);
                if exist(pf.orient_dir)==0
                    msgbox(strcat(['The file ',pf.orient_dir,' does not exist']),'Error','error');
                else
                    ext =  pf.orient_dir(end-3:end);
                    if strcmp(ext,'.dat')==1
                        pf.orient=load(pf.orient_dir);
                    elseif strcmp(ext,'.mat')==1
                        orient_mat=load(pf.orient_dir);        
                        pf.orient=getfield(vert_mat,cell2mat(fieldnames(orient_mat)));
                    end  
                end
            end
        end
        if isempty(pf.evol)==1
            if isempty(pf.evol_dir)==0
                set(pf_gui.str_evol,'visible','on','string',pf.evol_dir);
                if exist(pf.evol_dir)==0
                    msgbox(strcat(['The file ',pf.evol_dir,' does not exist']),'Error','error');
                else
                   ext = pf.evol_dir(end-3:end);
                   if strcmp(ext,'.dat')==1
                       pf.evol=load(pf.evol_dir);
                   elseif strcmp(ext,'.mat')==1
                       evol_mat=load(pf.evol_dir);
                       pf.evol=getfield(evol_mat,cell2mat(fieldnames(evol_mat)));
                   end
                end
            end
        else
            set(pf_gui.str_evol,'visible','on','string',pf.evol_dir);
        end
        if isempty(pf.cm_evol)==0
            if strcmp(pf.cm_evol,'')==0
                set(pf_gui.str_cm_evol,'visible','on','string',pf.cm_evol);
                set(pf_gui.flag_evol,'Value',1);
            end
        end
        if isempty(pf.ssp)==1
            if isempty(pf.ssp_dir)==0
                set(pf_gui.str_ssp,'visible','on','string',pf.ssp_dir);
                if exist(pf.ssp_dir)==0
                    msgbox(strcat(['The file ',pf.ssp_dir,' does not exist']),'Error','error');
                else
                ext = pf.ssp_dir(end-3:end);
                if strcmp(ext,'.dat')==1
                    pf.ssp=load(pf.ssp_dir);
                elseif strcmp(ext,'.mat')==1
                    ssp_mat=load(pf.ssp_dir);        
                    pf.ssp=getfield(ssp_mat,cell2mat(fieldnames(ssp_mat)));
                end 
                end
            end
        else
            set(pf_gui.str_ssp,'visible','on','string',pf.ssp_dir);
        end
        if isempty(pf.cov_matrix)==1
            if isempty(pf.cov_matrix_dir)==0
                set(pf_gui.str_cov_matrix,'visible','on','string',pf.cov_matrix);
                if exist(pf.cov_matrix_dir)==0
                    msgbox(strcat(['The file ',pf.cov_matrix_dir,' does not exist']),'Error','error');
                else
                    ext = pf.cov_matrix_dir(end-3:end);
                    if strcmp(ext,'.dat')==1
                        pf.cov_matrix=load(pf.cov_matrix_dir);
                    elseif strcmp(ext,'.mat')==1
                        cov_matrix_mat=load(pf.cov_matrix_dir);        
                        pf.cov_matrix=getfield(cov_matrix_mat,cell2mat(fieldnames(cov_matrix_mat)));
                    end 
                end
            end
        else
            set(pf_gui.str_cov_matrix,'visible','on','string',pf.cov_matrix);
        end
        if isempty(pf.np)==0
            set(pf_gui.str_np,'visible','on','string',pf.np);
        end
        if isempty(pf.data)==1
            if isempty(pf.data_dir)==0
                b.str_b=[];
                set(b.str_b,'string',num2str(pf.data_dir));
                set(pf_gui.str_boot_no,'visible','on','string',num2str(pf.data_dir));
                if exist(pf.data_dir)==0
                   msgbox(strcat(['The file ',pf.data_dir,' does not exist']),'Error','error');
                else
                    ext =  pf.data_dir(end-3:end);
                    if strcmp(ext,'.dat')==1 || strcmp(ext,'.mat')==1
                        if strcmp(ext,'.dat')==1
                            pf.data=load(pf.data_dir);
                            data=pf.data;
                        elseif strcmp(ext,'.mat')==1
                             data_mat=load(pf.data_dir);      
                             pf.data=getfield(data_mat,cell2mat(fieldnames(data_mat)));
                             data=pf.data;   
                        end
                        b.label_times_dat = [];
                        set(b.label_times_dat,'visible','on'); 
                        if isempty(pf.zero_time)==0
                          b.zero_time = [];
                            set(b.zero_time,'visible','on','string',pf.zero_time);
                        end
                        if isempty(pf.final_time)==0
                            set(b.final_time,'visible','on','string',pf.final_time);
                        end
                        if isempty(pf.freq)==0
                          b.freq = [];
                          set(b.freq,'visible','on','string',pf.freq);
                        end
                        b.dataset_number = [];
                        set(b.dataset_number,'visible','off');
                        b.label_ptimes_fif = [];
                        set(b.label_ptimes_fif,'visible','off');
                   elseif strcmp(ext,'.fif')==1
                       b.label_ptimes_fif=[];
                        set(b.label_ptimes_fif,'visible','on');
                        b.zero_time=[];
                        set(b.zero_time,'visible','off'); 
                        b.final_time=[];
                        set(b.final_time,'visible','off');
                        b.freq=[];
                        set(b.freq,'visible','off');  
                        if isempty(pf.dataset_number)==0
                            b.dataset_number=[];
                            set(b.dataset_number,'visible','on','string',pf.dataset_number);
                            b.cmd_dataset=[];
                            set(b.cmd_dataset,'enable','off');
                            aux_b_fif=fiff_read_evoked(pf.data_dir,pf.dataset_number);
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
                            pf.data=aux_b_fif.evoked.epochs(1:t,:)';
                            data=pf.data;
                        end
                    end
                end
            end
        else
            set(pf_gui.str_boot_no,'visible','on','string',pf.data_dir);
            data=pf.data;
        end
        if isempty(pf.sigma)==0
            b.sigma=[];
            set(b.sigma,'visible','on','string',pf.sigma);
        end
        if isempty(pf.sigma_par)==0
           set(pf_gui.str_sigma_par,'string',num2str(pf.sigma_par));
        end
        if isempty(pf.p1)==0
            b.p1=[];
            set(b.p1,'visible','on','string',num2str(pf.p1));
        end
        if isempty(pf.p2)==0
            b.p2=[];
            set(b.p2,'visible','on','string',num2str(pf.p2));
        end
        if isempty(pf.time_interval)==0
            if isempty(pf.t1)==0
                b.str_t1=[];
                set(b.str_t1,'visible','on','string',pf.time_interval(pf.t1));
            end
            if isempty(pf.t2)==0
                b.str_t2=[];
                set(b.str_t2,'visible','on','string',pf.time_interval(pf.t2));
            end 
            close
        end
    end
end
return;
