function [n,msg]=hades_check

%%% HADES_CHECK checks if all the input data required by the Particle Filter have
%%% been loaded

%   Copyright (C) 2010, Cristina Campi, Annalisa Pascarella, Michele Piana, Alberto
%   Sorrentino

n=0;
msg='';
global pf_gui;
global pf;


%%%%check if the autosave option is selected
if (get(pf_gui.flag_autosave,'Value')) == (get(pf_gui.flag_autosave,'Max'))
    pf.autosave =1;
end

%%%%check on the subjects dir
if isempty(pf.subjects_dir)==1 || strcmp(pf.subjects_dir,'0')==1
    msg=msgbox('You have to set the subjects directory','Error','error');
    return
else
    n=n+1;
end

%%%%set the results directory
if isempty(pf.results_dir)==1 || strcmp(pf.results_dir,'0')==1
        pf.results_dir=pf.subjects_dir;
end

%%%%check on the subject
if strcmp(get(pf_gui.str_subject,'string'),' ')==1 || strcmp(get(pf_gui.str_subject,'string'),'')==1 || isempty(get(pf_gui.str_subject,'string'))==1
   msg=msgbox('You have to set a subject/experiment','Error','error');
   return
else
    pf.subject = get(pf_gui.str_subject,'string');
    n=n+1;
end

%%%%set the lead field matrix
if isempty(pf.g_matrix)==0
    n=n+1;
else
    msg=msgbox( 'Load the lead field matrix, please','Error','error');
    return
end

%%%%set the noise covariance matrix
if isempty(pf.cov_matrix)==0
    n=n+1;
else
    if isempty(pf.g_matrix)==1
         msg=msgbox( 'Load the lead field matrix, please','Error','error');
         return
    else
        ns=size(pf.g_matrix,1);
        pf.cov_matrix = eye(ns);
        n=n+1;
    end
end

%%%%set the ssp matrix
if isempty(pf.ssp)==0
    n=n+1;
else
   if isempty(pf.g_matrix)==1
      msg=msgbox( 'Load the lead field matrix, please','Error','error');
      return
   else
        ns=size(pf.g_matrix,1);
        pf.ssp = eye(ns);
        n=n+1;
    end
end

%%%%set the pw matrix 
if isempty(pf.matrix_pw)==0
    n=n+1;
else
   if isempty(pf.g_matrix)==1
      msg=msgbox( 'Load the lead field matrix, please','Error','error');
      return
   else
        ns=size(pf.g_matrix,1);
        pf.matrix_pw = eye(ns);
        n=n+1;
    end
end

%%%%set the number of particles
if strcmp(get(pf_gui.str_np,'string'),'')==1
   msg=msgbox( 'You have to set the number of particles',...
        'Error','error');
    return
else
    pf.np=floor(str2num(get(pf_gui.str_np,'string')));
    n=n+1;
end

%%%% set the fitting parameter
if strcmp(get(pf_gui.str_sigma_par,'string'),'')==1
    pf.sigma_par=1;
else
    pf.sigma_par=(str2num(get(pf_gui.str_sigma_par,'string')));
end

%%%%check if the magnetic field has been loaded
if isempty(pf.data)==1 && isempty(pf.data_dir)
   msg=msgbox( 'You have to load a magnetic field',...
        'Error','error');
    return
elseif isempty(pf.data)==0
    [t_tot,ns_b] = size(pf.data);
    if isempty(pf.g_matrix)==1
      msg=msgbox( 'Load the lead field matrix, please','Error','error');
      return
    else
        if t_tot<pf.t2 && size(pf.g_matrix,1)==ns_b
            msg=msgbox('The choosen time interval does not fit the data','Error','error');
            return
        elseif ns_b~=size(pf.g_matrix,1) && pf.t2<=t_tot
            msg=msgbox('The number of sensors in the data does not correspond to the sensors number in the lead field matrix',...
            'Error','error');
            return
        elseif  pf.t2<=t_tot && size(pf.g_matrix,1)==ns_b
            n=n+1;
        end
    end
elseif isempty(pf.data_dir)==0
    name=pf.data_dir;
    ext = name(end-3:end);
    if strcmp(ext,'.dat')==1
        pf.data=load(name);
        n=n+1;
    elseif strcmp(ext,'.fif')==1
         if isempty(pf.dataset_number)==1
            msgbox( 'You have to set the dataset','Error', 'error');
            return
        else
            if isempty(pf.dataset_number)==1
               msgbox( 'You have to load the a correct dataset!','Error','error')
               return
            else
                aux_b_fif=fiff_read_evoked(name,pf.dataset_number);
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
                n=n+1;
                pf.zero_time=double(aux_b_fif.evoked.first)*10^-3;
                pf.final_time=double(aux_b_fif.evoked.last)*10^-3;
                pf.freq = aux_b_fif.info.sfreq;
                pf.time_interval=pf.zero_time:1/pf.freq:pf.final_time;
            end
        end
    end
end

%%%% set the evolution radius
cm_evol = str2num(get(pf_gui.str_cm_evol,'string'));
if isempty(pf.cm_evol)==1
    pf.cm_evol=-1000;
end

%%%%compute the evolution matrix
if (get(pf_gui.flag_evol,'Value')) == (get(pf_gui.flag_evol,'Max'))
    if strcmp(get(pf_gui.str_cm_evol,'string'),'')==1
            msg=msgbox( 'You have to set the radius for the neighbours', 'Error','error');
            return
        elseif isempty(cm_evol)==1
            msg=msgbox( 'You have to set a radius bigger than zero', 'Error','error');
            return
        else
            if isempty(pf.vertices)==1
                msg=msgbox('You have to load the vertices file','Error','error');
                return
            else
                if  pf.cm_evol~=cm_evol
                    pf.cm_evol=cm_evol;
                    pf.evol = hades_compute_neighbours(cm_evol,pf.vertices,pf.subjects_dir,pf.subject);
                    n=n+1;
                else
                    if isempty(pf.evol)==0
                        n=n+1;
                    else 
                        msg=msgbox( 'You have to load the neighbours', 'Error','error');
                        return        
                    end
                end
            end
    end
else
    if isempty(pf.evol)==0
        n=n+1;
    else 
        msg=msgbox( 'You have to load the neighbours', 'Error','error');
        return        
    end
end



