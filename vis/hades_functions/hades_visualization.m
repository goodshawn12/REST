%   HADES_VISUALIZATION opens the Results window and manages the functions
%   for visualizing the results.
% 
%   Copyright (C) 2010, Cristina Campi, Annalisa Pascarella, Michele Piana, Alberto
%   Sorrentino

function[]= hades_visualization(action,varargin)
if nargin<1,
   action='gui_visualization';
end;

feval(action,varargin{:});

return;

%%% this window contains the command buttons for saving the results in .stc
%%% and .w file, for running again the clustering procedure, for showing
%%%  the model selection, the recostructed dipoles and the clusters
function[]= gui_visualization()
global gui_vis;
global vis;
global pf;
global ris;
vis.mne=[];
vis.freesurfer=[];
scrsz = get(0,'ScreenSize');
gui_vis.fig=figure('HandleVisibility','Callback','Menubar','none',...
    'Name',strcat(['Results ',ris.subject]), 'NumberTitle','off', ...
        'Visible','on', 'BackingStore','off');
pos=get(gui_vis.fig,'pos');
    pos(1,3)=500;
    pos(1,4)=300;
    pos(1,1)=(scrsz(1,3)-pos(1,3))/2;
    pos(1,2)=(scrsz(1,4)-pos(1,4))/2;
set(gui_vis.fig,'pos',pos);
gui_vis.panel1=uipanel('parent',gui_vis.fig,'Title','Cluster the results','units','pixel',...
    'Position',[10 240 215 50],'FontWeight','bold');
%%% command button for runnig again the clustering procedure
gui_vis.cmd_clustering=uicontrol(gui_vis.panel1,'style', 'pushbutton','units', 'pixel',...
        'position',[5 5 200 30 ],...
        'enable','on',...
        'callback', 'hades_visualization(''clustering_show'');' ,...
        'string', 'cluster');
% %%% string for setting the number of components
gui_vis.label_comp = uicontrol(gui_vis.fig,'style', 'text',...
        'units', 'pixel',...
        'visible','off',...
        'string', 'clustering parameters',...
        'position', [250 253 150 20]);
gui_vis.popup_comp = uicontrol(gui_vis.fig,'style','popupmenu',...
          'String',{'location','orientation and location'},...
          'units', 'pixel',...
          'enable','on',...
          'visible','off',...
          'position', [410 225 90 50]);
% %%% string for setting the initial number of clusters       
gui_vis.label_NN = uicontrol(gui_vis.fig,'style', 'text',...
        'units', 'pixel',...
        'visible','off',...
        'string', 'number of clusters',...
        'position', [250 225 150 20]);      
gui_vis.str_NN = uicontrol(gui_vis.fig,'style','edit',...
          'String','',...
          'visible','off',...
          'units', 'pixel',...
          'enable','on',...
          'position', [410 225 50 20]);    
gui_vis.cmd_cluster=uicontrol(gui_vis.fig,'style','pushbutton',...
          'String','compute clusters',...
          'visible','off',...
          'units', 'pixel',...
          'enable','on',...
          'callback', 'hades_visualization(''clustering_ris'');',...
          'position', [250 200 150 20]);  
gui_vis.cmd_save_cluster=uicontrol(gui_vis.fig,'style','pushbutton',...
          'String','save',...
          'visible','off',...
          'units', 'pixel',...
          'enable','on',...
          'callback', 'hades_visualization(''clustering_save'');',...
          'position', [410 200 50 20]); 
%%%% for showing the results       
gui_vis.label_object = uicontrol(gui_vis.fig,'style', 'text',...
        'units', 'pixel',...
        'visible','off',...
        'string', 'select object ',...
        'position', [250 253 100 20]);
gui_vis.popup_object= uicontrol(gui_vis.fig,'style','popupmenu',...
          'String',{'estimated dipoles','clustered dipoles','estimated sources'},...
          'units', 'pixel',...
          'enable','on',...
          'visible','off',...
          'callback', 'hades_visualization(''set_object'');',...
          'position', [360 225 140 50]);      
gui_vis.label_vis = uicontrol(gui_vis.fig,'style', 'text',...
        'units', 'pixel',...
        'visible','off',...
        'string', '',...
        'position', [250 225 100 20]);     
gui_vis.str_vis = uicontrol(gui_vis.fig,'style', 'edit',...
        'units', 'pixel',...
        'visible','off',...
        'string', '',...
        'enable','on',...
        'position', [360 225 140 20]);       
gui_vis.label_split = uicontrol(gui_vis.fig,'style', 'text',...
        'units', 'pixel',...
        'visible','off',...
        'string', 'time split',...
        'position', [420 225 60 20]);  
gui_vis.flag_split = uicontrol(gui_vis.fig,'style', 'checkbox',...
        'units', 'pixel',...
        'visible','off',...
        'position', [480 225 30 20],...
        'callback','hades_visualization(''set_clusters_number'');');   
gui_vis.cmd_vis=uicontrol(gui_vis.fig,'style','pushbutton',...
          'String','',...
          'visible','off',...
          'units', 'pixel',...
          'enable','on',...
          'callback', 'hades_visualization(''vis_clusters_or_seeds'');',...
          'position', [250 200 150 20]); 
    
gui_vis.panel2=uipanel('parent',gui_vis.fig,'Title','Show the results','units','pixel',...
    'Position',[10 145 215 90],'FontWeight','bold');
%%% command button for showing the recosntructed dipoles or the clusters    
gui_vis.cmd_only_vertex=uicontrol(gui_vis.panel2,'style', 'pushbutton','units', 'pixel',...
        'position',[5 40 200 30],...
        'enable','on',...
        'callback', 'hades_visualization(''results_show'');' ,...
        'string', 'show the results');
%%% command button for showing the model selection
gui_vis.cmd_mod_sel=uicontrol(gui_vis.panel2,'style', 'pushbutton','units', 'pixel',...
        'position', [5 5 200 30],...
        'enable','on',...
        'callback', 'hades_visualization(''vis_mod_sel'');' ,...
        'string', 'show the model selection');
    
    
gui_vis.panel3=uipanel('parent',gui_vis.fig,'Title','Export the results','units','pixel',...
    'Position',[10 15 215 125],'FontWeight','bold');
%%% command button for saving .stc files
gui_vis.cmd_mne=uicontrol(gui_vis.panel3,'style', 'pushbutton','units', 'pixel',...
        'position', [5 75 200 30],...
        'enable','on',...
        'callback', 'hades_visualization(''vis_mne'');' ,...
        'string', 'write .stc files');
%%% command button for saving .w files    
gui_vis.cmd_freesurfer=uicontrol(gui_vis.panel3,'style', 'pushbutton','units', 'pixel',...
        'position', [5 40  200 30],...
        'enable','on',...
        'callback', 'hades_visualization(''vis_freesurfer'');' ,...
        'string', 'write .w files');
%%% command button for saving BS files    
gui_vis.cmd_freesurfer=uicontrol(gui_vis.panel3,'style', 'pushbutton','units', 'pixel',...
        'position', [5 5  200 30],...
        'enable','on',...
        'callback', 'hades_visualization(''vis_bs'');' ,...
        'string', 'write BrainStorm files');
%%% command button for closing this window    
gui_vis.cmd_exit=uicontrol(gui_vis.fig,'style', 'pushbutton',...
        'units', 'pixel',...
        'position', [350 20 150 30],...
        'enable','on',...
        'callback', 'close;',...
        'string', 'exit');  
return

%%% function for saving .stc files 
function[]=vis_mne
global fwd;
global ris;
global pf;
global gui_vis;
if isempty(ris.dip)==1
    msgbox('The PF algorithm did not reconstruct any source','Error','error');
    return
else
    if exist('mne_read_forward_solution.m')==0
        msgbox('The file mne_read_forward_solution.m is not in the path','Error','error');
           return
    else
        if isempty(fwd.location_fwd_fif)==1 && isempty(pf.g_matrix_dir)==1
            msgbox('You can write .stc file only if you have a lead field in .fif format. Please load it','Error','error');
            close(findobj('Name','Results','-and','-depth',6));
            return
        else
            if isempty(pf.g_matrix_dir)==0
                ext=pf.g_matrix_dir(end-3:end);
                if strcmp(ext,'.fif')==1
                    fwd_aux=mne_read_forward_solution(pf.g_matrix_dir);
                else
                     msgbox('You can write .stc file only if you have a lead field in .fif format. Please load it','Error','error');
                    close(findobj('Name','Results','-and','-depth',6));
                    return
                end
            else
                if fwd.location_fwd_fif==0
                    msgbox('You can write .stc file only if you have a lead field in .fif format. Please load it','Error','error');
                    close(findobj('Name','Results','-and','-depth',6));
                    return
                else
                    filename=fullfile(fwd.location_fwd_fif,fwd.name_fwd_fif);
                    if strcmp(filename(end-3:end),'.fif')==0 
                        msgbox('You can write .stc file only if you have a lead field in .fif format. Please load it','Error','error');
                        close(findobj('Name','Results','-and','-depth',6));
                        return
                    else
                        fwd_aux = mne_read_forward_solution(fullfile(fwd.location_fwd_fif,fwd.name_fwd_fif));
                    end
                end
            end
            vertices = [ fwd_aux.src(1,1).vertno'; fwd_aux.src(1,2).vertno'];   
            vert_coord=fwd_aux.source_rr;
            if isfield(ris,'results_dir')==0
                start_path=pwd;
                directory_name = uigetdir(start_path,'Select the directory where save the .stc files');
            else
                if exist(ris.results_dir)==0
                    start_path=pwd;
                    directory_name = uigetdir(start_path,'Select the directory where save the .stc files');
                else
                    directory_name = ris.results_dir;
                end
            end
            if directory_name==0
            else
               str_dip=get(gui_vis.str_vis,'string');
               time_dip=str2num(str_dip);
               if strcmp(str_dip,'')==0 && isempty(time_dip)==0
                    hades_write_stc(fwd_aux,vertices,vert_coord,ris.time_interval(ris.t1),double(1/ris.freq),ris.dip,ind_first,ind_fin,directory_name);
               else
                   hades_write_stc(fwd_aux,vertices,vert_coord,ris.time_interval(ris.t1),double(1/ris.freq),ris.dip,ris.t1,ris.t2,directory_name);
                 return
                end
            end
        end
    end
end

return



%%% function for saving .w files 
function[]=vis_freesurfer
global src;
global ris;
[src_aux.name,src_aux.location]=uigetfile('*-src.fif', 'Choose the <subject>-src.fif file');
if src_aux.name==0
else
     if exist('mne_read_source_spaces.m')==0
           msgbox('The file mne_read_source_spaces.m is not in the path','Error','error');
           return
     else
         if isempty(ris.dip)==1
             msgbox('The PF algorithm did not reconstruct any source','Error','error');
             return
         else
            src = mne_read_source_spaces(fullfile(src_aux.location,src_aux.name));
            if exist(ris.results_dir)==0
                start_path=pwd;
                directory_name = uigetdir(start_path,'Select the directory where save the .w files');
            else
                directory_name = ris.results_dir;
            end
            hades_write_w(src,ris.dip,'dip',directory_name);
         end
     end
end
return
function[]=vis_bs
global ris;
global gui_bs;
if isempty(ris.dip)==1
    msgbox('The PF algorithm did not reconstruct any source','Error','error');
    return
else

    scrsz = get(0,'ScreenSize');
    gui_bs.fig=figure('HandleVisibility','Callback','Menubar','none',...
    'Name','Export data to BrainStorm', 'NumberTitle','off', ...
        'Visible','on', 'BackingStore','off');
    pos=get(gui_bs.fig,'pos');
    pos(1,3)=500;
    pos(1,4)=200;
    pos(1,1)=(scrsz(1,3)-pos(1,3))/2;
    pos(1,2)=(scrsz(1,4)-pos(1,4))/2;
    set(gui_bs.fig,'pos',pos);
    gui_bs.label_DataFile = uicontrol(gui_bs.fig,'style', 'text',...
        'string', 'Data File',...
        'units', 'pixel',...
        'position', [10 160 150 30]);
    gui_bs.str_DataFile = uicontrol(gui_bs.fig,'style', 'edit',...
        'string', '',...
        'visible','on',...
        'enable','on',...
        'units', 'pixel',...
        'position', [170 160 320 30]);
    gui_bs.label_HeadModelFile = uicontrol(gui_bs.fig,'style', 'text',...
        'string', 'Head Model File',...
        'units', 'pixel',...
        'position', [10 130 150 30]);
    gui_bs.str_HeadModelFile = uicontrol(gui_bs.fig,'style', 'edit',...
        'string', '',...
        'visible','on',...
        'enable','on',...
        'units', 'pixel',...
        'position', [170 130 320 30]);
    gui_bs.label_SurfaceFile = uicontrol(gui_bs.fig,'style', 'text',...
        'string', 'Surface File',...
        'units', 'pixel',...
        'position', [10 100 150 30]);
    gui_bs.str_SurfaceFile = uicontrol(gui_bs.fig,'style', 'edit',...
        'string', '',...
        'visible','on',...
        'enable','on',...
        'units', 'pixel',...
        'position', [170 100 320 30]);
    %%% command button for closing this window    
    gui_bs.cmd_exit=uicontrol(gui_bs.fig,'style', 'pushbutton',...
        'units', 'pixel',...
        'position', [350 20 150 30],...
        'enable','on',...
        'callback', 'hades_write_BS',...
        'string', 'ok');  
end
return

function[]=clustering_show
global gui_vis;
set(gui_vis.label_comp,'visible','on');
set(gui_vis.popup_comp,'visible','on');
set(gui_vis.str_NN,'visible','on');
set(gui_vis.label_NN,'visible','on');
set(gui_vis.cmd_cluster,'visible','on');
set(gui_vis.cmd_save_cluster,'visible','on');
set(gui_vis.label_object,'visible','off');
set(gui_vis.popup_object,'visible','off');
set(gui_vis.cmd_vis,'visible','off');
set(gui_vis.str_vis,'visible','off');
set(gui_vis.label_vis,'visible','off');
set(gui_vis.label_split,'visible','off');
set(gui_vis.flag_split,'visible','off');
return
    
function[]=clustering_ris
global ris;
global gui_vis;
global pf;
global fwd;
ris.cluster=[];
ris.waveforms=[];
ris.waveforms_seed=[];
ris.seed=[];


if isempty(pf.orient)==1
   DIM = 3;  %%% if the dipoles orientation are not fixed
else
   DIM = 1;  %%% if the dipoles orientation are fixed
end

if isfield(pf, 'tol_pinv')==0
   pf.tol_pinv=max(size(pf.g_matrix)) * norm(pf.g_matrix) * eps(class(pf.g_matrix))
end

if isempty(ris.dip)==1
 msgbox('The PF algorithm did not reconstruct any source','Error','error');
 return
else
    if strcmp(get(gui_vis.str_NN,'string'),'')==1
        msgbox('You forgot to set the number of clusters!','Error','error')
        return
    else
        NN=str2num(get(gui_vis.str_NN,'string'));
        if isempty(NN)==1
            msgbox('Set the number of clusters','Error','error');
            return
        else
            if NN>size(ris.dip,1)
                msgbox(strcat(['The PF reconstructed ', num2str(size(ris.dip,1)),' dipoles. The chosen number of clusters is too large']),'Error','error');
                return
            elseif NN>size(ris.dip,1)/3 && size(ris.dip,1)>= NN
                 msgbox(strcat(['The chosen number of clusters could be too large']),'Caution','warn');
            end
            
            all_different=0;            
            i_wave = 1;
            NN0 = NN;
            dipoli = ris.dip;
            dip=ris.dip;
            %%% normalization of the values of the dipole moments 
            for i = 1:size(dipoli,1)
                dipoli(i,4:6) = dipoli(i,4:6)/norm(dipoli(i,4:6))/10;
            end
            if exist('kmeans')==0 && exist('ranksum')==2 
            %%% check if the employed functions are in the path
                msgbox('The file kmeams.m is not in the path','Error','error');
                return
            elseif exist('ranksum')==0 && exist('kmeans')==2
                msgbox('The file ranksum.m is not in the path','Error','error');
                return
            
            elseif exist('kmeans')==2 && exist('ranksum')==2 
                while all_different == 0  && NN > 0
                    hades_script_cluster;
                    hades_statistical_test_cluster;
                    NN = NN-1;
                end
                NN = NN+1;                
                clear waveforms seed;
                waveforms=zeros(size(struct_dipoli_cluster,2),size(ris.time_interval,2));

                seed=zeros(size(struct_dipoli_cluster,2),3);
                for k = 1:size(struct_dipoli_cluster,2)
                    for j = 1:size(struct_dipoli_cluster(k).dipoli_cluster,1)
                        waveforms(k,struct_dipoli_cluster(k).dipoli_cluster(j,end-1)) =  norm(struct_dipoli_cluster(k).dipoli_cluster(j,4:end-2));
                        seed(k,1:3)=struct_dipoli_cluster(k).centroide;
                    end
                end
                %%% store in the struct ris the values of the seed positions and waveforms
                ris.seed=seed;
                ris.waveforms = waveforms; 
                ris.cluster=struct_dipoli_cluster;
                ris.NoC=NN0;
                ris.waveforms_seed=zeros(size(ris.seed,1),size(ris.time_interval,2));
                for t=1:size(ris.seed,1)
                    selected_seed=ris.seed(t,1:3);
                    diff=(pf.vertices-repmat(selected_seed,size(pf.vertices,1),1)).^2;
                    somma = sum(diff');
                    somma = somma';
                    dist = sqrt(somma);
                    [dist_riord,ind_riord] = sort(dist,'ascend');
                    seed_ok(t) = ind_riord(1)%pf.vertices(ind_riord(1),1:3);
                end
                
                for t=ris.t1:ris.t2
                    vector_ind_t=zeros(1,size(ris.seed,1));
                    G = []; temp_mdd = [];
                    for n=1:size(ris.seed,1)
                        [val_t, ind_t]=find(ris.cluster(n).dipoli_cluster(:,end-1)==t);
                        if isempty(ind_t)==1
                        else
                            vector_ind_t(n)=ind_t(1);
                        end
                    end
                    for  n=1:size(ris.seed,1)
                        if vector_ind_t(n)~=0
                           G = [G pf.g_matrix(:,DIM*(seed_ok(n)-1)+1:DIM*seed_ok(n))];   
                        end
                    end
                    if isempty(G)==0

                        %% if orientations are not constrained to the
                        %% cortical normal direction, the lead fields of
                        %% individual dipoles are combined after removing
                        %% the "singular direction" to avoid source
                        %% waveform artefacts
                        
                        if DIM > 1
                            G_new = [];
                            for i = 1:size(G,2)/3
                                [u_base(:,:,i) sing_val(:,:,i) v_base(:,:,i)] = svd(G(:, 3*(i-1)+1:3*i));
%                                 sing_val
%                                 break;
                                G_new = [G_new reshape(u_base(1:size(G,1),1:size(G,1),i),size(G,1),size(G,1))*reshape(sing_val(1:size(G,1),1:2,i),size(G,1),2)*reshape(v_base(1:2,1:3,i),2,3)];
                            end
                            
                            size(G_new);
                            temp_mdd = pinv(G_new)*ris.data(t,:)';                                                            
                            
                        else                       
                            pf.tol_pinv_used = pf.tol_pinv;
                            temp_mdd = pinv(G,pf.tol_pinv)*ris.data(t,:)';
                        end
                            


                        %temp_mdd = temp_mdd';
                        aux_n=0;
                        for  n=1:size(ris.seed,1)
                            if vector_ind_t(n)~=0
                                aux_n=aux_n+1;
                                ris.waveforms_seed(n,t)=norm(temp_mdd((aux_n-1)*DIM+1: aux_n*DIM));
                            end
                        end
                    end
                end
               
                
                
                
                
                %%% store in the struct ris the struct containig the clusters
                %%% dipoles, the initial number of cluster, the number of
                %%% components
                str = get(gui_vis.popup_comp, 'String');
                val = get(gui_vis.popup_comp,'Value');
                 %%% setting of the number of components
                switch str{val};
                    case 'location' %  
                        ris.comp_cluster=3;
                    case 'orientation and location'
                        ris.comp_cluster=6;
                        aux_cluster=ris.cluster;
                        aux_waveforms=ris.waveforms;
                        aux_seed=ris.seed;
                        ris.cluster = []; ris.seed = []; ris.waveforms = []; ris.waveforms_seed = [];
                        [ris.cluster,ris.seed,ris.waveforms,ris.waveforms_seed]=hades_script_cluster_orientation(aux_cluster,aux_seed,aux_waveforms);
                end
                               
                %%% here the clusters are splitted in according with the
                %%% criteria presented in the manual
                ris.new_cluster = []; ris.new_seed = []; ris.new_waveforms = [];ris.new_waveforms_seed = [];
                count = 1;
                for n_c=1:size(ris.cluster,2)
                    dipoli_cluster = ris.cluster(n_c).dipoli_cluster;
                    t_amp = dipoli_cluster(:,7);
                    i = 1;
                    while i < length(t_amp)
                        t_i = i;
                        while t_amp(i+1)-t_amp(i) <= 10
                            i = i+1;
                            if i == length(t_amp)
                                break;
                            end
                        end
                        t_fin = i;
                        if size(dipoli_cluster(t_i:t_fin,:),1) > 1
                            ris.new_cluster(count).index = ris.cluster(n_c).index(t_i:t_fin,:);
                            ris.new_cluster(count).dipoli_cluster = dipoli_cluster(t_i:t_fin,:);
                            ris.new_cluster(count).centroide = mean(ris.new_cluster(count).dipoli_cluster(:,1: ris.comp_cluster));
                            ris.new_cluster(count).n_dip_cluster = size(ris.new_cluster(count).dipoli_cluster,1);
                            ris.new_seed(count,:) = ris.new_cluster(count).centroide;
                            count = count +1;
                        end
                        i = i+1;
                    end
                end
                t1_amp = zeros(1,size(ris.new_seed,1));
                for i=1:size(ris.new_seed,1)
                  t1_amp(i) = ris.new_cluster(i).dipoli_cluster(1,7);
                end
                [t_i i_t] = sort(t1_amp);
                ris.new_cluster = ris.new_cluster(i_t);
                ris.new_seed = ris.new_seed(i_t,:);
                waveforms = zeros(size(ris.new_seed,1), ris.t2);
                for k = 1:size(ris.new_seed,1)
                    for j = 1:size(ris.new_cluster(k).dipoli_cluster,1)
                        waveforms(k,ris.new_cluster(k).dipoli_cluster(j,end-1)) =  norm(ris.new_cluster(k).dipoli_cluster(j,4:end-2));
                    end
                end
                ris.new_waveforms = waveforms;  
                ris.new_waveforms_seed=zeros(size(ris.new_seed,1),size(ris.time_interval,2));
                for t=1:size(ris.new_seed,1)
                    selected_seed=ris.new_seed(t,1:3);
                    diff=(pf.vertices-repmat(selected_seed,size(pf.vertices,1),1)).^2;
                    somma = sum(diff');
                    somma = somma';
                    dist = sqrt(somma);
                    [dist_riord,ind_riord] = sort(dist,'ascend');
                    seed_ok(t) = ind_riord(1);
                end
                for t=ris.t1:ris.t2
                    vector_ind_t=zeros(1,size(ris.new_seed,1));
                    G = []; temp_mdd = [];
                    for n=1:size(ris.new_seed,1)
                        [val_t, ind_t]=find(ris.new_cluster(n).dipoli_cluster(:,end-1)==t);
                        if isempty(ind_t)==1
                        else
                            vector_ind_t(n)=ind_t(1);
                        end
                    end
                    for  n=1:size(ris.new_seed,1)
                        if vector_ind_t(n)~=0
                            G = [G pf.g_matrix(:,DIM*(seed_ok(n)-1)+1:DIM*seed_ok(n)+DIM)];   
                        end
                    end
                    if isempty(G)==0
                        temp_mdd = pinv(G,pf.tol_pinv)*ris.data(t,:)';
                        %temp_mdd = temp_mdd';
                        aux_n=0;
                        for  n=1:size(ris.new_seed,1)
                            if vector_ind_t(n)~=0
                                aux_n=aux_n+1;
                                ris.new_waveforms_seed(n,t)=norm(temp_mdd(aux_n,:));
                            end
                        end
                    end        
                end
                
                
                
                
                
                
                %%% the clusters and the relative waveforms are plotted
                f = figure;
                for t=1:size(ris.subject,2)
                    if ris.subject(t)=='_'
                        string_case(t)=' ';
                    else
                        string_case(t)=ris.subject(t);
                    end
                end
                set(f,'Name',strcat(['Case: ',string_case,' with  ', num2str(NN),' clusters']));
                scrsz = get(0,'ScreenSize');
                pos=get(f,'pos');
                pos(1,1)=scrsz(1,3)*.1;
                pos(1,2)=scrsz(1,4)*.1;
                pos(1,3)=scrsz(1,3)*.5;
                pos(1,4)=scrsz(1,4)*.5;
                set(f,'pos',pos);
                
                n_color=round(size(ris.cluster,2)/12)+1;
                paletta_base=[0 0 1;  1 0 0;  0 1 0;  1 0 1;  0 1 1;  1 1 0;   .7 .3 .3;  .3 .3 .8; 1 0.5 0; 0.5 0 1; 1 0 0.5;  0.5 0.5 1 ];
                paletta=repmat(paletta_base,n_color,1);
                simboli_base=[repmat('o ',12,1);repmat('x ',12,1);repmat('v ',12,1);repmat('p ',12,1)];
                simboli=repmat(simboli_base,n_color,1);
                subplot(2,2,1)
                for i = 1:size(ris.cluster,2)
                    A = [];
                    A = ris.cluster(i).dipoli_cluster;
                    plot3(A(:,1),A(:,2),A(:,3),simboli(i),'Linewidth',3,'Color',paletta(i,:));
                    hold on
                end
                plot3(pf.vertices(:,1),pf.vertices(:,2),pf.vertices(:,3),'k.','markersize',1)
                if isempty(fwd.name_sensors)==0 && isempty(fwd.name_versors)==0
                    xs=load(fullfile(fwd.location_sensors,fwd.name_sensors));
                    plot3(xs(:,1),xs(:,2),xs(:,3),'yo')
                end
                axis tight
                set(gca,'XTickLabel','');
                set(gca,'YTickLabel','');
                set(gca,'ZTickLabel','');
                view(0, 90)
                subplot(2,2,2)
                for i = 1:size(ris.cluster,2)
                    A = [];
                    A = ris.cluster(i).dipoli_cluster;
                    plot3(A(:,1),A(:,2),A(:,3),simboli(i),'Linewidth',3,'Color',paletta(i,:));
                    hold on
                end
                plot3(pf.vertices(:,1),pf.vertices(:,2),pf.vertices(:,3),'k.','markersize',1)
                if isempty(fwd.name_sensors)==0 && isempty(fwd.name_versors)==0
                    xs=load(fullfile(fwd.location_sensors,fwd.name_sensors));
                    plot3(xs(:,1),xs(:,2),xs(:,3),'yo')
                end
                view(90,0)
                axis tight
                set(gca,'XTickLabel','');
                set(gca,'YTickLabel','');
                set(gca,'ZTickLabel','');
                subplot(2,2,3)
                for i = 1:size(ris.cluster,2)
                    A = [];
                    A = ris.cluster(i).dipoli_cluster;
                    plot3(A(:,1),A(:,2),A(:,3),simboli(i),'Linewidth',3,'Color',paletta(i,:));
                    hold on
                end
                plot3(pf.vertices(:,1),pf.vertices(:,2),pf.vertices(:,3),'k.','markersize',1)
                if isempty(fwd.name_sensors)==0 && isempty(fwd.name_versors)==0
                    xs=load(fullfile(fwd.location_sensors,fwd.name_sensors));
                    plot3(xs(:,1),xs(:,2),xs(:,3),'yo')
                end
                view(0,0)
                axis tight
                set(gca,'XTickLabel','');
                set(gca,'YTickLabel','');
                set(gca,'ZTickLabel','');

    
                subplot(2,2,4)
                T=ris.time_interval;
                for i = 1:size(ris.cluster,2)
                    plot(1000*T(1,1:size(ris.waveforms,2)),10^9*ris.waveforms(i,:),'Marker',simboli(i),'Color',paletta(i,:));
                    hold on
                    axis([1000*T(ris.t1),1000*T(ris.t2),0,max(max(10^9*ris.waveforms))]);
                    xlabel('time [ms]');
                    ylabel('source waveforms [nAm]');
                end
            end
        end
    end
end
return

function[]=clustering_save
global ris;
[save_name_ris,save_path_ris] = uiputfile('*.mat','Save the result');
if save_path_ris == 0 
    return
end
save(fullfile(save_path_ris,save_name_ris),'ris','-mat');
          
return

%%% function for showing the results
function[]=results_show
global gui_vis;
global ris;
set(gui_vis.label_comp,'visible','off');
set(gui_vis.popup_comp,'visible','off');
set(gui_vis.str_NN,'visible','off');
set(gui_vis.label_NN,'visible','off');
set(gui_vis.cmd_cluster,'visible','off');
set(gui_vis.cmd_save_cluster,'visible','off');
set(gui_vis.label_object,'visible','on');
set(gui_vis.popup_object,'visible','on','value',1);
set(gui_vis.label_vis,'visible','off');
set(gui_vis.str_vis,'visible','off');
set(gui_vis.cmd_vis,'visible','off');
set(gui_vis.label_split,'visible','off');
set(gui_vis.flag_split,'visible','off');
set(gui_vis.label_vis,'visible','on','string',strcat([num2str(1000*ris.time_interval(ris.t1)),' to ', num2str(1000*ris.time_interval(ris.t2)),' [ms]']));
set(gui_vis.str_vis,'visible','on','string','');
set(gui_vis.cmd_vis,'string','show selected time points','visible','on')
return


%%% function for setting the selected object
function[]=set_object
global pf;
global gui_vis;
global ris;
if isempty(ris.dip)==1
   msgbox('The PF algorithm did not reconstruct any source','Error','error');
   return
else
    if isempty(pf.vertices)==1
        msgbox('You have to load the vertices file!!!','Error','error');
    else
        str = get(gui_vis.popup_object, 'String');
        val = get(gui_vis.popup_object,'Value');
        switch str{val};
        case 'estimated dipoles' %%% reconstructed dipoles      
            set(gui_vis.label_vis,'visible','on','string',strcat([num2str(1000*ris.time_interval(ris.t1)),' to ', num2str(1000*ris.time_interval(ris.t2)),' [ms]']));
            set(gui_vis.str_vis,'visible','on','string','','position',[360 225 140 20]);
            set(gui_vis.cmd_vis,'string','show selected time points','visible','on')
            set(gui_vis.label_split,'visible','off');
            set(gui_vis.flag_split,'visible','off');
        case 'clustered dipoles' %%% splitted clusters
            if isfield(ris,'new_cluster')==0 || isfield(ris,'new_seed')==0
                msgbox('There are not clustered dipoles','Error','error');
                set(gui_vis.label_vis,'visible','off');
                set(gui_vis.str_vis,'visible','off');
                set(gui_vis.cmd_vis,'visible','off');
                return   
            elseif isfield(ris,'new_cluster')==1 && isfield(ris,'new_seed')==1 && isempty(ris.dip)==0
                if (get(gui_vis.flag_split,'Value')) == (get(gui_vis.flag_split,'Max'))
                     n_c=size(ris.new_cluster,2);
                else
                    n_c=size(ris.cluster,2);
                end
                set(gui_vis.label_vis,'visible','on','string',strcat([num2str(n_c),' clusters']))
                set(gui_vis.str_vis,'visible','on','string','','position',[360 225 50 20]);
                set(gui_vis.cmd_vis,'string','Clustered dipoles','visible','on')
                set(gui_vis.label_split,'visible','on');
                set(gui_vis.flag_split,'visible','on');
            end
        case 'estimated sources'    %%% splitted seeds
            if isfield(ris,'new_cluster')==0 || isfield(ris,'new_seed')==0
                msgbox('There are not clustered dipoles','Error','error');
                set(gui_vis.label_vis,'visible','off');
                set(gui_vis.str_vis,'visible','off');
                set(gui_vis.cmd_vis,'visible','off');
                return   
            elseif isfield(ris,'new_cluster')==1 && isfield(ris,'new_seed')==1 && isempty(ris.dip)==0
                if (get(gui_vis.flag_split,'Value')) == (get(gui_vis.flag_split,'Max'))
                     n_c=size(ris.new_seed,1);
                else
                    n_c=size(ris.seed,1);
                end
                set(gui_vis.label_vis,'visible','on','string',strcat([num2str(n_c),' estimated sources']))
                set(gui_vis.str_vis,'visible','on','string','','position',[360 225 50 20]);
                set(gui_vis.cmd_vis,'string','Estimated sources','visible','on')
                set(gui_vis.label_split,'visible','on');
                set(gui_vis.flag_split,'visible','on');
            end
        end
    end
end
return

function[]=set_clusters_number
global ris;
global pf;
global gui_vis;
if isempty(ris.dip)==1
   msgbox('The PF algorithm did not reconstruct any source','Error','error');
   return
else
    if isempty(pf.vertices)==1
        msgbox('You have to load the vertices file!!!','Error','error');
    else
        str = get(gui_vis.popup_object, 'String');
        val = get(gui_vis.popup_object,'Value');
        switch str{val};
        case 'clustered dipoles' %%% splitted clusters
            if isfield(ris,'new_cluster')==0 || isfield(ris,'new_seed')==0
                msgbox('There are not clustered dipoles','Error','error');
                set(gui_vis.label_vis,'visible','off');
                set(gui_vis.str_vis,'visible','off');
                set(gui_vis.cmd_vis,'visible','off');
                return   
            elseif isfield(ris,'new_cluster')==1 && isfield(ris,'new_seed')==1 && isempty(ris.dip)==0
                if (get(gui_vis.flag_split,'Value')) == (get(gui_vis.flag_split,'Max'))
                     n_c=size(ris.new_cluster,2);
                else
                    n_c=size(ris.cluster,2);
                end
                set(gui_vis.label_vis,'visible','on','string',strcat([num2str(n_c),' clusters']));
            end
        case 'estimated sources'    %%% splitted seeds
            if isfield(ris,'new_cluster')==0 || isfield(ris,'new_seed')==0
                msgbox('There are not clustered dipoles','Error','error');
                set(gui_vis.label_vis,'visible','off');
                set(gui_vis.str_vis,'visible','off');
                set(gui_vis.cmd_vis,'visible','off');
                return   
            elseif isfield(ris,'new_cluster')==1 && isfield(ris,'new_seed')==1 && isempty(ris.dip)==0
                if (get(gui_vis.flag_split,'Value')) == (get(gui_vis.flag_split,'Max'))
                     n_c=size(ris.new_seed,1);
                else
                    n_c=size(ris.seed,1);
                end
                set(gui_vis.label_vis,'visible','on','string',strcat([num2str(n_c),' estimated sources']))
            end
        end
    end
end
return

%%% function for showing the selected clusters or seeds
function[]=vis_clusters_or_seeds
global pf;
global gui_vis;
global ris;
global fwd;
if isempty(ris.dip)==1
   msgbox('The PF algorithm did not reconstruct any source','Error','error');
   return
else
    %%% Set current data to the selected data set
    str = get(gui_vis.popup_object, 'String');
    val = get(gui_vis.popup_object,'Value');
    switch str{val};
    case 'estimated dipoles' %%% reconstructed dipoles
        str_dip=get(gui_vis.str_vis,'string');
        time_dip=str2num(str_dip);
        if strcmp(str_dip,'')==0 && isempty(time_dip)==0
            fig_dip=figure;
            set(fig_dip,'Name',strcat(['Estimated dipoles ',str_dip,' ms']));
            hold on
            ind_first=round((time_dip(1)-ris.zero_time)/1000*(ris.freq)+1);
            ind_fin=round((time_dip(end)-ris.zero_time)/1000*(ris.freq)+1);
            time_estimated_dip=ind_first:1:ind_fin;
            subplot(2,2,1)
            for t=1:size(time_estimated_dip,2)
                ind2plot=(find(ris.dip(:,end-1)==time_estimated_dip(t)));
                if isempty(ind2plot)==0
                   plot3(ris.dip(ind2plot,1),ris.dip(ind2plot,2),ris.dip(ind2plot,3),'ro','Linewidth',3);
                   hold on
                end
            end
            plot3(pf.vertices(:,1),pf.vertices(:,2),pf.vertices(:,3),'k.','markersize',1)
            if isempty(fwd.name_sensors)==0 && isempty(fwd.name_versors)==0
                xs=load(fullfile(fwd.location_sensors,fwd.name_sensors));
                plot3(xs(:,1),xs(:,2),xs(:,3),'yo')
            end
            view(0, 90)
            set(gca,'XTickLabel','');
            set(gca,'YTickLabel','');
            set(gca,'ZTickLabel','');
            axis tight
            subplot(2,2,2)
            for t=1:size(time_estimated_dip,2)
                ind2plot=(find(ris.dip(:,end-1)==time_estimated_dip(t)));
                if isempty(ind2plot)==0
                   plot3(ris.dip(ind2plot,1),ris.dip(ind2plot,2),ris.dip(ind2plot,3),'ro','Linewidth',3);
                   hold on
                end
            end
            plot3(pf.vertices(:,1),pf.vertices(:,2),pf.vertices(:,3),'k.','markersize',1)
            if isempty(fwd.name_sensors)==0 && isempty(fwd.name_versors)==0
                xs=load(fullfile(fwd.location_sensors,fwd.name_sensors));
                plot3(xs(:,1),xs(:,2),xs(:,3),'yo')
            end
            view(90,0)
            set(gca,'XTickLabel','');
            set(gca,'YTickLabel','');
            set(gca,'ZTickLabel','');
            axis tight
            subplot(2,2,3)
            for t=1:size(time_estimated_dip,2)
                ind2plot=(find(ris.dip(:,end-1)==time_estimated_dip(t)));
                if isempty(ind2plot)==0
                   plot3(ris.dip(ind2plot,1),ris.dip(ind2plot,2),ris.dip(ind2plot,3),'ro','Linewidth',3);
                   hold on
                end
            end
            plot3(pf.vertices(:,1),pf.vertices(:,2),pf.vertices(:,3),'k.','markersize',1)
            if isempty(fwd.name_sensors)==0 && isempty(fwd.name_versors)==0
                xs=load(fullfile(fwd.location_sensors,fwd.name_sensors));
                plot3(xs(:,1),xs(:,2),xs(:,3),'yo')
            end
            view(0,0)
            set(gca,'XTickLabel','');
            set(gca,'YTickLabel','');
            set(gca,'ZTickLabel','');
            axis tight
            subplot(2,2,4)
            for t=1:size(time_estimated_dip,2)
                ind2plot=(find(ris.dip(:,end-1)==time_estimated_dip(t)));
                if isempty(ind2plot)==0
                   plot3(ris.dip(ind2plot,1),ris.dip(ind2plot,2),ris.dip(ind2plot,3),'ro','Linewidth',3);
                   hold on
                end
            end
            plot3(pf.vertices(:,1),pf.vertices(:,2),pf.vertices(:,3),'k.','markersize',1)
            if isempty(fwd.name_sensors)==0 && isempty(fwd.name_versors)==0
                xs=load(fullfile(fwd.location_sensors,fwd.name_sensors));
                plot3(xs(:,1),xs(:,2),xs(:,3),'yo')
            end
            set(gca,'XTickLabel','');
            set(gca,'YTickLabel','');
            set(gca,'ZTickLabel','');
            axis tight
            
        else
            fig_dip=figure;
            set(fig_dip,'Name',strcat(['Estimated dipoles from ',num2str(1000*ris.time_interval(ris.t1)),' to ', num2str(1000*ris.time_interval(ris.t2)),' ms']));
            
            subplot(2,2,1)
            plot3(ris.dip(:,1),ris.dip(:,2),ris.dip(:,3),'ro','Linewidth',3);
            hold on
            plot3(pf.vertices(:,1),pf.vertices(:,2),pf.vertices(:,3),'k.','markersize',1)
            if isempty(fwd.name_sensors)==0 && isempty(fwd.name_versors)==0
                xs=load(fullfile(fwd.location_sensors,fwd.name_sensors));
                plot3(xs(:,1),xs(:,2),xs(:,3),'yo')
            end
            view(0, 90)
            set(gca,'XTickLabel','');
            set(gca,'YTickLabel','');
            set(gca,'ZTickLabel','');
            axis tight
            subplot(2,2,2)
            plot3(ris.dip(:,1),ris.dip(:,2),ris.dip(:,3),'ro','Linewidth',3);
            hold on
            plot3(pf.vertices(:,1),pf.vertices(:,2),pf.vertices(:,3),'k.','markersize',1)
            if isempty(fwd.name_sensors)==0 && isempty(fwd.name_versors)==0
                xs=load(fullfile(fwd.location_sensors,fwd.name_sensors));
                plot3(xs(:,1),xs(:,2),xs(:,3),'yo')
            end
            view(90,0)
            set(gca,'XTickLabel','');
            set(gca,'YTickLabel','');
            set(gca,'ZTickLabel','');
            axis tight
            subplot(2,2,3)
            plot3(ris.dip(:,1),ris.dip(:,2),ris.dip(:,3),'ro','Linewidth',3);
            hold on
            plot3(pf.vertices(:,1),pf.vertices(:,2),pf.vertices(:,3),'k.','markersize',1)
            if isempty(fwd.name_sensors)==0 && isempty(fwd.name_versors)==0
                xs=load(fullfile(fwd.location_sensors,fwd.name_sensors));
                plot3(xs(:,1),xs(:,2),xs(:,3),'yo')
            end
            view(0,0)
            set(gca,'XTickLabel','');
            set(gca,'YTickLabel','');
            set(gca,'ZTickLabel','');
            axis tight
            subplot(2,2,4)
            plot3(ris.dip(:,1),ris.dip(:,2),ris.dip(:,3),'ro','Linewidth',3);
            hold on
            plot3(pf.vertices(:,1),pf.vertices(:,2),pf.vertices(:,3),'k.','markersize',1)
            if isempty(fwd.name_sensors)==0 && isempty(fwd.name_versors)==0
                xs=load(fullfile(fwd.location_sensors,fwd.name_sensors));
                plot3(xs(:,1),xs(:,2),xs(:,3),'yo')
            end
            set(gca,'XTickLabel','');
            set(gca,'YTickLabel','');
            set(gca,'ZTickLabel','');  
            axis tight
        end
    case 'clustered dipoles' %%% splitted clusters
        if isfield(ris,'new_cluster')==0 || isfield(ris,'new_seed')==0
            msgbox('There are not clustered dipoles','Error','error');
            return
        elseif isfield(ris,'new_cluster')==1 && isfield(ris,'new_seed')==1 && isempty(ris.dip)==0
            if (get(gui_vis.flag_split,'Value')) == (get(gui_vis.flag_split,'Max'))
                hades_plot_clusters(ris.new_cluster,ris.new_waveforms);
            else
                hades_plot_clusters(ris.cluster,ris.waveforms);
            end
        end
    case 'estimated sources'   %%% splitted seeds 
        if isfield(ris,'new_cluster')==0 || isfield(ris,'new_seed')==0
            msgbox('There are not clustered dipoles','Error','error');
            return
        elseif isfield(ris,'new_cluster')==1 && isfield(ris,'new_seed')==1 && isempty(ris.dip)==0 
            if (get(gui_vis.flag_split,'Value')) == (get(gui_vis.flag_split,'Max'))
                hades_plot_seeds(ris.new_seed,ris.new_waveforms_seed);
            else
                hades_plot_seeds(ris.seed,ris.waveforms_seed);
            end
        end
    end
end
return

%%% function for showing the model selection
function[]=vis_mod_sel
global ris;
fig_mod=figure;
set(fig_mod,'Name','Model selection');
area(1000*ris.time_interval(ris.t1:1:ris.t1+size(ris.mod_sel,1)-1),ris.mod_sel);
axis([1000*ris.time_interval(ris.t1) 1000*ris.time_interval(ris.t1+size(ris.mod_sel,1)-1) 0 1]);
xlabel('time [ms]');
ylabel('model order probability');
legend('zero dipole','one dipole','two dipole','three dipole','four dipole','five dipole','location','best');
return


