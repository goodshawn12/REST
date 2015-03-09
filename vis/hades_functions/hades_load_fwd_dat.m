 function[]=hades_load_fwd_dat(loading_file)
%%% HADES_LOAD_FWD_DAT loads the lead field and the source space when a .dat file is selected
%%% in the load field window.
%   Copyright (C) 2010, Cristina Campi, Annalisa Pascarella, Michele Piana, Alberto Sorrentino
global pf;
global fwd;
if isempty(fwd.location_vertices)==1 || ischar(fwd.location_vertices)==0
    msgbox('You have to load the vertices matrix','Error','error');
else
    msgbox('Loading the lead field'); 
    pf.g_matrix_dir=loading_file;
    ext = loading_file(end-3:end);
    if strcmp(ext,'.dat')==1
        pf.g_matrix=load(loading_file);
        
    elseif strcmp(ext,'.mat')==1
        field_mat=load(loading_file);   
        pf.g_matrix=getfield(field_mat,cell2mat(fieldnames(field_mat)));
    end
    if isempty(pf.vertices_dir)==1
        pf.vertices_dir=fullfile(fwd.location_vertices,fwd.name_vertices);
    end
    ext =  pf.vertices_dir(end-3:end);
    if strcmp(ext,'.dat')==1
        pf.vertices=load(pf.vertices_dir);
    elseif strcmp(ext,'.mat')==1
        vert_mat=load(pf.vertices_dir);        
        pf.vertices=getfield(vert_mat,cell2mat(fieldnames(vert_mat)));
    end     
    if isempty(fwd.location_orient)==0
        if fwd.location_orient==0
            pf.orient=[];
            pf.orient_dir=[];
        else
            pf.orient_dir=fullfile(fwd.location_orient,fwd.name_orient);
            ext =  pf.orient_dir(end-3:end);
            if strcmp(ext,'.dat')==1
                pf.orient=load(pf.orient_dir);
            elseif strcmp(ext,'.mat')==1
                orient_mat=load(pf.orient_dir);        
                pf.orient=getfield(orient_mat,cell2mat(fieldnames(orient_mat)));
            end     
        end
    else
        if isempty(pf.orient_dir)==0
            ext =  pf.orient_dir(end-3:end);
            if strcmp(ext,'.dat')==1
                pf.orient=load(pf.orient_dir);
            elseif strcmp(ext,'.mat')==1
                orient_mat=load(pf.orient_dir);        
                pf.orient=getfield(orient_mat,cell2mat(fieldnames(orient_mat)));
            end  
        else
            pf.orient=[];
        end
    end
    close;
    close(fwd.fig);
    fwd.check=1;
    msgbox('The lead field was loaded');
end