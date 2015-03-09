function hades_load_fwd_fif(loading_file)
%%% HADES_LOAD_FWD_FIF loads the lead field and the source space when a .dat file is selected
%%% in the load field window.
%   Copyright (C) 2010, Cristina Campi, Annalisa Pascarella, Michele Piana, Alberto Sorrentino
global pf;
global fwd;
if exist('mne_read_forward_solution.m')==0
    close(fwd.fig)
    msg_err=msgbox('The file mne_read_forward_solution.m is not in the path','Error','error');
    figure(msg_err);
    return
else
    msgbox('Loading the lead field','replace'); 
    pf.g_matrix_dir=loading_file;
    fwd_aux = mne_read_forward_solution(loading_file);
    pf.g_matrix=fwd_aux.sol.data;
    pf.vertices=fwd_aux.source_rr;
    if size(fwd_aux.source_rr,1)==size(fwd_aux.source_nn,1)
         pf.orient=fwd_aux.source_nn;
    end
    pf.vertices_dir=[];
    pf.orient_dir=[];
    close;
    close(fwd.fig);
    fwd.check=1;
    msgbox('The lead field was loaded');
end 