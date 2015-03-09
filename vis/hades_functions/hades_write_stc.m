function[]=hades_write_stc(sub_fwd_sol,vertices,vert_coord,tmin,tstep,dip,t1,t2,directory_name)
%%% HAEDS_WRITE_STC write the particle filter results in .stc files.
%   Copyright (C) 2010, Cristina Campi, Annalisa Pascarella, Michele Piana, Alberto Sorrentino
global ris;
tot_dip = zeros(size(vertices,1),1);
dip_ist = zeros(size(vertices,1),t2-t1+1);

%%%% stc for the maximum of the blob
for t = 1:size(dip,1)
    if dip(t,end)>=t1 && dip(t,end)<=t2 
        tot_dip(dip(t,end),1)= tot_dip(dip(t,end),1)+1;
        dip_ist(dip(t,end),dip(t,end-1)-ris.t1+1) = dip_ist(dip(t,end),dip(t,end-1)-ris.t1+1)+1;
    end
end
dip_ist=[tot_dip,dip_ist(:,2:end)];
       
%%% filenames
if strcmp(ris.subject(end-3:end),'.mat')==1 
    filename_dip = strcat(strtrim(directory_name),'\',strtrim(ris.subject(1:end-4)),'_dipoles');
else
    filename_dip = strcat(strtrim(directory_name),'\',strtrim(ris.subject),'_dipoles');
end
tot_dip_stc_lh = struct('tmin',tmin,'tstep',tstep,'vertices',vertices(1:sub_fwd_sol.src(1,1).nuse),'data',dip_ist(1:sub_fwd_sol.src(1,1).nuse,:));
tot_dip_stc_rh = struct('tmin',tmin,'tstep',tstep,'vertices',vertices(sub_fwd_sol.src(1,1).nuse+1:end),'data',dip_ist(sub_fwd_sol.src(1,1).nuse+1:end,:));

%%%% write the tot_dip stc files
mne_write_stc_file(strcat(filename_dip,'-lh.stc'),tot_dip_stc_lh);
mne_write_stc_file(strcat(filename_dip,'-rh.stc'),tot_dip_stc_rh);
return