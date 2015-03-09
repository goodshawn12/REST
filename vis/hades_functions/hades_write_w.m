function[]=hades_write_w(src,obj,string,directory_name)
%%% HAEDS_WRITE_w write the particle filter results in .w files.
%   Copyright (C) 2010, Cristina Campi, Annalisa Pascarella, Michele Piana, Alberto Sorrentino

global ris;

filename_tot_dip = fullfile(directory_name,strcat(ris.subject,'_dipoles'));

n_sep = src(1,1).nuse;

obj_sort = sort(obj(:,end),'ascend');
obj_lh = [];
obj_rh = [];
i=1;j=1;
for t=1:size(obj_sort,1)
    if obj_sort(t)<=n_sep
        obj_lh(i)=obj_sort(t);
        i=i+1;
    else
        obj_rh(j)=obj_sort(t)-n_sep;
        j=j+1;
    end
end

i=1; j=1;
obj_lh=obj_lh';
obj_rh=obj_rh';
if isempty(obj_lh)==0
    for t=1:size(obj_lh,1)
        aux=obj_lh(1);
        da_elim = find(obj_lh==aux);
        vert_lh(i) = src(1,1).vertno(aux);
        data_lh(i) = size(da_elim,1);
        i = i + 1;
        obj_lh = obj_lh(da_elim(end)+1:end);
       if isempty(obj_lh)==1
           break
       end
    end
    vert_lh = uint32(vert_lh');
    data_lh = data_lh';
    struct_lh = struct('vertices',vert_lh,'data',data_lh);
    mne_write_w_file(strcat(filename_tot_dip,'-lh.w'),struct_lh );
end
if isempty(obj_rh)==0
    for t=1:size(obj_rh,1)
        aux=obj_rh(1);
        da_elim = find(obj_rh==aux);
        vert_rh(j) = src(1,2).vertno(aux);
        data_rh(j) = size(da_elim,1);
        j = j + 1;
        obj_rh = obj_rh(da_elim(end)+1:end);
       if isempty(obj_rh)==1
           break
       end
    end
    vert_rh = uint32(vert_rh');
    data_rh = data_rh';
    struct_rh = struct('vertices',vert_rh,'data',data_rh);
    mne_write_w_file(strcat(filename_tot_dip,'-rh.w'),struct_rh );
end


return