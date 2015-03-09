%%%%   HADES_SCRIPT_CLUSTER clusters the estimated dipoles
% 
%   Copyright (C) 2010, Cristina Campi, Annalisa Pascarella, Michele Piana, Alberto Sorrentino


%%% creation of the struct struct_dipoli_cluster containing the clustered
%%% dipoles
n_comp=3;
clear struct_dipoli_cluster 
n_ist = max(dipoli(:,end));
dips_sel = dipoli;
n_el = size(dips_sel,1);
dips_sel(:,1:n_comp);


[xt_cluster_kmeans centroidi_kmeans] = kmeans(dips_sel(:,1:n_comp),NN,'replicates',5,'emptyaction','drop');
clear xt_cluster n_dip_cluster centroidi
xt_cluster = xt_cluster_kmeans;   centroidi = centroidi_kmeans;  
ind=1;
NN_aux=NN;
for n_cluster=1:NN_aux
  clear index_cluster 
  index_cluster = find(xt_cluster == n_cluster);
  n_dip_cluster(n_cluster) = size(dipoli(index_cluster,1),1);
  
  if n_dip_cluster(n_cluster)==0
      NN=NN-1;
  else
      struct_dipoli_cluster(ind).index = index_cluster;
      struct_dipoli_cluster(ind).dipoli_cluster = dip(index_cluster,:);
      struct_dipoli_cluster(ind).n_dip_cluster = n_dip_cluster(n_cluster);
      struct_dipoli_cluster(ind).centroide = centroidi(n_cluster,:);  %AP 240408  
      ind=ind+1;
  end
end
   
