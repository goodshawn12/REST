%   Copyright (C) 2010, Cristina Campi, Annalisa Pascarella, Michele Piana, Alberto Sorrentino

function [new_clusters,new_seeds,new_waveforms,new_waveforms_seed]=hades_script_cluster_orientation(clusters,seeds,waveforms)
global ris;
global pf;
m=1;
new_seeds=[];
new_clusters=[];
new_waveforms=zeros(1,size(ris.time_interval,2));
new_waveforms_seed=zeros(1,size(ris.time_interval,2));
for t=1:size(clusters,2)
    all_different=0;
    MM=round(size(clusters(t).dipoli_cluster,1)/10);
    if MM>1 %%%the t-th cluster is clusterized taking account also of the orientations
        NN=MM;
          while all_different == 0  && NN > 0
            clear struct_dipoli_cluster_orient 
            dips_sel = clusters(t).dipoli_cluster;
            [xt_cluster_kmeans centroidi_kmeans] = kmeans(dips_sel(:,4:6),NN,'replicates',5,'emptyaction','drop');
            clear xt_cluster n_dip_cluster centroidi
            xt_cluster = xt_cluster_kmeans;   centroidi = centroidi_kmeans;   
            ind=1;
            NN_aux=NN;
            for n_cluster=1:NN_aux
              clear index_cluster 
              index_cluster = find(xt_cluster == n_cluster);
              n_dip_cluster(n_cluster) = size(dips_sel(index_cluster,1),1);
              if n_dip_cluster(n_cluster)==0
                 NN=NN-1;
              else
                struct_dipoli_cluster_orient(ind).index = index_cluster;
                struct_dipoli_cluster_orient(ind).dipoli_cluster = dips_sel(index_cluster,:);
                struct_dipoli_cluster_orient(ind).n_dip_cluster = n_dip_cluster(n_cluster);
                struct_dipoli_cluster_orient(ind).centroide = mean(dips_sel(index_cluster,1:6));  %AP 240408 
                ind=ind+1;
              end
            end
            for II = 1:size(struct_dipoli_cluster_orient,2)
                 for JJ = 1:size(struct_dipoli_cluster_orient,2)     
                    mean_cl1 = mean(struct_dipoli_cluster_orient(II).dipoli_cluster(:,4:6));
                    std_cl1 = std(struct_dipoli_cluster_orient(II).dipoli_cluster(:,4:6));
                    mean_cl2 = mean(struct_dipoli_cluster_orient(JJ).dipoli_cluster(:,4:6));
                    std_cl2 = std(struct_dipoli_cluster_orient(JJ).dipoli_cluster(:,4:6));
                    ctrl1 = abs(mean_cl1 - mean_cl2);
                    if II ~=JJ && (max(ctrl1-std_cl1) < 0 || max(ctrl1 - std_cl2) < 0)
                       all_different = 0;
                    end
                 end
            end
            NN=NN-1;
        end
        NN = NN+1;
        clear waveforms_orient seeds_orient
        seed_orient=zeros(size(struct_dipoli_cluster_orient,2),6);
        waveforms_orient=zeros(size(struct_dipoli_cluster_orient,2),size(ris.time_interval,2));
        for k = 1:size(struct_dipoli_cluster_orient,2)
            for j = 1:size(struct_dipoli_cluster_orient(k).dipoli_cluster,1)
                waveforms_orient(k,struct_dipoli_cluster_orient(k).dipoli_cluster(j,end-1)) =  norm(struct_dipoli_cluster_orient(k).dipoli_cluster(j,4:end-2));
                seeds_orient(k,1:6)=struct_dipoli_cluster_orient(k).centroide;
            end
        end
        %%% store in the struct ris the values of the seed positions and waveforms
         for nn=1:size(struct_dipoli_cluster_orient,2)
            wave=waveforms_orient(nn,:);     
            new_waveforms=[new_waveforms;wave];
      
            new_clusters(m).index = struct_dipoli_cluster_orient(nn).index;
            new_clusters(m).dipoli_cluster = struct_dipoli_cluster_orient(nn).dipoli_cluster;
            new_clusters(m).centroide = mean(struct_dipoli_cluster_orient(nn).dipoli_cluster(:,1:3));
            new_clusters(m).n_dip_cluster = size(struct_dipoli_cluster_orient(nn).dipoli_cluster,1);
            new_seeds(m,:) = seeds_orient(nn,:);
            m=m+1;
        end
    else
        wave=[waveforms(t,:), zeros(1,ris.t2-ris.t1+1-size(waveforms,2))];
        new_waveforms=[new_waveforms;wave];
        new_clusters(m).index = clusters(t).index;
        new_clusters(m).dipoli_cluster = clusters(t).dipoli_cluster;
        new_clusters(m).centroide = mean(clusters(t).dipoli_cluster(:,1:3));
        new_clusters(m).n_dip_cluster = size(clusters(t).dipoli_cluster,1);
        new_seeds(m,:) = [seeds(t,1:3),mean(clusters(t).dipoli_cluster(:,4:6),1)];
        m=m+1;
    end
end

new_waveforms=new_waveforms(2:end,:);

if isempty(pf.orient)==1
   DIM = 3;  %%% if the dipoles orientation are not fixed
else
   DIM = 1;  %%% if the dipoles orientation are fixed
end
new_waveforms_seed=zeros(size(new_seeds,1),size(ris.time_interval,2));

for t=1:size(new_seeds,1)
    selected_seed=new_seeds(t,1:3);
    diff=(pf.vertices-repmat(selected_seed,size(pf.vertices,1),1)).^2;
    somma = sum(diff');
    somma = somma';
    dist = sqrt(somma);
    [dist_riord,ind_riord] = sort(dist,'ascend');
    seed_ok(t) = ind_riord(1);
end
for t=ris.t1:ris.t2
    vector_ind_t=zeros(1,size(new_seeds,1));
    G = []; temp_mdd = [];
    for n=1:size(new_seeds,1)
        [val_t, ind_t]=find(new_clusters(n).dipoli_cluster(:,end-1)==t);
        if isempty(ind_t)==1
        else
          vector_ind_t(n)=ind_t(1);
        end
    end
    for  n=1:size(new_seeds,1)
      if vector_ind_t(n)~=0
        G = [G pf.g_matrix(:,DIM*(seed_ok(n)-1)+1:DIM*seed_ok(n)+DIM)];   
      end
    end
    if isempty(G)==0
        temp_mdd = pinv(G,pf.tol_pinv)*ris.data(t,:)';
        %temp_mdd = temp_mdd';
        aux_n=0;
        for  n=1:size(new_seeds,1)
            if vector_ind_t(n)~=0
                aux_n=aux_n+1;
                new_waveforms_seed(n,t)=norm(temp_mdd(aux_n,:));
            end
        end
    end
end