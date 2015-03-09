%%%%   HADES_STATISTICAL_TEST_CLUSTER performs a statistical test on the
%%%%   clusters computed in hades_script_cluster.
% 
%   Copyright (C) 2010, Cristina Campi, Annalisa Pascarella,Michele Piana,  Alberto
%   Sorrentino

%%% statistical test on the significance of the clusters

all_different = 1;
for II = 1:NN
    for JJ = 1:NN
        mean_cl1 = mean(struct_dipoli_cluster(II).dipoli_cluster(:,1:3));
        std_cl1 = std(struct_dipoli_cluster(II).dipoli_cluster(:,1:3));
        mean_cl2 = mean(struct_dipoli_cluster(JJ).dipoli_cluster(:,1:3));
        std_cl2 = std(struct_dipoli_cluster(JJ).dipoli_cluster(:,1:3));
        ctrl1 = abs(mean_cl1 - mean_cl2);
            if II ~=JJ && (max(ctrl1-std_cl1) < 0 || max(ctrl1 - std_cl2) < 0)
                all_different = 0;
            end
    end
end
