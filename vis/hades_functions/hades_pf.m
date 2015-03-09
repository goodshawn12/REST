%   HADES_PF runs the particle filter algorithm.
% 
%   Copyright (C) 2010, Cristina Campi, Annalisa Pascarella, Michele Piana, Alberto
%   Sorrentino



global pf;
global ris;
global string_ris;
global pf_gui;
global check_run;

set(pf_gui.cmd_GO,'enable','off');

ris.dip=[];
ris.mod_sel=[];
ris.cluster=[];
ris.new_cluster=[];
ris.waveforms=[];
ris.new_waveforms=[];
ris.seed=[];
ris.new_seed=[];
ris.comp_cluster=[];


fact=[1 1 2 6 24 120];

trial_1 = 1;
trial_2 = 1; 

%%% time interval for the analysis
if isempty(pf.t1)==1
   t1 = 1;
   t2 = size(pf.data,2);
else
   t1 = pf.t1;
   t2 = pf.t2; 
end

%%% check if the orientations are fixed

if isempty(pf.orient)==1
   DIM = 3;  %%% if the dipoles orientation are not fixed
else
   DIM = 1;  %%% if the dipoles orientation are fixed
end




NIST = t2-t1+1; %%% number of time points
NDIP = 5;       %%% maximun number of simultaneous dipoles
NVERT = size(pf.vertices,1);   %%% number of vertices 
NPART = pf.np;  %%% number of particles 
NSENS = size(pf.g_matrix,1);  %%% number of sensors
SIGMA_Q_IN = 10*(max(max(abs(pf.data)))/max(max(abs(pf.g_matrix)))); %%% value of the variance for the dipole moment initialization
SIGMA_Q = SIGMA_Q_IN/7;  %%% value of variance for the dipole moment evolution

%%% in the struct ris are copied some of the fields of the struct pf
ris.subjects_dir=pf.subjects_dir;
ris.results_dir=pf.results_dir;
ris.subject=pf.subject;
ris.data_dir=pf.data_dir;
ris.data=pf.data;
ris.evol=pf.evol_dir;
ris.zero_time=pf.zero_time;
ris.final_time=pf.final_time;
ris.freq=pf.freq;
ris.time_interval=pf.time_interval;
ris.t1=t1;
ris.t2=t2;
ris.particles=pf.np;
ris.ssp=pf.ssp_dir;
ris.cov_matrix=pf.cov_matrix_dir;
if ischar(pf.sigma)==1
    ris.sigma=str2num(pf.sigma);
else
    ris.sigma=pf.sigma;
end
ris.sigma_q_in = SIGMA_Q_IN;
ris.sigma_q = SIGMA_Q;

if isempty(pf.sigma_par)==1
    pf.sigma_par=1;
end
ris.sigma_par=pf.sigma_par;

if norm(pf.cov_matrix-eye(size(pf.g_matrix,1)))==0
    SIGMA = ris.sigma*ris.sigma_par  %%% mean noise variance
else
    SIGMA=ris.sigma_par
end

data = (pf.matrix_pw*pf.ssp*pf.data')';     %%% magnetic field
g_matrix = pf.matrix_pw*pf.ssp* pf.g_matrix; %%% lead field matrix
pf.tol_pinv=max(size(g_matrix)) * norm(g_matrix) * eps(class(g_matrix))


for tr = trial_1 : trial_2
   
    tic
    %%% allocation of the variables 
    mod_sel = zeros(NIST,NDIP+1);
    blob = zeros(NVERT,NIST);
    amp_blob = zeros(NVERT,NIST);

    ndip = zeros(NPART,1);
    pos = zeros(NPART,NDIP);
    amp = zeros(NPART,NDIP*DIM);
    cov_matrix = zeros(NDIP*DIM,NDIP*DIM,NPART);

    ndip_post = zeros(NPART,1);
    pos_post = zeros(NPART,NDIP);
    amp_post = zeros(NPART,NDIP*DIM);

  

    %%% initialization at the first time point
    ndip_post=floor((NDIP+1)*rand(NPART,1));  %%% number of dipoles in the particles
    for i = 1 : NPART
        %%% check that the number of dipoles is bigger than 0 e smaller
        %%% than NDIP
        if ndip_post(i) < 0
            ndip_post(i) = 0;
        elseif ndip_post(i) > NDIP
            ndip_post(i) = NDIP;
        end
        %%% position and dipole moment are associated to easch dipole in the particle 
        for j = 1 : ndip_post(i)
            pos_post(i,j) = floor((NVERT+1)*rand);
            amp_post(i,DIM*(j-1)+1:DIM*j) = SIGMA_Q_IN*randn(1,DIM);
            if pos_post(i,j) <= 0;
                pos_post(i,j) = 1;
            elseif pos_post(i,j)> NVERT
                pos_post(i,j) = NVERT;
            end
        end
    end

 
    bar_pf=waitbar(0,'please wait','Name','the particle filter is running');
    for t = t1 : t2
        waitbar((t-t1+1)/NIST);
        pos = []; amp = []; ndip = []; cov_matrix = [];

        pos = pos_post;
        amp = amp_post;
        ndip = ndip_post;

        num_dip(:,t-t1+1) = ndip;
        for i = 1 : NPART
                n = ndip(i); %%% number of dipole of the i-th particle
                if n>0
                    G = []; Gq = []; GtS = []; GtG = []; K = [];
                    pos_aux1 = []; pos_aux2 = []; pos_aux3 = [];pos_tot = [];
                    if DIM == 3                        
                        pos_aux1 = 3*(pos(i,1:n)-1)+1;
                        pos_aux2 = 3*(pos(i,1:n)-1)+2;
                        pos_aux3 = 3*(pos(i,1:n)-1)+3;
                        pos_tot = [pos_aux1' pos_aux2' pos_aux3'];
                        pos_tot = reshape(pos_tot',1,[]);
                    elseif DIM == 1
                        pos_tot = pos(i,1:n);
                    end

                    G = g_matrix(:,pos_tot); %%% lead field related to the dipoles of the i-th particles
                    Gq = G*amp(i,1:n*DIM)';	% NSENS x 1
                    w(i) = (1/fact(n+1))*exp((-norm(data(t,:)'-Gq)^2)/(2*NSENS*SIGMA*SIGMA)); %%% computation of the weight
                    blob(pos(i,1:n),t-t1+1) = blob(pos(i,1:n),t-t1+1)+w(i); %%% setting of the weight in the blob
                else
                    w(i) = (1/fact(n+1))*exp((-norm(data(t,:)')^2)/(2*NSENS*SIGMA*SIGMA));
                end
                mod_sel(t-t1+1,n+1)=mod_sel(t-t1+1,n+1)+w(i); %%% seeting of the weight in the model selection
        end
        
        %%% normalization
        peso_tot = sum(w);
        w = w/peso_tot;
        w_tot(:,t-t1+1)=w;
        blob(:,t-t1+1) = blob(:,t-t1+1)/peso_tot;
        mod_sel(t-t1+1,:) = mod_sel(t-t1+1,:)/peso_tot;
        
        %%% resampling
        c = 0; 
        c(1) = w(1); 

        for k = 2 : NPART 
            c(k) = c(k-1) + w(k); 
        end
        c(NPART + 1) = 1;
        u = 0; 
        u(1) = 1./NPART .* rand(1);
        y=1;
        for j = 1 : NPART 
            u(j) = u(1) + (j-1) ./ NPART; 
            while u(j) >= c(y) 
                y = y + 1;       
            end 
            if y == NPART +1
                pos_post(j,:) = pos(NPART,:); 
                amp_post(j,:)= amp(NPART,:);
                ndip_post(j)= ndip(NPART);
            else
                pos_post(j,:) = pos(y,:);
                amp_post(j,:)= amp(y,:);
                ndip_post(j) = ndip(y);
            end
        end    
        %%% end resampling


        %%% particles evolution
        for i = 1 : NPART
            n = ndip_post(i);
            if n > 0 
                for j = 1 : n  
                    vicini = find(pf.evol(pos_post(i,j),:)>0);
                    el = floor((size(vicini,2)+1)*rand);
                    if el <= 0
                        el = 1;
                    elseif el > size(vicini,2)
                        el = size(vicini,2);
                    end
                    pos_post(i,j) = pf.evol(pos_post(i,j),el);
                    amp_post(i,DIM*(j-1)+1:DIM*j) = amp_post(i,DIM*(j-1)+1:DIM*j)+SIGMA_Q * randn(1,DIM);
                end
            end
        %%% end evolution
        
        %%% creation and death of dipoles
            if n == 0 
                if rand > 2/3  %%% a dipole is created
                    ndip_post(i) = 1;
                    pos_post(i,1) = floor((NVERT+1)*rand);
                    if pos_post(i,1) <= 0;
                        pos_post(i,1) = 1;
                    elseif pos_post(i,1)> NVERT
                        pos_post(i,1) = NVERT;
                    end
                    amp_post(i,1:DIM) = SIGMA_Q_IN*randn(1,DIM);
                end
            elseif n == NDIP
                if rand > 2/3 %%% a dipole dies
                    ind_morto = floor((NDIP+1)*rand);
                    if ind_morto < 1
                        ind_morto = 1;
                    elseif ind_morto > NDIP
                        ind_morto = NDIP;
                    end
                    ndip_post(i) = NDIP-1;
                    if ind_morto == 1
                        pos_post(i,:) = [pos_post(i,ind_morto+1:end) 0];
                        amp_post(i,:) = [amp_post(i,ind_morto*DIM+1:end) , zeros(1,DIM)];
                    elseif ind_morto == NDIP
                        pos_post(i,:) = [pos_post(i,1:ind_morto-1), 0];
                        amp_post(i,:) = [amp_post(i,1:(ind_morto-1)*DIM),  zeros(1,DIM)];
                    else
                        pos_post(i,:) = [pos_post(i,1:ind_morto-1),pos_post(i,ind_morto+1:end) 0];
                        amp_post(i,:) = [amp_post(i,1:(ind_morto-1)*DIM),amp_post(i,ind_morto*DIM+1:end) , zeros(1,DIM)];
                    end
                end
            else
                var_al = rand;
                if var_al < 1/3   %%% a dipole is created
                    ndip_post(i) = ndip_post(i)+1;
                    pos_post(i,ndip_post(i) ) = floor((NVERT+1)*rand);
                    if pos_post(i,ndip_post(i) ) <= 0;
                        pos_post(i,ndip_post(i) ) = 1;
                    elseif pos_post(i,ndip_post(i) )> NVERT
                        pos_post(i,ndip_post(i) ) = NVERT;
                    end
                    amp_post(i,DIM*(ndip_post(i)-1)+1:DIM*ndip_post(i)) = SIGMA_Q_IN * randn(1,DIM);
                elseif var_al> 2/3 %%% a dipole dies
                    ind_morto = floor((ndip_post(i)+1)*rand);
                    if ind_morto < 1
                        ind_morto = 1;
                    elseif ind_morto > ndip_post(i)
                        ind_morto = ndip_post(i);
                    end
                    ndip_post(i) = ndip_post(i)-1;
                    if ind_morto == 1
                        pos_post(i,:) = [pos_post(i,ind_morto+1:end) 0];
                        amp_post(i,:) = [amp_post(i,ind_morto*DIM+1:end) , zeros(1,DIM)];
                    else
                        pos_post(i,:) = [pos_post(i,1:ind_morto-1),pos_post(i,ind_morto+1:end) 0];
                        amp_post(i,:) = [amp_post(i,1:(ind_morto-1)*DIM),amp_post(i,ind_morto*DIM+1:end) , zeros(1,DIM)];
                    end
                end
            end
        end

    end
    close(bar_pf);
    toc
        
    ris.mod_sel=mod_sel;
    ris.blob=blob;
    %%% searching for the maximum of the blob
    bar_post=waitbar(0,'please wait','Name','the post processing procedure is running');
    dip = [];
    k = 1;
    tic;
    for t = t1:t2
       waitbar((t-t1+1)/NIST);
        [max_mod, ind_mod] = max(mod_sel(t-t1+1,:)); %selection of the higher probability space
        if ind_mod > 1 
            max_ind_g = ind_mod-1;  G = []; temp_mdd = [];
            j=1;
            b_sel = blob(:,t-t1+1);
            clear aux1 vicini aux2 val ind lista_massimi val_massimi val_max ind_max f aux3 val_f ind_f dip_in index_vertices_in
            for i = 1 : size(b_sel,1)
                aux1 = pf.evol(i,:);
                vicini = aux1(aux1>0);
                aux2 = b_sel(vicini);
                f(i) = sum(aux2);  
                [val ind] = max(aux2);
                if ind == 1
                    lista_massimi(j)=i; 
                    val_massimi(j)=val;
                    j=j+1;
                end
            end
            [val_max,ind_max] = sort(val_massimi,'descend');
            for i=1:max_ind_g
                dip_in(i,:) = pf.vertices(lista_massimi(ind_max(i)),:);
                index_vertices_in(i) = lista_massimi(ind_max(i));
            end
            for i=1:max_ind_g
                G = [G pf.g_matrix(:,DIM*(index_vertices_in(i)-1)+1:DIM*index_vertices_in(i))];
            end
            temp_mdd = pinv(G,pf.tol_pinv)*data(t,:)';
            temp_mdd = temp_mdd';
            for i = 1:size(dip_in,1)
                if DIM ==3
                    dip(k,:)=[dip_in(i,:),temp_mdd((i-1)*DIM+1:i*DIM),t,index_vertices_in(i)];
                elseif DIM ==1
                    dip(k,:)=[dip_in(i,:),temp_mdd((i-1)*DIM+1:i*DIM)*pf.orient(index_vertices_in(i),:),t,index_vertices_in(i)];
                end
              k=k+1;
            end
        end
    end
    toc;
    close(bar_post);
    ris.dip=dip;


end
set(pf_gui.cmd_GO,'enable','on');
clear data
clear g_matrix
if pf.autosave==1
   save(string_ris,'ris','-mat');
   %%%%% clustering and visualization of the results 
   hades_visualization
else
   hades_visualization
end


    



