function[evol]=hades_compute_neighbours(radius,vertices,dir,subject)

%%% HADES_COMPUTE_NEIGHBOURS computes the neighbours for the points of the source space. 
%%% It can be employed for pre-computing the neighbours matrix, running in the Matlab Command Window 
%%% the following command:
%%% [neighbours]=hades_compute_neighbours(radius,vertices,dir,subject)
%%% where
%%% neighbours   output neighbours matrix;
%%% radius       radius [cm] within compute the neighbours;
%%% vertices     the coordinates of the points of the source space, stored in a \textit{nvert}$\times$3 matrix, \textit{nvert} number of points of the source space;
%%% dir          the directory where save the matrix
%%% subject 	 name of the subject.
%%% The output neighbours matrix is automatically saved in a .dat file with name 
%%% `<subject>_neigbours_<radius>.dat

%   Copyright (C) 2010, Cristina Campi, Annalisa Pascarella, Michele Piana, Alberto Sorrentino

evol=[];
radius = radius*10^-2;
bar_evol=waitbar(0,'please wait','Name','Computing the evolution');
for i = 1 : size(vertices,1);
    waitbar(i/size(vertices,1));
    clear aux3
    aux1 = vertices(i,:);
    aux2 = repmat(aux1,size(vertices,1),1);
    diff = vertices-aux2;
    diff = diff.^2;
    somma = sum(diff');
    somma = somma';
    dist = sqrt(somma);
    [dist_riord,ind_riord] = sort(dist,'ascend');
    aux3 = find(dist_riord<radius);
    evol(i,1:size(aux3,1)) = ind_riord(aux3)';
end

check_neighbours=zeros(size(evol,1),1);
n=1;
a=evol(n,1);
check_neighbours(a,1)=1;
a_next=evol(a,:);
a_next=a_next(find(a_next>0));
aux_norm=-100;
check_neighbours(a_next',1)=1;
while(n<=size(evol,1)) 
    a=evol(n,1);
    if check_neighbours(a,1) ==0
        n=n+1;
    else
        check_neighbours(a,1)=1;
        a_next=evol(a,:);
        a_next=a_next(find(a_next>0));
        check_neighbours(a_next',1)=1;
        n=n+1;
    end
    if n==size(evol,1)
        if (norm(check_neighbours-ones(size(check_neighbours))))==0
            break
        elseif aux_norm==norm(check_neighbours-ones(size(check_neighbours)))
            break
        else
            aux_norm=norm(check_neighbours-ones(size(check_neighbours)));
           n=1; 
        end
    end
end

if (norm(check_neighbours-ones(size(check_neighbours))))==0
    save(strcat(fullfile(strtrim(dir),strtrim(subject)),'_neighbours_',num2str(radius*10^2),'.mat'),'evol','-mat');
    close(bar_evol);
    else
    choice_save_neighbours = questdlg('The source space is not globally conected by the neighbours matrix. Do you want to save the matrix?', ...
            'Saving neighbours matrix...', 'yes','no','no');
        switch choice_save_neighbours
        case 'yes'
            save(strcat(fullfile(strtrim(dir),strtrim(subject)),'_neighbours_',num2str(radius*10^2),'.mat'),'evol','-mat');
        case 'no'
            close all
        end
end