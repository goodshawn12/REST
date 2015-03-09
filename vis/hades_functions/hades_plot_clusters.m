function hades_plot_clusters(struct_dipoli_cluster,waveforms)

%%% HADES_PLOT_CLUSTERS function for plotting the clusters

% Copyright (C) 2010, Cristina Campi, Annalisa Pascarella, Michele Piana, Alberto Sorrentino
global ris;
global gui_vis;
global pf;
global fwd;

n_c=size(struct_dipoli_cluster,2);
if isfield(ris,'time_interval')==0
    T=ris.times;
else
    T=ris.time_interval;
end
str_cluster=get(gui_vis.str_vis,'string');
if strcmp(str_cluster,'')==1 
    c1=1:1:n_c;
                min_c=min(c1);
            max_c=max(c1);
else
    %%% vector with the indices of clusters to plot
    c1=str2num(str_cluster);
    if 0<norm(floor(c1)-c1)  
        msgbox('You have to insert integer number','Error','error');
        return
    else
        if isempty(c1)==1
            min_c=1;
            max_c=n_c;
            c1=min_c:1:max_c;
        else
            min_c=min(c1);
            max_c=max(c1);
        end
    end
end
if min_c<=0 && max_c<=n_c
     msgbox('The clusters enumeration starts from 1!!','Error','error');
elseif 0<min_c && n_c<max_c
     msgbox('You have chosen too many cluster!!','Error','error');
elseif min_c<=0 && n_c<max_c
     msgbox('You have selected wrong clusters!!','Error','error');
elseif 0<min_c && max_c<=n_c
    n_color=round(n_c/12)+1;
    paletta_base=[0 0 1;  1 0 0;  0 1 0;  1 0 1;  0 1 1;  1 1 0;   .7 .3 .3;  .3 .3 .8; 1 0.5 0; 0.5 0 1; 1 0 0.5;  0.5 0.5 1 ];
    paletta=repmat(paletta_base,n_color,1);
    simboli_base=[repmat('o ',12,1);repmat('x ',12,1);repmat('v ',12,1);repmat('p ',12,1)];
    simboli=repmat(simboli_base,n_color,1);
    fig_cluster=figure;
    for t=1:size(ris.subject,2)
        if ris.subject(t)=='_'
            string_case(t)=' ';
        else
         string_case(t)=ris.subject(t);
        end
    end
    set(fig_cluster,'Name',strcat(['Case: ',string_case,', clustered dipoles ',num2str(c1)]));
    subplot(2,2,1)
    %%title(['Case: ',string_case,', clusters ',num2str(c1)])
    for t = 1:size(c1,2)
        i=c1(t);
        A = [];
        A = struct_dipoli_cluster(i).dipoli_cluster;
        hold on
        plot3(A(:,1),A(:,2),A(:,3),simboli(i),'Linewidth',5,'Color',paletta(i,:));
    end
    plot3(pf.vertices(:,1),pf.vertices(:,2),pf.vertices(:,3),'k.','markersize',1)
    if isempty(fwd.name_sensors)==0 && isempty(fwd.name_versors)==0
        xs=load(fullfile(fwd.location_sensors,fwd.name_sensors));
        plot3(xs(:,1),xs(:,2),xs(:,3),'bo')
    end
    view(0, 90)
    set(gca,'XTickLabel','');
    set(gca,'YTickLabel','');
    set(gca,'ZTickLabel','');
    axis tight
    subplot(2,2,2)
    for t = 1:size(c1,2)
        i=c1(t);
        A = [];
        A = struct_dipoli_cluster(i).dipoli_cluster;
        hold on
        plot3(A(:,1),A(:,2),A(:,3),simboli(i),'Linewidth',5,'Color',paletta(i,:));
    end
    plot3(pf.vertices(:,1),pf.vertices(:,2),pf.vertices(:,3),'k.','markersize',1)
    if isempty(fwd.name_sensors)==0 && isempty(fwd.name_versors)==0
        xs=load(fullfile(fwd.location_sensors,fwd.name_sensors));
        plot3(xs(:,1),xs(:,2),xs(:,3),'bo')
    end
    view(90,0)
    set(gca,'XTickLabel','');
    set(gca,'YTickLabel','');
    set(gca,'ZTickLabel','');
    axis tight
    subplot(2,2,3)
    for t = 1:size(c1,2)
        i=c1(t);
        A = [];
        A = struct_dipoli_cluster(i).dipoli_cluster;
        hold on
        plot3(A(:,1),A(:,2),A(:,3),simboli(i),'Linewidth',5,'Color',paletta(i,:));
    end
    plot3(pf.vertices(:,1),pf.vertices(:,2),pf.vertices(:,3),'k.','markersize',1)
    if isempty(fwd.name_sensors)==0 && isempty(fwd.name_versors)==0
        xs=load(fullfile(fwd.location_sensors,fwd.name_sensors));
        plot3(xs(:,1),xs(:,2),xs(:,3),'bo')
    end
    view(0,0)
    set(gca,'XTickLabel','');
    set(gca,'YTickLabel','');
    set(gca,'ZTickLabel','');
    axis tight
    subplot(2,2,4)
    for t = 1:size(c1,2)
        i=c1(t);
        plot(1000*T(1:size(waveforms,2)),10^9*waveforms(i,:),'Marker',simboli(i),'Color',paletta(i,:));
        hold on
    end  
    hold on
    xlabel('time [ms]');
    ylabel('source waveforms [nAm]');
    axis tight
    axis([1000*T(ris.t1),1000*T(ris.t2),0,max(max(10^9*waveforms))]);
end

