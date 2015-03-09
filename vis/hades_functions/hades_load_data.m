%   Copyright (C) 2010, Cristina Campi, Annalisa Pascarella, Michele Piana, Alberto Sorrentino

function[]=hades_load_data
global b_fwd;
global b;
global pf;
scrsz = get(0,'ScreenSize');
b_fwd.fig=figure('HandleVisibility','Callback','Menubar','none',...
    'Name','Loading the dataset ', 'NumberTitle','off','Visible','on', 'BackingStore','off');
    pos=get(b_fwd.fig,'pos');
    pos(1,3)=300;
    pos(1,4)=80;
    pos(1,1)=(scrsz(1,3)-pos(1,3))/2;
    pos(1,2)=(scrsz(1,4)-pos(1,4))/2;
    set(b_fwd.fig,'pos',pos); 

    clear comment
    comment=[];
    for t=1:b.dataset_tot
        if t==1
            comment=[comment,b.fif_tot.evoked(1,t).comment];
        else
            comment=[comment,'|',b.fif_tot.evoked(1,t).comment];
        end
    end
    b_fwd.label= uicontrol(b_fwd.fig,'style', 'text',...
        'units', 'pixel',...
        'string', 'Select the dataset',...
         'position', [10 40 100 20]);
     if isempty(pf.dataset_number)==1
         val=1;
     else
         val=pf.dataset_number;
     end
    b_fwd.popup=uicontrol(b_fwd.fig,'style', 'popup',...
        'units', 'pixel',...
        'String', comment,...
        'value',val,...
        'enable','on',...
         'position', [110 30 190 30],...
        'callback','hades_main(''set_dataset_fif'');');
end