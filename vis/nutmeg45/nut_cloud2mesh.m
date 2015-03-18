function mesh=nut_cloud2mesh(position)
% MESH1 = NUT_CLOUD2MESH(POSITION)
%
% This function is used to covert the cloud to the mesh.
%


pp=pi*[-40:40]/40; 
dbg=0;
im=0;

%Generating mesh slice by slice
for k=min(position(:,3)):(max(position(:,3)))
    %get the index for the non zero entries in the given slice
    in0=find(position(:,3)==k);
    if(length(in0)>0)        
        if(im==0) % upper cap
            im=im+1;
            %find the mesh value for the first slice
            x0=mean(position(in0,1));
            y0=mean(position(in0,2));
            mesh(:,im,1)=x0*ones(size(pp));
            mesh(:,im,2)=y0*ones(size(pp));
            mesh(:,im,3)=k*ones(size(pp));
        end
        
        if(dbg)
            clf;
            plot(position(in0,1),position(in0,2),'r*');
            hold on;
            plot(x0,y0,'ko');
        end
        
        [r,p]=nut_phs(position(in0,1)-x0,position(in0,2)-y0);
        
        if(dbg)
            plot(r.*cos(p)+x0,r.*sin(p)+y0,'b.');
        end
        
        [p,in]=sort(p); r=r(in);
        
        p=[p;p(1)];
        r=[r;r(1)];
        in=find(p(2:end)~=p(1:end-1));
        p=p(in);
        r=r(in);  
        if length(p)==0
            continue
        end
        p=[p(end)-2*pi; p; p(1)+2*pi;];
        r=[r(end); r; r(1); ];
        
        rr=interp1(p,r,pp);
        x=rr.*cos(pp)+x0;
        y=rr.*sin(pp)+y0;
        z=k*ones(1,length(x));
        im=im+1;
        
        if(dbg)
            plot(x,y,'g');
            pause;
        end
        
        mesh(:,im,1)=x;
        mesh(:,im,2)=y;
        mesh(:,im,3)=z;        
    end
end

% lower cap
% find the mesh value for the last slice
im=im+1;
mesh(:,im,1)=x0*ones(size(pp));
mesh(:,im,2)=y0*ones(size(pp));
mesh(:,im,3)=mesh(1,im-1,3)*ones(size(pp));
return;


