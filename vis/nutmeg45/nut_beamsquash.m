function nut_beamsquash
% squash s_beam into two dimensions

%% coordinates onto right-handed Cartesian plane!!! (this is a neuroaxial view)

monochrome = true;
swapxy = true
plotrect = true
timeidx = 201:601;
pointsup=true

plotbounds = [-50 50 -60 60];

% rectbounds: [-x +x; -y +y]
rectbounds = [-40 40
           10 50];

% src = [24 -30
%        -6 -38
%        -20 -24];
%    
%        obf = [-24 22 40
%                -18 38 40
%                16 28 40]

src = [0 -300000000000000
       -0 -30
       -200000 -240000];

obf = [0 30
      -180000 300008
       160000 280000];


% rect = [-40 10
%         -40 50
%         40  50
%         40 10];
    
rect = [rectbounds(1,1) rectbounds(2,1)
        rectbounds(1,1) rectbounds(2,2)
        rectbounds(1,2)    rectbounds(2,2)
        rectbounds(1,2)   rectbounds(2,1)];


if(swapxy)
    %% CUIDADO!!!! MUY IMPORTANTE!!!!
    %% from this point on, x->y and y->x due to squashing down left-handed MEG head
    src = src(:,[2 1]);
    obf = obf(:,[2 1]);
    rect = rect(:,[2 1]);
    plotbounds = plotbounds(:,[3 4 1 2]);
end

global beam rivets
s_beam = rivets.s_beam;
s_II = rivets.s_II;
% s_beam = abs(s_II);
t = beam.timewindow;
voxels = beam.voxels(:,[2 1]);  % squash down z-axis, swap x and y
[projvoxels,i,projmapping]=unique(voxels,'rows');
projbeam=zeros(size(projvoxels,1),size(s_beam,2));
for i=1:length(projmapping)
    projbeam(projmapping(i),:) = projbeam(projmapping(i),:) + s_beam(i,:);
end

[crap, maxvoxelindex] = max(max(projbeam'));
[crap, maxtimeindex] = max(max(projbeam));

% [x,y,z]=nut_meshthis(projvoxels,projbeam(:,maxtimeindex),2);
[x,y,z]=nut_meshthis(projvoxels,mean(projbeam,2),2);
figure;
subplot(1,2,1);
if(monochrome)
    contour(x,y,z,'k');
else
    contourf(x,y,z); colormap(hot);
end
if(swapxy)
    xobj = xlabel('MEG y-axis (mm)');
    yobj = ylabel('MEG x-axis (mm)');
    axis(plotbounds);
    a=plotbounds;
    text(a(2)+16,a(3),'P','FontName','Times');
    text(a(2)+16,a(4),'A','FontName','Times');
    text(a(1),a(3)-16,'R','FontName','Times');
    text(a(2),a(3)-16,'L','FontName','Times');

else
    xobj = xlabel('MEG x-axis (mm)');
    yobj = ylabel('MEG y-axis (mm)');
end
% axis equal; axis tight;
hold on;
markerobj(1) = plot(src(1,1),src(1,2),'o');
% if(exist('src2','var'))
    markerobj(2) = plot(src(2,1),src(2,2),'s');
%     if(exist('src3','var'))
        markerobj(3)= plot(src(3,1),src(3,2),'d');
%         if(exist('obf1','var'))
            markerobj(4)= plot(obf(1,1),obf(1,2),'o');
%             if(exist('obf2','var'))
                markerobj(5)= plot(obf(2,1),obf(2,2),'o');
%                 if(exist('obf3','var'))
                    markerobj(6)= plot(obf(3,1),obf(3,2),'o');
%                 end
%             end
%         end
%     end
% end

if(plotrect)
    patchobj = patch(rect(:,1),rect(:,2),'white');
    set(patchobj,'FaceColor','none','LineStyle','-.');
end
hold off;
% title(['t=' num2str(beam.timewindow(maxtimeindex))]);

if(monochrome)
    set([gca xobj yobj],'FontName','Times');
else
    %poster style
    set(gca,'LineWidth',2,'FontWeight','bold','FontSize',14);
    set([xobj yobj],'FontWeight','bold','FontSize',14);
end

if(swapxy)
    set(gca,'XDir','reverse');
end

set(gca,'DataAspectRatio',[1 1 1])

% markerobj = get(gca,'Children');
if(monochrome)
%     set(markerobj,'LineWidth',1.5,'MarkerSize',20,'MarkerFaceColor',[0 0 0],'Color',[0 0 0]);
    set(markerobj,'MarkerFaceColor',[0 0 0],'Color',[0 0 0]);
else
    set(markerobj,'LineWidth',1.5,'MarkerSize',20,'MarkerFaceColor',[0 0.7 1],'Color',[0 0 0]);
end


if(pointsup)
    subplot(1,2,2);
%     voxelidx = dsearchn(beam.voxels,[src(:,[2 1]) 40]);
 voxelidx = dsearchn(beam.voxels,[projvoxels(maxvoxelindex,:) 40]);
    plot(t(timeidx),s_II(voxelidx,timeidx),'k'); xobj = xlabel('Time (ms)'); yobj = ylabel('Amplitude');
    if(monochrome)
        set([gca xobj yobj],'FontName','Times');
    else
        set(gca,'LineWidth',2,'FontWeight','bold','FontSize',14);
        set([xobj yobj],'FontWeight','bold','FontSize',14);
        set(get(gca,'Children'),'LineWidth',2);
    end
else
    for ii=1:size(src,1)
        subplot(3,2,2*ii);
        voxelidx = dsearchn(beam.voxels,[src(ii,[2 1]) 40]);
    % maxvoxelindex2 = dsearchn(beam.voxels,[projvoxels(maxvoxelindex,:) 40]);
        plot(t(timeidx),s_II(voxelidx,timeidx),'k'); xobj = xlabel('Time (ms)'); yobj = ylabel('Amplitude');
        if(monochrome)
            set([gca xobj yobj],'FontName','Times');
        else
            set(gca,'LineWidth',2,'FontWeight','bold','FontSize',14);
            set([xobj yobj],'FontWeight','bold','FontSize',14);
            set(get(gca,'Children'),'LineWidth',2);
        end
    end
end