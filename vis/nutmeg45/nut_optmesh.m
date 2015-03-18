function mesh=nut_optmesh(mesh)
% MESH1 = NUT_OPTMESH(MESH1)
%
% Optimizes the already generated mesh
% SYNTAX: mesh=nut_optmesh(mesh)
%

DEBUG = false;

% reduction in z, caps excluded
ttt=(mesh(:,4:end-1,1:2)+mesh(:,2:end-3,1:2))*0.5 - mesh(:,3:end-2,1:2);
ttt=sqrt(sum(sum(ttt.^2,1),3));
ttt=conv2(ttt,[1 2 1]/4,'same');
thld=0.5*max(ttt);
in=find(ttt<thld);
itmp=[1:length(ttt)];
itmp(in(1:2:end))=0;
in=find(itmp);

if(DEBUG)
    figure;
    subplot(2,1,1);
    plot(ttt);
    hold on;
    plot(in,thld*ones(size(in)),'r*'); 
    title(['z, retained: ' num2str(length(in)/length(ttt)) ]);
end
mesh=mesh(:,[1 2+in end],:);

% reduction in xy
ttt=(mesh(3:end,:,1:2)+mesh(1:end-2,:,1:2))*0.5 - mesh(2:end-1,:,1:2);
ttt=sqrt(sum(sum(ttt.^2,2),3));
ttt=conv2(ttt,[1 2 1]/4,'same');
thld=0.25*max(ttt);
in=find(ttt<thld);
itmp=[1:length(ttt)];
itmp(in(1:2:end))=0;
in=find(itmp);
if(DEBUG)
    subplot(2,1,2);
    plot(ttt);
    hold on;
    plot(in,thld*ones(size(in)),'r*');
    title(['xy, retained: ' num2str(length(in)/length(ttt)) ]);
end

mesh=mesh([1 1+in end],:,:);

return
