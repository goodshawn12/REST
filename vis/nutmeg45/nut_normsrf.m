function n0=nut_normsrf(mesh)
% N0 = NUT_NORMSRF(MESH1)
% Finds the normal to the surface
%

mxy=[mesh(end,:,:);mesh;mesh(1,:,:)];
mz=[mesh(:,1,:) mesh mesh(:,end,:)];
n0=cross(mxy(3:end,:,:)-mxy(1:end-2,:,:),mz(:,3:end,:)-mz(:,1:end-2,:) );
n0(:,end,3)=1;
n0(:,1,3)=-1;

rr=sqrt(sum(n0.^2,3));

n0(:,:,1)=n0(:,:,1)./rr;
n0(:,:,2)=n0(:,:,2)./rr;
n0(:,:,3)=n0(:,:,3)./rr;

return;
