function nut_show_head(mesh)
% NUT_SHOW_HEAD(MESH1)
%
% Show the mesh in 3D.

% light=[258 256 256]; 
% xs=light(1);ys=light(2);zs=light(3);
% r0(:,:,1)=mesh(:,:,1)-xs;
% r0(:,:,2)=mesh(:,:,2)-ys;
% r0(:,:,3)=mesh(:,:,3)-zs;
% % find the shading for the surface lighting
% meshtemp = reshape(mesh,size(mesh,1)*size(mesh,2),3);
% meshtemp = nut_voxels2mm(meshtemp);
% mesh = reshape(meshtemp,size(mesh));
% 
% shd=dot(nut_normsrf(mesh),r0,3)./sqrt(sum(r0.^2,3));
%surf(mesh(:,:,1),mesh(:,:,2),mesh(:,:,3),shd);
surfl(mesh(:,:,1),mesh(:,:,2),mesh(:,:,3),[10 -60]);  %fourth arg here is light source direction in degrees pointing down on head
shading interp; 
% colormap(gray);
axis equal; 
view([232.5,10]); 
grid off;
sepia = gray;
sepia = [sepia(:,1) .9*sepia(:,2) .8*sepia(:,3)];
colormap(sepia);
return
