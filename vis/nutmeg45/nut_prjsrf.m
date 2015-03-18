function [prjpts_mm]=nut_prjsrf(mesh,hsCoords_mri_mm)
% [PRJPTS_MM, PRJPTS] = 
%	NUT_PRJSRF(MESH1, hsCoords_mri, ORIGIN, VOXELSIZE)
% This function project the points to the Head surface.
%


nn=nut_normsrf(mesh);
[n1,n2,n3]=size(mesh);
mesh=reshape(mesh,n1*n2,n3);
nn=reshape(nn,n1*n2,3);
%h = waitbar(0,'Please wait, Projecting to surface...');
for k=1:size(hsCoords_mri_mm,1)    
	r0=hsCoords_mri_mm(k,:);
	[tmp,in]=min((mesh(:,1)-r0(1)).^2 + (mesh(:,2)-r0(2)).^2 + (mesh(:,3)-r0(3)).^2);

	rn=mesh(in,:);
	n=nn(in,:);
	
	prjpts_mm(k,:)=(n*(rn-r0)')*n+r0;
	%waitbar(k/(size(hsCoords_mri_mm,1)),h);
end
%converting voxels to mm
%prjpts_mm = nut_voxels2mm(prjpts);

% close(h);
return
