function mesh=nut_mesh_smoothing(mesh)
% MESH1 = NUT_MESH_SMOOTHING(MESH1)
%
% This function is used to smooth the mesh.
%

mesh(3:end-2,:,1:2)= mesh(3:end-2,:,1:2)/2 + ...
	              mesh(2:end-3,:,1:2)/4 + ...
                      mesh(4:end-1,:,1:2)/4;  

mesh(:,3:end-2,1:2)= mesh(:,3:end-2,1:2)/2 + ...
                      mesh(:,2:end-3,1:2)/4 + ...
                      mesh(:,4:end-1,1:2)/4;
return
