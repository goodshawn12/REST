function [xyz_o,y_o,z_o]=nut_mm2voxels(xyz_i,y_i,z_i)
% [XYZ_O,Y_O,Z_O]=NUT_MM2VOXELS(X_I,Y_I,Z_I)
% [XYZ_O]=NUT_MM2VOXELS(XYZ_I)
% Takes SPM mm coordinates and converts them to voxel coordinates
%
% Uses global st.

global st

switch nargin
	case 3  % if given x,y,z in separate N x 1 vectors
		xyz=[xyz_i y_i z_i ones(length(z_i),1)]*inv(st.vols{1}.mat)';
		xyz_o=xyz(:,1);  % really just x in this case
		y_o=xyz(:,2);
		z_o=xyz(:,3);
	case 1  % if given N x 3 matrix with x,y,z coords
       	xyz_o=[xyz_i ones(size(xyz_i,1),1)]*inv(st.vols{1}.mat)';
        xyz_o(:,4)=[];   % destroy last column of ones
	otherwise
        error('You''re doing something wrong.');
end
