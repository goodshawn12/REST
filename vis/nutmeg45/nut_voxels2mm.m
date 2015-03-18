function [xyz_o,y_o,z_o]=nut_voxels2mm(xyz_i,y_i,z_i,tfmat)
% [X_O,Y_O,Z_O]=NUT_VOXELS2MM(X,Y,Z,TFMAT)
% [XYZ_O]=NUT_VOXELS2MM(XYZ,TFMAT)
%
% Takes SPM voxel coordinates and converts them to mm coordinates
%
% X,Y,Z,XYZ  SPM voxel coordinates
% TFMAT      (optional) transformation matrix. If this input is not given,
%            data needs to be displayed in nut_results viewer.

if rem(nargin,2)
    global st
    tfmat = st.vols{1}.mat;
elseif nargin==2
    tfmat = y_i; 
end

if nargin>2  % if given x,y,z in separate N x 1 vectors
    xyz=[xyz_i y_i z_i ones(length(z_i),1)] * tfmat';
    xyz_o=xyz(:,1);  % really just x in this case
    y_o=xyz(:,2);
    z_o=xyz(:,3);
else  % if given N x 3 matrix with x,y,z coords
    xyz_o=[xyz_i ones(size(xyz_i,1),1)] * tfmat';
    xyz_o(:,4)=[];   % destroy last column of ones
end

