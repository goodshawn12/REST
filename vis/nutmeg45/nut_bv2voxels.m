function [xyz_o,y_o,z_o]=nut_bv2voxels(xyz_i,y_i,z_i)
% [XYZ_O,Y_O,Z_O]=NUT_BV2VOXELS(X_I,Y_I,Z_I)
% [XYZ_O]=NUT_BV2VOXELS(XYZ_I)
% Takes BrainVisa/Anatomist mm voxel coordinates and converts them to SPM voxel coordinates
%
% Uses global st from SPM.

error('this program is obsolete and does not work correctly!')

global st

mat = st.vols{1}.mat;
R=spm_imatrix(st.vols{1}.mat);
R = spm_matrix([0 0 0 R(4:6)]);
R = R(1:3,1:3);
dim = st.vols{1}.dim*inv(R);
dim_mm = abs(dim .* diag(st.vols{1}.mat(1:3,1:3)*inv(R))');
voxsize_mm = abs(diag(st.vols{1}.mat(1:3,1:3)*inv(R))');

% construct transform matrix
% BrainVisa coordinates begin at opposite corner of volume, using mm
% rather than voxel indices
% note this is voxels2bv tfm, so we apply inv(tfm) later for bv2voxels
tfm = -diag(voxsize_mm);
tfm(1:4,4) = [st.vols{1}.dim.*voxsize_mm 1];

switch nargin
	case 3  % if given x,y,z in separate N x 1 vectors
        error('TODO: x,y,z in separate vectors not supported yet. use [x y z] in Nx3 matrix for now...')
	case 1  % if given N x 3 matrix with x,y,z coords
        xyz_o = nut_coordtfm(xyz_i,inv(tfm));
    otherwise
        error('You''re doing something wrong.');
end

