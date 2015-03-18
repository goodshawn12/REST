function [xyz_o,y_o,z_o]=nut_voxels2bv(xyz_i,y_i,z_i)
% [XYZ_O,Y_O,Z_O]=NUT_VOXELS2BV(X_I,Y_I,Z_I)
% [XYZ_O]=NUT_VOXELS2BV(XYZ_I)
% Takes SPM voxel coordinates and converts them to BrainVisa/Anatomist mm coordinates
%
% Uses global st from SPM.


error('this program is obsolete and does not work correctly!')


global st


% determine MRI dimensions in voxels and mm
R=spm_imatrix(st.vols{1}.mat);
R = spm_matrix([0 0 0 R(4:6)]);
R = R(1:3,1:3);
dim = st.vols{1}.dim*inv(R);
dim_mm = abs(dim .* diag(st.vols{1}.mat(1:3,1:3)*inv(R))');
voxsize_mm = abs(diag(st.vols{1}.mat(1:3,1:3)*inv(R))');

% construct transform matrix
% BrainVisa coordinates begin at opposite corner of volume, using mm
% rather than voxel indices
tfm = -diag(voxsize_mm);
tfm(1:4,4) = [st.vols{1}.dim.*voxsize_mm 1];

switch nargin
	case 3  % if given x,y,z in separate N x 1 vectors
        error('TODO: x,y,z in separate vectors not supported yet. use [x y z] in Nx3 matrix for now...')
% 		xyz=[xyz_i y_i z_i ones(length(z_i),1)]*inv(st.vols{1}.mat)';
% 		xyz_o=xyz(:,1);  % really just x in this case
% 		y_o=xyz(:,2);
% 		z_o=xyz(:,3);
	case 1  % if given N x 3 matrix with x,y,z coords
%        	xyz_o = (st.vols{1}.dim - xyz_i) .* voxsize_mm;
%       	xyz_o_old = (repmat(st.vols{1}.dim,[size(xyz_i,1) 1]) - xyz_i) .* repmat(voxsize_mm,[size(xyz_i,1) 1]);
        xyz_o = nut_coordtfm(xyz_i,tfm);

%         % SANITY CHECK -- make sure new matrix transform is equivalent to
%         % old method
%         errcheck = xyz_o - xyz_o_old;
%         if(any(abs(errcheck(:))>1e-10)) % test will have non-zero entries if check fails
%             warning('brainvisa transform may have failed! double-check results...');
%             errcheck
%         end

    otherwise
        error('You''re doing something wrong.');
end
