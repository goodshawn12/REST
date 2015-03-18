function [xyz_o,y_o,z_o]=nut_nmag2meg(xyz_i,y_i,z_i)
% [XYZ_O]=NUT_NMAG2MEG(XYZ_I)
% [X_O,Y_O,Z_O]=NUT_NMAG2MEG(X_I,Y_I,Z_I)
% Takes NeuroMag coords (mm) and converts them to CTF/BTi MEG head space coordinates (mm)

switch nargin
	case 3  % if given x,y,z in separate N x 1 vectors
		% [xyz_o,y_o,z_o]=nut_coordtfm(xyz_i,y_i,z_i,inv(coreg.meg2mri_tfm));
        xyz_o = y_i;
        y_o = -x_i;
        z_o = z_i;
    case 1  % if given N x 3 matrix with x,y,z coords
        xyz_o = [xyz_i(:,2) -xyz_i(:,1) xyz_i(:,3)];
    otherwise
        error('Incorrect number of input arguments.');
end
