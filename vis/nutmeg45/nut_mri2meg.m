function [xyz_o,y_o,z_o]=nut_mri2meg(xyz_i,y_i,z_i)
% [XYZ_O]=NUT_MRI2MEG(XYZ_I)
% [X_O,Y_O,Z_O]=NUT_MRI2MEG(X_I,Y_I,Z_I)
% Takes MRI coords (mm) and converts them to MEG head space coordinates (mm)

global beam

if(isfield(beam,'coreg'))
    coreg=beam.coreg;
else
    global nuts
    if(isfield(nuts,'coreg'))
        coreg = nuts.coreg;
    else
        error('no coregistration info found in beam or nuts');
    end
end

switch nargin
	case 3  % if given x,y,z in separate N x 1 vectors
		[xyz_o,y_o,z_o]=nut_coordtfm(xyz_i,y_i,z_i,inv(coreg.meg2mri_tfm));
	case 1  % if given N x 3 matrix with x,y,z coords
		xyz_o=nut_coordtfm(xyz_i,inv(coreg.meg2mri_tfm));
	otherwise
        error('Incorrect number of input arguments.');
end
