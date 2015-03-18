function [xyz_o,y_o,z_o]=nut_meg2mri(xyz_i,y_i,z_i,coreg)
% [XYZ_O]=NUT_MEG2MRI(XYZ,COREG)
% [X_O,Y_O,Z_O]=NUT_MEG2MRI(X,Y,Z,COREG)
%
% Takes MEG head space coordinates (mm) and converts them to MRI coords (mm)
%
% X,Y,Z,XYZ     MEG voxel coordinates
% COREG         (optional) coregistration information. 

if rem(nargin,2)  % if 1 or 3 input params
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
elseif nargin==2
    coreg = y_i; 
end

if nargin>2  % if given x,y,z in separate N x 1 vectors
    [xyz_o,y_o,z_o]=nut_coordtfm(xyz_i,y_i,z_i,coreg.meg2mri_tfm);
else  % if given N x 3 matrix with x,y,z coords
    xyz_o=nut_coordtfm(xyz_i,coreg.meg2mri_tfm);
end