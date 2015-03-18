function [xyz_o,y_o,z_o] = nut_coordtfm(varargin)
% [XYZ_O]=NUT_COORDTFM(XYZ_I,A)
% [X_O,Y_O,Z_O]=NUT_COORDTFM(X_I,Y_I,Z_I,A)
% Performs coordinate transformation using transform matrix A.
% Supply either xyz_i as an N x 3 matrix containing cartesian coordinates
% or x_i, y_i, and z_i in separate N x 1 vectors.
% A should be 4 x 4.


switch nargin
	case 4  % if given x,y,z in separate N x 1 vectors
        [x_i,y_i,z_i,A] = deal(varargin{:});

        xyz=[x_i y_i z_i ones(length(z_i),1)]*A';
		xyz_o=xyz(:,1);  % really just x in this case
		y_o=xyz(:,2);
		z_o=xyz(:,3);
        if(~isequal(xyz(:,4),ones(size(xyz(:,4)))))% if fourth column is not a bunch of ones, something went wrong
            warning('something may have gone awry with this coordinate transform.');
        end
    case 2  % if given N x 3 matrix with x,y,z coords
       	[xyz_i,A] = deal(varargin{:});
        xyz_o=[xyz_i ones(size(xyz_i,1),1)]*A';
%         if(~isequal(xyz_o(:,4),ones(size(xyz_o(:,4)))))% if fourth column is not a bunch of ones, something went wrong
        if(norm(xyz_o(:,4)-1,'fro')>1e-13)% if fourth column is not a bunch of ones, something went wrong
            warning('something may have gone awry with this coordinate transform.');
        end
        xyz_o(:,4)=[];   % destroy last column of ones
	otherwise
        error('You''re doing something wrong.');
end
