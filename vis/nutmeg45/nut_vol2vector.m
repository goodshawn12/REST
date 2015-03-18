function [vector,voxelsblob]=nut_vol2vector(vol,beam)
% NUT_VOL2VECTOR  transforms a 3D data volume to a vector.
%
%  [vector,coordinates]=nut_vol2vector(vol,¦beam¦)
%
% vol           3D data volume
% beam          (optional) beam structure with voxel coordinates or Nx3 matrix
%               with blob coordinates
% vector        Nx1 data vector
% coordinates   Nx3 blob coordinates

if nargin==1 && nargout==1
    vector = vol(:);
    return
end

if nargin<2
    voxelsblob = nut_coordgrid(1:size(vol,1),1:size(vol,2),1:size(vol,3));
else
    if isstruct(beam)
        voxelsblob = nut_voxels2blob(beam);  
    else
        voxelsblob = beam;
    end
end
voxelsblob_round=round(voxelsblob);  %this does seem necessary..once in awhile, something is not perfect integer

if ( max(max(abs(voxelsblob_round-voxelsblob))) > 1e-3 )
    error('Non integer voxel coordinates are currently not supported.')
end
numvox = size(voxelsblob_round,1);

if nargin<2 && ndims(vol)==3
    vector = vol(:); % this is much faster
else                 % but in this case we are not sure that the coordinates are in good order so we have to do it slow.
    vector = zeros(numvox,size(vol,4),size(vol,5));
    for vv=1:numvox 
        vector(vv,:,:) = vol(voxelsblob_round(vv,1),voxelsblob_round(vv,2),voxelsblob_round(vv,3),:,:);
    end
end
