function vol=nut_vector2vol(vector,beam,dimxyz)
% NUT_VECTOR2VOL  transforms data vector from voxel coordinates to a 3D volume.
%
%  vol=nut_vector2vol(vector,beam)
%
% vector    Nx1 data vector
% beam      beam structure with voxel coordinates and coreg info, or Nx3
%           matrix with blob coordinates
% vol       3D volume

if isstruct(beam)
    voxelsblob = nut_voxels2blob(beam);  
else
    voxelsblob = beam;
end
voxelsblob_round=round(voxelsblob);  %this does seem necessary..once in awhile, something is not perfect integer

isround = ( max(max(abs(voxelsblob_round-voxelsblob))) <= 1e-3 );
if ~isround
    if nargin<3, disp('The blob coordinates are not close to integer. Using a linear interpolation.'), end
    % Create an integer image grid, that matches the non-integer coords without leaving
    % holes, and without adding far away voxels.
    gridlow=floor(voxelsblob);
    gridhigh=ceil(voxelsblob);
    gridmid=round(voxelsblob);
    voxelsblob_round=union(gridlow,gridhigh,'rows');
    voxelsblob_round=union(voxelsblob_round,gridmid,'rows');
    voxelsblob_round=voxelsblob_round( all(voxelsblob_round,2) , : );   
    clear grid*
    
    % Linear interpolation
    F = TriScatteredInterp(double(voxelsblob), vector);
    vector = F(double(voxelsblob_round));
end

if nargin<3
    dimxyz = max(voxelsblob_round);
end
vol = nan([dimxyz size(vector,2) size(vector,3)]);
for vv=1:size(voxelsblob_round,1)  % vv num of voxels
    vol(voxelsblob_round(vv,1),voxelsblob_round(vv,2),voxelsblob_round(vv,3),:,:)=vector(vv,:,:);
end


%     % Nearest neighbor interpolation
%     nn=dsearchn(gridall,voxelsblob);
%     for vv=1:size(voxelsblob_round,1)
%         vol(gridall(nn(vv),1),gridall(nn(vv),2),gridall(nn(vv),3))=vector(vv);
%     end
