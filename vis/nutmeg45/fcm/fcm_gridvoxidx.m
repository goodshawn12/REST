function grid=fcm_gridvoxidx(nuts,jump,doalign)
% FCM_GRIDVOXIDX  returns the indices of grid voxels.
%
%  grid = fcm_gridvoxidx(nuts,grid_spacing)
%
% Each "grid_spacing" voxel is used as grid voxel.

if nargin<3, doalign=0; end

dointerp = any(any( round(nuts.voxels(1:5,:)) - nuts.voxels(1:5,:) > 0.01 ));
if dointerp   % if voxel coordinates are not close to integer, we have to interpolate
    % Make ideal grid
    mx = ceil(1.3 .* max(nuts.voxels));
    mn = floor(1.3 .* min(nuts.voxels));
    mx(1) = max(mx(1),-mn(1));
    mn(1) = -mx(1);
    x = [fliplr(-nuts.voxelsize(1)*jump:-nuts.voxelsize(1)*jump:mn(1)) 0:nuts.voxelsize(1)*jump:mx(1)];
    y = [fliplr(-nuts.voxelsize(2)*jump:-nuts.voxelsize(2)*jump:mn(2)) 0:nuts.voxelsize(2)*jump:mx(2)];
    z = [fliplr(-nuts.voxelsize(3)*jump:-nuts.voxelsize(3)*jump:mn(3)) 0:nuts.voxelsize(3)*jump:mx(3)];
    voxgrid = nut_coordgrid(x,y,z);
    
    % Optional: correct for excentric MRI
    if doalign~=0
        %voxgrid=alignvoxels(voxgrid);
        voxgrid(:,1)=voxgrid(:,1)+doalign;
    end

    % Find voxels that are closest to the ideal grid.
    [near,d] = dsearchn(nuts.voxels,voxgrid);
    f = (d < max(nuts.voxelsize)); %*floor(jump/2)
    grid = near(f)';
else    
    voxels = round(nuts.voxels);
    for k=1:3
        dum{k}=unique(voxels(:,k));
        beg(k)=rem(length(dum{k}),jump); if beg(k)==0, beg(k)=jump; end
        beg(k)=1+floor(beg(k)/2);
    end
    grid=find(ismember(voxels(:,1),[dum{1}(beg(1):jump:end)]) & ismember(voxels(:,2),[dum{2}(beg(2):jump:end)]) & ismember(voxels(:,3),[dum{3}(beg(3):jump:end)]))';
    clear dum beg
end
