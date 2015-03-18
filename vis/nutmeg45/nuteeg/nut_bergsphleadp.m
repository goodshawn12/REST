%%
% Calculates the lead potentials for each voxel using the spherical head model.
% voxels - a list of voxel coordinates (Nvoxels x 3 matrix) [mm]
% sensorcoord - the coordinates of the sensor (1x3 matrix) [mm]
% lsc - the coordinates of the local sphere center (1x3 matrix) [mm]
% @return lx1 lead field for unit dipole in x direction. (1 x Nvoxels) [V]
% @return ly1 lead field for unit dipole in y direction. (1 x Nvoxels) [V]
% @return lz1 lead field for unit dipole in z direction. (1 x Nvoxels) [V]
% @author Daniel D.E. Wong
%%
function [lx1,ly1,lz1] = nut_bergsphleadp(voxels,sensorcoord,lsc,vol)

if nargin == 3
    %r(1) = 95.2;
    %r(2) = 88.6;
    %r(3) = 76.2;
    r(1) = norm(sensorcoord-lsc);
    r(2) = 0.93*r(1);
    r(3) = 0.85*r(1);
    ci = [0.33 0.0132 0.33];
else    % 4
    if isfield(vol,'r')
        r = vol.r;
    else
        % Assume surfaces are fairly concentric locally
        r = zeros(1,length(vol.bnd));
        for i = 1:length(r)
            r(i) = norm(vol.bnd(i).vertices(dsearchn(vol.bnd(i).vertices,sensorcoord),:)-lsc);
        end
    end
    ci=vol.cond;
end

% Project the sensor coordinate to the sphere surface
rs = sensorcoord-lsc; rs = rs/norm(rs);
sensorcoord = lsc + rs*r(1);

[mu,lambda] = nut_findbergparams(r,ci);
lx1 = zeros(size(voxels,1),1); ly1 = lx1; lz1 = lx1;
rq = voxels-repmat(lsc,size(voxels,1),1);
for k = 1:length(mu)
    if exist('vol','var')
        [lx1_,ly1_,lz1_] = nut_1lyrsphleadp(rq*mu(k),sensorcoord,[0 0 0],vol);
    else
        [lx1_,ly1_,lz1_] = nut_1lyrsphleadp(rq*mu(k),sensorcoord,[0 0 0]);
    end
    lx1 = lx1 + lambda(k)*lx1_;
    ly1 = ly1 + lambda(k)*ly1_;
    lz1 = lz1 + lambda(k)*lz1_;
end