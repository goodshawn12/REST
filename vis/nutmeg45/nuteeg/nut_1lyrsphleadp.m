%%
% Calculates the lead potentials for each voxel using the spherical head
% model. The sphere radius is taken to be the distance between the sensor
% and the sphere center. Conductivity is taken as the conductivity of the
% innermost layer provided by the vol structure.
% @param voxels - a list of voxel coordinates (Nvoxels x 3 matrix) [mm]
% @param sensorcoord - the coordinates of the sensor (1x3 matrix) [mm]
% @param lsc - the coordinates of the local sphere center (1x3 matrix) [mm]
% @param vol - the head volume structure.
% @return lx1 lead field for unit dipole in x direction (1 x Nvoxels) [T]
%%
function [lx1,ly1,lz1] = nut_1lyrsphleadp(voxels,sensorcoord,lsc,vol)

rq = norm(sensorcoord-lsc);
if nargin == 3
    ci = [0.33 0.0132 0.33];
else    % 4
    ci=vol.cond;
end

lx1 = zeros(1,size(voxels,1)); ly1 = lx1; lz1 = lx1;

% Convert coordinates with respect to lsc
r = (sensorcoord-lsc)/1000; %[m]
rn = norm(r);
r = repmat(r,size(voxels,1),1);
rq = (voxels-repmat(lsc,size(voxels,1),1))/1000;    %[m]
rqn = nut_rownorm(rq);
q = repmat(1e-9,size(voxels,1),1);

% Method of images
im = find(rqn > rn);
if length(im)
    rqno = rqn(im);             % Backup original rqn
    rqn(im) = rn^2./rqn(im);    % Scale source radius
    rq(im,:) = rq(im,:)./repmat(rqno,1,3).*repmat(rqn(im),1,3); % Update actual source radial vector
    q(im) = q(im).*rn./rqno; % Scale source moment
end

d = r-rq;
dn = nut_rownorm(d);
F = dn.*(rn.*dn + rn.^2 - dot(rq,r,2));

c1 = 1/4/pi/ci(end)./rqn.^2 .* (2*dot(d,rq,2)./dn.^3+1./dn-1./rn);
c2 = 1/4/pi/ci(end)./rqn.^2 .* (2./dn.^3 + (dn+rn)./rn./F);

q = 1e-9;   % We would like a 1 nAm dipole
v = (repmat(c1-c2.*dot(r,rq,2),1,3).*rq + repmat(rqn.^2,1,3).*repmat(c2,1,3).*r)*q;

lx1 = v(:,1); ly1 = v(:,2); lz1 = v(:,3);