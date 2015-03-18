function O = nut_smoothvol(V,beam,vFWHM)
% NUT_SMOOTHVOL  performs 3D smoothing on functional data.
%
%  Usage:
%     O = nut_smoothvol(V,beam,vFWHM)
%
%  V     functional data (voxels*time*frequency)
%  beam  beam structure with voxel coordinates and voxelsize
%  fFHWM width of Gaussian smoothing kernel.
%
% Calls spm_smooth.

if nargin<3 || isempty(vFWHM)
    vFWHM = [20 20 20];
elseif isscalar(vFWHM)
    vFWHM = vFWHM .* ones(1,3);
end

[nv,nt,nf] = size(V);
O = zeros(nv,nt,nf);
for t=1:nt
    for f=1:nf
        vol = nut_vector2vol(V(:,t,f),beam);
        volM = double(isfinite(vol) & vol>0);
        dimxyz = size(vol);            
        outvol = zeros(dimxyz);
        outvolM= outvol;
        spm_smooth(volM,outvolM,vFWHM./beam.voxelsize);
        % for some reason, SPM changes outvol AND outvolM, so we have to
        % save it temporarily!
        save temp outvolM   
        spm_smooth(vol,outvol,vFWHM./beam.voxelsize);
        load temp, delete temp.mat
        O(:,t,f) = nut_vol2vector(outvol./outvolM,beam);
    end
end