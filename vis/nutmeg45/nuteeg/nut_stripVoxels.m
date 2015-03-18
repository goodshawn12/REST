%%
% Removes voxels that are outside the brain volume.
% @author Daniel D.E. Wong
%%
function nut_stripVoxels

global nuts

[file,path] = uigetfile('*_vol.mat','Open Head Surface Mesh File');
if isequal(file,0); return; end;
load(strcat(path,file));

prompt = 'Minimum distance to interior surface (mm):';
title = 'Strip voxels';
lines = 1;
def{1} = '5';
mindist = inputdlg(prompt,title,lines,def);
if isempty(mindist); answer = def; end;
mindist = str2num(cell2mat(mindist));

% Convert to MEG coordinates
for i = 1:length(vol.bnd)
    vol.bnd(i).vertices = nut_mri2meg(vol.bnd(i).vertices);
end

brainmesh=PrepareTriangleMesh(vol.bnd(end).vertices,vol.bnd(end).faces(:,[1 3 2]));

goodvoxels = zeros(size(nuts.voxels,1),1);
count = 0;
for v = 1:size(nuts.voxels,1)
    tri = dsearchn(brainmesh.mp,nuts.voxels(v,:));
    trivox_uv = nuts.voxels(v,:) - brainmesh.mp(tri,:); trivox_uv = trivox_uv/norm(trivox_uv);
    dist = norm(brainmesh.n(tri,:).*dot(brainmesh.n(tri,:),nuts.voxels(v,:)-brainmesh.mp(tri,:))/sum(brainmesh.n(tri,:).^2));
    if dot(brainmesh.un(tri,:),trivox_uv) < 0 & dist > mindist  % voxel is inside the boundary
        count = count+1;
        goodvoxels(count) = v;
    end
end
nuts.voxels = nuts.voxels(goodvoxels(1:count),:);
if isfield(nuts,'Lp')
    nuts.Lp = nuts.Lp(:,:,goodvoxels(1:count));
end