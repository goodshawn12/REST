%%
% Removes voxels that are at least a certain distance away from the
% cortical surface.
% @author Daniel D.E. Wong
%%
function nut_stripcsvoxels
global nuts

[file,path] = uigetfile('*_cs.mat','Open Cortical Surface Mesh File');
if isequal(file,0); return; end;
load(strcat(path,file));

prompt = 'Strip voxels outside this range (mm):';
title = 'Strip voxels';
lines = 1;
def{1} = '5';
mindist = inputdlg(prompt,title,lines,def);
if isempty(mindist); answer = def; end;
mindist = str2num(cell2mat(mindist));

brainmesh=PrepareTriangleMesh(nut_mri2meg(cs.vertices),cs.faces(:,[1 3 2]));

goodvoxels = zeros(size(nuts.voxels,1),1);
count = 0;
for v = 1:size(nuts.voxels,1)
    tri = dsearchn(brainmesh.mp,nuts.voxels(v,:));
    if norm(brainmesh.mp(tri,:)-nuts.voxels(v,:)) < mindist
        count = count+1;
        goodvoxels(count) = v;
    end
end

nuts.voxels = nuts.voxels(goodvoxels(1:count),:);
nuts.Lp = nuts.Lp(:,:,goodvoxels(1:count));