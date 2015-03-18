%%
% Obtains voxel orientations from the unit normal of the nearest point on a
% cortical surface mesh.
% @author Daniel D.E. Wong
%%
function nut_findvoxelorientations
global nuts

[file,path] = uigetfile('*_cs.mat','Open Cortical Surface Mesh File');
load(strcat(path,file));

mesh = PrepareTriangleMesh(nut_mri2meg(cs.vertices),cs.faces(:,[1 3 2]));

ori = zeros(size(nuts.voxels));
for vox = 1:size(nuts.voxels,1)
    ori(vox,:) = mesh.un(dsearchn(mesh.mp,nuts.voxels(vox,:)),:);
end

[file,path] = uiputfile('_ori.mat','Save to MAT File');
if ~isequal(file,0) & ~isequal(path,0)
    save(strcat(path,file),'ori');
end