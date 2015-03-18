%%
% Exports lead potential.
% @author Daniel D.E. Wong
%%
function nut_exportLp

global nuts
LpData.Lp = nuts.Lp;
LpData.voxels = nuts.voxels;
LpData.voxelsize = nuts.voxelsize;

[file,path] = uiputfile('Lp_*.mat','Export Lead Potential');
if ~isequal(file,0) & ~isequal(path,0)
    save(strcat(path,file),'LpData');
end