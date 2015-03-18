%%
% Imports lead potential data exported from NUTEEG.
% @author Daniel D.E. Wong
%%
function nut_importLp

global nuts
[file,path] = uigetfile('Lp_*.mat','Open Lead Potential File');
if isequal(file,0); return; end;
load(strcat(path,file));

nuts.Lp = LpData.Lp;
nuts.voxels = LpData.voxels;
nuts.voxelsize = LpData.voxelsize;

nut_enabler;