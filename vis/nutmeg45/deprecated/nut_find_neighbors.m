function [neighbors,searchindices,distances_out]=nut_find_neighbors(voxels,coord,inradius,exradius)
% function [neighbors,searchindices,distances_out]=nut_find_neighbors(voxels,coord,inradius,exradius)

[MEGvoxelindex,distances] = dsearchn(voxels,coord);
distances = nut_rownorm(nut_coord_diff(voxels,voxels(MEGvoxelindex,:)));
searchindices1 = find(distances >= exradius);
searchindices2 = find(distances <= inradius);
searchindices = intersect(searchindices1,searchindices2);
neighbors=voxels(searchindices,:);
distances_out=distances(searchindices);




