%%
% Generates a tri-mesh sphere and expands/contracts it to fit a surface.
% This function was created as a solution to BrainStorm meshes having
% overlapping faces after reducepatch is applied. The sphere winding will
% be clockwise.
% @param nfv the surface to be fitted.
% @param r the radius of the initial sphere. 300 [mm] usually works well.
% @param f the number of faces the sphere should be reduced to. Should be
% less than 5120.
% @return sph the shrink-fitted sphere.
% @author Daniel D.E. Wong
%%
function sph = shrinkfitsphere(nfv,r,f)
warning('off')
if f > 5120; f = 5120; end;
sph = reducepatch(sphere_tri('ico',4,r),f);   % reducepatch seems to work fine with spheres

for i = 1:size(sph.vertices,1)
    for j = 1:size(nfv.faces,1)
        isect = triangle_intersect(nfv.vertices(nfv.faces(j,:),:),sph.vertices(i,:));
        if ~isempty(isect)
            sph.vertices(i,:) = isect;
            break;
        end
    end
end


%%
% Determines whether the line between a point and the origin intersects the
% specified triangle.
% @param tri_pnts the points of the triangle.
% @param pnt the point used for the line-plane intersection.
% @return isect the coordinates of the intersection if one exists.
% Otherwise an empty matrix is returned.
%%
function isect = triangle_intersect(tri_pnts, pnt)

tuv = inv([pnt' tri_pnts(2,:)'-tri_pnts(1,:)' tri_pnts(3,:)'-tri_pnts(1,:)']) * (pnt'-tri_pnts(1,:)');
if tuv>= 0 & tuv<=1 & tuv(2)+tuv(3)<=1
    isect = pnt*(1-tuv(1));
else
    isect = [];
end