%%
% Determines whether the line between a point and the origin intersects the
% specified triangle.
% @param tri_pnts the points of the triangle.
% @param pnt the point used for the line-plane intersection.
% @return isect the coordinates of the intersection if one exists.
% Otherwise an empty matrix is returned.
% @author Daniel D.E Wong
%%
function isect = nut_triangle_intersect(tri_pnts, pnt)

tuv = inv([pnt' tri_pnts(2,:)'-tri_pnts(1,:)' tri_pnts(3,:)'-tri_pnts(1,:)']) * (pnt'-tri_pnts(1,:)');
if tuv>= 0 & tuv<=1 & tuv(2)+tuv(3)<=1
    isect = pnt*(1-tuv(1));
else
    isect = [];
end