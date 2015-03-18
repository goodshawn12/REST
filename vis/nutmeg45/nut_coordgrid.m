function xyz=nut_coordgrid(x,y,z)
%  xyz = nut_coordgrid([x_vector],[y_vector],[z_vector])
%  xy = nut_coordgrid([x_vector],[y_vector])
%
% Form Nx2 or Nx3 list of coordinates forming a grid from given vectors.
 

if(exist('z','var'))
    [meshx,meshy,meshz]=ndgrid(x,y,z);
    xyz=single([meshx(:) meshy(:) meshz(:)]);
else
    [meshx,meshy]=ndgrid(x,y);
    xyz=single([meshx(:) meshy(:)]);
end
        
