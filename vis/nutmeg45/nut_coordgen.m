function XYZ = nut_coordgen(lims, stepsize)
% nut_coordgen([xmin xmax ymin ymax zmin zmax])
% generates Nx3 list of coordinates from specified bounds

if(nargin==1)
    stepsize =1;
end

[meshx,meshy,meshz]=ndgrid(lims(1):stepsize:lims(2),lims(3):stepsize:lims(4),lims(5):stepsize:lims(6));
XYZ = [meshx(:) meshy(:) meshz(:)];
