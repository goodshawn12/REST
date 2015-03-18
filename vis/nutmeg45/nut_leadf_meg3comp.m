function [lead_x, lead_y, lead_z]=nut_leadf_xyz(r_source,r_det,normal_v,lsc)
% [LEAD_F, LEAD_G]=NUT_LEADF(R_SOURCE,R_DET,NORMAL_V,LSC)
%
%  calculation of lead field for f and g component 
%  in the spherical homogeneous conductor
%  using Sarvas formula 
%  appeared in Phys. Med. Biol. Vol.32, 11-22, 1987
%
%  R_SOURCE : source coordinate   (mm)         (raw vector)
%  R_DET    : detector coordinate (mm)       (raw vector)
%  NORMAL_V : detector normal direction     (raw vector)
%  LSC      : local sphere center  (mm)      (raw vector)

% const = 1e-7 % in SI units (input = meters, output = Tesla), simulating a 1 Am dipole strength
const = 1e-10;  % input = millimeters, output = Tesla, simulating a 1 nAm dipole strength

numvoxels = size(r_source,1);
numsensors = size(r_det,1);

r_det = nut_coord_diff(r_det,lsc);
r_source = nut_coord_diff(r_source,lsc);

Qx = repmat([1 0 0],[numvoxels 1]);
aa1x = cross(Qx,r_source,2);
clear Qx
temp = repmat(aa1x,[1 1 numsensors]);
temp = shiftdim(temp,2);
aa1x = reshape(temp,[numsensors*numvoxels 3]);

Qy = repmat([0 1 0],[numvoxels 1]);
aa1y = cross(Qy,r_source,2);
clear Qy
temp = repmat(aa1y,[1 1 numsensors]);
temp = shiftdim(temp,2);
aa1y = reshape(temp,[numsensors*numvoxels 3]);


Qz = repmat([0 0 1],[numvoxels 1]);
aa1z = cross(Qz,r_source,2);
clear Qz
temp = repmat(aa1z,[1 1 numsensors]);
temp = shiftdim(temp,2);
aa1z = reshape(temp,[numsensors*numvoxels 3]);


r_det = repmat(r_det,numvoxels,1);
normal_v = repmat(normal_v,numvoxels,1);
temp = repmat(r_source,[1 1 numsensors]);
temp = shiftdim(temp,2);
r_source = reshape(temp,[numsensors*numvoxels 3]);
clear temp;

a_vec = r_det - r_source;
a_mag = nut_rownorm(a_vec);
rd_mag = nut_rownorm(r_det);
%
   F = a_mag.*(rd_mag.*a_mag + rd_mag.^2 - nut_dotfast(r_source,r_det,2));
   constInvFsq = const./F.^2;
%    coef1 = (1./rd_mag).*a_mag.^2 + (1./a_mag).*nut_dotfast(a_vec,r_det,2)+2.*a_mag+2.*rd_mag;
%    coef2 = a_mag + 2*rd_mag +  (1./a_mag).*nut_dotfast(a_vec,r_det,2);
   coef1 = a_mag.^2./rd_mag + nut_dotfast(a_vec,r_det,2)./a_mag + 2.*a_mag + 2.*rd_mag;
   coef2 = a_mag + 2*rd_mag +  nut_dotfast(a_vec,r_det,2)./a_mag;
   clear a_vec a_mag rd_mag
   gradF= nut_rowmult(coef1,r_det) - nut_rowmult(coef2,r_source);
   clear r_source coef1 coef2
   
   % derive leads_f by assuming Qefg=[0 1 0]   
   B = nut_rowmult(constInvFsq,nut_rowmult(F,aa1x) - nut_rowmult(nut_dotfast(r_det,aa1x,2),gradF));
   clear aa1x
   lead_x = nut_dotfast(B,normal_v,2);
   % clear B

   %  derive leads_g by assuming Qefg=[0 0 1]
   B = nut_rowmult(constInvFsq,nut_rowmult(F,aa1y) - nut_rowmult(nut_dotfast(r_det,aa1y,2),gradF));
   clear aa1y
   lead_y = nut_dotfast(B,normal_v,2);

      %  derive leads_g by assuming Qefg=[0 0 1]
   B = nut_rowmult(const./F.^2,nut_rowmult(F,aa1z) - nut_rowmult(nut_dotfast(r_det,aa1z,2),gradF));
   clear r_det aa1z
   lead_z = nut_dotfast(B,normal_v,2);
   clear gradF

lead_x =  reshape(lead_x,numsensors,numvoxels);
lead_y =  reshape(lead_y,numsensors,numvoxels);
lead_z =  reshape(lead_z,numsensors,numvoxels);
