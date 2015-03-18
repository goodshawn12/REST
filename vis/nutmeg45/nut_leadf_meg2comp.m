function [lead_f, lead_g]=nut_leadf(r_source,r_det,normal_v,lsc)
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

% Qefg2 is a unit vector in z-direction; TODO: change variable names to make sense!
Qefg2 = repmat([0 0 1],[numvoxels 1]); %  derive leads_f by computing orthogonal to source-z plane  -> latitude
Qg = cross(Qefg2,r_source,2);
clear Qefg2;
Qg = Qg./repmat(nut_rownorm(Qg),[1 3]);  % normalize
Qg(isnan(Qg))=0;  % normalization of [0 0 0] should remain 0 instead of NaN

Qf = cross(Qg,r_source,2);   % compute orthorgonal to plane formed by source-Qg plane -> longitude
Qf = Qf./repmat(nut_rownorm(Qf),[1 3]); % normalize
Qf(isnan(Qf))=0;  % normalization of [0 0 0] should remain 0 instead of NaN

%%%% no longer necessary... replaced with above cross products!)
% Qefg1=[0 1 0]; %  derive leads_f by assuming Qefg=[0 1 0]
% Qefg2=[0 0 1]; %  derive leads_g by assuming Qefg=[0 0 1]
% 
% Qf = zeros(size(r_source));
% Qg = zeros(size(r_source));
% for i=1:length(r_source);
%    Qf(i,:)=efg2xyz(Qefg1, r_source(i,:), lsc);
%    Qg(i,:)=efg2xyz(Qefg2, r_source(i,:), lsc);
% end
%%%%

aa1f = cross(Qf,r_source,2);
temp = repmat(aa1f,[1 1 numsensors]);
temp = shiftdim(temp,2);
aa1f = reshape(temp,[numsensors*numvoxels 3]);

aa1g = cross(Qg,r_source,2);
temp = repmat(aa1g,[1 1 numsensors]);
temp = shiftdim(temp,2);
aa1g = reshape(temp,[numsensors*numvoxels 3]);

r_det = repmat(r_det,numvoxels,1);
normal_v = repmat(normal_v,numvoxels,1);
temp = repmat(r_source,[1 1 numsensors]);
temp = shiftdim(temp,2);
r_source = reshape(temp,[numsensors*numvoxels 3]);
clear Qf Qg temp;

a_vec = r_det - r_source;
a_mag = nut_rownorm(a_vec);
rd_mag = nut_rownorm(r_det);
%
   F = a_mag.*(rd_mag.*a_mag + rd_mag.^2 - nut_dotfast(r_source,r_det,2));
   coef1 = (1./rd_mag).*a_mag.^2 + (1./a_mag).*nut_dotfast(a_vec,r_det,2)+2.*a_mag+2.*rd_mag;
   coef2 = a_mag + 2*rd_mag +  (1./a_mag).*nut_dotfast(a_vec,r_det,2);
   clear a_vec a_mag rd_mag
   vec_F= nut_rowmult(coef1,r_det) - nut_rowmult(coef2,r_source);
   clear r_source coef1 coef2
   
   % derive leads_f by assuming Qefg=[0 1 0]   
   B = nut_rowmult(const./F.^2,nut_rowmult(F,aa1f) - nut_rowmult(nut_dotfast(r_det,aa1f,2),vec_F));
   clear aa1f
   lead_f = nut_dotfast(B,normal_v,2);
   % clear B

   %  derive leads_g by assuming Qefg=[0 0 1]
   B = nut_rowmult(const./F.^2,nut_rowmult(F,aa1g) - nut_rowmult(nut_dotfast(r_det,aa1g,2),vec_F));
   clear r_det aa1g
   lead_g = nut_dotfast(B,normal_v,2);
   clear vec_F

lead_f =  reshape(lead_f,numsensors,numvoxels);
lead_g =  reshape(lead_g,numsensors,numvoxels);
