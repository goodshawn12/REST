function [lead_f, lead_g]=nut_leadf_slow(r_source,r_det,normal_v,lsc)
% [LEAD_F, LEAD_G]=NUT_LEADF_SLOW(R_SOURCE,R_DET,NORMAL_V,LSC)
%
%  Original lead field computation -- needs to be wrapped in nested for
%  loops since it can take only one r_source and r_det at a time... slooow
%
%  calculation of lead field for f and g component 
%  in the spherical homogeneous conductor
%  using Sarvas formula 
%  appeared in Phys. Med. Biol. Vol.32, 11-22, 1987
%
%  R_SOURCE : source coordinate   (cm)         (raw vector)
%  R_DET    : detector coordinate (cm)       (raw vector)
%  NORMAL_V : detector normal dorection     (raw vector)
%  LSC      : local sphere center  (cm)      (raw vector)
%



%    const=10^-7;  % for r (cm), Q (nAm), B (pT) 
const = 1e-10;  % input = millimeters, output = Tesla, simulating a 1 nAm dipole strength 
%
   r_source = (r_source -lsc);
   r_det    = (r_det    -lsc);
   lsc      = (lsc      -lsc);  % origin set at local sphere center
%
   a_vec=r_det - r_source;
   a_mag=norm(a_vec);
   rd_mag=norm(r_det);
%
   F = a_mag*(rd_mag*a_mag + rd_mag^2 - r_source*r_det');
   coef1 = (1/rd_mag)*a_mag^2 + (1/a_mag)*a_vec*r_det'+2*a_mag+2*rd_mag;
   coef2 = a_mag + 2*rd_mag +  (1/a_mag)*a_vec*r_det';
   vec_F= coef1*r_det - coef2*r_source;
   
%
%  derive leads_g

   % compute orthogonal to source-z plane  -> latitude (east/west on globe)
   Qlat = cross([0 0 1],r_source); % [0 0 1] in cartesian (i.e., direction of z-axis)
   Qlat = Qlat/norm(Qlat);  % normalize to unit vector
   aa1=cross(Qlat,r_source);
   B = (const/F^2) *(F*aa1 - (aa1*r_det')*vec_F);
   lead_g = B*normal_v';

%
%  derive leads_f
% 
   % compute orthogonal to plane formed by source-Qlat plane -> longitude (north/south)
   Qlon = cross(Qlat,r_source);
   Qlon = Qlon/norm(Qlon); % normalize to unit vector
   aa1=cross(Qlon,r_source);
   B = (const/F^2) *(F*aa1 - (aa1*r_det')*vec_F);
   lead_f = B*normal_v';
