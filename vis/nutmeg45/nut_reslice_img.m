function [z,newxyz]=nut_reslice_img(y,xyz,Vin)
% resamples [y,xyz] of old analyze image to new one, using spm2's V structure
% exact numbers of resampling is hardcoded for now for a specific purpose,
% but can be modified easily for any resampling size

smooth=1;

for ii=1:size(y,3)
    ztmp=resample(y(:,:,ii),3,4,smooth);
    z(:,:,ii)=resample(ztmp',3,4,smooth)';
end
clear ztmp
ztmp=permute(z,[3 1 2]);
for ii=1:size(ztmp,3)
    ztmptmp(:,:,ii)=resample(ztmp(:,:,ii),6,5,smooth);
end
clear z
z=permute(ztmptmp,[2 3 1]);

newmat=[5 0 0 Vin.mat(1,4);...
    0 5 0 Vin.mat(2,4);...
    0 0 5 Vin.mat(3,4);...
    0 0 0 1];

newxyz=nut_coordtfm(nut_coordtfm(xyz',inv(Vin.mat)),newmat);

Vinmod=Vin;
Vinmod.mat=newmat;
Vinmod.dim=[size(z,1) size(z,2) size(z,3) Vin.dim(4)];
Vinmod.fname=[Vin.fname(1:end-4) 'megresamp.img'];
Vinmod.private.hdr.hist.origin(1:3)=nut_coordtfm(nut_coordtfm(Vin.private.hdr.hist.origin(1:3),Vin.mat),inv(newmat));
Vinmod.private.hdr.dime.dim=[4 size(z,1) size(z,2) size(z,3) 1 0 0 0];
Vinmod.private.hdr.dime.pixdim=[0 5 5 5 1 0 0 0];

Vout=spm_write_vol(Vinmod,z);
