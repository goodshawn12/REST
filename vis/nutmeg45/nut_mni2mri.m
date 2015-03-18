function [xyz_o]=nut_mni2mri(xyz_i,coreg)
% [XYZ_O]=NUT_MNI2MRI(XYZ_I,coreg)
% Takes MNI coords (mm) and converts them to MRI coordinates (mm)
%
% XYZ_I is Nx3 list of coords
%
% coreg:  optional input (from either beam.coreg or nuts.coreg)
%         else we will find it from global beam or global nuts
%         Either way, it must contain coreg.norm_mripath and coreg.mripath
%
% this is modified from code originally written by John Ashburner:
% http://www.sph.umich.edu/~nichols/JG2/get_orig_coord2.m


if nargin<2 || isempty(coreg)
    global beam
    if(isempty(beam))
        global nuts
        if(isfield(nuts,'coreg'))
            coreg = nuts.coreg;
        else
            error('no coregistration info found in global beam or nuts');
        end
    else
        coreg = beam.coreg;
    end
end
    
if isstruct(coreg)
    if ~isfield(coreg,'norm_mripath')
        [dum,mrifile,dum] = fileparts(coreg.mripath);
        if ( mrifile(1)=='w' || strcmp(mrifile,'T1') || strncmp(mrifile,'avg',3) )  % if struct MRI is already normalized.     
            norm_mripath=coreg.mripath;       
        else
            error('must have a normed MRI loaded');
        end
    else
        norm_mripath=coreg.norm_mripath;
    end
    mripath=coreg.mripath;
else
    error('second input, if specified, must be either blank or coreg structure')
end


if size(xyz_i,2)~=3, error('coord must be an N x 3 matrix'); end;

if strcmp(mripath,norm_mripath)
    % This may happen if the displayed data is already spatially
    % normalized.
    xyz_o = xyz_i;
    return
end

[path,name,ext]=fileparts(norm_mripath);
Vorig=spm_vol(mripath);
norminput=fullfile(path,[name(2:end) ext]);
Vnorminput=spm_vol(norminput);
snmat=fullfile(path,[name(2:end) '_sn.mat']);

t  = load(snmat);
coord = xyz_i';

Mat = inv(t.VG.mat);
xyz = Mat(1:3,:)*[coord ; ones(1,size(coord,2))];
Tr  = t.Tr;
Affine = t.Affine;
d   = t.VG.dim(1:3);
Mult = t.VF.mat*Affine;

if (prod(size(Tr)) == 0),
    affine_only = 1;
    basX = 0; tx = 0;
    basY = 0; ty = 0;
    basZ = 0; tz = 0;
else,
    affine_only = 0;
    basX = spm_dctmtx(d(1),size(Tr,1),xyz(1,:)-1);
    basY = spm_dctmtx(d(2),size(Tr,2),xyz(2,:)-1);
    basZ = spm_dctmtx(d(3),size(Tr,3),xyz(3,:)-1);
end;

if affine_only,
    xyz2 = Mult(1:3,:)*[xyz ; ones(1,size(xyz,2))];
else,
    Tr1 = reshape(Tr(:,:,:,1),size(Tr,1)*size(Tr,2),size(Tr,3));
    Tr2 = reshape(Tr(:,:,:,2),size(Tr,1)*size(Tr,2),size(Tr,3));
    Tr3 = reshape(Tr(:,:,:,3),size(Tr,1)*size(Tr,2),size(Tr,3));
    
    xyz2 = zeros(size(xyz));
    xyztmp = zeros(size(xyz));
    
    for i=1:size(xyz,2),
        bx = basX(i,:);
        by = basY(i,:);
        bz = basZ(i,:);
        % 		tx = reshape(...
        % 			reshape(Tr(:,:,:,1),size(Tr,1)*size(Tr,2),size(Tr,3))...
        % 			*bz', size(Tr,1), size(Tr,2) );
        % 		ty = reshape(...
        % 			reshape(Tr(:,:,:,2),size(Tr,1)*size(Tr,2),size(Tr,3))...
        % 			*bz', size(Tr,1), size(Tr,2) );
        % 		tz =  reshape(...
        % 			reshape(Tr(:,:,:,3),size(Tr,1)*size(Tr,2),size(Tr,3))...
        % 			*bz', size(Tr,1), size(Tr,2) );
        tx = reshape(Tr1*bz', size(Tr,1), size(Tr,2) );
        ty = reshape(Tr2*bz', size(Tr,1), size(Tr,2) );
        tz = reshape(Tr3*bz', size(Tr,1), size(Tr,2) );
        % 		xyz2(:,i) = Mult(1:3,:)*[xyz(:,i) + [bx*tx*by' ; bx*ty*by' ; bx*tz*by']; 1];
        xyztmp(:,i) = [bx*tx*by' ; bx*ty*by' ; bx*tz*by'];
    end;
    
    xyz2 = Mult(1:3,:)*[ xyz + xyztmp ; ones(1,size(xyz,2))];
    
end;
orig_coord = xyz2';

xyz_o=nut_coordtfm(orig_coord,Vorig.mat*inv(Vnorminput.mat));





