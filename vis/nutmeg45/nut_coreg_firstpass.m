function nut_coreg_firstpass

global coreg st
%[filename_ahs, pathname1]=uigetfile('*.ahs','Please select the head shape.');
% if isequal(filename_ahs,0)|isequal(pathname1,0)
%     % screw this, user hit cancel
%     return;
% else
%     % Read the Head Shape file, get the headshape points and the fiducials       
%     [a,b,c,d,e,f,g,h,i,j,k,l,m]=textread([pathname1 filename_ahs],'%f%f%f%n%n%n%n%n%n%n%n%s%s',2, 'headerlines',1);
%     [a(3),b(3),c(3),d,e,f,g,h,i,j,k,l]=textread([pathname1 filename_ahs],'%f%f%f%n%n%n%n%n%n%n%n%s',1, 'headerlines',3);
%     fiducials_meg=[a,b,c]*10;
%     fiducials_meg(:,4)=1;
%     fiducials_meg=fiducials_meg';
%     [a1,b1,c1,d,e,f,g,h,i,j,k,l]=textread([pathname1 filename_ahs],'%f%f%f%n%n%n%n%n%n%n%n%s','headerlines',6);
%     datapoints_meg=[a1,b1,c1]*10;
%     % datapoints_meg(:,4)=1;
%     datapoints_meg=[datapoints_meg;fiducials_meg'];     
%     clear a1 b1 c1 a b c d e f g h i j k l m;
%     nuts.headshapefile=[pathname1 filename_ahs];
% end


% hsCoord_mri=nut_meg2mri(coreg.hsCoord);
hsCoord_mri = nut_coordtfm(coreg.hsCoord,coreg.meg2mri_tfm);

% project the points to Head surface
[prjpts_mm]=nut_prjsrf(coreg.mesh,hsCoord_mri); %,ORIGIN,VOX);

% Find the final stage transformation matrix
[coreg.hsCoord_mrimm_ls, coreg.mean_sq_err]= nut_coreg_leastsquares(coreg.hsCoord,prjpts_mm); %,VOX,ORIGIN);


% following jazz is needed because SPM's blob freaks out when given
% negative coordinates, and needs a dilation to know these are big voxels
% (chose 4 mm^3 voxel size)
voxelsize = [4 4 4];
if min(coreg.hsCoord_mrimm_ls(:))<0
    shift=min(coreg.hsCoord_mrimm_ls(:));
else
    shift=0;
end
translation_tfm = [ voxelsize(1)            0            0 shift
                               0 voxelsize(2)            0 shift
                               0            0 voxelsize(3) shift
                               0            0            0 1 ];

hs_coord = nut_coordtfm(coreg.hsCoord_mrimm_ls,inv(translation_tfm));

keep=find(prod(double(hs_coord > 0.5),2)); % discard coords with nonpositive voxels
hs_blobs  = hs_coord(keep,:);
spm_orthviews('rmblobs',1);
spm_orthviews('addblobs',1,hs_blobs',zeros(size(keep)),translation_tfm);
if(strcmp(spm('ver'),'SPM2'))
    % only for SPM2; SPM8 freaks out if we delete the colorbar
    delete(st.vols{1}.blobs{1}.cbar);
end
spm_orthviews('redraw',1);
nut_coreg_enabler