function [hsCoord_mrimm_ls, mean_sq_err] = nut_coreg_leastsquares(hsCoord_meg,mri_prjpts)
% [hsCoord_FINAL, MEAN_SQ_ERR] =NUT_COREG_LEASTSQUARES(hsCoord_MEG, mri_prjpts)
%
% Uses ST globals.
%
% hsCoord_meg in MEG mm space and mri_prjpts in MRI mm space
% finds estimation of transformation between the two.
	
global st coreg
%normalize MEG and MRI data
meg_origin=[mean(hsCoord_meg(:,1)) mean(hsCoord_meg(:,2)) mean(hsCoord_meg(:,3))]; 
normal_meg=[hsCoord_meg(:,1)-meg_origin(1) hsCoord_meg(:,2)-meg_origin(2) hsCoord_meg(:,3)-meg_origin(3)];
mri_origin=[mean(mri_prjpts(:,1)) mean(mri_prjpts(:,2)) mean(mri_prjpts(:,3))];
normal_mri=[mri_prjpts(:,1)-mri_origin(1) mri_prjpts(:,2)-mri_origin(2) mri_prjpts(:,3)- mri_origin(3)];
% meg_origin=[mean(nuts.datapoints_meg(:,1)) mean(nuts.datapoints_meg(:,2)) mean(nuts.datapoints_meg(:,3))]; 
% normal_meg=[nuts.datapoints_meg(:,1)-meg_origin(1) nuts.datapoints_meg(:,2)-meg_origin(2) nuts.datapoints_meg(:,3)-meg_origin(3)];
% mri_origin=[mean(mri_prjpts(:,1)) mean(mri_prjpts(:,2)) mean(mri_prjpts(:,3))];
% normal_mri=[mri_prjpts(:,1)-mri_origin(1) mri_prjpts(:,2)-mri_origin(2) mri_prjpts(:,3)- mri_origin(3)];

%find covariance
for i=1:3
    for j=1:3
        covariance(i,j)=mean(normal_meg(:,i).*normal_mri(:,j));
    end
end
%find rotation matrix
[U,S,V] = svd(covariance,0);
rotation=V*U';
%find translation
translation= mri_origin' - rotation * meg_origin';
%form transformation matrix
transform_matrix_final=[rotation translation];
transform_matrix_final(4,:)=0;
transform_matrix_final(4,4)=1;
%bring MEG points to MRI space
coreg.meg2mri_tfm=transform_matrix_final;

%hsCoord_mrimm_ls=nut_meg2mri(hsCoord_meg);
%hsCoord_mrimm_ls=nut_meg2mri(hsCoord_meg);
hsCoord_mrimm_ls=nut_coordtfm(hsCoord_meg,coreg.meg2mri_tfm);
hsCoord_mrimm_ls=nut_coordtfm(hsCoord_meg,coreg.meg2mri_tfm);
%estimation=estimation';
%find error
error=mri_prjpts-hsCoord_mrimm_ls(:,1:3);
%convert Transformation matrix such that when applied to MEG points in mm
%yields points in MRI space in voxels
%transform_matrix_final_vx = inv(st.vols{1}.mat)*transform_matrix_final;
%nuts.coreg.meg2mri_vx_tfm = transform_matrix_final_vx;

%hsCoord_final=nut_meg2mri_vx(hsCoord_meg);
% datapoints_final=transform_matrix_final_vx*datapoints_meg';
% datapoints_final=datapoints_final';
% datapoints_final(:,4)=[];

%find mean square error
mean_sq_err=sqrt(mean(error(:,1).*error(:,1))+mean(error(:,2).*error(:,2))+mean(error(:,3).*error(:,3)));
%save('-APPEND',workfile,'transform_matrix_final','transform_matrix_final_vx','datapoints_final','mean_sq_err');
return
