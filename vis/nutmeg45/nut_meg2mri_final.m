function [transform_matrix_final_vx, datapoints_final, mean_sq_err] = nut_coreg_leastsquares(datapoints_meg,mri_data)
% [TRANSFORM_MATRIX_FINAL_VX, DATAPOINTS_FINAL, MEAN_SQ_ERR] = ...
%	NUT_MEG2MRI_FINAL(DATAPOINTS_MEG, MRI_DATA)
%
% This function will allow you to do estimation between the MEG points and
% MRI points.
%
% Uses NUTS and ST globals.
	
global nuts st;
%normalize MEG and MRI data
meg_origin=[mean(nuts.datapoints_meg(:,1)) mean(nuts.datapoints_meg(:,2)) mean(nuts.datapoints_meg(:,3))]; 
normal_meg=[nuts.datapoints_meg(:,1)-meg_origin(1) nuts.datapoints_meg(:,2)-meg_origin(2) nuts.datapoints_meg(:,3)-meg_origin(3)];
mri_origin=[mean(mri_data(:,1)) mean(mri_data(:,2)) mean(mri_data(:,3))];
normal_mri=[mri_data(:,1)-mri_origin(1) mri_data(:,2)-mri_origin(2) mri_data(:,3)- mri_origin(3)];
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
estimation=transform_matrix_final*nuts.datapoints_meg';
estimation=estimation';
%find error
error=mri_data-estimation(:,1:3);
%convert Transformation matrix such that when applied to MEG points in mm
%yields points in MRI space in voxels
transform_matrix_final_vx = inv(st.vols{1}.mat)*transform_matrix_final;
datapoints_final=transform_matrix_final_vx*nuts.datapoints_meg';
datapoints_final=datapoints_final';
datapoints_final(:,4)=[];
%find mean square error
mean_sq_err=sqrt(mean(error(:,1).*error(:,1))+mean(error(:,2).*error(:,2))+mean(error(:,3).*error(:,3)));
nuts.coreg.meg2mri_tfm=transform_matrix_final;
%save('-APPEND',workfile,'transform_matrix_final','transform_matrix_final_vx','datapoints_final','mean_sq_err');
return
