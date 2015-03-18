function nut_adjust_coregistry
% function nut_adjust_coregistry
% This function would allow the user to do the iteration
% over the MRI data to reduce the error between the estimated
% MRI headshape data and the points projected to the head surface

global coreg st
if(~isfield(coreg,'mesh'))
             msgbox('Click Mesh Generation -> nut_cloud2mesh first','Mesh not found','warn')
             return            
end
%check if hsCoord exist in the data file
if(~isfield(coreg,'hsCoord'))  %what about hsCoord_mrimm_ls
        msgbox('Go to Mark/Edit Fiducials -> Save fiducials and show coregistry in 2D first',...
               'Could not find data in structure NUTS','warn')
        return;
else
        %Generate a prompt for user to respond
        def={'3'};
        prompt   = {'Select the number of iterations'};
        title    = 'Iteration Value ';
        lines = 1;
        answer   = inputdlg(prompt,title,lines,def);
        if (isempty(answer)) 
            % msgbox('You pressed "Cancel"');
            return;
        end
        %get the value of iteration from the PROMPT line to the workspace
        iteration=['iteration=[' answer{1} '];'];
        eval(iteration);
        if (isempty(iteration)==1)
            msgbox('Not a valid iteration value', 'Unable to do iterations','warn')
            return
        end
        %calling the function
        [transform_matrix_final_iter,coreg.hsCoord_iter,coreg.mean_sq_err_iter]=nut_coregistry(coreg.mesh,iteration,coreg.hsCoord,coreg.hsCoord_mrimm_ls,coreg.meg2mri_tfm);
        coreg.meg2mri_tfm_lsq=coreg.meg2mri_tfm;
        coreg.meg2mri_tfm=transform_matrix_final_iter;
        coreg.iteration=iteration;
        
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
        hs_coord = nut_coordtfm(coreg.hsCoord_iter,inv(translation_tfm));
		keep=find(prod(double(hs_coord > 0.5),2)); % discard coords with nonpositive voxels
		hs_blobs  = hs_coord(keep,:);
		spm_orthviews('rmblobs',1);
		spm_orthviews('addblobs',1,hs_blobs',zeros(size(keep)),translation_tfm);
        if(strcmp(spm('ver'),'SPM2'))
            % only for SPM2; SPM8 freaks out if we delete the colorbar
            delete(st.vols{1}.blobs{1}.cbar);
        end
		spm_orthviews('redraw',1);


end


function [transform_matrix_final_iter,hsCoord_iter,mean_sq_err_iter]=nut_coregistry(mesh,iteration,hsCoord_orig,hsCoord_ls,transform_matrix_final)  
%get the values from the edit text box of the GUI nut_meg_mri_toolbox
% fov=str2num(get(findobj('Tag','nut_field_of_view'),'String'));
% width=str2num(get(findobj('Tag','nut_slice_width'),'String'));
global st
% DIM = st.vols{1}.dim(1:3);
% fov=st.vols{1}.private.hdr.dime.dim(2)*st.vols{1}.private.hdr.dime.pixdim(2);
% width=st.vols{1}.private.hdr.dime.pixdim(4);

transform=eye(4);
h = waitbar(0,'Please wait while we iterate...');
%find the transformation matrix after the iteration 
for k=1:iteration    
%    [h,res]=prjsrf1(mesh,hsCoord_ls,iteration,k,h);
    res=nut_prjsrf(mesh,hsCoord_ls);
%    [w,r0,err]=nut_cnv3(res,hsCoord_ls,1);
    [w,r0,err_mm]=nut_cnv3(res,hsCoord_ls,1);
    hsCoord_ls=hsCoord_ls*w+nut_spread(hsCoord_ls,r0);        
    transform=[[w,r0'];[0,0,0,1]]*transform;  
    waitbar(k/iteration,h);
end
%changing voxels to mm
% % err_mm(:,1)=err(:,1)*fov/DIM(1);
% % err_mm(:,2)=err(:,2)*fov/DIM(2);
% % err_mm(:,3)=err(:,3)*width;
%err_mm=nut_voxels2mm(err);
mean_sq_err_iter=sqrt(mean(err_mm(:,1).*err_mm(:,1))+mean(err_mm(:,2).*err_mm(:,2))+mean(err_mm(:,3).*err_mm(:,3)));
close(h)
%getting the MEG points in MRI space by applying the transformation matrix
%found above
transform_matrix_final_iter=transform*transform_matrix_final;
hsCoord_iter=nut_coordtfm(hsCoord_orig,transform_matrix_final_iter);
% hsCoord_iter=hsCoord_iter';
% hsCoord_iter(:,4)=[];

% function p=spread(head,r0)
% p=r0(:)*ones(1,size(head,1),1);
% p=p';

% function [h,res]=prjsrf1(mesh,cloud,iteration,i,h)
% %project the points to the head surface
% nn=nut_normsrf(mesh);
% [n1,n2,n3]=size(mesh);
% mesh=reshape(mesh,n1*n2,3);
% nn=reshape(nn,n1*n2,3);
% for k=1:size(cloud,1)
%  r0=cloud(k,:);
%  [tmp,in]=min((mesh(:,1)-r0(1)).^2 + (mesh(:,2)-r0(2)).^2 + (mesh(:,3)-r0(3)).^2);
%  rn=mesh(in,:);
%  n=nn(in,:);
%  res(k,:)=(n*(rn-r0)')*n+r0;
%  waitbar((k*i)/(iteration*size(cloud,1)),h);
% end

