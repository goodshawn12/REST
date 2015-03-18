function nut_verify_coregistry_in_3d
% NUT_VERIFY_COREGISTRY_IN_3D
%
% Plots the the nut_coregistered data for all the stages
%
% Uses NUTS global.
errordlg('I don''t think so. This function needs a lot of work!');
return
        global nuts;
            %check if required file and variables exist
            if  isempty(nuts)
                msgbox('Load the image','Volume not found','warn')
                return;
            end    
			
%%%% EXPERIMENTAL HACK -- do everything in physical mm rather than voxels
%%%% (SPM2 only!)
% if(spm('ver')=='SPM2')
            if ~isfield(nuts.coreg,'fiducials_mri')
                msgbox('Go to Mark/Edit Fiducials -> Save fiducials and show coregistry in 2D first',...
                'Could not find data in NUTS','warn')
                return;            
            end
            if ~isfield(nuts.coreg,'mesh')
                msgbox('Click Mesh Generation -> nut_cloud2mesh first','Mesh not found','warn')
                return            
            end

    global st
    VOX = nut_pixdim(st.vols{1});
%     ORIGIN = st.vols{1}.private.hdr.hist.origin(1:3);
    ORIGIN = st.vols{1}.dim(1:3)/2;


	if(0)
            for i=1:3
                        nuts.datapoints_mri_vx(:,i) = (nuts.datapoints_mri_vx(:,i)-ORIGIN(i))*VOX(i);
                        nuts.datapoints_final(:,i) = (nuts.datapoints_final(:,i)-ORIGIN(i))*VOX(i);
                        if isfield(nuts.coreg,'datapoints_final_iter')
                            nuts.datapoints_final_iter(:,i) = (nuts.datapoints_final_iter(:,i)-ORIGIN(i))*VOX(i);
                        end
            end
    else
        nuts.datapoints_mri_vx = nut_voxels2mm(nuts.datapoints_mri_vx);
        nuts.datapoints_final = nut_voxels2mm(nuts.datapoints_final);
        if isfield(nuts.coreg,'datapoints_final_iter')
            nuts.datapoints_final_iter = nut_voxels2mm(nuts.datapoints_final_iter);
        end
    end
    
            figure
            subplot(221)            
                nut_show_head(nuts.coreg.mesh);
                hold on;
                plot3(nuts.datapoints_mri_vx(:,1),nuts.datapoints_mri_vx(:,2),nuts.datapoints_mri_vx(:,3),'.');
                hold off;
                title('First Stage Co-registry')
           subplot(222)           
                nut_show_head(nuts.coreg.mesh);
                hold;
                plot3(nuts.prjpts_mm(:,1),nuts.prjpts_mm(:,2),nuts.prjpts_mm(:,3),'.');
                hold off;
                title('Projected Points')
           subplot(223)           
                nut_show_head(nuts.coreg.mesh);
                hold on;
                plot3(nuts.datapoints_final(:,1),nuts.datapoints_final(:,2),nuts.datapoints_final(:,3),'.');
                hold off
                title('Final Stage Co-registry')
                xlabel(strcat(['The error is ', num2str(nuts.coreg.mean_sq_err),' mm']));
          %plot only if "datapoints_final_iter" exists
          if isfield(nuts.coreg,'datapoints_final_iter')
                subplot(224)
                nut_show_head(nuts.coreg.mesh);
                hold on;
                plot3(nuts.datapoints_final_iter(:,1),nuts.datapoints_final_iter(:,2),nuts.datapoints_final_iter(:,3),'.');
                hold off
                title(strcat(['Final Stage Co-registry after ', num2str(nuts.iteration),' iteration']))
                xlabel(strcat(['The error is ', num2str(nuts.coreg.mean_sq_err_iter),' mm']));
           end
% else
%                 if exist('fiducials_mri_vx')==0
%                 msgbox('Go to Mark/Edit Fiducials -> Save fiducials and show coregistry in 2D first',...
%                 'Could not find data file data_meg_mri.mat','warn')
%                 return;            
%             end
%             if exist('mesh')==0
%                 msgbox('Click Mesh Generation -> nut_cloud2mesh first','Mesh not found','warn')
%                 return            
%             end
%             figure
%             subplot(221)            
%                 nut_show_head(mesh);
%                 hold on;
%                 plot3(datapoints_mri_vx(:,1),datapoints_mri_vx(:,2),datapoints_mri_vx(:,3),'.');
%                 hold off;
%                 title('First Stage Co-registry')
%            subplot(222)           
%                 nut_show_head(mesh);
%                 hold;
%                 plot3(prjpts(:,1),prjpts(:,2),prjpts(:,3),'.');
%                 hold off;
%                 title('Projected Points')
%            subplot(223)           
%                 nut_show_head(mesh);
%                 hold on;
%                 plot3(datapoints_final(:,1),datapoints_final(:,2),datapoints_final(:,3),'.');
%                 hold off
%                 title('Final Stage Co-registry')
%                 xlabel(strcat(['The error is ', num2str(mean_sq_err),' mm']));
%           %plot only if "datapoints_final_iter" exists
%           if exist('datapoints_final_iter')==1
%                 subplot(224)
%                 nut_show_head(mesh);
%                 hold on;
%                 plot3(datapoints_final_iter(:,1),datapoints_final_iter(:,2),datapoints_final_iter(:,3),'.');
%                 hold off
%                 title(strcat(['Final Stage Co-registry after ', num2str(iteration),' iteration']))
%                 xlabel(strcat(['The error is ', num2str(mean_sq_err_iter),' mm']));
%            end
%        end
    
       
    
