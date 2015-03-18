function nut_edit_fiducials
% NUT_EDIT_FIDUCIALS
%
% This function is actuated when edit fiducials in 3d is pressed.
error('I am not going to let you run this file. It needs a lot of work!');
global nuts st
switch get(gcbo,'Tag')
    case {'nut_edit_fiducials_3d'}
        %check if "fiducials_mri_vx" exists in the data file data_meg_mri.mat
        if ~isfield(nuts.coreg,'fiducials_mri')
            msgbox('Go to Mark/Edit Fiducials -> Save fiducials and show coregistry in 2D first','Could not find data file data_meg_mri.mat','warn')
            return;
        end
        %check if "mesh" exists in the structure nuts
        if ~isfield(nuts.coreg,'mesh')
            msgbox('Click Mesh Generation -> nut_cloud2mesh first','Mesh not found','warn')
            return            
        end
        %check if "datapoints_final_iter" exists in the data file data_meg_mri.mat
        if isfield(nuts.coreg,'datapoints_final_iter')
            %clear "datapoints_final_iter" from the data_meg_mri.mat
            nuts.coreg = rmfield(nuts.coreg,'datapoints_final_iter');
        end
        %getting values from GUI edit test boxes
%         fov=str2num(get(findobj('Tag','nut_field_of_view'),'String'));
%         width=str2num(get(findobj('Tag','nut_slice_width'),'String'));
        DIM = st.vols{1}.dim(1:3);
        %DIM=sscanf(get(findobj('Tag','nut_voxel_dimension'),'String'), '%g %g %g');
        %checking for invalid value
%         if (isempty(fov)==1 | isempty(width)==1 | isempty(DIM)==1)
%             msgbox('Please give field of view, voxel dimension, and slice width','Unable to display fiducials','warn')
%             return
%         end
        %displaying the fiducials in 3D for editting
        %%%% EXPERIMENTAL
        fiducials=nuts.coreg.fiducials_mri_mm;
        figure;
        nut_show_head(nuts.coreg.mesh);
        hold on;      
        plot3(fiducials(:,1),fiducials(:,2),fiducials(:,3),'.');
        hold off
        title('Please hit Enter each time you select a point on Figure from mouse')
        fiducials=fiducials';
        % getting user input for each fiducial
        buttonName=questdlg('Do you want to change LPA','Fiducials Edition') ;              
        if strcmp(buttonName,'Cancel')==1
            return
        elseif strcmp(buttonName,'Yes')==1                   
            %waiting for user to hit ENTER
            pause
            %getting the selected point from the graph
            p=select3d(gca);
            if isempty(p)~=1
                fiducials(:,1)=p;
                hold on
                %plotting the marked fiducial
                plot3(fiducials(1,1),fiducials(2,1),fiducials(3,1),'*');  
                hold off
            end
        else
        end
        buttonName=questdlg('Do you want to change RPA','Fiducials Edition') ;              
        if strcmp(buttonName,'Cancel')==1
            return
        elseif strcmp(buttonName,'Yes')==1                   
            %waiting for user to hit ENTER
            pause
            %getting the selected point from the graph
            p=select3d(gca);
            if isempty(p)~=1
                fiducials(:,2)=p;
                hold on
                %plotting the marked fiducial
                plot3(fiducials(1,2),fiducials(2,2),fiducials(3,2),'*');  
                hold off
            end
        else
        end
        buttonName=questdlg('Do you want to change Nasion','Fiducials Edition') ;              
        if strcmp(buttonName,'Cancel')==1
            return
        elseif strcmp(buttonName,'Yes')==1                   
            %waiting for user to hit ENTER
            pause
            %getting the selected point from the graph
            p=select3d(gca);
            if isempty(p)~=1
                fiducials(:,3)=p;
                hold on
                %plotting the marked fiducial
                plot3(fiducials(1,3),fiducials(2,3),fiducials(3,3),'*');  
                hold off
            end
        else
        end                
        fiducials=fiducials';                
        origin_vx=DIM/2;
%         voxelwidth=[fov/DIM(1) fov/DIM(2) width];
        % EXPERIMENTAL -- with SPM2, we are displaying in "real" dimensions, so
        % this affects how we record new fiducials....
            %Changing fiducials from voxel to mm to store in a file
            fiducials_mri_mm(:,1)=fiducials(:,1);
            fiducials_mri_mm(:,2)=fiducials(:,2);
            fiducials_mri_mm(:,3)=fiducials(:,3);
        [filename pathname] = uiputfile('*.txt','Save the fiducials to desired file (OPTIONAL)...');
        %save the marked fiducials to the desired file (OPTIONAL)
        if ischar(filename)
            fid=fopen([pathname filename],'w');                
            fprintf(fid,'%g %g %g\n %g %g %g\n %g %g %g',fiducials_mri_mm(1,1),fiducials_mri_mm(1,2),fiducials_mri_mm(1,3),fiducials_mri_mm(2,1),fiducials_mri_mm(2,2),fiducials_mri_mm(2,3),fiducials_mri_mm(3,1),fiducials_mri_mm(3,2),fiducials_mri_mm(3,3));
            fclose(fid);		
        end
        % making axis system in MRI space same as MEG space 
        origin_mri=(fiducials_mri_mm(1,:)+fiducials_mri_mm(2,:))/2;
        x_mri=fiducials_mri_mm(3,:)-origin_mri;
        x_mri=x_mri/sqrt(x_mri*x_mri');
        z_mri=cross(x_mri,fiducials_mri_mm(2,:)-origin_mri);
        z_mri=z_mri/sqrt(z_mri*z_mri');
        y_mri=cross(x_mri,z_mri);
        transformation_matrix=[x_mri' y_mri' z_mri' origin_mri'];
        transformation_matrix(4,:)=0;
        transformation_matrix(4,4)=1;
        %Converting transformation matrix such that when applied to
        %MEG points in mm gives result in voxels
        transformation_matrix_vx = inv(st.vols{1}.mat)*transformation_matrix;
        % Bringing MEG points to MRI space 
        datapoints_mri_vx=transformation_matrix_vx*(nuts.datapoints_meg)';
        datapoints_mri_vx=datapoints_mri_vx';
        datapoints_mri_vx(:,4)=[];    
        fiducials_mri_vx=fiducials;
        % projecting the MRI points to Head Surface
        [prjpts_mm,prjpts]=nut_prjsrf(nuts.coreg.mesh,datapoints_mri_vx); %,ORIGIN,VOX);
        % Finding final coregistry
        [transform_matrix_final_vx, datapoints_final, mean_sq_err]= nut_meg2mri_final(nuts.datapoints_meg,prjpts_mm); %,VOX,ORIGIN);
        nuts.datapoints_mri_vx=datapoints_mri_vx;
        nuts.transformation_matrix_vx=transformation_matrix_vx;
        nuts.coreg.meg2mri_tfm_vx=transform_matrix_final_vx;
        nuts.coreg.fiducials_mri_mm=fiducials_mri_mm;
        nuts.fiducials_mri_vx=fiducials_mri_vx;
        nuts.datapoints_final=datapoints_final;
        nuts.prjpts_mm=prjpts_mm;
        return            
    otherwise
        warndlg('Callback not yet defined for this button',...
            'UNDEFINED CALLBACK');
        return
end


