% function transformation_matrix=nut_coreg_fiducials(headshape,fiducials_mri_mm)
function coreg=nut_coreg_fiducials(headshape,coreg)
global st;
% global coreg st;

if isfield(coreg,'fiducials_mri_mm')
    fiducials_mri_mm=coreg.fiducials_mri_mm;
end

DIM = st.vols{1}.dim(1:3);
VOX = nut_pixdim(st.vols{1});

if ~exist('fiducials_mri_mm','var')
    l = sscanf(get(findobj('Tag','nut_left_text'),'String'), 'MRI: %g %g %g');
    r = sscanf(get(findobj('Tag','nut_right_text'),'String'), 'MRI: %g %g %g');
    n = sscanf(get(findobj('Tag','nut_nose_text'),'String'), 'MRI: %g %g %g');

    if(isempty(l) | isempty(r) | isempty (n))
        return  % not ready yet if all three fiducials aren't set
    end
    fiducials_mri_mm=[l';r';n'];

else
    l=fiducials_mri_mm(1,:)';
    r=fiducials_mri_mm(2,:)';
    n=fiducials_mri_mm(3,:)';
end

if strcmp(l,'')==0 & strcmp(r,'')==0 & strcmp(n,'')==0 
    % If head shape file is already loaded once, don't ask user to load it again.
    %% TODO:
    %% Figure out what to do if headshape file already loaded (i.e.,
    %% datapoints_meg already exists)    
    if(headshape) % && ~isfield(nuts,'datapoints_meg'))
        disp('this shouldnt be happening')
%         [filename_ahs, pathname1]=uigetfile('*.ahs','Please select the head shape.');
%         if isequal(filename_ahs,0)|isequal(pathname1,0)
%             % screw this, user hit cancel
%             return;
%         else
%             % Read the Head Shape file, get the headshape points and the fiducials       
%             [a,b,c,d,e,f,g,h,i,j,k,l,m]=textread([pathname1 filename_ahs],'%f%f%f%n%n%n%n%n%n%n%n%s%s',2, 'headerlines',1);
%             [a(3),b(3),c(3),d,e,f,g,h,i,j,k,l]=textread([pathname1 filename_ahs],'%f%f%f%n%n%n%n%n%n%n%n%s',1, 'headerlines',3);
%             fiducials_meg=[a,b,c]*10;
%             fiducials_meg(:,4)=1;
%             fiducials_meg=fiducials_meg';
%             [a1,b1,c1,d,e,f,g,h,i,j,k,l]=textread([pathname1 filename_ahs],'%f%f%f%n%n%n%n%n%n%n%n%s','headerlines',6);
%             datapoints_meg=[a1,b1,c1]*10;
%             % datapoints_meg(:,4)=1;
%             datapoints_meg=[datapoints_meg;fiducials_meg'];     
%             clear a1 b1 c1 a b c d e f g h i j k l m;
%             nuts.headshapefile=[pathname1 filename_ahs];
%         end
    end
    % Make same axis system in MRI same as that of MEG
    % (this assumes BTi's convention)
    % image orientation affects cross product sign, so check for it here:
    % value of 1 = neurological (L is L); 2 = radiological (L is R)
    
    %coreg.orientation = get(findobj('Tag','nut_orientation_menu'),'Value');
    
    origin_mri=(fiducials_mri_mm(1,:)+fiducials_mri_mm(2,:))/2;
    x_mri=fiducials_mri_mm(3,:)-origin_mri;
    x_mri=x_mri/sqrt(x_mri*x_mri');

    coreg.orientation = 1;
    
    z_mri=cross(x_mri,fiducials_mri_mm(1,:)-origin_mri);
    if(coreg.orientation == 2) z_mri = -z_mri; end;
    z_mri=z_mri/sqrt(z_mri*z_mri');
    y_mri=cross(z_mri,x_mri);
    if(coreg.orientation == 2) y_mri = -y_mri; end;
    transformation_matrix=[x_mri' y_mri' z_mri' origin_mri'];
    transformation_matrix(4,:)=0;
    transformation_matrix(4,4)=1;
    transformation_matrix_vx=transformation_matrix;    
    
    % Change the transformation matrix such that it yields points in MRI
    % space in voxels when MEG points are given in mm
    transformation_matrix_vx = inv(st.vols{1}.mat)*transformation_matrix;
    fiducials_mri_vx = [fiducials_mri_mm ones(3,1)]*inv(st.vols{1}.mat)';
    fiducials_mri_vx(:,4) = [];
    
    if(headshape)
        disp('this shouldnt be happening either!')
%         datapoints_mri_vx=transformation_matrix_vx*datapoints_meg';
%         datapoints_mri_vx=datapoints_mri_vx';    
%         datapoints_mri_vx(:,4)=[];        
%         
%         % project the points to Head surface
%         [prjpts_mm,prjpts]=nut_prjsrf(coreg.mesh,datapoints_mri_vx); %,ORIGIN,VOX);
%         % Find the final stage transformation matrix
%         %%%% FIXME: Ugly code -- nut_meg2mri_final is poorly named
%         %%%% (conflicts with nut_meg2mri) and sneakily inserts coreg.meg2mri_tfm
%         %%%% should return it as function output instead.
%         [meg2mri_tfm_vx, datapoints_final, coreg.mean_sq_err]= nut_meg2mri_final(datapoints_meg,prjpts_mm); %,VOX,ORIGIN);
%         
%         
%         % following jazz is needed because SPM's blob freaks out when given
%         % negative coordinates, and needs a dilation to know these are big voxels
%         % (chose 4 mm^3 voxel size)
%         voxelsize = [4 4 4];
%         translation_tfm = [ voxelsize(1)            0            0 0
%                                        0 voxelsize(2)            0 0
%                                        0            0 voxelsize(3) 0
%                                        0            0            0 1 ];
% 
%         hs_coord = nut_coordtfm(datapoints_final,inv(translation_tfm));
%         keep=find(prod(double(hs_coord > 0.5),2)); % discard coords with nonpositive voxels
%         hs_blobs  = hs_coord(keep,:);
%         spm_orthviews('rmblobs',1);
%         spm_orthviews('addblobs',1,hs_blobs',zeros(size(keep)),st.vols{1}.mat*translation_tfm);
%         delete(st.vols{1}.blobs{1}.cbar);
%         spm_orthviews('redraw',1);
    else
        %%% FIXME: functionize following transformation
        % following line commented cuz datapoints_meg isn't used after this
        % anyway...right??
        %datapoints_meg = (inv(transformation_matrix) * [fiducials_mri_mm [1;1;1]]')';
        if(0)  %testing if these are needed here anymore
            datapoints_mri_vx=fiducials_mri_vx;
            datapoints_mri_mm=fiducials_mri_mm;
        end
        % project the points to Head surface
        % following three lines commented, cuz we don't want to project to
        % surface at this time.
%         if(isfield(coreg,'mesh'))  % do projection only if we have optional mesh
%             [nuts.prjpts_mm,prjpts]=nut_prjsrf(coreg.mesh,datapoints_mri_vx); %,ORIGIN,VOX);
%         end
        % Find the final stage transformation matrix
%         if(0)         %  BUG: nut_meg2mri_final is munging the axes orientation... 
%             [transform_matrix_final_vx, datapoints_final, mean_sq_err]= nut_meg2mri_final(datapoints_meg,prjpts_mm);  %,VOX,ORIGIN);
%         else  % pretend like we really did nut_meg2mri_final
            coreg.mean_sq_err = NaN;  % dummy value
            coreg.meg2mri_tfm = double(transformation_matrix);
            %coreg.meg2mri_vx_tfm=transformation_matrix_vx;
%                 meg2mri_tfm_vx = transformation_matrix_vx;
%                 datapoints_final = fiducials_mri_vx;
%         end
    end
    % nuts.datapoints_mri_vx=datapoints_mri_vx;
    coreg.fiducials_mri_mm=fiducials_mri_mm;
else
    msgbox('Fiducials not yet defined',...
        'UNABLE TO LOAD DATA TO FILE data_meg_mri.mat','warn');
end
