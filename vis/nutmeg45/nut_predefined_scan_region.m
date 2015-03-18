function nut_predefined_scan_region% function NUT_PREDEFINED_SCAN_REGION%
% Loads in a pre-defined mask file (generally an cortical white and gray% matter mask) for restricting the beamforming algorithm% 01/28/04 : WADE : Wrote it.%
% See also nut_view_scan_region%
% Uses ST and WORKFILE globals.
global st;global workfile; % Check for existing anatomy image data
 if (isempty(st)==1 | isempty(st.vols{1})==1)
   msgbox('Please load the image','Image not loaded','error')
   return
 end 
% get volume from the 3D image
 V=spm_vol(st.vols{1}.fname);
 Y = spm_read_vols(V);
% We read in the 3d volume at this point. Is this necessary? [classFileName, classPathName] = uigetfile( ...       {'*.img'}, ...        'Pick an analyze class file');     % Analyze class files are binary analyze format files. They must be the    % same size as the current image and contain a mask of 1s and 0s that    % define where the beamforming is to take place.     [classFilePath,classFileNameBase,ext,dummy2]=fileparts(fullfile(classPathName,classFileName));   % Load it in using spm_volfn=fullfile(classPathName,classFileName) [classV]=spm_vol(fn); fprintf('\nReading class file\n');  [classImage]=spm_read_vols(classV);  % This is a binary image. Immediately binarize it... classImage(classImage~=0)=1; % Mask Y ( remember 'Y'? ) with classImage % Do some bounds checking here first...  Y=Y.*classImage; % Generate a header V : just like the original image but with a different % name outVol=spm_vol(st.vols{1}.fname); [imFilePath,imFileName,dummy1,dummy2]=fileparts(st.vols{1}.fname);  outVol.fname=fullfile(imFilePath,'rawdata'); spm_write_vol(outVol,Y);
% Update the viewer window
spm_image('init',fullfile(imFilePath,'rawdata.img'));
%spm_orthviews('Reposition',MRIpos);
%set(findobj('Tag','nut_view_VOI_button'),'Enable','on');
set(findobj('Tag','nut_lead_field_button'),'Enable','on');VOI=(classImage==1);
save('-APPEND',workfile,'VOI');