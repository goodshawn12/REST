function nut_fiducial_callback
% NUT_FIDUCIAL_CALLBACK
%
% This callback routine is called when any button on 'Mark_Fiducials' GUI 
% is pressed. The file contain a SWITCH construction and implementations of
% a great score of functions.
%
% Uses a lot of globals.

global coreg;
% if isempty(nuts)
%     msgbox('Load the image','Volume not found','warn')
%     return;
% end
global st;
switch get(gcbo,'Tag')
    case {'nut_left'}
        nut_leftcallback;
%         nut_coreg_fiducials(false);
        coreg=nut_coreg_fiducials(false,coreg);
    case {'nut_right'}
        nut_rightcallback;
%         nut_coreg_fiducials(false);
        coreg=nut_coreg_fiducials(false,coreg);
    case {'nut_nose'}
        nut_nosecallback;
%         nut_coreg_fiducials(false,coreg);
        coreg=nut_coreg_fiducials(false,coreg);
    case {'nut_lsc'}
        nut_lsccallback;
%         nut_coreg_fiducials(false);
        coreg=nut_coreg_fiducials(false,coreg);
    case {'nut_saveFiducials'}
        nut_saveFiducialscallback;
    case {'nut_showleft'}
        nut_showleftcallback;
    case {'nut_showright'}
        nut_showrightcallback;
    case {'nut_shownose'}
        nut_shownosecallback;
    case {'nut_showlsc'}
        nut_showlsccallback;
%     case {'nut_load_lsc'}   %maybe this could be implemented in the future?
%         nut_load_lsc;
    case {'nut_load_fiducials'}
        nut_load_fiducials;
        %nut_coreg_fiducials(false);
%     case {'nut_loadHeadshape'}
%         nut_loadHeadshape(true)
    otherwise
        warndlg('Callback not yet defined for this button',...
            'UNDEFINED CALLBACK');
end
nut_coreg_enabler


%%--------------------------------------------
function nut_leftcallback
global screws
%get the cursor value
u=spm_orthviews('Pos');
set(screws.coreg.handles.nut_left_text,'String',sprintf('MRI: %.1f %.1f %.1f',u));
screws.l=u';
if (~isempty(get(screws.coreg.handles.nut_nose_text,'String')) && ~isempty(get(screws.coreg.handles.nut_right_text,'String')))
    global coreg
    coreg.fiducials_mri_mm(1,:) = u';
    coreg.fiducials_mri_mm(2,:) = sscanf(get(screws.coreg.handles.nut_right_text,'String'), 'MRI: %g %g %g');
    coreg.fiducials_mri_mm(3,:) = sscanf(get(screws.coreg.handles.nut_nose_text,'String'), 'MRI: %g %g %g');
end
%nut_coreg_enabler


%%--------------------------------------------
function nut_rightcallback
global screws
%get the cursor value
u=spm_orthviews('Pos');
set(screws.coreg.handles.nut_right_text,'String',sprintf('MRI: %.1f %.1f %.1f',u));
screws.r=u';
if (~isempty(get(screws.coreg.handles.nut_left_text,'String')) && ~isempty(get(screws.coreg.handles.nut_nose_text,'String')))
    global coreg
    coreg.fiducials_mri_mm(1,:) = sscanf(get(screws.coreg.handles.nut_left_text,'String'), 'MRI: %g %g %g');
    coreg.fiducials_mri_mm(2,:) = u';
    coreg.fiducials_mri_mm(3,:) = sscanf(get(screws.coreg.handles.nut_nose_text,'String'), 'MRI: %g %g %g');
end
%nut_coreg_enabler


%%--------------------------------------------
function nut_nosecallback
global screws
%get the cursor value
u=spm_orthviews('Pos');
set(screws.coreg.handles.nut_nose_text,'String',sprintf('MRI: %.1f %.1f %.1f',u));
screws.n=u';
if (~isempty(get(screws.coreg.handles.nut_left_text,'String')) && ~isempty(get(screws.coreg.handles.nut_right_text,'String')))
    global coreg
    coreg.fiducials_mri_mm(1,:) = sscanf(get(screws.coreg.handles.nut_left_text,'String'), 'MRI: %g %g %g');
    coreg.fiducials_mri_mm(2,:) = sscanf(get(screws.coreg.handles.nut_right_text,'String'), 'MRI: %g %g %g');
    coreg.fiducials_mri_mm(3,:) = u';
end
%nut_coreg_enabler

%%--------------------------------------------
function nut_lsccallback
%get the cursor value
global coreg
%u=nut_mri2meg(spm_orthviews('Pos')');
u=nut_coordtfm(spm_orthviews('Pos')',inv(coreg.meg2mri_tfm));
if (u(1)==0 & u(2)==0 & u(3)==0)
    msgbox('Please load the image','Image not loaded','error')
    return
else
%     hx_obj = findobj('Tag','nut_lsc_text');
%     %set the static text box with the above value
%     set(hx_obj,'String',sprintf('MEG: %.1f %.1f %.1f',u));
    % Switch the Show button ON
    global nuts
    nuts.meg.lsc = u;
    if(isfield(nuts.meg,'lsc_sensor_labels'))
        nuts.meg=rmfield(nuts.meg,'lsc_sensor_labels');
    end
    clear nuts;

    nut_coreg_enabler
end


%%--------------------------------------------
function nut_saveFiducialscallback
%Load the "data_meg_mri.mat" file in the present function
global coreg;
if isfield(coreg,'datapoints_final_iter')
    %delete "datapoints_final_iter" from "data_meg_mri.mat" file
    %clear datapoints_final_iter;
    coreg = rmfield(coreg,'datapoints_final_iter');       %%%is this necessary?
end
global st;
%check if Image is loaded or not
if (isempty(st)==1 | isempty(st.vols{1})==1)
    msgbox('Please load the image','Image not loaded','error')
    return
end
%Get the fiducials from Mark/Edit fiducials GUI text box
l = sscanf(get(findobj('Tag','nut_left_text'),'String'), 'MRI: %g %g %g');
r = sscanf(get(findobj('Tag','nut_right_text'),'String'), 'MRI: %g %g %g');
n = sscanf(get(findobj('Tag','nut_nose_text'),'String'), 'MRI: %g %g %g');
%check if the fiducials text box is empty
if strcmp(l,'')==0 & strcmp(r,'')==0 & strcmp(n,'')==0
    fiducials_mri_mm=[l';r';n'];
    [filename1 pathname1] = uiputfile('*.txt','Save the fiducials to desired file (OPTIONAL)...');
    if ischar(filename1)
        %Save the fiducials to desired file (OPTIONAL)
        fid=fopen([pathname1 filename1],'w');
        fprintf(fid,'%g %g %g\n %g %g %g\n %g %g %g',fiducials_mri_mm(1,1),fiducials_mri_mm(1,2),...
            fiducials_mri_mm(1,3),fiducials_mri_mm(2,1),fiducials_mri_mm(2,2),fiducials_mri_mm(2,3),...
            fiducials_mri_mm(3,1),fiducials_mri_mm(3,2),fiducials_mri_mm(3,3));
        fclose(fid);
    end
end


%%--------------------------------------------


%%--------------------------------------------
function nut_showleftcallback
% take cross hair to the specified position
l = sscanf(get(findobj('Tag','nut_left_text'),'String'), 'MRI: %g %g %g');
spm_orthviews('Reposition',l);

%%--------------------------------------------
function nut_showrightcallback
% take cross hair to the specified position
l = sscanf(get(findobj('Tag','nut_right_text'),'String'), 'MRI: %g %g %g');
spm_orthviews('Reposition',l);

%%--------------------------------------------
function nut_shownosecallback
% take cross hair to the specified position
l = sscanf(get(findobj('Tag','nut_nose_text'),'String'), 'MRI: %g %g %g');
spm_orthviews('Reposition',l);

%%--------------------------------------------
function nut_showlsccallback
% take cross hair to the specified position
%l = transpose(sscanf(get(findobj('Tag','nut_lsc_text'),'String'), 'MEG: %g %g %g'));

% in the case of multiple spheres, shows average location for sphere center
global nuts
l = mean(nuts.meg.lsc,1);

% nut_meg2mri doesn't work here, since we haven't yet inserted meg2mri_tfm
% into nuts.coreg...
% spm_orthviews('Reposition',nut_meg2mri(l));

% so this is a silly workaround, but i'm not losing any sleep over it
global coreg
spm_orthviews('Reposition',nut_coordtfm(l,coreg.meg2mri_tfm));

%%--------------------------------------------
