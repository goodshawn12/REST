function nut_coreg_enabler
% figures out what information is loaded into NUTMEG and
% enables/disables GUI buttons/menus as necessary

global coreg screws
handles = screws.coreg.handles;

if isfield(coreg,'orientation')
    set(handles.nut_orientation_menu,'Value',coreg.orientation);
end

% get handlenames of all buttons and menus, place in enableditems
% i.e., begin with everything enabled
handlenames = fieldnames(handles);
enableditems = {};
if ~isfield(coreg,'mripath')
    coreg.mripath = which('blank.img');
    %errordlg('yikes...what happened to your MRI?  you should at least have the original blank.img somewhere.  you must lose your keys often')
end
[mrifilepath,mrifilename]=fileparts(coreg.mripath);
% textitems={};

for ii = 2:(length(handlenames)-1)
    % overly complicated way of testing if the handlename contains "button" or "menu"
    if(~isempty(regexpi(get(handles.(handlenames{ii}),'Style'),'button')) || ~isempty(regexpi(get(handles.(handlenames{ii}),'Style'),'menu')))
        enableditems{end+1} = handlenames{ii};
    end
%     if (~isempty(regexpi(get(handles.(handlenames{ii}),'Style'),'text')))
%         textitems{end+1}=handlenames{ii};
%     end

end

% now selectively disable:
disableditems = {};
if(~isfield(coreg,'fiducials_mri_mm'))
    %disableditems=union(disableditems,{'nut_shownose','nut_showright','nut_showleft','nut_saveFiducials'});
%     disableditems=union(disableditems,{'nut_adjust_coreg','nut_edit_fiducials_3d','nut_show_fiducials_3d','nut_verify_coregistry_in_3d'});
    disableditems=union(disableditems,{'nut_edit_fiducials_3d','nut_show_fiducials_3d','nut_verify_coregistry_in_3d'});
    if ~isfield(screws,'l'), set(handles.nut_left_text,'String',[]);end
    if ~isfield(screws,'r'), set(handles.nut_right_text,'String',[]);end
    if ~isfield(screws,'n'), set(handles.nut_nose_text,'String',[]);end
else
    set(handles.nut_left_text,'String',['MRI: ' num2str(coreg.fiducials_mri_mm(1,:),'%.1f ')]);
    set(handles.nut_right_text,'String',['MRI: ' num2str(coreg.fiducials_mri_mm(2,:),'%.1f ')]);
    set(handles.nut_nose_text,'String',['MRI: ' num2str(coreg.fiducials_mri_mm(3,:),'%.1f ')]);
end

if regexpi(mrifilename,[filesep 'blank.img'])
    disableditems=union(disableditems,{'nut_genmesh','nut_load_fiducials','nut_nose','nut_left','nut_right'});
    %disableditems=union(disableditems,{'nut_loadHeadshape'});
end

if isempty(get(handles.nut_nose_text,'String'))
    disableditems=union(disableditems,{'nut_shownose','nut_saveFiducials'});
end
if isempty(get(handles.nut_right_text,'String'))
    disableditems=union(disableditems,{'nut_showright','nut_saveFiducials'});
end
if isempty(get(handles.nut_left_text,'String'))
    disableditems=union(disableditems,{'nut_showleft','nut_saveFiducials'});
end

global nuts
if(isfield(nuts,'meg') && isfield(nuts.meg,'lsc'))
    if(size(nuts.meg.lsc,1)>1)
        set(handles.nut_lsc_text,'String','multiple spheres');
    else
        set(handles.nut_lsc_text,'String',['MEG: ' num2str(nuts.meg.lsc,' %.1f ')]);
    end
else
    set(handles.nut_lsc_text,'String',[]);
end
clear nuts % does NOT erase global nuts, just prevents inadvertent mishaps after this point
if isempty(get(handles.nut_lsc_text,'String'))
    %disableditems=union(disableditems,{'nut_showlsc','nut_saveFiducials'});
    disableditems=union(disableditems,{'nut_showlsc'});
end

if ~isfield(coreg,'mesh')
    disableditems=union(disableditems,{'nut_edit_fiducials_3d','nut_show_fiducials_3d','nut_verify_coregistry_in_3d','nut_adjust_coreg','nut_coreg_firstpass'});
    %disableditems=union(disableditems,{'nut_loadHeadshape'});
end

if ~isfield(coreg,'hsCoord')
    disableditems=union(disableditems,{'nut_adjust_coreg','nut_coreg_firstpass'});
end

if isfield(coreg,'orientation')
    if get(handles.nut_orientation_menu,'Value')==2
        disableditems=union(disableditems,{'nut_normMRI'});
        if isfield(coreg,'norm_mripath')
           coreg = rmfield(coreg,'norm_mripath');
        end
        set(handles.nut_normimage_nametext,'String',[]);
    end
end
    

% remove explicitly disabled items from enableditems
enableditems = setdiff(enableditems,disableditems);


% enable enableditems
for i=1:length(enableditems)
    set(handles.(enableditems{i}),'Enable','On')
end

% disable disableditems
for i=1:length(disableditems)
    set(handles.(disableditems{i}),'Enable','Off')
end
