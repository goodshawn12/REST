function roi = fcm_paintroi
% FCM_PAINTROI  allows painting seed region of interest 
%
%   roi = fcm_paintroi
%

global nuts st

if isempty(findobj('tag','nutmegfig'))
    if isempty(nuts)
        error('You must load session file first.')
    end
    nutmeg(nuts);
end

figure(st.fig)
texthandle=uicontrol('Style','text','BackgroundColor',[1 1 1],'String','Set crosshair to center of seed region (right click when done)','Units','normalized','Position',[.6 .55 .3 .1],'FontWeight','bold','FontSize',16);
set(st.fig,'SelectionType','normal','currentcharacter',' ')
while ~strcmp(get(st.fig, 'SelectionType'),'alt') && double(get(st.fig,'currentcharacter'))~=13  % right mouse button or return key
    pause(.2)
end

voxels=nut_select_VOI(false);  

if isequal(nuts.coreg.meg2mri_tfm,eye(4))
    roi.MRIvoxels=voxels; % coordinates are in MRI space
else
    roi.MEGvoxels=voxels; % coordinates are in MEG space
end

if isfield(nuts,'coreg')
    fcm_roiidx(roi);
end

[filename, pathname] = uiputfile('*.mat', 'Save painted ROI? (Cancel otherwise)'); 
if ~(isequal(filename,0) && ~isequal(pathname,0))
    save(fullfile(pathname, filename),'-struct','roi');
end