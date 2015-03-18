function iserr=nut_autonorm(newvoxsize)
%NUT_AUTONORM  spatially normalizes (functional) beamformer images of all subjects
%              and conditions.
% nut_autonorm(new_voxel_size)


if nargin==0, newvoxsize=[]; end

global rivets
% handles=guihandles(rivets.fig);
% viewer=get(rivets.fig,'Tag');

% set(handles.nut_threshold_menu,'Value',1)
% set(handles.nut_thresholdpos_text,'String','0'); %we want to renorm all voxels.
% set(handles.nut_thresholdneg_text,'String','0'); %we want to renorm all voxels.

voi = rivets.voi;
num2norm = length(voi.subjnr);

hw = waitbar(0,'Spatial normalization...');
iserr=true;
for k=1:num2norm
    if(~exist(voi.pathnames{k},'file') && ~exist([voi.pathnames{k} '.mat'],'file'))
        delete(hw)
        iserr=true;
        errordlg(['The file ' voi.pathnames{k} ' does not exist.']);
        return
    end
    [fpath,ffile,fext]=fileparts(voi.pathnames{k});
    currfile=fullfile(fpath,[ffile '_spatnorm.mat']);
    if ~exist(currfile,'file')
%                  set(handles.nut_getset_patientnr,'Value',find(subjlist==s))
%                  set(handles.nut_getset_conditionnr,'Value',find(condlist==c))
%                 feval(viewer,voi.pathnames{k})
%                 iserr=nut_normalize_beam(rivets.beampath,[],newvoxsize);
        iserr=nut_normalize_beam(voi.pathnames{k},[],newvoxsize);
        if iserr, delete(hw), return, end      % error test, false if normalization worked correctly
    end
    waitbar(k/num2norm,hw);
end
delete(hw)