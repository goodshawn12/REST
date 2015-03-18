function VOIvoxels=fcm_selectvoxels(nuts)
% FCM_SELECTVOXELS  lets you mark seed voxels in the MRI.
% Usages:
%  VOIvoxels = fcm_selectvoxels(nuts);
%  VOIvoxels = fcm_selectvoxels;


if isempty(findobj('tag','nutmegfig'))
    if nargin>0
        nutmeg(nuts);
    else
        nutmeg
        global nuts
        nuts.coreg=nut_CoregistrationTool(nuts.coreg);
    end
end

% [fHDR, pHDR] = uigetfile('*.hdr','Select T2 weighted MRI');
% if isequal(fHDR,0), return, end
% fullHDR=fullfile(pHDR,fHDR);
% 
% dos(['"C:\Program Files\MRIcro\mricro.exe" ' fullHDR])    %open T2 image in MRIcro

input('Set green crosshair to center of seed region. [Press return when done].','s');
VOIvoxels=nut_select_VOI;

