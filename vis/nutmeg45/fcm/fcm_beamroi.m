function beam = fcm_beamroi(beam,rfile)
% FCM_BEAMROI  transforms voxel data to ROIs.
%
% Usages:
%  beam = fcm_beamroi(beam,roideffile,method)
%  beam = fcm_beamroi(beam,subjectroifile)
%  beam = fcm_beamroi(beam);
%
% roideffile        ROI definition file containing the MNI coordinates of
%                   each ROI to be analyzed.
% subjectroifile    File containing the voxel coordinates of each ROI in
%                   the given individual subject (created with
%                   FCM_VOXEL2ROI).
% If neither roideffile nor subjectroifile is given, the default ROI
% definition file AAL_ROI is used.


if nargin<2
    if ~isfield(beam,'R')
        roideffile = [fileparts(which('fcm_gui')) filesep 'templates' filesep 'AAL_ROI.mat'];
        beam.R = fcm_voxel2roi(beam.voxels,beam.coreg,'one',roideffile);
    end
else
    load(rfile)
    if exist('R','var')
        beam.R=R;
    elseif exist('ROI','var')
        clear ROI
        beam.R = fcm_voxel2roi(beam.voxels,beam.coreg,'one',rfile);
    else
        error('Second input parameter must be a ROI definition file or the subject''s ROI file.')
    end
end
  
for k=1:length(beam.s)
    if size(beam.s{k},1)>length(beam.R.goodvoxels)
        beam.s{k}=beam.s{k}(beam.R.goodvoxels,:,:);
        if k==1
            beam.voxels = beam.voxels(beam.R.goodvoxels,:);
        end
    end
    for k3=1:size(beam.s{1},3)
        beam.rois{k}(:,:,k3) = beam.R.voxel2roi_tfm.' * beam.s{k}(:,:,k3);
    end
end
 