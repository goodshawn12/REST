function [roi,beam] = fcm_functroi(radius,posneg)
% roi = fcm_functroi(radius,posneg)
%
% nut_results_viewer must be open, threshold selected, and relevant time
% points/freqs selected.
%
% radius        in mm from crosshair
% posneg        set to 'pos' or 'neg' to include positive or negative voxels only.

global rivets beam nuts
posneg = strcmpi(posneg,'neg') +1;

dist = nut_rownorm(nut_coord_diff(rivets.voxelsMRI,rivets.voxelsMRI(rivets.currMEGvoxelindex,:)));
isgood = ( (dist <= radius) & rivets.threshold(:,rivets.timeselect,posneg) );

roi.MNIvoxels  = nut_mri2mni(nut_meg2mri(beam.voxels(isgood,:)));
roi.timeselect = rivets.timeselect;
roi.freqselect = rivets.freqselect;

if all(isfield(nuts,{'voxels' 'coreg'}))
    fcm_roiidx(roi);
end

if nargout>1
    beam=struct('s',{{ones(size(roi.MNIvoxels,1),1)}},'voxels',beam.voxels(isgood,:),'voxelsize',beam.voxelsize, ...
        'timepts',1,'timewindow',[0 2],'bands',[0 1],'srate',1,'coreg',beam.coreg);
    %beam= beam2;
    save s_beam_roitemplate beam
    tv s_beam_roitemplate
end