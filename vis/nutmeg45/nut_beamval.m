function [beam_val,peak_coord] = nut_beamval(coord)
% get current s_beam magnitude

global beam

% [a,maxvoxels]=max(beam.s_beam(:,1:dsearchn(beam.timewindow,200)));
[a,maxvoxels]=max(beam.s_beam);
[b,maxtime]=max(a);

% i know, this is stupid. but it gets around a freaky MATLAB bug that
% causes an "undefined command" error when beam.voxelsMRI is manipulated
voxelsMRI = beam.voxelsMRI;

MEGvoxelindex=maxvoxels(maxtime);
voxelsMRI(MEGvoxelindex,:);

timept = maxtime;
timept2 = maxtime;
time = beam.timewindow(timept);
time2 = beam.timewindow(timept2);
peak_coord = mean(beam.s_beam(MEGvoxelindex,timept:timept2),2);

% if coord not specified, use coordinates currently selected in SPM viewer
if(~exist('coord','var'))
    coord = transpose(spm_orthviews('Pos'));
end

MEGvoxelindex = dsearchn(voxelsMRI,coord);
beam_val = beam.s_beam(MEGvoxelindex,timept);

