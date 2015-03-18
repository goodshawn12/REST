function nut_save_likelihood(lmap,algorithm,outname,beam);
% beam input optional, but if you do input it, it should have all the fields specified below.  if you don't input beam, then nuts should exist
% and contain the necessary fields (not the case for saketini).
global bolts nuts
beam.timepts=1;
% beam.s_ph=0;
% beam.s_z=0;
beam.params = bolts.params;
beam.params.algorithm = algorithm;
beam.params.beamformertype = algorithm;
if ~exist('beam')
beam.params.filename = nuts.meg.filename;
beam.params.goodchannels = nuts.meg.goodchannels;
beam.coreg = nuts.coreg;
beam.srate = nuts.meg.srate;
beam.voxels = nuts.voxels;
beam.voxelsize = nuts.voxelsize;
end

beam.params.preprocessing = nuts.preprocessing;
beam.params.invtype = nuts.preprocessing.invtype;
beam.sa{1}=sqrt(lmap(:)-min(lmap(:)));

save([outname(1:end-4) '_lqv' algorithm(4:end) '.mat'],'-struct','beam');

%     beam.s_th=sqrt(lqvs_0(:)-min(lqvs_0(:)));
%     save(['s_beam_lqvs_0_' params.thisrun '.mat'],'beam')
%     beam.s_th=sqrt(lratmap(:)-min(lratmap(:)));
%     save(['s_beam_lrats_' params.thisrun '.mat'],'beam')
