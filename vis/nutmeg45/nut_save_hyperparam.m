function nut_save_hyperparam(champ_out,outname,beam)

%% Set up variables
% beam input optional, but if you do input it, it should have all the fields specified below.  if you don't input beam, then nuts should exist
% and contain the necessary fields
global bolts nuts
algorithm='nut_Champagne';


if ~exist('beam','var')
    beam.params.filename = nuts.meg.filename;
    beam.params.goodchannels = nuts.meg.goodchannels;
    beam.coreg = nuts.coreg;
    beam.srate = nuts.meg.srate;
    beam.voxels = nuts.voxels;
    beam.voxelsize = nuts.voxelsize;
end

beam.params = bolts.params;
beam.params.algorithm = algorithm;
beam.params.beamformertype = algorithm;
beam.params.preprocessing = nuts.preprocessing;
beam.params.invtype = nuts.preprocessing.invtype;

%% Save time course data

beam.timepts=champ_out.timepts;
beam.sa{1}=[];
%%beam.sa{1}(:,:,1,1)=champ_out.s_th;
%%beam.sa{1}(:,:,1,2)=champ_out.s_ph;

for i=1:size(champ_out.sources,1)
    beam.sa{1}(:,:,1,i)=squeeze(champ_out.sources(i,:,:));
end
save([outname(1:end-4) '_time.mat'],'-struct','beam');

%% Save hyperparam data

beam.timepts=1;
beam.sa{1}=[];
beam.sa{1}(:,1,1,1)=champ_out.hyper1;
save([outname(1:end-4) '_hyperparam.mat'],'-struct','beam');
