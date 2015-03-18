function [lfp, voxels] = nut_ftopenmeeg(voxelsize,varargin)
global nuts

% Generate voxels in MRI coords
voxels = nut_make_voxels(voxelsize);
voxels_mrim = nut_meg2mri(voxels)/1000; % voxels in spm-m

if isfield(nuts.coreg,'volfile') && exist(nuts.coreg.volfile,'file')
    load(nuts.coreg.volfile);
    [path,file,~] = fileparts(nuts.coreg.volfile);
else
    [file,path] = uigetfile('*_vol.mat','Open Head Surface Mesh File');
    load(fullfile(path,file));
    [path,file,~] = fileparts([path file]);
end
if length(file) > 4 && strcmp(file(end-3:end),'_vol'); file = file(1:end-4); end;
vol_tmp = vol; clear vol;

% Convert nutmeg vol (already in mri coords) to fieldtrip style
for ii = 1:length(vol_tmp.bnd)
    if isfield(vol_tmp.bnd(ii),'tri'); vol_tmp.bnd(ii).faces = vol_tmp.bnd(ii).tri; end;
    if isfield(vol_tmp.bnd(ii),'pnt'); vol_tmp.bnd(ii).vertices = vol_tmp.bnd(ii).pnt; end;
    vol.bnd(ii).tri = vol_tmp.bnd(ii).faces;
    vol.bnd(ii).pnt = vol_tmp.bnd(ii).vertices;
end
vol.cond = vol_tmp.cond;
vol.type = 'openmeeg';
vol.basefile = ft_getopt(vol_tmp,'basefile',file);
vol.path = ft_getopt(vol_tmp,'path',path);
vol.unit = ft_getopt(vol_tmp,'unit','mm');

if ~nuts.meg.eegflag
    sens.type = 'meg';
    sens.coilpos = nuts.meg.sensorCoord;
    sens.coilori = nuts.meg.sensorOrient;
else
    sens.type = 'eeg';
    sens.elecpos = nuts.meg.sensorCoord;
    sens.chanpos = nuts.meg.sensorCoord;
end
if isfield(nuts.meg,'chanmixMtx'); sens.tra = nuts.meg.chanmixMtx{1}; end;
sens.label = nuts.meg.sensor_labels;
sens.unit = 'mm';
sens = ft_transform_sens(nuts.coreg.meg2mri_tfm, sens); % Convert sens to spm coords

vol = ft_convert_units(vol,'m');
sens = ft_convert_units(sens,'m');

lfp = ft_leadfield_openmeeg (voxels_mrim, vol, sens);
lfp = reshape(lfp,[size(lfp,1), 3, size(lfp,2)/3]);