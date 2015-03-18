function nuts = nut_spm2nuts(matfile)
% converts SPM's MEEG data format into nutmeg session
% nuts = nut_fcdc2nuts(data,grid);
% 'data' is the fieldtrip structure containing raw data as well as sensor info
% 'grid' is the fieldtrip structure containing lead field and voxel coordinates
%
% NOTE: probably not finished yet, just a template!   jz

nuts.meg.filename = 'FTdata'
nuts.meg.srate = data.hdr.Fs;
% data.hdr.nChans
% data.hdr.nSamples
% data.hdr.nSamplesPre
% data.hdr.nTrials


nuts.meg.sensorCoord = data.grad.pnt*1000;
nuts.meg.Gcoef = 0;
% BTi data in meters, but other data formats might have different units???
% this could be serious problem with fieldtrip's leadfield computation
% probably *does not* scale!

% data.grad.tra contains synthetic third gradiometer coefficients
nuts.meg.grad_order = 0;  % for BTi data
% nuts.meg.refSensorCoord;
% nuts.meg.refSensorOrient;
% nuts.meg.ref_sensor_labels;



nuts.meg.sensorOrient = data.grad.ori;

if(length(data.trial) > 1)
    nuts.meg.data = zeros(size(data.trial{1},2), size(data.trial{1},1), length(data.trial));
    for ii=1:length(data.trial)
        nuts.meg.data(:,:,ii) = data.trial{ii}';
    end
else
    nuts.meg.data = (single(cell2mat(data.trial)))';
    data.trial = [];
end

% nuts.meg.data = cat(3,data.trial{:});
% nuts.meg.data = permute(nuts.meg.data,[2 1 3]);
nuts.meg.latency = data.time{1}'*1000; 

nuts.meg.goodchannels = [1:size(nuts.meg.data,2)];

nuts.coreg.mripath = which('blank.img');
nuts.coreg.meg2mri_tfm = eye(4);


if(exist('grid','var'))
    nuts.voxels = 10*grid.pos(grid.inside,:);
    nuts.voxelsize = [10 10 10]
    warning('SSD: doesn''t really matter, but we should calculate proper nuts.voxelsize');
    tempLp = grid.leadfield(grid.inside);
    nuts.Lp = cat(3,tempLp{:});
end

if(0)  % assign original names
    nuts.meg.sensor_labels = data.grad.label;
    % data.grad.label has chosen sensors as ordered in raw file
    % data.label has chosen sensors, user-specified order
    % data.hdr.label contains all channels from raw file
else   % assign labels based on left/right position

    leftchans = find(nuts.meg.sensorCoord(:,2) >= 0);
    rightchans = find(nuts.meg.sensorCoord(:,2) < 0);
    for i = 1:length(leftchans)
        channel_string = int2str(i);
        nuts.meg.sensor_labels{leftchans(i)} = ['L' '0'*ones(1,4-length(channel_string)) channel_string ' (' data.grad.label{leftchans(i)} ')'];
    end
    for i = 1:length(rightchans)
        channel_string = int2str(i);
        nuts.meg.sensor_labels{rightchans(i)} = ['R' '0'*ones(1,4-length(channel_string)) channel_string ' (' data.grad.label{rightchans(i)} ')'];
    end
end

% selected channels



% lead field
% nuts.meg.lsc
% nuts.meg.lsc_sensor_labels
