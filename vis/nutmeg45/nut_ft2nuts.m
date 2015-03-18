function nuts = nut_ft2nuts(data,grid,mripath,mristruct,fid,mmvoxflag)
% nuts = nut_ft2nuts(data,grid,mripath,mristruct,fid,mmvoxflag)
%
% converts FieldTrip's data format into nutmeg session
%
% nuts = nut_ft2nuts(data,grid);
% 'data' is the fieldtrip structure containing raw data as well as sensor info
% 'data' can be 'raw' or 'timelock'
% 'grid' is the fieldtrip structure containing lead field and voxel coordinates
%
% mripath, mristruct, fid: (see nut_ftmriN2coreg) (may be empty)
% mmvoxflag (see nut_ftmriN2coreg) (may be empty)

% warning('hard coded for neuromag import; disable above line for import of ctf etc.')

if ft_datatype(data,'raw')
  nuts.meg.filename = 'FTraw';
elseif ft_datatype(data,'timelock')
  nuts.meg.filename = 'FTtimelock';
elseif ft_datatype(data,'source')
  % FIXME: beam=nut_ft2beam(data);
  error('call nut_ft2beam instead');
elseif any(strcmp(ft_datatype(data),{'freq' 'spike' 'comp' 'volume' 'dip' 'mvar' 'freqmvar'}))
  error(['sorry cannot handle ft_datatype ' ft_datatype(data) ' yet']);
end
if iscell(data.time)
  nuts.meg.srate=1/(data.time{1}(2)-data.time{1}(1));
else
  nuts.meg.srate=1/(data.time(2)-data.time(1));
end

if isfield(data,'fsample')
  nuts.meg.srate=data.fsample;
elseif isfield(data,'hdr') & isfield(data.hdr,'Fs')
  nuts.meg.srate = data.hdr.Fs;
elseif isfield(data,'time')
  if ft_datatype(data,'timelock')
    nuts.meg.srate=1/(data.time(2)-data.time(1))
  elseif ft_datatype(data,'raw')
    nuts.meg.srate=1/(data.time{1}(2)-data.time{1}(1))
  end
end
% data.hdr.nChans
% data.hdr.nSamples
% data.hdr.nSamplesPre
% data.hdr.nTrials


% this should take care of different units across different systems,e.g.BTi
% FIXME: keep track of FT's leadfield and other computation units as well!

if isfield(data,'grad')
  chanmm=ft_convert_units(data.grad,'mm');
  % data.grad.tra contains synthetic third gradiometer coefficients
  if isfield(data.grad,'balance')
    if strcmp(data.grad.balance.current,'none')
      nuts.meg.Gcoef = 0;
    elseif strcmp(data.grad.balance.current,'G3BR')
      nuts.meg.Gcoef = 3;
    elseif strcmp(data.grad.balance.current,'G1BR')
      nuts.meg.Gcoef = 1;
    elseif strcmp(data.grad.balance.current,'G2BR')
      nuts.meg.Gcoef = 2;
    else
      error('unexpected gradient correction')
    end
  else
    nuts.meg.Gcoef = 0;
  end
%   nuts.meg.sensorOrient = data.grad.ori; % rownorm is 1, so unitless
elseif isfield(data,'elec')
  chanmm=ft_convert_units(data.elec,'mm');
else
  warning('no channel information. include .grad or .elec?')
end
chanmm=ft_datatype_sens(chanmm);
if isfield(chanmm,'tra')
  nuts.meg.chanmixMtx{1} = chanmm.tra(match_str(data.grad.label,data.label),:)
end

nuts.meg.grad_order = 0;  % for BTi data
% FIXME: is the above line still needed?

% in FT, ref sensors are treated like regular sensors.
if any(ft_chantype(data.label,'MEGREF'))
  % FIXME: do something here
  % DO NOT USE:
  % nuts.meg.refSensorCoord;
  % nuts.meg.refSensorOrient;
  % nuts.meg.ref_sensor_labels;
  % THESE ARE NOW OBSOLETE!
  % Any reference scheme should be reflected by nuts.meg.chanmixMtx
end


if strcmp(nuts.meg.filename,'FTraw')
  % FIXME check if trials are same or variable length
  if(length(data.trial) > 1)
    nuts.meg.data = zeros(size(data.trial{1},2), size(data.trial{1},1), length(data.trial));
    for ii=1:length(data.trial)
      try
        nuts.meg.data(:,:,ii) = data.trial{ii}';
      catch
        error('you must input trials of equal length')
      end
    end
  else
    nuts.meg.data = (single(cell2mat(data.trial)))';
    data.trial = [];
  end
  nuts.meg.latency = data.time{1}'*1000;
elseif strcmp(nuts.meg.filename,'FTtimelock')
  if strcmp(data.dimord,'rpt_chan_time')
    nuts.meg.data=permute(data.trial,[3 2 1]);
  else
    error('why is timelock data not rpt_chan_time?')
  end
  nuts.meg.latency = data.time'*1000;
end


data.trial=[]; % clear for memory


nuts.coreg.mripath = which('blank.img');
nuts.coreg.meg2mri_tfm = eye(4);

if(exist('grid','var'))
  %     nuts.voxels = 10*grid.pos(grid.inside,:);
  %     nuts.voxelsize = [10 10 10]
  %     warning('SSD: doesn''t really matter, but we should calculate proper nuts.voxelsize');
  %     tempLp = grid.leadfield(grid.inside);
  %     nuts.Lp = cat(3,tempLp{:});
  [nuts.Lp,nuts.voxels,nuts.voxelsize]=nut_ftgrid2nutsLpvox(grid);
end
if isfield(nuts,'Lp')
    if isempty(nuts.Lp)
        nuts=rmfield(nuts,'Lp');
    end
end

global coreg
coreg=nuts.coreg;
if(exist('mripath','var')) & size(mripath,1)~=0
  if ~exist('fid','var')
    warning('you must also give fiducial file along with MRI file')
  end
  if ~exist('mristruct','var')
    error('you must also give mristruct along with MRI file; may be empty!')
  end
  if ~exist('mmvoxflag','var')
    error('you must also give mmvoxflag along with MRI file')
  end
  nut_ftmriN2coreg(mripath,mristruct,fid,mmvoxflag)
elseif (exist('mristruct','var'))
  nut_ftmriN2coreg([],mristruct)
else
  nut_ftmriN2coreg([]);
end
nuts.coreg=coreg;
clear global coreg
nuts.coreg.orientation=1;
nuts.meg.system=ft_senstype(data);

if(1)  % assign original names
  nuts.meg.sensor_labels=data.label';
  nuts.meg.ft_gradeleclabels = chanmm.label';
    % data.grad.label has chosen sensors as ordered in raw file
    % data.label has chosen sensors, user-specified order
    % data.hdr.label contains all channels from raw file
  
%  nuts.meg.sensorCoord = chanmm.coilpos(match_str(chanmm.label,ft_channelselection(chanmm.label,data.label)),:);
 if isfield(chanmm,'coilpos')
    nuts.meg.sensorCoord = chanmm.coilpos;
    nuts.meg.sensorOrient = chanmm.coilori;
    nuts.meg.eegflag = 0;
 end
%   nuts.meg.sensorCoord = chanmm.pnt;
  
  % ensuring order of channels is kept consistent
  nuts.meg.data=nuts.meg.data(:,match_str(nuts.meg.sensor_labels,ft_channelselection(chanmm.label,data.label)),:);
  nuts.meg.sensor_labels=nuts.meg.sensor_labels(match_str(nuts.meg.sensor_labels,ft_channelselection(chanmm.label,data.label)));
  
  switch(ft_senstype(data))
    % coordinate system may need to be transformed depending on MEG system
    case 'neuromag306'
      nuts.meg.sensorCoord(:,1) = -nuts.meg.sensorCoord(:,1);
      nuts.meg.sensorCoord = nuts.meg.sensorCoord(:,[2 1 3]);
      
      nuts.meg.sensorOrient(:,1) = -nuts.meg.sensorOrient(:,1);
      nuts.meg.sensorOrient = nuts.meg.sensorOrient(:,[2 1 3]);
    case {'ctf151' 'ctf275'}
      senres=data.hdr.orig.res4.senres;
      % next few lines copied from nut_import_CTF
      channeltypes = [senres.sensorTypeIndex];
      coilselect = find((channeltypes == 0) | (channeltypes == 1) | (channeltypes == 5));
      jj=1;
      for ii=1:length(coilselect)
        % note: parse_sensor_label is a subfunction also used in nut_import_CTF
        sensor_labels=parse_sensor_label(data.hdr.orig.res4.chanNames(coilselect(ii),:));
        if(size(senres(coilselect(ii)).pos,2) == 2) % if it's a gradiometer pair, create 2nd entry for outer coil
          nuts.meg.rawsensor_labels{jj} = [sensor_labels 'A'];
          jj = jj+1;
          nuts.meg.rawsensor_labels{jj} = [sensor_labels 'B'];
        else
          % create label for simple magnetometer
          nuts.meg.rawsensor_labels{jj} = sensor_labels;
        end
        jj = jj+1;
      end
      
  end
  if isfield(chanmm,'elecpos')
     nuts.meg.rawsensor_labels = nuts.meg.sensor_labels;
     posidx = zeros(1,length(nuts.meg.sensor_labels));
     for ii = 1:length(posidx)
         posidx(ii) = find(strcmpi(nuts.meg.sensor_labels{ii},chanmm.label));
     end
     nuts.meg.sensorCoord = chanmm.elecpos(posidx,:);
     nuts.meg.eegflag = 1;
  end
else   % assign labels based on left/right position
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % neurospin neuromag addition by baptiste
  tmp = 1;
  for j = 5:5:510
    NeuromagSensorCoord(tmp,:) = data.grad.pnt(j,:);
    NeuromagSensorCoord(tmp+1,:) = data.grad.pnt(j,:);
    NeuromagSensorCoord(tmp+2,:) = data.grad.pnt(j,:);
    tmp = tmp + 3;
  end
  
  leftchans = find(NeuromagSensorCoord(:,2) >= 0);
  rightchans = find(NeuromagSensorCoord(:,2) < 0);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % original code
  %     leftchans = find(nuts.meg.sensorCoord(:,2) >= 0);
  %     rightchans = find(nuts.meg.sensorCoord(:,2) < 0);
  for i = 1:length(leftchans)
    channel_string = int2str(i);
    nuts.meg.sensor_labels{leftchans(i)} = ['L' '0'*ones(1,4-length(channel_string)) channel_string ' (' data.grad.label{leftchans(i)} ')'];
  end
  for i = 1:length(rightchans)
    channel_string = int2str(i);
    nuts.meg.sensor_labels{rightchans(i)} = ['R' '0'*ones(1,4-length(channel_string)) channel_string ' (' data.grad.label{rightchans(i)} ')'];
  end
end

nuts.meg.goodchannels = [1:size(nuts.meg.data,2)];

% selected channels



% lead field
% nuts.meg.lsc
% nuts.meg.lsc_sensor_labels

function sensorName = parse_sensor_label(temp)

% sensorName = parse_sensor_label(temp)
% parse sensor label names

temp = deblank(temp');
temp(temp>127) = 0;
temp(temp<0) = 0;
temp = strtok(temp,char(0));
temp = strtok(temp,'-');
sensorName = char(temp)';

return


