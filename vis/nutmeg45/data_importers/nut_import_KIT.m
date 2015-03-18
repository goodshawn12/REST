function nut_import_kit(megdata_filename);

% this might be done better by Fieldtrip...to investigate further

global nuts ndefaults

meg=load('-ASCII',megdata_filename);

% drop last 10 characters from selected filename, so we can construct filenames of parameter files
kit_basename = megdata_filename(1:(end-10));

sensor_filename = [kit_basename 'sns.txt'];

sensor2head = kit_matrix([kit_basename 'matrix.txt']);   % read sensor2head transform matrix
[lsc_kit r] = kit_sphere([kit_basename 'sphere.txt']);

sensor2head(4,:) = [0 0 0 1]; % augment transformation matrix to make 4x4


sensor_filename

fid=fopen(sensor_filename);

lineL = fgetl(fid);
lineL = fgetl(fid);
lineL = fgetl(fid);

%  amdmax=192;
% [channelnum,channeltype,coil_th,coil_ph,coilsize,baseline]=deal(zeros(amdmax,1));
% coil_coord = zeros(amdmax,3);

ichannel = 0;
while 1
    lineL = fgetl(fid);
    if(~ischar(lineL)),break,end;
    [crap,channeltype] = strread(lineL,'%f%s',1);
    if(strcmp(channeltype,'AxialGradioMeter'))
        ichannel = ichannel+1;
        [channelnum(ichannel,1),channeltype(ichannel,1),coil_coord(ichannel,1),coil_coord(ichannel,2),coil_coord(ichannel,3),coil_th(ichannel,1),coil_ph(ichannel,1),coilsize(ichannel,1),baseline(ichannel,1)] = strread(lineL,'%f%s%f%f%f%f%f%f%f',1);
    end
end

fclose(fid);

meg = meg(:,channelnum+1);  % discard bogus channels

baseline=baseline(1);  % all baseline values should be the same; just keep first one

coil_th=coil_th*pi/180;  % convert to radians
coil_ph=coil_ph*pi/180;  % convert to radians
normal_direction = [sin(coil_th).*cos(coil_ph) sin(coil_th).*sin(coil_ph) cos(coil_th)];
coil_coord2 = coil_coord + normal_direction.*baseline;

%%%% convert from KIT sensor coords to head coords -- should already be in millimeters
coil_coord = nut_coordtfm(coil_coord,sensor2head);
coil_coord2 = nut_coordtfm(coil_coord2,sensor2head);

% convert KIT coord convention to BTi/CTF
% (x_ctf = -y_KIT; y_ctf = x_KIT; z_ctf = z_KIT)
lsc = [-lsc_kit(2) lsc_kit(1) lsc_kit(3)];
nuts.meg.lsc=lsc;
coil_coord = [-coil_coord(:,2) coil_coord(:,1) coil_coord(:,3)];
coil_coord2 = [-coil_coord2(:,2) coil_coord2(:,1) coil_coord2(:,3)];

normal_direction = coil_coord2-coil_coord;
normal_direction = normal_direction./repmat(nut_rownorm(normal_direction),[1 3]);
% normal_direction3 = nut_coordtfm(normal_direction,[sensor2head(:,1:3) [0;0;0;1]]);


% WE NEED TO PROMPT USER FOR LATENCY INFO
prompt   = 'Enter starting latency and sampling rate (e.g., [-100 256] for 256 Hz starting at -100 ms)';
title    = 'Input latency info';
lines = 1;
def{1}   = num2str('-100 1000');

answer   = inputdlg(prompt,title,lines,def);
if (isempty(answer)) errordlg('MEG loading cancelled. Freak.');return; end;
timeinfo = str2num(answer{1});

start_time = timeinfo(1);
srate = timeinfo(2);

no_points = size(meg,1);
latency = (1000*(0:(no_points-1))/srate) + start_time;  % in ms
if ndefaults.meg.single
    nuts.meg.data = meg;
else
    %     nuts.meg.data=double(nuts.meg.data);
    nuts.meg.data=double(meg);
end

%%% lsc now specified in coregistration for CTF, we still do this line up
%%% in BTi and KIT code though.
% nuts.meg.lsc = lsc;
nuts.meg.sensorCoord(:,:,1) = coil_coord;
if(exist('coil_coord2','var'))
    nuts.meg.sensorCoord(:,:,2) = coil_coord2;
end
if ~exist('sensor_labels','var')
    leftchans = find(coil_coord(:,2) >= 0);
    rightchans = find(coil_coord(:,2) < 0);
    for i = 1:length(leftchans)
        channel_string = int2str(i);
        sensor_labels{leftchans(i)} = ['L' '0'*ones(1,4-length(channel_string)) channel_string];
    end
    for i = 1:length(rightchans)
        channel_string = int2str(i);
        sensor_labels{rightchans(i)} = ['R' '0'*ones(1,4-length(channel_string)) channel_string];
    end
end
if size(sensor_labels,1) > size(sensor_labels,2)
    sensor_labels = sensor_labels'; % must be a row vector
end
nuts.meg.sensor_labels = sensor_labels;

nuts.meg.sensorOrient = normal_direction;

if(exist('Gcoef','var'))  % store the good bits for synthetic gradient correction
    nuts.meg.Gcoef = Gcoef;
    nuts.meg.grad_order = grad_order;
    nuts.meg.refSensorCoord(:,:,1) = ref_coord;
    nuts.meg.refSensorCoord(:,:,2) = ref_coord2;

    % refGradiometerChans = find(all(ref_coord2,2)); % find reference coils that are gradiometers, i.e., ref_coord2(chan,:) ~= [0 0 0])
    refMagnetometerChans = find(~any(ref_coord2,2)); % find reference coils that are magnetometers, i.e., ref_coord2(chan,:) == [0 0 0])
    nuts.meg.refSensorCoord(refMagnetometerChans,:,2) = NaN;  % magnetometers don't have second coil, so set to NaN

    nuts.meg.refSensorOrient = ref_direction;
    nuts.meg.ref_sensor_labels = ref_sensor_labels;
else
    nuts.meg.Gcoef = 0;
    nuts.meg.grad_order = 0;
end

nuts.meg.type=repmat({'meg'},size(nuts.meg.data,2),1);

nuts.meg.goodchannels = 1:size(nuts.meg.data,2);
nuts.meg.srate = srate;
nuts.meg.latency = latency';   % transpose for ease of use later

return;


function meg2norm=kit_matrix(matrix_name)
%read KIT sensor coord to head coord transformation matrix
disp('point conversion info file')
%      [filename, pathname, filterindex] = uigetfile('*matrix.txt', 'conversion info file');
%       matrix_name=[pathname,filename];
fid=fopen(matrix_name);

str = fgetl(fid);

while  length(str)~=length('[Transform Matrix]')|str~='[Transform Matrix]'
    str = fgetl(fid);
end
str1 = fgetl(fid);
str2 = fgetl(fid);
str3 = fgetl(fid);

ku1=findstr(str1,'=')+1;
ku246=findstr(str1,'*')-1;, ku2=ku246(1);, ku4=ku246(2);, ku6=ku246(3);
ku357=findstr(str1,'+')+1;, ku3=ku357(1);, ku5=ku357(2);, ku7=ku357(3);
ku8=findstr(str1,'[')-1;

mat11=str2num( str1(ku1:ku2) );, mat12=str2num( str1(ku3:ku4) );,
mat13=str2num( str1(ku5:ku6) );, mat14=str2num( str1(ku7:ku8) );

ku1=findstr(str2,'=')+1;
ku246=findstr(str2,'*')-1;, ku2=ku246(1);, ku4=ku246(2);, ku6=ku246(3);
ku357=findstr(str2,'+')+1;, ku3=ku357(1);, ku5=ku357(2);, ku7=ku357(3);
ku8=findstr(str2,'[')-1;

mat21=str2num( str2(ku1:ku2) );, mat22=str2num( str2(ku3:ku4) );,
mat23=str2num( str2(ku5:ku6) );, mat24=str2num( str2(ku7:ku8) );


ku1=findstr(str3,'=')+1;
ku246=findstr(str3,'*')-1;, ku2=ku246(1);, ku4=ku246(2);, ku6=ku246(3);
ku357=findstr(str3,'+')+1;, ku3=ku357(1);, ku5=ku357(2);, ku7=ku357(3);
ku8=findstr(str3,'[')-1;

mat31=str2num( str3(ku1:ku2) );, mat32=str2num( str3(ku3:ku4) );,
mat33=str2num( str3(ku5:ku6) );, mat34=str2num( str3(ku7:ku8) );

meg2norm=[mat11 mat12 mat13 mat14;mat21 mat22 mat23 mat24;mat31 mat32 mat33 mat34];

fclose(fid);


function [lsc,r]= kit_sphere(spferef_name)
% read KIT sphere origin
disp('point sphere info file');
%     [filename, pathname, filterindex] = uigetfile('*sphere.txt', 'sphere info file');
%     spferef_name=[pathname,filename];
fid=fopen(spferef_name);
%
for i=1:4
    str = fgetl(fid);
end
str = fgetl(fid);
idx1 = findstr(str,'center');
idx2 = findstr(str,',');
idx3 = findstr(str,')');
idx4 = findstr(str,'radius');
idx5 = findstr(str,'[');

center_x = str2num(str(idx1(1,1)+8:idx2(1,2)-1));
center_y = str2num(str(idx2(1,2)+1:idx2(1,3)-1));
center_z = str2num(str(idx2(1,3)+1:idx3(1,1)-1));
r = str2num(str(idx4(1,1)+7:idx5(1,1)-1));

% fprintf('lsc=[%f,%f,%f] \n head_radius=%f \n',center_x,center_y,center_z,radius);
fclose(fid);
lsc=[center_x center_y center_z];   % should be in mm



