function [Lp,voxels]=nut_compute_lead_field(voxelsize,lfcomp,chanmixMtx);
% NUT_COMPUTE_LEAD_FIELD
%
% note: 'voxelsize' is either scalar (isotropic voxels)
%       ELSE it will be interpreted as a voxel list 'voxels'
%
% computes a vector lead field (forward solution)
% lfcomp = # of desired components (2 or 3 -- single sphere MEG only)
% 

global nuts ndefaults

if ~exist('lfcomp')
    lfcomp = 3
end

disp('Computing lead field...');

if(all(size(voxelsize)==[1 1])) % if voxel list supplied as argument
    voxels=nut_make_voxels(voxelsize);
else
    warning('this is a total hack; ideally we should do nut_make_voxels before nut_compute_lead_field call')
    warning('the ''voxelsize'' input is being used as your voxel list')
    voxels = voxelsize;
end

lsc = nuts.meg.lsc;

if ((size(lsc,1)>1) & (lfcomp==2))
    lfcomp=3;
    disp('local multispheres detected, thus changing lfcomp to 3.  is this what you want?')
end

tic;

% if 2 component
switch lfcomp
    case 23
        % three-component intracranial EEG lead field
        [lx1, ly1, lz1] = nut_leadf_iEEG(voxels,nuts.meg.sensorCoord);
        [lx2, ly2, lz2] = deal(zeros(size(nuts.meg.sensorCoord,1),size(voxels,1)));
        [ref_lx1,ref_lx2,ref_ly1,ref_ly2,ref_lz1,ref_lz2] = deal(0);

    case 1  % scalar lead field, to implement someday
        error('Not implemented yet.')
    case 2  % 2-component theta-phi (spherical w/o radius) coordinate system; only makes sense for MEG single sphere
      % map lsc entries to channels...
        if(size(lsc,1) == 1) % single sphere means all channels use same sphere center (i.e., the only entry in nuts.meg.lsc)
            indexmap = ones(size(nuts.meg.sensorCoord,1),1);
        else % multisphere means each channel has its own sphere center (resulting lead field is necessarily rank 3!)
            indexmap=zeros(size(nuts.meg.sensorCoord,1),1);
            for i=1:size(nuts.meg.sensorCoord,1)
                idxtmp=strmatch(nuts.meg.rawsensor_labels{i},nuts.meg.lsc_sensor_labels);
                if(isempty(idxtmp)) % in case of gradiometer channel with extra character
                    indexmap(i)=strmatch(nuts.meg.rawsensor_labels{i}(1:(end-1)),nuts.meg.lsc_sensor_labels);
                else
                    indexmap(i)=idxtmp;
                end
            end
        end

        Nchannels = size(nuts.meg.sensorCoord,1);
        for i=1:Nchannels
            nut_progress(i,Nchannels,2);
            [lx(i,:),ly(i,:)] = nut_leadf_meg2comp(voxels,nuts.meg.sensorCoord(i,:,1),nuts.meg.sensorOrient(i,:),lsc(indexmap(i),:));
        end
    case 3  % 3-component Cartesian coordinate system
        % map lsc entries to channels...
        if(size(lsc,1) == 1) % single sphere means all channels use same sphere center (i.e., the only entry in nuts.meg.lsc)
            indexmap = ones(size(nuts.meg.sensorCoord,1),1);
        else % multisphere means each channel has its own sphere center (resulting lead field is necessarily rank 3!)
            indexmap=zeros(size(nuts.meg.sensorCoord,1),1);
            for i=1:size(nuts.meg.sensorCoord,1)
                idxtmp=strmatch(nuts.meg.rawsensor_labels{i},nuts.meg.lsc_sensor_labels);
                if(isempty(idxtmp)) % in case of gradiometer channel with extra character
                    indexmap(i)=strmatch(nuts.meg.rawsensor_labels{i}(1:(end-1)),nuts.meg.lsc_sensor_labels);
                else
                    indexmap(i)=idxtmp;
                end
            end
        end

        Nchannels = size(nuts.meg.sensorCoord,1);
        for i=1:Nchannels
            nut_progress(i,Nchannels,2);
            [lx(i,:),ly(i,:),lz(i,:)] = nut_leadf_meg3comp(voxels,nuts.meg.sensorCoord(i,:,1),nuts.meg.sensorOrient(i,:),lsc(indexmap(i),:));
        end
    otherwise
        error('Lead fields cannot have more than 3 components. NUTMEG does not support parallel universes.')
end

% multiplying by chanmixMtx forms gradiometers and, if "synthetic gradient correction"
% is present, also takes into account reference sensor contributions
if(exist('lz','var'))
    Lp(:,3,:) = chanmixMtx * lz;
end
Lp(:,2,:) = chanmixMtx * ly;
Lp(:,1,:) = chanmixMtx * lx;

watchoff;



function nut_compute_meglf_multsph_sensorloop
% multiple spheres... gotta do it sensor by sensor, since each sensor has its own lsc...
[lx1,lx2,ly1,ly2,lz1,lz2]=deal(zeros(size(nuts.meg.sensorCoord,1),size(voxels,1)));

% map lsc entries to channels...
for i=1:size(nuts.meg.sensorCoord,1)
    indexmap(i)=strmatch(nuts.meg.rawsensor_labels{i},nuts.meg.lsc_sensor_labels);
end

% one voxel at a time
Nvoxels = size(voxels,1);
Nchannels = size(nuts.meg.sensorCoord,1);

Nchannels = size(nuts.meg.sensorCoord,1);
for i=1:Nchannels
    %           waitbar(i/Nchannels,h);
    [lx1(i,:),ly1(i,:),lz1(i,:)] = nut_leadf_meg3comp(voxels,nuts.meg.sensorCoord(i,:,1),nuts.meg.sensorOrient(i,:),lsc(indexmap(i),:));
    [lx2(i,:),ly2(i,:),lz2(i,:)] = nut_leadf_meg3comp(voxels,nuts.meg.sensorCoord(i,:,2),nuts.meg.sensorOrient(i,:),lsc(indexmap(i),:));
end

lx1 = lx1-lx2;
clear lx2
ly1 = ly1-ly2;
clear ly2
lz1 = lz1-lz2;
clear lz2


function nut_compute_meglf_multsph_sensorvoxloop
% multiple spheres... sensor by sensor AND voxel by voxel (sloooooooow)
[lx1,lx2,ly1,ly2,lz1,lz2]=deal(zeros(size(nuts.meg.sensorCoord,1),size(voxels,1)));

% map lsc entries to channels...
for i=1:size(nuts.meg.sensorCoord,1)
    indexmap(i)=strmatch(nuts.meg.rawsensor_labels{i},nuts.meg.lsc_sensor_labels);
end

% one voxel at a time
Nvoxels = size(voxels,1);
Nchannels = size(nuts.meg.sensorCoord,1);

for i=1:Nchannels
    waitbar(i/Nchannels,h);
    for j=1:Nvoxels
        [lx1(i,j),ly1(i,j),lz1(i,j)] = nut_leadf_xyz_slow(voxels(j,:),nuts.meg.sensorCoord(i,:,1),nuts.meg.sensorOrient(i,:),lsc(indexmap(i),:));
        [lx2(i,j),ly2(i,j),lz2(i,j)] = nut_leadf_xyz_slow(voxels(j,:),nuts.meg.sensorCoord(i,:,2),nuts.meg.sensorOrient(i,:),lsc(indexmap(i),:));;
    end
    fprintf('.');
end

lx1 = lx1-lx2;
clear lx2
ly1 = ly1-ly2;
clear ly2
lz1 = lz1-lz2;
clear lz2



function nut_compute_meglf_1sph2comp_sensorloop
% single sphere MEG lead field (2-component)
% computed sensor by sensor
[lx1,lx2,ly1,ly2]=deal(zeros(size(nuts.meg.sensorCoord,1),size(voxels,1)));

% one sensor at a time
Nchannels = size(nuts.meg.sensorCoord,1);
innercoords = nuts.meg.sensorCoord(:,:,1);
outercoords = nuts.meg.sensorCoord(:,:,2);
for i=1:Nchannels
    %            waitbar(i/Nchannels,h);
    [lx1(i,:),ly1(i,:)] = nut_leadf_meg2comp(voxels,innercoords(i,:),nuts.meg.sensorOrient(i,:),lsc);
    [lx2(i,:),ly2(i,:)] = nut_leadf_meg2comp(voxels,outercoords(i,:),nuts.meg.sensorOrient(i,:),lsc);
    lx1(i,:) = lx1(i,:)-lx2(i,:);
    ly1(i,:) = ly1(i,:)-ly2(i,:);
end


function nut_compute_meglf_1sph2comp_voxloop
% single sphere MEG lead field (2-component)
% computed voxel by voxel (slooooooow)
[lx1,lx2,ly1,ly2]=deal(zeros(size(nuts.meg.sensorCoord,1),size(voxels,1)));

% one voxel at a time
Nvoxels = size(voxels,1);
innercoords = nuts.meg.sensorCoord(:,:,1);
outercoords = nuts.meg.sensorCoord(:,:,2);
for i=1:Nvoxels
    %            waitbar(i/Nvoxels,h);
    [lx1(:,i),ly1(:,i)] = nut_leadf_meg2comp(voxels(i,:),innercoords,nuts.meg.sensorOrient,lsc);
    [lx2(:,i),ly2(:,i)] = nut_leadf_meg2comp(voxels(i,:),outercoords,nuts.meg.sensorOrient,lsc);
    lx1(:,i) = lx1(:,i)-lx2(:,i);
    ly1(:,i) = ly1(:,i)-ly2(:,i);
end


function nut_compute_meglf_1sph2comp_voxsensorloop
% single sphere MEG lead field (2-component)
% computed sensor by sensor AND voxel by voxel (even slooooooower)
% but easier to understand
[lx1,lx2,ly1,ly2]=deal(zeros(size(nuts.meg.sensorCoord,1),size(voxels,1)));

% one voxel at a time
Nvoxels = size(voxels,1);
Nchannels = size(nuts.meg.sensorCoord,1);
innercoords = nuts.meg.sensorCoord(:,:,1);
outercoords = nuts.meg.sensorCoord(:,:,2);

for i=1:Nchannels
    %            waitbar(i/Nchannels,h);
    for j=1:Nvoxels
        [lx1(i,j),ly1(i,j)] = nut_leadf_slow(voxels(j,:),innercoords(i,:),nuts.meg.sensorOrient(i,:),lsc);
        [lx2(i,j),ly2(i,j)] = nut_leadf_slow(voxels(j,:),outercoords(i,:),nuts.meg.sensorOrient(i,:),lsc);
        lx1(i,j) = lx1(i,j)-lx2(i,j);
        ly1(i,j) = ly1(i,j)-ly2(i,j);
    end
end

function nut_compute_meglf_1sph2comp_voxsensorloop_somewhatfaster
% single sphere MEG lead field (2-component)
% computed sensor by sensor AND voxel by voxel (somewhat faster -- partially vectorized)
% but easier to understand

[lx1, ly1] = nut_leadf_less_slow(voxels,nuts.meg.sensorCoord(:,:,1),nuts.meg.sensorOrient,lsc);
% waitbar(0.5,h);
[lx2, ly2] = nut_leadf_less_slow(voxels,nuts.meg.sensorCoord(:,:,2),nuts.meg.sensorOrient,lsc);
lx1 = lx1-lx2;
ly1 = ly1-ly2;

function nut_compute_meglf_1sph3comp_voxsensorloop
% single sphere MEG lead field (3-component)
% computed sensor by sensor AND voxel by voxel (somewhat faster -- partially vectorized)
% but easier to understand
[lx1,lx2,ly1,ly2,lz1,lz2]=deal(zeros(size(nuts.meg.sensorCoord,1),size(voxels,1)));

% one voxel at a time
Nvoxels = size(voxels,1);
Nchannels = size(nuts.meg.sensorCoord,1);
innercoords = nuts.meg.sensorCoord(:,:,1);
outercoords = nuts.meg.sensorCoord(:,:,2);

for i=1:Nchannels
    %            waitbar(i/Nchannels,h);
    for j=1:Nvoxels
        [lx1(i,j),ly1(i,j),lz1(i,j)] = nut_leadf_xyz_slow(voxels(j,:),innercoords(i,:),nuts.meg.sensorOrient(i,:),lsc);
        [lx2(i,j),ly2(i,j),lz2(i,j)] = nut_leadf_xyz_slow(voxels(j,:),outercoords(i,:),nuts.meg.sensorOrient(i,:),lsc);
    end
    fprintf('.');
end
lx1 = lx1-lx2;
ly1 = ly1-ly2;
lz1 = lz1-lz2;
return


function nut_compute_meglf_1sph3comp_noloop
% fast vectorized version (but requires lots of RAM)

[lx1, ly1, lz1] = nut_leadf_meg3comp(voxels,nuts.meg.sensorCoord(:,:,1),nuts.meg.sensorOrient,lsc);
%       waitbar(0.5,h);

if(size(nuts.meg.sensorCoord,3)==2)
    [lx2, ly2, lz2] = nut_leadf_meg3comp(voxels,nuts.meg.sensorCoord(:,:,2),nuts.meg.sensorOrient,lsc);
else  % magnetometers
    [lx2, ly2, lz2]=deal(zeros(size(nuts.meg.sensorCoord,1),size(voxels,1)));
end

lx1 = lx1-lx2;
clear lx2
ly1 = ly1-ly2;
clear ly2
lz1 = lz1-lz2;
clear lz2

function nut_compute_meglf_1sph3comp_sensorloop
% sensor by sensor -- a little slower, but a lot less RAM

[lx1,lx2,ly1,ly2,lz1,lz2]=deal(zeros(size(nuts.meg.sensorCoord,1),size(voxels,1)));

% one sensor at a time
Nchannels = size(nuts.meg.sensorCoord,1);
for i=1:Nchannels
    %            waitbar(i/Nchannels,h);
    [lx1(i,:),ly1(i,:),lz1(i,:)] = nut_leadf_meg3comp(voxels,nuts.meg.sensorCoord(i,:,1),nuts.meg.sensorOrient(i,:),lsc);
    [lx2(i,:),ly2(i,:),lz2(i,:)] = nut_leadf_meg3comp(voxels,nuts.meg.sensorCoord(i,:,2),nuts.meg.sensorOrient(i,:),lsc);
end
lx1 = lx1-lx2;
clear lx2
ly1 = ly1-ly2;
clear ly2
lz1 = lz1-lz2;
clear lz2

function  nut_compute_meglf_1sph2comp_noloop
% fast vectorized version (but requires lots of RAM)

[lx1, ly1] = nut_leadf_meg2comp(voxels,nuts.meg.sensorCoord(:,:,1),nuts.meg.sensorOrient,lsc);
%        waitbar(0.5,h);
%     pack;
if(size(nuts.meg.sensorCoord,3)==2)
    [lx2, ly2, lz2] = nut_leadf_meg2comp(voxels,nuts.meg.sensorCoord(:,:,2),nuts.meg.sensorOrient,lsc);
else  % magnetometers
    [lx2, ly2, lz2]=deal(zeros(size(nuts.meg.sensorCoord,1),size(voxels,1)));
end


lx1 = lx1-lx2;
ly1 = ly1-ly2;
%%%% sensor by sensor...
%     [lx1,lx2,ly1,ly2]=deal(zeros(size(nuts.meg.sensorCoord,1),size(voxels,1)));
%     Nchannels = size(nuts.meg.sensorCoord,1);
%     for ii=1:Nchannels
%          waitbar(ii/Nchannels,h);
%         [lx1(ii,:),ly1(ii,:)] = nut_leadf_meg2comp(voxels,nuts.meg.sensorCoord(ii,:,1),nuts.meg.sensorOrient(ii,:),lsc);
%         [lx2(ii,:),ly2(ii,:)] = nut_leadf_meg2comp(voxels,nuts.meg.sensorCoord(ii,:,2),nuts.meg.sensorOrient(ii,:),lsc);
%     end
