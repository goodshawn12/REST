function voxels = nut_makeMNIvoi(voxelsize,labelselect,VOIvoxels)
% nut_makeMNIvoi(voxelsize,'wholehead') to select head VOI
% nut_makeMNIvoi(voxelsize,'wholebrain') to select only brainy brain (no
% eyes, skull, etc.)
% nut_makeMNIvoi(voxelsize,labelselect) to select specific MNI labels
% nut_makeMNIvoi(voxelsize,[],VOIvoxels) to add on additional VOI
% Creates a VOI for a given subject using MNI headshape
% Requires that both subject MRI and normalized MRI have been loaded.

global ndefaults

% we start with untrimmed MEG voxel grid at desired resolution,
% then convert to MNI coords and keep points that are in common with
% MNI template VOI

% create empty "row" if VOIvoxels doesn't exist, to obviate matlab
% acrobatics further down


switch(labelselect)
    case 'wholehead'
        %this is found in nutmeg/templates
        load('MNIvoxels');
    case 'wholebrain' % just the brain; no eyes, skull, etc.
%        load('ihb_DataBase.cdb','-mat'); % load MNI labels
%        [meshx,meshy,meshz]=ndgrid(MNIdb.minX:MNIdb.voxX:MNIdb.maxX,MNIdb.minY:MNIdb.voxY:MNIdb.maxY,MNIdb.minZ:MNIdb.voxZ:MNIdb.maxZ);
%         load('MNIdaemon');

        MNIdm=nut_loadMNI;
        MNIvoxels = MNIdm.coords(find(MNIdm.data),:);
        clear MNIdm
%         [meshx,meshy,meshz]=ndgrid(MNIdb.minX:MNIdb.maxX,MNIdb.minY:MNIdb.maxY,MNIdb.minZ:MNIdb.maxZ);
%         MNIvoxels = [meshx(:) meshy(:) meshz(:)];
%         MNIvoxels(find(MNIdb.data(:,1)==1),:)=[]; % remove "unidentified" voxels
%         clear MNIdb meshx meshy meshz
    otherwise % we have a specific MNI label
%         load('MNIdaemon');
%         for ii=1:5
%             labelidx=strmatch(lower(labelselect),lower(MNIdm.labels{ii}),'exact');
%             if(~isempty(labelidx))
%                 break
%             end
%         end

%         if(isempty(labelidx))
%             error('label not found');
%         end
%         
%         switch(ndefaults.mni2tal)
%             case 'icbm'
%                 coordidx=find(ismember(MNIdm.data_icbm,labelidx));
%             case 'brett'
%                 coordidx=find(ismember(MNIdm.data_brett,labelidx));
%         end
%         MNIvoxels=MNIdm.coords(coordidx,:);

%        MNIdm=nut_loadMNI;
        MNIvoxels = nut_MNIlabelselect(labelselect);

        clear MNIdm coordidx labelindex
end

if(0)
        load('ihb_DataBase.cdb','-mat'); % load MNI labels
        [meshx,meshy,meshz]=ndgrid(MNIdb.minX:MNIdb.voxX:MNIdb.maxX,MNIdb.minY:MNIdb.voxY:MNIdb.maxY,MNIdb.minZ:MNIdb.voxZ:MNIdb.maxZ);
        MNIdb.coords = single([meshx(:) meshy(:) meshz(:)]);
        clear meshx meshy meshz

        for ii=1:5
            labelindextmp=strmatch(lower(labelselect),lower(MNIdb.cNames{ii}),'exact');
            if(~isempty(labelindextmp))
                labelindex = [ii labelindextmp];
            end
        end
        
        if(~exist('labelindex','var'))
            error('label not found');
        end
        
        coordidx=find(MNIdb.data(:,labelindex(1)) == labelindex(2));
        MNIvoxels=MNIdb.coords(coordidx,:);
        clear MNIdb coordidx labelindex labelindextmp
end
    
    
if(~exist('voxelsize','var'))
    voxelsize = 5;
end

if(~exist('VOIvoxels','var'))
    %VOIvoxels = zeros(0,3);
    voxellist = MNIvoxels;
else
    voxellist = union(MNIvoxels,round(VOIvoxels),'rows');
    clear VOIvoxels
end

% start with standard bounds for MEG space
meggrid = nut_coordgrid(single([-100:voxelsize:100]),single([-80:voxelsize:80]),single([-50:voxelsize:140]));

% silly workaround for case where several crazy outside-the-head voxels somehow map
% to MNI of (0,0,0)
meggridmni=nut_mri2mni(nut_meg2mri(meggrid),[],0); % this calls global nuts inside
% search and destroy! take out NaNs while we're at it.
killvoxels=find(all(meggridmni==0 | isnan(meggridmni), 2));
meggrid(killvoxels,:) = [];
meggridmni(killvoxels,:) = [];
clear MNIvoxels killvoxels

meggridmni = round(meggridmni);

% this stuff would make it faster, but have to keep track of removed voxels
% or you really screw yourself.
% remove NaN voxels from MNI space and corresponding MEG space
if(0)
    killvoxels = any(isnan(meggridmni),2);
    meggrid(killvoxels,:)=[];
    meggridmni(killvoxels,:)=[];
    meggridmni=unique(meggridmni,'rows');
end

% MATLAB R2007b introduced the new try-catch syntax, so for previous
% versions we have to use lasterror instead.
%
% http://www.mathworks.com/help/techdoc/rn/bq08o1n-1.html#bq24_2g-1


try
    [selvoxels, select] = intersect(meggridmni,voxellist,'rows');
catch % out of memory
    if strcmp(lasterror.identifier,'MATLAB:nomem')
        nummni = size(voxellist,1);
        [dum, select1] = intersect(meggridmni,voxellist(1:ceil(nummni/2),:),'rows');
        [dum, select2] = intersect(meggridmni,voxellist(ceil(nummni/2)+1:end,:),'rows');
        select = union(select1,select2); 
        clear select1 select2 dum
    else
        rethrow(lasterror)
    end
end

voxels = meggrid(select,:);
