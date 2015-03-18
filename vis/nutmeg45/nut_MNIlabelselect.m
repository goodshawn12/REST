function [MNIvoxels,selectidx] = nut_MNIlabelselect(labelselect,hemisphere,MNIvoxels)
% [MNIvoxels,selectidx] = nut_MNIlabelselect(labelselect,hemisphere,MNIvoxels)
% given a list of MNI coordinates, returns indices of that list
% corresponding to specified Talaraich label.
% 'hemisphere' is optional; 'right' or 'left' can be specified
% --
% Note that you may need to convert voxels from MEG space to MNI space,
% like so:
% MNIvoxels = nut_mri2mni(nut_meg2mri(beam.voxels,beam.coreg,0));

% TalDB = nut_loadMNI;
TalDB = load('talairachDB');

% convert to MNI
TalDB.coords = round(nut_tal2mni(TalDB.coords));

if(exist('MNIvoxels','var'))
    MNIvoxels = round(MNIvoxels);

    [ispresent,map2db] = ismember(MNIvoxels,TalDB.coords,'rows');
    map2db(~ispresent) = 1; % map voxels outside brain to first voxel of TalDB.coords (which is also outside brain)
else  % if MNIvoxels not supplied, we use the database list as-is
    MNIvoxels = TalDB.coords;
    map2db = 1:length(TalDB.data);
end
    

% find MNI "scores" containing selected label
[crap,scoreidx]=ind2sub(size(TalDB.labels),strmatch(lower(labelselect),lower([TalDB.labels{:}]),'exact'));

selectidx = find(ismember(TalDB.data(map2db),scoreidx));
MNIvoxels = sort(MNIvoxels(selectidx,:));

if(exist('hemisphere','var'))
    switch(hemisphere)
        case 'left'
            MNIvoxels(find(MNIvoxels(:,1)>0),:) = [];
        case 'right'
            MNIvoxels(find(MNIvoxels(:,1)<0),:) = [];
    end
end