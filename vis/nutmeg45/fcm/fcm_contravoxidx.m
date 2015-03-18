function nuts=fcm_contravoxidx(nuts)

selvoxcoord=nuts.voxels(nuts.selvox.ipsi,:);

% find corresponding indices of selected voxels
dointerp = any(any( round(nuts.voxels(1:5,:)) - nuts.voxels(1:5,:) > 0.01 ));
if dointerp
    selvoxcoord(:,1) = -selvoxcoord(:,1);      % flip left-right, we assume that non-integer voxels are in MRI space !
    numsel=length(nuts.selvox.ipsi);
    coll=zeros(round(numsel*1.5),1);   % for faster performance, reserve some workspace
    count=0;
    if isscalar(nuts.selvox.ipsi)
        nuts.selvox.contra = dsearchn(nuts.voxels,selvoxcoord);
    else % ROI
        for k=1:numsel
            D=nut_rownorm(nut_coord_diff(nuts.voxels,selvoxcoord(k,:)));
            f=find(D<nuts.voxelsize(1));
            numf=length(f);
            coll(count+1:count+numf)=f;
            count=count+numf;
        end
        coll = coll(1:count);           % remove zeros
        nuts.selvox.contra = unique(coll);
    end
else
    selvoxcoord(:,2) = -selvoxcoord(:,2);      % flip left-right, we assume that integer voxels are in MEG space !    
    nuts.selvox.contra = find(ismember(nuts.voxels,selvoxcoord,'rows'));
end