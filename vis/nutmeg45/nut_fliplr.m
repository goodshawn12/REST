function newstru=nut_fliplr(stru,subj2flip)
% NUT_FLIPLR performs left-right flip of voxel coordinates. Voxel
% coordinates must be integer!
%
% Usages:
%  beam = nut_fliplr(beam);
%  pop  = nut_fliplr(pop,subj2flip);
 

if any(stru.voxels(1:5,:)-round(stru.voxels(1:5,:))>1e-2)
    error('Voxel coordinates are not close to integer.')
end

left = find(stru.voxels(:,1)<0);
leftvox = stru.voxels(left,:);
lookfor = [-leftvox(:,1) leftvox(:,2:3)];

mid = find(stru.voxels(:,1)==0);
midvox = stru.voxels(mid,:);

[good,flipidx] = ismember(lookfor,stru.voxels,'rows');
bad = ~ismember(stru.voxels,[leftvox;midvox;lookfor],'rows');

if( length(good)/length(left) < 0.8 )
    error('Many voxel coordinates are not symmetrical.')
end

newstru = stru;

if nargin<2    
    if iscell(stru.s)    % beam structure
        for jj=1:length(stru.s)
            newstru.s{jj}(flipidx(good),:,:) = stru.s{jj}(left(good),:,:);
            newstru.s{jj}(left(good),:,:)  = stru.s{jj}(flipidx(good),:,:);
            newstru.s{jj}(bad,:,:) = [];
        end
    else        % ?? just in case
        newstru.s(flipidx(good),:,:) = stru.s(left(good),:,:);
        newstru.s(left(good),:,:) = stru.s(flipidx(good),:,:);
        newstru.s(bad,:,:) = [];
    end
else    % pop structure
    subj2flip=subj2flip(:)'      % need to make sure numbers are in row because of for...end loop
    if max(subj2flip)>size(stru.s,1), error('Indices in subj2flip do not match with pop structure.'), end
    for k=subj2flip
        newstru.s(k,flipidx(good),:,:) = stru.s(k,left(good),:,:);
        newstru.s(k,left(good),:,:) = stru.s(k,flipidx(good),:,:);
    end
    newstru.s(:,bad,:,:) = [];
end

newstru.voxels(bad,:) = [];

    
    

