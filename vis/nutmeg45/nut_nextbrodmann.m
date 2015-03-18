function nut_nextbrodmann(posmni,num)
% nut_nextbrodmann(MNIcoord)
%    displays the next MNI Brodmann label that is not unknown.
% nut_nextbrodmann(MNIcoord,number)
%    displays the next NUMBER Brodmann voxels that are not unknown.

if nargin<2, num=1; end

load('ihb_DataBase.cdb','-mat'); % load MNI labels
[meshx,meshy,meshz]=ndgrid(MNIdb.minX:MNIdb.voxX:MNIdb.maxX,MNIdb.minY:MNIdb.voxY:MNIdb.maxY,MNIdb.minZ:MNIdb.voxZ:MNIdb.maxZ);
MNIdb.coords = [meshx(:) meshy(:) meshz(:)];

distances = nut_rownorm(nut_coord_diff(MNIdb.coords,posmni));
f = find(distances <= 30);
ind2 = MNIdb.data(f,5);
good = find(ind2>1);
if num==1
    [d,idx]=min(distances(f(good)));
else
    [d,idx]=sort(distances(f(good)));
end

for k=1:num
    ind=MNIdb.data(f(good(idx(k))),:);

    for ii=1:5
        labels(ii,1) = MNIdb.cNames{ii}(ind(ii));
    end

    fprintf('\nFound at distance %1.1f mm\n',d(k))
    labels
end

