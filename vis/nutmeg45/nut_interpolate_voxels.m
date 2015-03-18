function beam=nut_interpolate_voxels(beam)
% NUT_INTERPOLATE_VOXELS  estimates activation values of missing voxels by
% using a linear interpolation of all existing neighbouring voxels.
%
%   beam = nut_interpolate_voxels(beam)
%
% beam   structure containing SPATIALLY NORMALIZED (!!!) beamforming data

load MNIvoxels
v=find(~rem(MNIvoxels(:,1),beam.voxelsize(1)) & ~rem(MNIvoxels(:,2),beam.voxelsize(2)) & ~rem(MNIvoxels(:,3),beam.voxelsize(3)));
MNIvoxels=MNIvoxels(v,:); 
clear v

%MNItfm = [beam.voxelsize(1) 0 0 min(MNIvoxels(:,1)); 0 beam.voxelsize(2) 0 min(MNIvoxels(:,2)); 0 0 beam.voxelsize(3) min(MNIvoxels(:,3)); 0 0 0 1];
%voxelgrid = nut_coordtfm(beam.voxels,inv(MNItfm));

missing=setdiff(MNIvoxels,beam.voxels,'rows');
num = size(missing,1);
nuv = size(beam.voxels,1);

sinterp = nan(num,length(beam.timepts),size(beam.bands,1));
for k=1:num
    lim = beam.voxelsize*1.1;
    vot = abs(beam.voxels-repmat(missing(k,:),[nuv 1]));
    v = find(vot(:,1)<lim(1) & vot(:,2)<lim(2) & vot(:,3)<lim(3));      % find neighbouring voxels
    if ~isempty(v)
        dist = (1./sqrt(sum(vot(v,:).^2,2)))';                           % calculate Euclidian distances to missing voxel
        totdist = sum(dist);                                            % shorter distances weigh more
        for t=1:length(beam.timepts)
            for f=1:size(beam.bands,1)
                sinterp(k,t,f) = (dist./totdist) * beam.s{1}(v,t,f);     % linear interpolation
            end
        end
    end
end
v = find(~isnan(sinterp(:,1,1)));
beam.s{1}=cat(1,beam.s{1},sinterp(v,:,:));
beam.voxels=cat(1,beam.voxels,missing(v,:));
