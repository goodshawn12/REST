function rawz=nut_zscore(raw)
% raw is data to be transformed to z-distribtuion (mean = 0, std = 1)
% raw should be size: time x channels_or_voxels x trials
% zvalues are computed per channel/voxel independently, using time points from all trials

[nt,ns,ni]=size(raw);
rawz=zeros([nt ns ni]);
for jj=1:ns
    locmean=mean(reshape(squeeze(raw(:,jj,:)),nt*ni,1));
    locstd=std(reshape(squeeze(raw(:,jj,:)),nt*ni,1));
    rawz(:,jj,:)=(raw(:,jj,:)-locmean)/locstd;
end