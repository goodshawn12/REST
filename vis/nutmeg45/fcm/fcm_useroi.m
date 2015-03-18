function fcm_useroi(roifile)
% fcm_useroi(roifile)

global fuse
load(roifile);
fuse.roi=length(R.goodroi);
