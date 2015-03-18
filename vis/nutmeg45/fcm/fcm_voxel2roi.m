function R = fcm_voxel2roi(voxels,coreg,method,roideffile,badroiidx)
% FCM_VOXEL2ROI  calculates the transformation matrices for switching
%                beween voxels and ROIs
% R = fcm_voxel2roi(voxels,coreg,¦method¦,¦roideffile¦,¦rois2exclude¦)
%                  (parameters in ¦¦ are optional)
%   
% R             (output) structure containing the transformation matrices
% voxels        voxel coordinates from session file
% coreg         coregistered MRI info from session file
% method        'one' ROI time series are equal to the time series of a single
%                     center voxel of each ROI (default).
%               'mean' ROI time series are the average of the time series
%                      of all voxels of a given ROI (this is usually
%                      problematic because voxel dipoles of a ROI can have flipped
%                      polarities making the mean meaningless. If you
%                      choose this function, make sure to fix the dipole
%                      orienations approriately).
% roideffile    optional matlab file containing ROI definitions (AAL definitions by default)
% rois2exclude  optional index numbers of ROIs to be excluded


global fuse

if nargin<3 || isempty(method)
    method='one';
end
if nargin<4 || isempty(roideffile)
    if isstruct(fuse) && isfield(fuse,'roidef')
        roideffile = fuse.roidef;
    else
        roideffile = [fileparts(which('fcm_gui')) filesep 'templates' filesep 'AAL_ROI.mat'];
    end
end    

usefastalgo = true; % try to set to false if no voxels found for many ROIs.

load(roideffile);

if nargin>4 && ~isempty(badroiidx)
    fprintf('\tExcluding ROI(s): '), fprintf('%s ',ROI.label{badroiidx}), fprintf('\n')
    bad=ismember(ROI.nr,badroiidx);
    ROI.nr(bad)=[];
    ROI.voxels(bad,:)=[];
    clear bad
end

roilist = unique(ROI.nr);
roilist = roilist(:).';  % make row vector
numroi=length(ROI.label);

SMV=nut_mri2mni(nut_meg2mri(voxels,coreg),coreg);

if usefastalgo   % This is much much much much faster, but rejects more voxels as not beloning to any ROI.
    [vb,tfm] = mm2voxels(ROI.voxels);
    RV = nut_vector2vol(ROI.nr,vb);
    SMVB = nut_coordtfm(SMV,inv(tfm));
    newnr = interp3(RV,SMVB(:,2),SMVB(:,1),SMVB(:,3),'nearest',NaN);
%     F = TriScatteredInterp(double(ROI.voxels),double(ROI.nr),'nearest');
%     newnr=F(double(SMV));
    good = isfinite(newnr);
    T = newnr(good);  % roi number of each voxel, voxels without ROI are removed
else
    voxelsize=10;
    [idx,d] = dsearchn(ROI.voxels,SMV);
    good=(d<voxelsize/2); clear d
    T=ROI.nr(idx(good));  % roi number of each voxel, voxels without ROI are removed
end
numgood = sum(good);

if strcmpi(method,'one')  % find center voxel of each ROI
    SMV = SMV(good,:);
    idx = nan(numroi,1);
    for k=1:numroi 
        curroi=find(T==roilist(k));
        if ~isempty(curroi)
            M=mean(SMV(curroi,:),1);
            idx(k) = curroi(dsearchn(SMV(curroi,:),M));
        end
    end
    T2 = nan(numgood,1);
    goodroi = isfinite(idx);
    T2(idx(goodroi)) = roilist(goodroi);
else
    T2=T;
end

R.voxel2roi_tfm = zeros(numgood,numroi,'single');
R.roi2voxel_tfm = zeros(numgood,numroi,'single');
R.goodvoxels = find(good);
R.goodroi = [];
R.numvoxinroi = zeros(1,numroi);
R.roilabel = ROI.label;
R.roidef = roideffile;

for k=roilist
    icurr= (T==k);
    icurr2=(T2==k);
    R.numvoxinroi(k)=sum(icurr);
    if R.numvoxinroi(k)<1
        warning('NUTMEG:fcm:ROIwithoutVoxel','ROI "%s" has no voxels.',ROI.label{k})
    else
        R.goodroi = [R.goodroi k];
        R.voxel2roi_tfm(icurr2,k) = 1/sum(icurr2);
        R.roi2voxel_tfm(icurr,k) = 1;
    end
end

if length(R.goodroi)~=82
    R.voxel2roi_tfm = R.voxel2roi_tfm(:,R.goodroi);
    R.roi2voxel_tfm = R.roi2voxel_tfm(:,R.goodroi);
end

if ~isempty(fuse)
    fuse.roi=length(R.goodroi);
else
    clear global fuse
end

%===========
function [voxelsblob,blob2mri_tfm] = mm2voxels(voxelsMNI)

voxelsize=[2 2 2];   % true for AAL

blob2mri_tfm = [ voxelsize(1)        0                 0             min(voxelsMNI(:,1))-voxelsize(1)
                            0          voxelsize(2)        0             min(voxelsMNI(:,2))-voxelsize(2)
                            0                 0         voxelsize(3)     min(voxelsMNI(:,3))-voxelsize(3)
                            0                 0                 0                                1 ];
voxelsblob = nut_coordtfm(voxelsMNI,inv(blob2mri_tfm));