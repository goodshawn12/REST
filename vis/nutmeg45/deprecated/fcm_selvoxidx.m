function nuts=fcm_selvoxidx(VOIvoxels)
% FCM_SELVOXIDX  allows selecting seed voxels of interest and returns their
%       indices as well as the indices of homologous contralateral voxels.
%       The output is stored in the fields nuts.selvox.ipsi and
%       nuts.selvox.contra.
%
% Usage:
%
%   ¦nuts¦ = fcm_selvoxidx(¦VOIvoxels¦)
%     (parameters in ¦¦ are optional)
%
% nuts         NUTMEG session structure.
% VOIvoxels    structural coordinates of seed region of interest (optional).

global nuts

if nargin<1 
    if ~isfield(nuts,'VOIvoxels')
        VOIvoxels=fcm_selectvoxels(nuts);
    else
        VOIvoxels=nuts.VOIvoxels;
    end
end

% Try to make integer if necessary
% if ~(isempty(find(rem(VOIvoxels,1))))   
%     realshift = rem(VOIvoxels(1,:),1);
%     VOIvoxels = VOIvoxels - repmat(realshift,[size(VOIvoxels,1) 1]);
%     if ~(isempty(find(rem(VOIvoxels,1))))
%         error('Does not work.')
%     end
%     %voxels = voxels - repmat(realshift,[size(voxels,1) 1]);
% end

% Resample VOIvoxels to the same voxelsize as voxels
VVresampled =  zeros(size(VOIvoxels));
for ii = 1:3
    VVresampled(:,ii) = round(VOIvoxels(:,ii) ./ nuts.voxelsize(ii)) .* nuts.voxelsize(ii);
end
VVresampled = unique(VVresampled,'rows');

% Add integer correction back to resampled VOIvoxels
% if exist('realshift','var')
%     VVresampled = VVresampled + repmat(realshift,[size(VVresampled,1) 1]);
% end

% find corresponding indices of selected voxels
dointerp = any(any( round(nuts.voxels(1:5,:)) - nuts.voxels(1:5,:) > 0.01 ));
if dointerp
    numVVresampled=size(VVresampled,1);
    coll=zeros(round(numVVresampled*1.5),1);   % for faster performance, reserve some workspace
    count=0;
    for k=1:numVVresampled
        D=nut_rownorm(nut_coord_diff(nuts.voxels,VVresampled(k,:)));
        f=find(D<nuts.voxelsize(1));
        numf=length(f);
        coll(count+1:count+numf)=f;
        count=count+numf;
    end
    coll = coll(1:count);           % remove zeros
    svox = unique(coll);
else
    svox = find(ismember(nuts.voxels,VVresampled,'rows'));
end

% find contralateral voxels

% Make sure that brain is symmetric around middle
%     if ~isaligned
%         voxels=alignvoxels(nuts.voxels);
%         VVresampled=alignvoxels(VVresampled);
%     end
if dointerp
    VVresampled(:,1) = -VVresampled(:,1);      % flip left-right, we assume that non-integer voxels are in MRI space !
    coll=zeros(round(numVVresampled*1.5),1);   % for faster performance, reserve some workspace
    count=0;
    for k=1:numVVresampled
        D=nut_rownorm(nut_coord_diff(nuts.voxels,VVresampled(k,:)));
        f=find(D<nuts.voxelsize(1));
        numf=length(f);
        coll(count+1:count+numf)=f;
        count=count+numf;
    end
    coll = coll(1:count);           % remove zeros
    cvox = unique(coll);
else
    VVresampled(:,2) = -VVresampled(:,2);      % flip left-right, we assume that integer voxels are in MEG space !
    cvox = find(ismember(nuts.voxels,VVresampled,'rows'));
end

nuts.selvox.ipsi = svox;
nuts.selvox.contra = cvox;
    



