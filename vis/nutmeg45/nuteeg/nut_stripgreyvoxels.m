function nut_stripgreyvoxels(greyfile)

global nuts st

if nargin<1
	[file,path] = uigetfile('*.img;*.nii','Open Grey Matter MRI');
	if isequal(file,0); return; end;
	greyfile=fullfile(path,file);
end
V = spm_vol(greyfile);
pixdim = diag(V.mat);
pixdim = abs(pixdim(1:3)');

voxelsMRI = double(nut_meg2mri(nuts.voxels));

blob2mri_tfm = [ pixdim(1)                0                  0 pixdim(1)*floor(st.bb(1,1)/pixdim(1))
                    0               pixdim(2)                0 pixdim(2)*floor(st.bb(1,2)/pixdim(2))
                    0                 0                 pixdim(3) pixdim(3)*floor(st.bb(1,3)/pixdim(3))
                    0                 0                     0                                                     1 ];
voxelsblob = nut_coordtfm(voxelsMRI,inv(blob2mri_tfm));
clear voxelsMRI
voxelsblob_round=round(voxelsblob);  
clear voxelsblob

Y = spm_read_vols(V);

%%% SOLUTION 1: forget it, this takes for ever
% [x y z] = ind2sub(size(Y),find(Y));
% greyvoxels = [x y z]; 
% clear x y z Y
% numvox = size(voxelsblob,1);
% numgrey = size(greyvoxels,1);
% dist = zeros(numvox,1);
% for k = 1:numvox
%     vbr = repmat(voxelsblob(k,:),numgrey,1);
%     dist(k) = min(sum((greyvoxels-vbr).^2,2));
% end
% good = ( dist < min(nuts.voxelsize)/2 );

%%% SOLUTION 2: Works, but has some holes in grey matter voxels.
% numvox = size(voxelsblob_round,1);
% good = true(numvox,1);
% for k=1:numvox
%     good(k) = ( Y(voxelsblob_round(k,1),voxelsblob_round(k,2),voxelsblob_round(k,3)) > 0 );
% end

%%% SOLUTION 3: 
h = waitbar(0,'Extracting grey matter voxels...');
numvox = size(voxelsblob_round,1);
vr = [1 1 1];
good = true(numvox,1);
for k=1:numvox
    x = [voxelsblob_round(k,1)-vr(1):voxelsblob_round(k,1)+vr(1)];
    y = [voxelsblob_round(k,2)-vr(2):voxelsblob_round(k,2)+vr(2)];
    z = [voxelsblob_round(k,3)-vr(3):voxelsblob_round(k,3)+vr(3)];
    %lims = [voxelsblob_round(k,1)-vr(1) voxelsblob_round(k,1)+vr(1) voxelsblob_round(k,2)-vr(2) voxelsblob_round(k,2)+vr(2) voxelsblob_round(k,3)-vr(3) voxelsblob_round(k,3)+vr(3)];
    %NV = subvolume(Y,lims);  % Subvolume does the same thing but rotates  the MRI!!!!!
    NV = Y(x,y,z);
    good(k) = any(NV(:));
    waitbar(k/numvox,h);
end
delete(h);

nuts.voxels = nuts.voxels(good,:);
if isfield(nuts,'Lp')
    nuts.Lp =  nuts.Lp(:,:,good);
end


    
    