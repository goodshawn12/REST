function vox_trans = nut_import_smacsolutionpoints(fSPI,fHDR)
% vox_trans = nut_import_smacsolutionpoints(solpointfile,mrifile)

% Read solutionpoint file
voxels = dlmread(fSPI, ' ');

% Get MRI info and SPM transformation matrix
V = spm_vol(fHDR);
mridim = V.dim;
mrimat = V.mat;
if mrimat(1,1)<0
    mrimat(1,:)=-mrimat(1,:);  % make sure nothing is flipped in SPM8.
end
%%% old solution which failed with some MRIs
% hdr=spm_read_hdr(fHDR);
% mridim = hdr.dime.dim(2:4);
% mripixdim = hdr.dime.pixdim(2:4);
% if any(hdr.hist.origin(1:3))
%     mriorigin_pix = hdr.hist.origin(1:3);
% else
%     mriorigin_pix = (mridim+1)/2;
% end
% mriorigin_mm = mriorigin_pix .* mripixdim;
% mrimat = [mripixdim(1) 0 0 -mriorigin_mm(1) ; 0 mripixdim(2) 0 -mriorigin_mm(2) ; 0 0 mripixdim(3) -mriorigin_mm(3) ; 0 0 0 1];

% Rotate and transform coordinates
vox_trans=zeros(size(voxels));
vox_trans(:,1) = voxels(:,3);
vox_trans(:,2) = mridim(2) - voxels(:,1);
vox_trans(:,3) = mridim(3) - voxels(:,2);

% Bring voxels to mm data according to SPM convention
% This corresponds to: vox_trans = nut_voxels2mm(vox_trans);
vox_trans = [vox_trans ones(size(vox_trans,1),1)] * mrimat';
vox_trans(:,4)=[];   % destroy last column of ones


%======================
% This is an old solution that worked but not always.

% vox_trans = nut_import_smacsolutionpoints(solpointfile.spi,mrifile.hdr)

% if nargin<3, otype='SPM'; end
% 
% voxels = dlmread(fSPI, ' ');
% origin_mm = nut_read_mriorigin(fHDR,otype);
% 
% voxels = voxels+1.5;    % Revert correction made in SMAC

% Rotate and transform functional data
% vox_trans=zeros(size(voxels));
% vox_trans(:,1)=voxels(:,3)-origin_mm(3);
% vox_trans(:,2)=origin_mm(1)-voxels(:,1);
% vox_trans(:,3)=origin_mm(2)-voxels(:,2);
