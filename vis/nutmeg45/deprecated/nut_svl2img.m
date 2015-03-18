function nut_svl2img(SVLfile)
% nut_svl2img(SVLfile)
% converts CTF SAM output (*.svl) into analyze/spm image (*.hdr + *.img pair)
% NOTE: currently leaves data in MEG head coordinates, and does not
%       provide transformation to subject's MRI
% TODO: option to normalize/transform to MNI coordinates

[SAMimage,coords,SAMparams] = nut_read_svl(SVLfile);

beam.voxels = coords*1000; % convert from m to mm
beam.voxelsize = repmat(SAMparams.stepsize*1000,[1 3]);

global defaults
defaults.analyze.flip=0; % needed to placate SPM2's trigger-happy flipping

[outputdir,outputname]=fileparts(SVLfile);
SN.fname=fullfile(outputdir,[outputname '.img']);

SN.mat=diag([beam.voxelsize 1]);
SN.mat(1:3,4)=min(beam.voxels) - beam.voxelsize;

dimxyz=round(range(beam.voxels)./beam.voxelsize + 1);
SN.dim=[dimxyz 64];  % 16 means float, 64 means double precision
SN.pinfo=[1;0;0];

% rearrange order of voxels for MRI volume convention
[beam.voxels,sortidx]=sortrows(beam.voxels,[3 2 1]);
SAMimage = SAMimage(sortidx);

SAMimageout = reshape(10*log10(weirdFtostandardF(SAMimage)),dimxyz);
SN = spm_write_vol(SN,SAMimageout);



function F=weirdFtostandardF(Fweird)
% converts CTF's weird F value to a standard F ratio
F = zeros(size(Fweird));
F(Fweird>=0) = Fweird(Fweird>=0)+1;
F(Fweird<0) = 1./(1-Fweird(Fweird<0));
