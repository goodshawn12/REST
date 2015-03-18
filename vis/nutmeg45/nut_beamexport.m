function nut_beamexport(outname)
% nut_beamexport('outputname.img')
% exports displayed activation blob as an SPM/Analyze format image
% Set threshold to 0 if entire activation map is desired; also, uncheck
% "normalize intensities" box if raw values are desired.
% ORIENTATION OF RESULTING IMAGE NOT THOROUGHLY TESTED!
% In particular, exported SAM images may need to be spatially normalized and/or
% resliced to match coordinates of corresponding MRI.
% Also note that SPM may freak out if, e.g., "s_beam_abc.mat" is exported to
% "s_beam_abc.img" -- be sure to name the image file something different so
% SPM doesn't try to mine transformation matrices from the mat file.

% fuck spm -- if true, removes NaN values for better compatibility with external image viewers
fuckspm = false;

global beam st
SN.fname=outname;
SN.mat=st.vols{1}.blobs{1}.mat;
[dimx,dimy,dimz]=size(st.vols{1}.blobs{1}.vol);
SN.dim=[dimx dimy dimz];  % 16 means float, 64 means double precision
% SN.dim=[dimx dimy dimz 64];  % 16 means float, 64 means double precision
SN.pinfo=[1;0;0];

% following two lines dealing with beam.* probably no longer necessary...
beam.voxelsblob=nut_coordtfm(nut_meg2mri(beam.voxels),inv(st.vols{1}.blobs{1}.mat));
beam.blobtransform = [ beam.voxelsize(1)                0                  0 beam.voxelsize(1)*floor(st.bb(1,1)/beam.voxelsize(1))
                                       0 beam.voxelsize(2)                 0 beam.voxelsize(2)*floor(st.bb(1,2)/beam.voxelsize(2))
                                       0                 0 beam.voxelsize(3) beam.voxelsize(3)*floor(st.bb(1,3)/beam.voxelsize(3))
                                       0                 0                 0                                                     1 ];

data = st.vols{1}.blobs{1}.vol;
if(fuckspm)
    data(isnan(data))=0;   % for compatibility with external programs
end
SN = spm_write_vol(SN,data);
