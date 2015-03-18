function bvtfm = nut_read_bvtfm(bvfile)
% reads BrainVisa transformation (BrainVisa mm -> Nifti mm) from
% associated minf header file

minffid = fopen([bvfile '.minf']);
hdr=fgetl(minffid);
tfm_idx = strfind(hdr,'''transformations'':') + 21;
bvtfm=sscanf(hdr(tfm_idx:end),'%f,',[4 4])';
