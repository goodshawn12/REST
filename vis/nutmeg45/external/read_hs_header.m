function hdr = read_hs_header(fid)

hdr.version = fread(fid, 1, 'uint32=>uint32');
hdr.timestamp = fread(fid, 1, 'int32=>int32');
hdr.checksum = fread(fid, 1, 'int32=>int32');
hdr.npoints = fread(fid, 1, 'int32=>int32');