function index = read_hs_index(fid)

index.lpa = fread(fid, 3, 'double');
index.rpa = fread(fid, 3, 'double');
index.nasion = fread(fid, 3, 'double');
index.cz = fread(fid, 3, 'double');
index.inion = fread(fid, 3, 'double');