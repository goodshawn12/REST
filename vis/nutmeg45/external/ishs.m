function test = ishs(filename)

if nargin < 1
    error('Too few input argumets');
end

if ~exist(filename, 'file')
    test = false;
    return
end

fid = fopen(filename, 'r', 'b');

if fid == -1
    error('Cannot open file %s', filename);
end

hdr = read_hs_header(fid);

fclose(fid);

if isempty(hdr.version) || ...        
    isempty(hdr.timestamp) || ...
     isempty(hdr.checksum) || ...
      isempty(hdr.npoints)
  
      test = false;
      return
end

%looks as hs_file
test = true;