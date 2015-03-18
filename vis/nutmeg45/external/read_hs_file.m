function hs = read_hs_file(filename)

%hs = read_hs_file(filename)
%
%Reads MSI hs_file
%
%Converts C-style strings into MATLAB string by removing zeros
%Otherwise there is no data conversion

%Revision 1.0  09/11/06  eugene.kronberg@uchsc.edu

if nargin < 1
    error('Too few input argumets');
end

if ~ishs(filename)
    error('File %s is not head shape file', filename);
end

fid = fopen(filename, 'r', 'b');

if fid == -1
    error('Cannot open file %s', filename);
end

hs.hdr = read_hs_header(fid);
hs.index = read_hs_index(fid);
hs.point = fread(fid, [3 double(hs.hdr.npoints)], 'double');

fclose(fid);