function [NFV,hdr]=readdfs(filename)
% READDFS Writes a Duff Surface file (dfs).
% [NFV,hdr]=READDFS(FILENAME) reads the file specified by FILENAME string.
%
% DFS has the following structure:
% NFV.faces : the face data,
% NFV.vertices : the face data,
% Examples
%
% surface=readdfs('brain.dfs');
%
% Author : David Shattuck (shattuck@loni.ucla.edu)

% This is the header information from C++
% char magic[8]; // Magic number (DUFFSURF on little-endian machines or byte-swapped equivalent on other architectures)
% char version[4]; // A number in the format 1.1.1.1
% int32 hdrsize; // Size of complete header (i.e., offset of first data element)
% int32 mdoffset; // Start of metadata.
% int32 pdoffset; // Start of patient data header.
% int32 nTriangles; // Number of triangles
% int32 nVertices; // Number of vertices
% int32 nStrips; // Number of triangle strips
% int32 stripSize; // size of strip data
% int32 normals; // 4 Int32 <normals> Start of vertex normal data (0 if not in file)
% int32 uvStart; // Start of surface parameterization data (0 if not in file)
% int32 vcoffset; // vertex color
% uint8 precision; // Vertex Precision -- usually float32 or float64
% uint8 pad[3]; // padding
% float64 orientation[4][4]; //4x4 matrix, affine transformation to world coordinates

fid=fopen(filename,'rb','ieee-le');
if (fid<0) error('unable to open file'); end;
%hdr.magic = ['D' 'U' 'F' 'F' 'S' 'U' 'R' 'F']';
hdr.magic=char(fread(fid,8,'char'));
hdr.version=fread(fid,4,'char');
hdr.hdrsize=fread(fid,1,'int32');
hdr.mdoffset=fread(fid,1,'int32');
hdr.pdoffset=fread(fid,1,'int32');
hdr.nTriangles=fread(fid,1,'int32');
hdr.nVertices=fread(fid,1,'int32');
hdr.nStrips=fread(fid,1,'int32');
hdr.stripSize=fread(fid,1,'int32');
hdr.normals=fread(fid,1,'int32');
hdr.uvStart=fread(fid,1,'int32');
hdr.vcoffset=fread(fid,1,'int32');
hdr.precision=fread(fid,1,'int32');
hdr.orientation=fread(fid,[4 4],'float64');

fseek(fid,hdr.hdrsize,-1);
NFV.faces = fread(fid,[3 hdr.nTriangles],'int32')+1;
NFV.vertices=fread(fid,[3 hdr.nVertices],'float32');
NFV.faces=NFV.faces';
NFV.vertices=NFV.vertices';
if (hdr.normals>0)
display('reading normals.');
fseek(fid,hdr.normals,-1);
NFV.normals = fread(fid,[3 hdr.nVertices],'float32')';
end;
if (hdr.vcoffset>0)
display('reading color.');
fseek(fid,hdr.vcoffset,-1);
NFV.vcolor = fread(fid,[3 hdr.nVertices],'float32')';
end;
if (hdr.uvStart>0)
display('reading uv.');
fseek(fid,hdr.uvStart,-1);
uv = fread(fid,[2 hdr.nVertices],'float32');
NFV.u = uv(1,:);
NFV.v = uv(2,:);
end;
fclose(fid);