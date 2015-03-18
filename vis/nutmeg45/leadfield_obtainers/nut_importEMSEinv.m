function [L,voxels] = importEMSEinv(invfile)

% read INV file
fid=fopen(invfile,'r','a');
magicnum=fscanf(fid,'%s',1);
if(~strcmp(magicnum,'454d5345'))
    error('This doesn''t look like an EMSE forward model.');
end
majorrev = fscanf(fid,'%d',1);
minorrev = fscanf(fid,'%d',1);
emsetype = fscanf(fid,'%s',1);
emsemode = fscanf(fid,'%s',1);
emsestate = fscanf(fid,'%s',1);
emseoption = fscanf(fid,'%s',1);

tanspacedim = fscanf(fid,'%d',1);
checksum = fscanf(fid,'%d',1);
numrows = fscanf(fid,'%d',1);
numchans = fscanf(fid,'%d',1);
fread(fid,1); % read line feed

L2 = fread(fid,numrows*numchans,'double');

for ii=1:6
    fgetl(fid);
end
numvoxels = fscanf(fid,'%d',3);
startcoords = fscanf(fid,'%f',3);
voxelsize = fscanf(fid,'%f',3);
fgetl(fid);fgetl(fid);

coordflag = zeros(numvoxels(1),numvoxels(2),numvoxels(3));
coordflag2 = zeros(prod(numvoxels),1);

iter=0;

for kk=1:numvoxels(3)
    fscanf(fid,'%d',3);
    for jj=1:numvoxels(2)
        fscanf(fid,'%d',2);
        for ii=1:numvoxels(1)
            iter=iter+1;
            temp=fscanf(fid,'%f',7);
            if(temp(4)==5)
                coordflag(ii,jj,kk)=1;
                coordflag2(iter)=1;
            end
        end
    end
end

fclose(fid);
