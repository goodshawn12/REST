function [L,voxels] = importEMSElf(fwdfile,smpfile)

% read FWD file
fid=fopen(fwdfile,'r','a');
magicnum=fscanf(fid,'%s',1);
if(~strcmp(magicnum,'454d5345'))
    error('This doesn''t look like an EMSE forward model.');
end
majorrev = fscanf(fid,'%d',1);
minorrev = fscanf(fid,'%d',1);
emsetype = fscanf(fid,'%s',1);
emsemode = fscanf(fid,'%s',1);
numrows = fscanf(fid,'%d',1);
numchans = fscanf(fid,'%d',1);
tanspacedim = fscanf(fid,'%d',1);
fread(fid,1,'char'); % read line feed


L = fread(fid,inf,'double');
fclose(fid);

% read SMP file
fid=fopen(smpfile,'r','a');
% option=fscanf(fid,'%d',1);
% state=fscanf(fid,'%d',1);
% condnum=fscanf(fid,'%d',1);
% headradius=fscanf(fid,'%f',2)
for ii=1:7
    fgetl(fid);
end
numvoxels = fscanf(fid,'%d',3);
startcoords = fscanf(fid,'%f',3);
voxelsize = fscanf(fid,'%f',3);
fgetl(fid);fgetl(fid);

coordflag = zeros(numvoxels(1),numvoxels(2),numvoxels(3));

for kk=1:numvoxels(3)
    fscanf(fid,'%d',3);
    for jj=1:numvoxels(2)
        fscanf(fid,'%d',2);
        for ii=1:numvoxels(1)
            temp=fscanf(fid,'%f',7);
            if(temp(4)==5)
                coordflag(ii,jj,kk)=1;
            end
        end
    end
end

fclose(fid);

keep=find(coordflag(:));
L=reshape(L,3,length(keep),304);
% L=reshape(L2,3,length(keep),271);

L=permute(L,[3 1 2]);

for ii=1:3
    endcoords(ii) = startcoords(ii)+(numvoxels(ii)-1)*voxelsize(ii);
end

voxelsizemm = voxelsize*1000;
startcoordsmm = startcoords*1000;
endcoordsmm=endcoords*1000;
coords=nut_coordgrid(startcoordsmm(1):voxelsizemm(1):endcoordsmm(1),startcoordsmm(2):voxelsizemm(2):endcoordsmm(2),startcoordsmm(3):voxelsizemm(3):endcoordsmm(3));
voxels=coords(keep,:);