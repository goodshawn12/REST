function W=nut_import_smacis(filename)
% W = nut_import_cartoolis(filename)
% reads SMAC inverse solution files (*.is) containing weight matrices.

fid=fopen(filename,'r');
if fid<0, error('Invalid filename.'), end

% Read Header
prectxt=fread(fid,4,'*char')';
switch prectxt
    case 'IS01'
        prec = 'float32';
    case 'IS02'
        prec = 'double';
    otherwise
        error('Invalid IS file.')
end

numelec=fread(fid,1,'int32');
numsolpts=fread(fid,1,'int32');
isinversescalar= fread(fid,1,'char');
if isinversescalar
    numvec = 1;
else
    numvec = 3;
end

tmp = fread(fid,Inf,prec);
tmp = reshape(tmp,[numsolpts*numvec numelec]);

if ~isinversescalar
    W=zeros(numsolpts,3,numelec);
    for k=1:3
        W(:,k,:)=tmp(k:3:end,:);
    end
else
    W=tmp;
end
    

