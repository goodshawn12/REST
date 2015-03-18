function R = nut_import_cartoolris(risfilename)
% reads "result of inverse solution" (RIS) files from Cartool.
%   R = nut_import_cartoolris(risfilename)


fid = fopen(risfilename,'r');

% RIS Header
test=fread(fid,4,'*char')';
if ~strcmp(test,'RI01'), error('Not a valid RIS file.'), end
numsolutionpoints=fread(fid,1,'int32');
numtimeframes=fread(fid,1,'int32');
samplingfrequency=fread(fid,1,'float32');
isinversescalar=fread(fid,1,'char');

% Data
if isinversescalar
    si = [numsolutionpoints numtimeframes];
else
    si = [numsolutionpoints*3 numtimeframes];
end
data=fread(fid,si,'float32');     % float32 = single!

% Close file
fclose(fid);

% Bring to NUTMEG format
if isinversescalar
    R = data.';
else
    R=zeros(numtimeframes,3,numsolutionpoints);
    for k=1:3
        R(:,k,:)=data(k:3:end,:).';
    end
end