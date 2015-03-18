function virtu(filtdatafile,weightfile,currjob,numjobs)
% VIRTU  Creating of virtual SAM channels. For compilation and batch use.
%
% Usage:
%   virtu filtdatafile weightfile currentjob

if nargin<3
    disp(' ')
    disp('Usages:')
    disp('virtu filtdatafile weightfile currentjob')
    disp('virtu filtdatafile weightfile currentjob totaljobs')
    disp(' ')
    return
end

if nargin<4, numjobs=5; 
elseif ischar(numjobs), numjobs=str2double(numjobs);
end

if ischar(currjob), currjob=str2double(currjob); end

% Load filtereddata
load(filtdatafile);
if ~exist('F','var'), error('Not a valid filtdata file.'), end
if ~isstruct(F), error('filtdata file does not contain F structure.'), end

% Load weights file
load(weightfile)
if ~exist('W','var'), error('Not a valid weights file.'), end

% Calculate and write to file
[lentrial,dum,numtrial] = size(F.data);
numvirt = size(W,2);
junks = ceil(numvirt/numjobs);

fid = fopen(sprintf('virtfile.%03d',currjob),'w');
if currjob==1
    fwrite(fid,'NUTMEG VCH v2','uint8');
    fwrite(fid,[lentrial;numtrial;numvirt;F.srate],'uint16');
    fwrite(fid,F.latency,'double');
end
curv=(currjob-1)*junks + [1:junks];
if any(curv>numvirt)
    if currjob<numjobs
        error('Strange things happen here')
    end
    curv=[curv(1):numvirt];
end
for cc=curv
    for kk=1:numtrial
        D(:,kk) = F.data(:,:,kk) * W(:,cc);
    end
    fwrite(fid,D(:),'double');     % 64 bits = 8 bytes 
end
fclose(fid); 

