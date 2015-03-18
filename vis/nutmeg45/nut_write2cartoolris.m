function nut_write2cartoolris(risfilename,data,srate)
% NUT_WRITE2CARTOOLRIS  writes "result of inverse solution" (RIS) files for
%                       display in Cartool
%
%   nut_write2cartoolris(risfilename,data,samplingrate)
%
% risfilename       Name of new RIS file
% data              Inverse solution data (2D matrix, voxels*time)
% samplingrate      in Hz


% Prepare Data
numsolutionpoints = size(data,1);
numtimeframes     = size(data,2);
samplingfrequency = srate;
isinversescalar   = 1;              % Set to 0 here if you want to export vectorial solutions.

fid = fopen(risfilename,'w');

% RIS Header
fwrite(fid,'RI01','char');
fwrite(fid,numsolutionpoints,'int32');
fwrite(fid,numtimeframes,'int32');
fwrite(fid,samplingfrequency,'float32');
fwrite(fid,isinversescalar,'char');

% Data
fwrite(fid,data,'float32');     % float = single!

fclose(fid); 