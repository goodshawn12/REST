function [active,control]=nut_read_SAMparam(paramfile)

fid=fopen(paramfile);

fscanf(fid,'%d',1); % skip first line
fscanf(fid,'%s',1); % skip marker string
active = fscanf(fid,'%f',2); % read active window (two doubles)
fscanf(fid,'%d',1); % skip first line
fscanf(fid,'%s',1); % skip marker string
control = fscanf(fid,'%f',2); % read control window (two doubles)
fclose(fid);

% convert to ms and transpose
active = 1000*active';
control = 1000*control';