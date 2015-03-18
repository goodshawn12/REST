function nut_add_channel2vch(virtfile,data,newvirtfile)
% nut_add_channel2vch(virtfile,data,¦newvirtfile¦)

if nargin<2, help nut_add_channel2vch, return, end

if ~exist(virtfile,'file'), error('NUT_ADD_CHANNEL2VCH: virtual channel file not found.'), end

if nargin>2,
    copyfile(virtfile,newvirtfile);
    virtfile=newvirtfile;
end

fid=fopen(virtfile,'r+'); 

% Get version
test = fread(fid,13,'uint8=>char')';
isv2 = strcmp(test,'NUTMEG VCH v2');
if ~isv2, fseek(fid,0,'bof'); end

% Read Header
datadims=fread(fid,3,'uint16')';    % [lentrial numtrials numvirt]

if ~all(datadims(1:2) == [size(data,1),size(data,3)])
    fclose(fid);
    error('NUT_ADD_CHANNEL2VCH: data to add must have the same dimensions as the virtual channels.')
end

% Update with new channel count
fseek(fid,-2,0);
fwrite(fid,datadims(3)+size(data,2),'uint16');

% Add new channels
fseek(fid,0,1);
for k=1:size(data,2)
    D=data(:,k,:);
    fwrite(fid,D(:),'double');
end

fclose(fid);
