function [X,tim,fs] = nut_read_vch(virtfile,chan,trials)
% NUT_READ_VCH  reads nutmeg virtual channels in .vch format
%   [data,timevector,srate] = nut_read_vch(virtfile,channels,trials)
%   [data,timevector,srate] = nut_read_vch(virtfile,'all',trials)
%   [header,timevector,srate] = nut_read_vch(virtfile,0)
%
% virtfile      filename to read (*.vch)
% channels      indices of virtual channels to be read. Set to 0 to obtain header 
%               data in 'data' output [num_samples num_trials num_virtualchannels].
%               Set to 'all' to read all virtual channels (but you may run
%               out of memory).
% trials        (optional) indices of trials to be read. Default is all trials.

if ~exist(virtfile,'file'), error('NUT_READ_VCH: file not found.'), end

% Initialize
X  = [];
fs = [];
tim= [];

% Open virtual channels file
fid=fopen(virtfile,'r'); 

% Get version
test = fread(fid,13,'uint8=>char')';
isv2 = strcmp(test,'NUTMEG VCH v2');
if ~isv2, fseek(fid,0,'bof'); end

% Read Header
datadims=fread(fid,3,'uint16')';    % [lentrial numtrials numvirtualchannels]
if any(chan > datadims(3))
    fclose(fid);
    error(sprintf('NUT_READ_VCH: the file only contains %d virtual channels.',datadims(3)))
end
if ischar(chan), chan=1:datadims(3); end
if nargin<3 || ischar(trials), trials=1:datadims(2); end

if isv2
    fs  = fread(fid,1,'uint16');
    tim = fread(fid,datadims(1),'double');
end

if (chan==0)
    %fprintf('Samples:        \t%d\n',datadims(1))
    %fprintf('Trials:         \t%d\n',datadims(2))
    %fprintf('Virtual Channels\t%d\n',datadims(3))
    X = datadims;
    fclose(fid);
    return
end

% Offsets
oHeader = ftell(fid);               % header offset
sVirt   = prod(datadims(1:2));      % data size of 1 virtual channel
oData   = sVirt * 8;                % offset per virtual channel: double is 64 bits = 8 bytes;
numchan = length(chan);
numtri  = length(trials);

% Read data
X = zeros(datadims(1),numchan,numtri);
for cc=1:numchan
    fseek(fid , oHeader + (chan(cc)-1)*oData , 'bof');
    x = fread(fid,sVirt,'double');
    temp = reshape(x,[datadims(1) 1 datadims(2)]);
    X(:,cc,:) = temp(:,1,trials);
end

fclose(fid);