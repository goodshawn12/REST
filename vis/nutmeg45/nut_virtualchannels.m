function nut_virtualchannels(virtfile,filtdatafile,weightfile)
% NUT_VIRTUALCHANNELS  Calculates and saves virtual channels. 
% The created file is in binary format with the extension .vch
%
% Usages:
%   All inputs are optional.
%   nut_virtualchannels(virtfile,filtdatafile,weightfile)
%   nut_virtualchannels(virtfile,filtdata,W)



% Load filtereddata
if (nargin<2 || isempty(filtdatafile))
    [fS, pS] = uigetfile('filtdata_*.mat','Select filtered data file');
    if isequal(fS,0), return, end
    filtdatafile=fullfile(pS,fS);
end
if ischar(filtdatafile)
    load(filtdatafile);
	if exist('meg','var'), F=meg; clear meg, end
    if ~exist('F','var'), error('Not a valid filtdata file.'), end
else
    F=filtdatafile;
end
clear filtdatafile

% Load weights file
if (nargin<3 || isempty(weightfile))
    [fW, pW] = uigetfile('weights_*.mat;W_*.mat;*.is','Select weights file');
    if isequal(fW,0), return, end
    weightfile=fullfile(pW,fW);
end    
if ischar(weightfile)
    [pW,fW,eW]=fileparts(weightfile);
    switch eW
        case {'.mat' ''}
            load(weightfile)
            if ~exist('W','var') 
                if exist('Wact','var')
                    W = Wact; clear Wact
                else
                    error('Not a valid weights file.')
                end
            end
        case '.is'
            W=nut_import_smacis(weightfile);
            if ndims(W)>2
                W = permute(W,[3 2 1]);
            else
                W = W';
            end
        otherwise
            error('unknown weight file format.')
    end
            
else
    W = weightfile;
end

if ndims(W)>2
    [nsens,nor,nvox] = size(W);
    if nor>1
        [fO, pO] = uigetfile('*_ori.mat','Select voxel orientation file');
        if isequal(fO,0), return, end
        load(fullfile(pO,fO));

        Wn = zeros(nsens,nvox);
        for k=1:nsens
            Wn(k,:)=dot(squeeze(W(k,:,:)),ori');
        end
        W = Wn;
        clear Wn
    else
        W = squeeze(W);
    end
end

% Prompt for virtual channels file name if necessary
if (nargin<1 || isempty(virtfile))
    if exist('fW','var')
        [fV, pV] = uiputfile('*.vch', 'Save as...',[strrep(fW,'W_','virtual_') '.vch']);
    else
        [fV, pV] = uiputfile('*.vch', 'Save as...');
    end
    if isequal(fV,0), return, end
    virtfile=fullfile(pV,fV);
elseif ~strcmpi(virtfile(end-3:end),'.vch')
    virtfile = [virtfile '.vch'];
end

% Legacy compatibility
isnewver=isstruct(F);
if ~isnewver
    tmp=F;
    F=struct('data',tmp);  clear tmp
end

% Calculate and write to file
fid = fopen(virtfile,'w');

[lentrial,dum,numtrial] = size(F.data);
numvirt = size(W,2);
if ~isnewver
    F.srate=0;
    F.latency=zeros(lentrial,1);
end

fprintf('Writing: ')
perc10=floor(numvirt/10);
if perc10>0, step10=[perc10:perc10:numvirt];
else, step10=numvirt;
end, clear perc10

fwrite(fid,'NUTMEG VCH v2','uint8');
fwrite(fid,[lentrial;numtrial;numvirt;F.srate],'uint16');
fwrite(fid,F.latency,'double');

for cc=1:numvirt
    for kk=1:numtrial
        D(:,kk) = F.data(:,:,kk) * W(:,cc);
    end
    fwrite(fid,D(:),'double');     % 64 bits = 8 bytes
	if any(step10==cc)
        fprintf('...%d%%',ceil((cc/numvirt)*100))
    end    
end
fprintf('\n')

fclose(fid); 

