function CCO=fcm_ccoherevirt(virtfile,chfile,ftype,lofrq,hifrq,nfft)
%FCM_CCOHEREVIRT  Calculates complex coherency of virtual channel files.
%
% ¦CC¦=fcm_ccoherevirt(virtfile,connections,ftype,lofrq,hifrq,¦nfft¦)
%   ( parameters in ¦¦ are optional )
%  
%   CC          output structure containing the complex valued functional connectivity
%               data.
%   virtfile    name of virtual channel file.
%   connections must be a Nx2 matrix indicating the N connections between 2 channels 
%               to be calculated, or the name of a text file containing the connection
%               definitions.
%   ftype       'a'=[a]ll frequency bins, 'b'=frequency [b]and.
%   lofrq       lower limit of the frequency range to be analyzed. 
%   hifrq       upper limit of the frequency range to be analyzed
%   nfft        indicates the number of frequency bins used for FFT. Default is the next
%               power of 2 of the trial length.

CCO=[];
if nargin<5
    help fcm_ccoherevirt
    return
end

tic;

if ischar(chfile)
    % Open connection file
    [chpath,chfile,chext]=fileparts(chfile);
    if isempty(chext), chext='.txt'; end
    ch = load(fullfile(chpath,[chfile chext]));
else    % numeric array
    if nargout==0, error('You must specify an output variable when using a matlab variable for connections definition.'), end
    ch = chfile;
    chfile = 'matvar';
end
    
if ischar(lofrq), lofrq=str2num(lofrq); end
if ischar(hifrq), hifrq=str2num(hifrq); end
if nargin<6 || isempty(nfft)
    nfft=2^nextpow2(fs);
elseif ischar(nfft)
    nfft=str2double(nfft);
end

% Get header of virtual channels file
fid=fopen(virtfile,'r'); 
if fid==-1, error('Invalid virtual channel filename.'), end
test = fread(fid,13,'uint8=>char')';
isv2 = strcmp(test,'NUTMEG VCH v2');
if ~isv2, fclose(fid); error('Not a valid virtual channel file.'), end

datadims=fread(fid,3,'uint16')';    % [lentrial numtrials numvirtualchannels]
fs  = fread(fid,1,'uint16');
tim = fread(fid,datadims(1),'double');
tim = [tim(1) tim(end)];

if any(ch(:)>datadims(3))
    error('At least one of the connection indices does not exist in virtual channel file.')
end
if ( length(unique(ch(:))) < datadims(3) )
    disp('You seem to be analyzing only a subset of voxels. Is this what you want?')
end

% Info Display
fprintf('Sampling Rate: %d\n',fs)

% Freq vectors
frq=[max(lofrq,1/(datadims(1)/fs)) min(hifrq,fs/2)];    % check input freq limits
ffft=[0:fs/nfft:fs/2]';
select = find(ffft>=frq(1) & ffft<=frq(2));
if strcmpi(ftype,'b'), f = [ffft(select(1)) ffft(select(end))];
else,   f  = ffft(select);
end
fd = (fs/nfft)/2;

% Prepare
wdw=repmat(hanning(datadims(1)),[1 datadims(2)]);
ncomp=size(ch,1);
if strcmpi(ftype,'b')
    Coh=zeros(1,ncomp);
else
    Coh=zeros(length(select),ncomp);
end
oHeader = ftell(fid);               % header offset
sVirt   = prod(datadims(1:2));      % data size of 1 virtual channel
oData   = sVirt * 8;                % offset per virtual channel: double is 64 bits = 8 bytes;

perc10=floor(ncomp/10);
if perc10>0, step10=[perc10:perc10:ncomp];
else, step10=ncomp;
end, clear perc10

% Main Part
fprintf('Calculating complex coherency: ')
for c=1:ncomp
    % Find and load current virtual channels
    fseek(fid , oHeader + (ch(c,1)-1)*oData , 'bof');
    x = fread(fid,sVirt,'double');
    x = reshape(x,datadims(1:2));
    fseek(fid , oHeader + (ch(c,2)-1)*oData , 'bof');
    y = fread(fid,sVirt,'double');
    y = reshape(y,datadims(1:2));

    %Pxx = zeros(length(select),1); 
    %Pyy = zeros(length(select),1); 
    %Pxy = zeros(length(select),1); 

    %for k=1:datadims(2)    % numtrials
    %xw = x(:,k); 
    %yw = y(:,k);             
    x = wdw.*detrend(x);  % Windowing and linear detrend
    y = wdw.*detrend(y);
    Xx = fft(x,nfft,1);
    Yy = fft(y,nfft,1);
    Xy = Yy.*conj(Xx);
    Xx = sum(abs(Xx(select,:)).^2,2);
    Yy = sum(abs(Yy(select,:)).^2,2);
    Xy = sum(Xy(select,:),2);   
    %end        
    switch lower(ftype)
    case 'b'
        Xx = sum(Xx);
        Yy = sum(Yy);
        Xy = sum(Xy); 
        Coh(c) = Xy/sqrt(Xx*Yy);             % complex coherence 
    case 'a'
        Coh(:,c) = Xy./sqrt(Xx.*Yy);         % complex coherence 
    end
    if any(step10==c)
        fprintf('...%d%%',ceil((c/ncomp)*100))
    end    
end
fclose(fid);
fprintf('\n')

% Output
CC.coh=Coh;
CC.comps=ch;
CC.frq=[f(:,1)-fd f(:,end)+fd];
CC.time=tim;
CC.N=datadims(2);

if nargout==0
    if ~exist('comcoh','dir'), mkdir comcoh, end
    save(fullfile('comcoh',['CC' chfile(2:end)]),'CC')
else
    CCO=CC;
end

fprintf('Done (fcm_ccoherevirt).\n')
tempus=toc;
fprintf('Elapsed time is %d minutes and %d seconds.\n',round(tempus/60),round(rem(tempus,60)))
