function spli1(filtfile,wfile,chfile,n,stp,ftype,varargin)
%SPLI1  calculates source phase lag index for continuous single-trial data.
%
%  spli1(filtfile,weightfile,connections,len,step,ftype,¦lofrq¦,¦hifrq¦,¦nfft¦,¦options¦)
%        ( parameters in ¦¦ are optional )
% 
% Input:
%  filtfile     file containing filtered sensor data (filtdata_*.mat).
%  weightfile   file containing weight matrix from inverse solution (W_*.mat)
%  connections  can be a Nx2 matrix containing the indices of voxel pairs
%               to be analyzed (created with fcm_conndef) or the name of a text 
%               file containing these indices (created with fcm_preparejobs).
%  len          the length of one data segment in sec. 
%  step         the time step between data segments in sec. 
%  ftype        's'=frequency [s]pectrum, 'b'=frequency [b]and.
%  lofrq        lower bound of frequencies to be analyzed (in Hz).
%  hifrq        upper bound of frequencies to be analyzed (in Hz).
%  nfft         Indicates the number of frequency bins created by FFT. 
%               Should be a power of 2.
%  options      parameter/value pairs to specify additional properties:
%               For multi-tapering with Slepian sequences: 
%                   spli1(...,'taper',freqwidth_in_Hz)
%               For calculating PLI with extracerebral channels
%               (e.g., EMG for cortico-muscular PLI)
%                   spli1(...,'ext',extdata), with extdata being the non-
%                   event-locked multi-trial extracerebral data (samples X
%                   channels X trials) or a corresponding matlab file containing 
%                   variable "ext".
%               If ftype is 'b', you can specify the method to calculate the PLI:
%                   spli31(...,'bandmethod','hilbert') is the official but very slow way
%                   spli31(...,'bandmethod','fasthilbert') also uses a Hilbert transform,  
%                       but then estimates the phase difference in each data segment rather 
%                       than at each sample (default).
%                   spli31(...,'bandmethod','fft') used a FFT for calculating phase 
%                       differences (use at own risk).
%               To calculate FC between anatomical ROIs:
%                   spli1(...,'roi',roifile), where roifile contains a structure 
%                   R which can be created with fcm_voxel2roi.
%               To use independent components:
%                   spli1(...,'ica',icadata), where icadata contains
%                   variable IC, the non-artifact independent components and 
%                   variable A, the mixing matrix.
%
% Output:
%  results are saved in a directory "pli".
%
% References:
%   - Stam et al. Hum Brain Mapp 2007; 28: 1178–1193


if nargin<6
    fprintf('\nusage: spli1 filtfile wfile chfile len step ftype [lofrq] [hifrq] [nfft] [options]\n\n')
    return
end

tic;

if ischar(n),     n=str2double(n); end
if ischar(stp),   stp=str2double(stp); end

% Open all files
if ischar(chfile)
    [chpath,chfile,chext]=fileparts(chfile);
    if isempty(chext), chext='.txt'; end
    ch = load(fullfile(chpath,[chfile chext]));
else    % if comps given as matrix
    ch=chfile;
    chfile='comps';
end
ncomp=size(ch,1);

load(filtfile);
if exist('F','var'), meg=F; clear F, end    % legacy compatibility
[dum,filtfile,dum]=fileparts(filtfile);
[lentrial,numsens,numtrial] = size(meg.data);
frq = sscanf(fliplr(strtok(fliplr(filtfile),'_')),'%fto%fHz')';

[dum,dum,eW]=fileparts(wfile);
switch eW
    case {'.mat' ''}
        load(wfile)
        if ~exist('W','var') 
            if exist('Wact','var')
                W = Wact; clear Wact
            else
                error('Not a valid weights file.')
            end
        end
    case '.is'
        W=nut_import_smacis(wfile);
        if ndims(W)>2
            W = permute(W,[3 2 1]);
        else
            W = W';
        end
    otherwise
        error('unknown weight file format.')
end
if ndims(W)>2
    [nsens,nor,nvox] = size(W);
    if nor>1
        if ~exist('ori','var')
            [fO, pO] = uigetfile('*_ori.mat','Select voxel orientation file');
            if isequal(fO,0), return, end
            load(fullfile(pO,fO));
        end

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
[dum,wfile]=fileparts(wfile); clear dum

% Deal with variable inputs.
mtaper=false; isext=false; forcehilbert=(ftype=='b'); useroi=false; ds=1; bandmethod='fasthilbert'; %defaults
lofrq=[]; hifrq=[]; nfft=[];
if ~isempty(varargin)
    ii=find(strcmpi('taper',varargin),1); 
    mtaper = ~isempty(ii);
    if mtaper
        taper = varargin{ii+1};
        if ischar(taper), taper=str2double(taper); end
        varargin(ii:ii+1)=[];
    end
    
    ii=find(strcmpi('ext',varargin),1); 
    isext = ~isempty(ii);
    if isext
        ext = varargin{ii+1};
        if ischar(ext)
            load(ext);  %must contain variable ext
        end
        if ~isequal([lentrial numtrial],[size(ext,1) size(ext,3)])
            error('Extracerebral data must have same trial length and same number of trials as cerebral data.')
        end
        fprintf('Using extracerebral channels.\n')
        varargin(ii:ii+1)=[];
    end
    
    ii=find(strcmpi('roi',varargin),1); 
    useroi=~isempty(ii);
    if useroi
        fprintf('Using ROIs.\n')
        roifile=varargin{ii+1};
        load(roifile,'R');   % must contain structure R
        roidef = R.roidef;
        W = W(:,R.goodvoxels);
        W = W * R.voxel2roi_tfm; clear R
        varargin(ii:ii+1)=[];
    end
    
    ii=find(strcmpi('ica',varargin),1); 
    if ~isempty(ii)
        fprintf('Using ICA data.\n')
        load(varargin{ii+1},'IC','A');   % must contain variables IC and A
        if ~isequal([size(meg.data,1) size(meg.data,3)],[size(IC,1) size(IC,3)])
            error('Independent components must have same trial length and same number of epochs as original data.')
        end
        meg.data=IC; clear IC
        numsens = size(meg.data,2);
        gain=max(abs(A));
        A=A./repmat(gain,[size(A,1) 1]);
        W = A' * W;
        W = W.*repmat(gain',[1 size(W,2)]);
        varargin(ii:ii+1)=[];
    end
    
    ii=find(strcmpi('bandmethod',varargin),1);
    if ~isempty(ii) 
        if ftype=='b'
            bandmethod=lower(varargin{ii+1}); 
            forcehilbert=~isempty(strfind(bandmethod,'hilbert'));
        end
        varargin(ii:ii+1)=[];
    end   
    
    ii=find(strcmpi('filt',varargin),1); 
    if ~isempty(ii)
        if forcehilbert
            load(varargin{ii+1});   % must contain filt structure with filter settings
            if ~exist('filt','var'), error('Invalid filter configuration file.'), end
        end
        varargin(ii:ii+1)=[];        
    end    
    
    ii=find(strcmpi('ds',varargin),1); 
    if ~isempty(ii)
        if strcmp(bandmethod,'hilbert')
            ds = varargin{ii+1};
            if ischar(ds), ds=str2num(ds); end        
            fprintf('Looking only at 1/%d samples.\n',ds)
        end
        varargin(ii:ii+1)=[];
    end    
    
    switch length(varargin)
    case 3
        lofrq=varargin{1}; hifrq=varargin{2}; nfft=varargin{3};
    case 2
        lofrq=varargin{1}; hifrq=varargin{2};
    case 1
        nfft=varargin{1};
    end
    if ischar(lofrq), lofrq=str2num(lofrq); end
    if ischar(hifrq), hifrq=str2num(hifrq); end
    if ischar(nfft), nfft=str2num(nfft); end
        
end
if mtaper && ftype=='b', error('Multi tapering does not make sense when assessing freq bands.'), end
clear varargin

% Info Display
fprintf('\nspli1 %s %s %s %3.1f %3.1f %s %3.1f %3.1f %d\n\n',filtfile,wfile,chfile,n,stp,ftype,lofrq,hifrq,nfft)

% Prepare
fs = meg.srate;
tim = [meg.latency(1) meg.latency(end)];
ns = floor(n*fs);
stp= floor(stp*fs);
nbt= floor((lentrial-ns)/stp) + 1;    % number of time points
all2all = (~isext && ncomp>nbt)
numvirt=size(W,2);
if forcehilbert
    if all2all, numh=numvirt; else numh=numtrial; end
    if ~rem(ns,2)    % even 
        h    = 2*ones(ns/2+1,numh);
        h([1 end],:) = 1;
        neg0 = zeros(ns/2-1,numh);
        select = [1:ns/2+1];
    else    % odd 
        h    = 2*ones((ns+1)/2,numh);
        h(1,:) = 1;
        neg0 = zeros((ns-1)/2,numh);
        select = [1:(ns+1)/2];
    end
    numfrq = length(select);
    nfft = ns;
    if ~isempty(lofrq) && ( lofrq-frq(1)>1 || hifrq-frq(2)<-1 )
        frq=[lofrq hifrq];
        if ~exist('filt','var')  % default filter parameters
            load ellipauto
        end
        fprintf('Filtering sensor data: ')
        meg.data = nut_filter2(meg.data,filt.class,filt.type,filt.order,lofrq,hifrq,meg.srate,1);
        if isext
            ext = nut_filter2(ext,filt.class,filt.type,filt.order,lofrq,hifrq,meg.srate,1);
        end
        fprintf('done.\n')        
    end
else
    if isempty(lofrq), lofrq=frq(1); hifrq=frq(2); end
    frq=[lofrq min(hifrq,fs/2)];    % check input freq limits
    ffft=[0:fs/nfft:fs/2]';
    select = find(ffft>=frq(1) & ffft<=frq(2));
    numfrq = length(select);
    fv  = ffft(select);
    fd = (fs/nfft)/2;
end

if mtaper
    taper = [ns/fs taper];
    numtaper = floor( 2*taper(2)*taper(1) - 1);
    NW = taper(1)*taper(2);
    N  = taper(1)*fs;     
    wdw = dpss(N, NW, 'calc');
    wdw = wdw(:,1:numtaper);
    wdw = permute( repmat(wdw,[1 1 numsens]), [1 3 2]);  
    fprintf('Using %d tapers.\n',numtaper)
else
    if ~forcehilbert
        wdw=repmat(hanning(ns),[1 numsens]); 
        if isext, wdwe = wdw(:,size(ext,2),:); end
    end
    numtaper=1;
end

% Sensor Fourier Transform
%-------------------------
fprintf('Calculating sensor Fourier transform: ')

try
    FF = complex(zeros(numfrq,numsens,numtaper,nbt));
catch ME
    if strcmp(ME.identifier,'MATLAB:nomem');
        fprintf('\nOut of memory, trying to use single precision... ')
        FF=complex(zeros(numfrq,numsens,numtaper,nbt,'single'));
        fprintf('ok.\n')
    else
        rethrow(ME)
    end
end    
index = 1:ns;
for k=1:nbt
    tmp = detrend(meg.data(index,:));  
    if ~forcehilbert, tmp = wdw .* repmat(tmp,[1 1 numtaper]); end
    tmp = fft(tmp, nfft, 1); 
    FF(:,:,:,k) = tmp(select,:,:);
    index = index + stp; 
end
clear meg
if isext 
    EE = complex(zeros(numfrq,size(ext,2),numtaper,nbt,class(FF)));
    index = 1:ns;
    for k=1:nbt
        tmp = detrend(ext(index,:));  
        if ~forcehilbert, tmp = wdwe .* repmat(tmp,[1 1 numtaper]); end
        tmp = fft(tmp,nfft,1);
        EE(:,:,:,k) = tmp(select,:,:);
        index = index + stp; 
    end
end
clear tmp ext
fprintf('done.\n')

if forcehilbert
    % Try to speedup some more
    % We can remove the filtered frequencies for now and then add them back as 0 later.
    D=mean(mean(abs(diff(log10(abs(FF(:,1:3:end,1:3:end))+1e-12))),3),2);     % find bins where spectrogram approaches 0
    f=(D>0.02); clear D
    if (numfrq-sum(f)>10) && (sum(abs(diff(f)))<3)  % only possible if filtered data present
        select=[1;find(f)+1]; clear f
        fprintf('Speedup on: working only on %d of %d freq bins.\n',length(select),ns)    
        FF = FF(select,:,:,:);
        if isext, EE = EE(select,:,:,:); end
    end
else
    select = 1:numfrq;
end

% Source connectivity
% -------------------
fprintf('Calculating source PLI: ')
if all2all
    % For large numbers of connections, this might be faster

    if ftype=='b'
        Coh = zeros(ncomp,1);
        perc10=floor(nbt/10);
        if perc10>=1, numstep=10; step10=round(linspace(perc10,nbt,numstep));
        else step10=numfrq; numstep=1;
        end, clear perc10
        perccount=0;        
        XY = zeros(numvirt);
        for k=1:nbt
            x = zeros(numfrq,numvirt); 
            x(select,:) = FF(:,:,1,k) * W;
            
            switch bandmethod
                case 'fasthilbert'
                    % this is the same as x=hilbert(data)
                    x = ifft(cat(1,h.*x,neg0)); 
                    % get the PLI of the entire matrix real quick.
                    XY = XY + sign(angle(x' * real(x))); % we could window the signal here (.* wdw) but hilbert does not do it either.

                case 'hilbert' 
                    % this is the same as x=angle(hilbert(data))
                    x = angle(ifft(cat(1,h.*x,neg0))); 
                    for t=1:ds:ns
                        X = repmat(x(t,:),[numvirt 1]);
                        XY = XY + sign(sin(X-X.'));
                    end         

                case 'fft'
                    warning('Calculating PLI with FFT is inofficial and highly experimental')
                    x = mean(x,1);
                    x = angle(x);
                    X = repmat(x,[numvirt 1]);
                    XY = XY + sign(sin(X-X.'));
            end
            
            if any( step10==k )
                perccount=perccount+100/numstep;
                fprintf('...%d%%',perccount)
            end               
        end
        
        switch bandmethod
            case 'hilbert'
                XY = XY ./ (nbt*length(1:ds:ns));
            otherwise
                XY =XY ./ nbt;
        end

        % We bring the matrix to vector form to save more than half of
        % the memory        
        for c=1:ncomp
            Coh(c) = XY(ch(c,1),ch(c,2)); 
        end
        
    else
        Coh=zeros(ncomp,1,numfrq);
        perc10=floor(numfrq/10);
        if perc10>=1, numstep=10; step10=round(linspace(perc10,numfrq,numstep));
        else step10=numfrq; numstep=1;
        end, clear perc10
        perccount=0;

        for f=1:numfrq
            x = zeros(nbt,numvirt,numtaper); 
            for ee=1:numtaper 
                x(:,:,ee) = permute(FF(f,:,ee,:),[4 2 1 3]) * W;
            end       
            if mtaper
                x=mean(x,3);
            end
            
            x = angle(x); 
            
            XY = zeros(numvirt);
            for t=1:nbt
                X = repmat(x(t,:),[numvirt 1]);
                XY = XY + sign(sin(X-X.'));
            end    
            
            XY = XY ./ nbt;  % mean    

            % We bring the matrix to vector form to save more than half of
            % the memory
            for c=1:ncomp
                Coh(c,1,f) = XY(ch(c,1),ch(c,2)); % complex coherence 
            end   

            if any( step10==f )
                perccount=perccount+100/numstep;
                fprintf('...%d%%',perccount)
            end   
        end % for f=1:numfrq
    end

    clear X XY
    fprintf('\n')
    
else  % if only few connections need to be calculated, it's faster this way
    
    if ftype=='b', Coh=zeros(ncomp,1);
    else Coh=zeros(ncomp,numfrq);
    end
    perc10=floor(ncomp/10);
    if perc10>10, step10=[perc10:perc10:ncomp]; numstep=10;
    else step10=ncomp; numstep=1;
    end, clear perc10
    perccount=0;
    for c=1:ncomp

        x = zeros(numfrq,nbt,numtaper); y=x;
        for k=1:nbt      
            for ee=1:numtaper 
                x(select,k,ee) = FF(:,:,ee,k) * W(:,ch(c,1));
                if ~isext, y(select,k,ee) = FF(:,:,ee,k) * W(:,ch(c,2)); end
            end 
        end
        if isext            
            y = permute(EE(:,ch(c,2),:,:),[1 4 3 2]);  
        end
        
        if ftype=='b'
            % this is the same as X=mean(angle(hilbert(data)))
            x = angle(ifft(cat(1,h.*x,neg0))); 
            y = angle(ifft(cat(1,h.*y,neg0))); 
            Coh(c) = mean(mean(sign(sin(x-y)),1),2);
        else
            if mtaper
                x = mean(x,3); y=mean(y,3);
            end
            x = angle(x).';
            y = angle(y).';
            Coh(c,:) = mean(sign(sin(x-y)),1);
        end
               
        if any(step10==c)
            perccount=perccount+100/numstep;
            fprintf('...%d%%',perccount)
        end    
    end
    if ftype~='b', Coh=reshape(Coh,[ncomp 1 numfrq]); end

    fprintf('\n')

end  % if all2all

% Output
CC.coh=Coh;
CC.method = 'pli';
CC.diminfo = 'connection * time * frequency';
CC.comps=ch;
if ftype=='b'
    CC.frq=frq;
else
    CC.frq=[fv(:,1)-fd fv(:,end)+fd];
end
CC.time=tim;
CC.N=numtrial;
CC.len=lentrial;
if mtaper
    CC.taperwidth=taper(2);
    CC.numtaper=numtaper;
end
if useroi
    CC.roidef = roidef;
    CC.roifile = roifile;
end

if ~exist('pli','dir'), mkdir pli, end
save(fullfile('pli',['CC' chfile(2:end)]),'CC')

fprintf('Done (spli1).\n')
tempus=toc;
fprintf('Elapsed time is %d minutes and %d seconds.\n',round(tempus/60),round(rem(tempus,60)))
