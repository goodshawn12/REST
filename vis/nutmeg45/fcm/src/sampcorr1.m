function sampcorr1(filtfile,wfile,chfile,n,stp,ftype,varargin)
%SAMPCORR1  calculates source amplitude envelope correlations for continuous single-trial data.
%           An estimate of  mean signal amplitude is obtained for each data segment 
%           and correlated between voxels. 
%
%  sampcorr1(filtfile,weightfile,connections,len,step,ftype,¦lofrq¦,¦hifrq¦,¦nfft¦,¦options¦)
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
%                   sampcorr1(...,'taper',freqwidth_in_Hz)
%               For calculating amplitude correlation with extracerebral channels
%               (e.g., EMG for cortico-muscular correlation)
%                   sampcorr1(...,'ext',extdata), with extdata being the non-
%                   event-locked multi-trial extracerebral data (samples X
%                   channels X trials) or a corresponding matlab file containing 
%                   variable "ext".
%               To calculate FC between anatomical ROIs:
%                   sampcorr1(...,'roi',roifile), where roifile contains a structure 
%                   R which can be created with fcm_voxel2roi.
%               To use independent components:
%                   sampcorr1(...,'ica',icadata), where icadata contains
%                   variable IC, the non-artifact independent components and 
%                   variable A, the mixing matrix.
%               If ftype='b' you can specify the method for obtaining band amplitudes:
%                   sampcorr1(...,'bandmethod',method), where method is either
%                   'hilbert' (default) or 'FFT'. FFT is faster but tends to loose
%                   some information due to the required window taper. If ftype='s',
%                   method is always 'FFT'.

%
% Output:
%  results are saved in a directory "ampcorr".
%
% References:
%   - Bruns et al. NeuroReport 2000; 11: 1509-1514
%   - Brookes et al. PNAS 2011; 108: 16783-16788


if nargin<6
    fprintf('\nusage: sampcorr1 filtfile wfile chfile len step ftype [lofrq] [hifrq] [nfft] [options]\n\n')
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
mtaper=false; isext=false; forcehilbert=(ftype=='b'); useroi=false; 
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
        if ftype=='b', forcehilbert=strcmpi(varargin{ii+1},'hilbert'); end
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
fprintf('\nsampcorr1 %s %s %s %3.1f %3.1f %s %3.1f %3.1f %d\n\n',filtfile,wfile,chfile,n,stp,ftype,lofrq,hifrq,nfft)

% Prepare
fs = meg.srate;
tim = [meg.latency(1) meg.latency(end)];
ns = floor(n*fs);
stp= floor(stp*fs);
nbt= floor((lentrial-ns)/stp) + 1;    % number of time points
all2all = (~isext && ncomp>nbt);
numvirt=size(W,2);
if ftype=='b' && forcehilbert
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
fprintf('Calculating source amplitude correlation: ')
if all2all
    % For large numbers of connections, this is much much much much .... faster
    if ftype=='b'
        Coh = zeros(ncomp,1);
        x = zeros(nbt,numvirt);
        for k=1:nbt
            s = complex(zeros(numfrq,numvirt)); 
            s(select,:) = FF(:,:,1,k) * W;
            
            if forcehilbert
                % this is the same as X=mean(abs(hilbert(data)))
                x(k,:) = mean(abs(ifft(cat(1,h.*s,neg0))),1); 
            else
                x(k,:) = abs(mean(s,1));
            end
        end
        clear s
        x = zscore(x);

        XY = x.' * x;  clear x
        XY = XY ./ sqrt( diag(XY) * diag(XY).' );
        %R(c) = (x0./norm(x0))' * (y0 ./ norm(y0));

        % We bring the matrix to vector form to save more than half of
        % the memory
        for c=1:ncomp
            Coh(c) = XY(ch(c,1),ch(c,2)); 
        end      
        fprintf('...100%%')
        
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
            x = zscore(abs(x)); % get normalized amplitude
            XY = x.' * x;  clear x
            XY = XY ./ sqrt( diag(XY) * diag(XY).' );  % amplitude correlation

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
            % this is the same as X=mean(abs(hilbert(data)))
            x = mean(abs(ifft(cat(1,h.*x,neg0))),1); 
            y = mean(abs(ifft(cat(1,h.*y,neg0))),1); 
        else
            if mtaper
                x = mean(x,3); y=mean(y,3);
            end
            x = abs(x);
            y = abs(y);
        end
        x = zscore(x.');
        y = zscore(y.');
        
        Coh(c,:) = sum(x.*y,1) ./ sqrt( sum(x.^2,1) .* sum(y.^2,1) );
        
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
CC.method = 'ampcorr';
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

if ~exist('ampcorr','dir'), mkdir ampcorr, end
save(fullfile('ampcorr',['CC' chfile(2:end)]),'CC')

fprintf('Done (sampcorr1).\n')
tempus=toc;
fprintf('Elapsed time is %d minutes and %d seconds.\n',round(tempus/60),round(rem(tempus,60)))
