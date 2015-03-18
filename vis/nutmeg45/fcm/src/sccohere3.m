function sccohere3(filtfile,wfile,chfile,n,stp,ftype,lofrq,hifrq,nfft,varargin)
%SCCOHERE3  Calculates source complex coherency for event-locked multi-trial data.  
%
%  sccohere3(filtfile,weightfile,connections,len,overlap, ...
%        ftype,lofrq,hifrq,nfft,¦options¦)
%          (parameters in ¦¦ are optional)
% 
% Input:
%  filtfile     file containing filtered sensor data (filtdata_*.mat).
%  weightfile   file containing weight matrix from inverse solution (W_*.mat)
%  connections  can be a Nx2 matrix containing the indices of voxel pairs
%               to be analyzed (created with fcm_conndef) or the name of a text 
%               file containing these indices (created with fcm_preparejobs).
%  len          the length of one coherence time-window in sec. 
%  overlap      the overlap between coherence time-windows in sec. 
%  ftype        's'=frequency [s]pectrum, 'b'=frequency [b]and.
%  lofrq        lower bound of frequencies to be analyzed (in Hz).
%  hifrq        upper bound of frequencies to be analyzed (in Hz).
%  nfft         indicates the number of frequency bins created by FFT. 
%               Should be a power of 2.
%  options      parameter/value pairs to specify additional properties:
%               For multi-tapering with Slepian sequences: 
%                   sccohere3(...,'taper',freqwidth_in_Hz)
%               For calculating complex coherence with extracerebral channels
%               (e.g., EMG for cortico-muscular coherence)
%                   sccohere3(...,'ext',extdata), with extdata being the event-locked
%                   multi-trial extracerebral data (samples X channels X trials) 
%                   or a corresponding matlab file containing variable "ext".
%               To calculate FC between anatomical ROIs:
%                   sccohere3(...,'roi',roifile), where roifile contains a structure 
%                   R which can be created with fcm_voxel2roi
%               To use independent components:
%                   sccohere3(...,'ica',icadata), where icadata contains
%                   variable IC, the non-artifact independent components and 
%                   variable A, the mixing matrix.
%               To obtain additional baseline corrected outputs:
%                   sccohere3(...,'baseline',basestart,basestop), where
%                   basestart and basestop are seconds relative to the
%                   beginning of each trial.
%
% Output:
%  results are saved in a directory "ccohere".
%
% Use:  SCCOHERE1 for continuous single-trial data.
%       SCCOHERE31 for non-event-locked multi-trial data.
%
% Reference:
%  - Nolte et al. Clin Neurophysiol 2004; 115: 2292-2307

if nargin<9
    fprintf('\nusage: sccohere3 filtfile wfile chfile ftype lofrq hifrq nfft len overlap basestart basestop [options]\n\n')
    return
end

tic;

if ischar(lofrq), lofrq=str2num(lofrq); end
if ischar(hifrq), hifrq=str2num(hifrq); end
if ischar(nfft),  nfft=str2double(nfft); end
if ischar(n),     n=str2double(n); end
if ischar(stp),   stp=str2double(stp); end
if ( ~strcmpi(ftype,'b') && ( ~isscalar(lofrq) || ~isscalar(hifrq) ) )
    error('Cannot calculate spectrum for multiple frequency bands.')
end

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

% Info Display
% fprintf('\nsccohere3 %s %s %s %s %3.1f %3.1f %d %3.1f %3.1f\n\n',filtfile,wfile,chfile,ftype,lofrq,hifrq,nfft,n,stp)

% Deal with additional parameters
mtaper=false; isext=false; dosavevirt=false; useroi=false;
if nargin>9
    ii=find(strcmpi('taper',varargin));
    mtaper = ~isempty(ii);
    if mtaper
        taper = varargin{ii+1};
        if ischar(taper), taper=str2double(taper); end
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
    end
    
    ii=find(strcmpi('savevirt',varargin),1);
    dosavevirt=~isempty(ii);
    if dosavevirt
        tempfile=varargin{ii+1};
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
    end
    
    ii=find(strcmpi('baseline',varargin));
    if ~isempty(ii)
        basestart = varargin{ii+1};
        basestop  = varargin{ii+2};
        if ischar(basestart)
            basestart=str2double(basestart);
        end
        if ischar(basestop)
            basestop=str2double(basestop);
        end
        baseline=[basestart basestop].*1000; clear basestart basestop
    end

    ii=find(strcmpi('tmp',varargin),1);
    if ~isempty(ii)
        virtfile = varargin{ii+1};
    end
end
clear varargin

% Prepare
numvirt= size(W,2);
fs  = meg.srate;
lat = meg.latency;
ns  = floor(n*fs);
stps= floor(stp*fs);
nbt = floor((lentrial-ns)/stps) + 1;    % number of time points

if any(hifrq>fs/2), error('Your high frequency cutoff is greater than the Nyquist frequency.'), end
frq=[lofrq(:) hifrq(:)]; clear lofrq hifrq
numband= size(frq,1);
ffft=[0:fs/nfft:fs/2]';
for b=1:numband
    bandsel{b}=find(ffft>=frq(b,1) & ffft<=frq(b,2));
end
select = unique(vertcat(bandsel{:}));
numfrq = length(select);
if strcmpi(ftype,'b')
    fv = frq;
    for b=1:numband
        bandsel{b}=find(ismember(select,bandsel{b}));
    end
else
    fv  = ffft(select);
    fd = (fs/nfft)/2;
    fv = [fv-fd fv+fd];
    bandsel={1:length(select)};
end
clear frq ffft fd

all2all = (~isext && ncomp>40 && ncomp>numfrq*nbt);

if mtaper
    numtaper = floor( 2*taper*n - 1);
    NW = n*taper;
    if NW < 1
        error('Multi Taper: the product LEN * FREQWIDTH_IN_HZ must be >1');
    end
    wdw = dpss(ns, NW, 'calc');
    wdw = wdw(:,1:numtaper);
    wdw = permute( repmat(wdw,[1 1 numsens]), [1 3 2]);  
    fprintf('Using %d tapers.\n',numtaper)
else
    wdw=repmat(hanning(ns),[1 numsens]);
    numtaper=1;
end
if isext, wdwe = wdw(:,size(ext,2),:); end

% Sensor Fourier Transform
%-------------------------
fprintf('Calculating sensor Fourier transform: ')
try
    try
        FF = complex(zeros(numfrq,numsens,numtaper,nbt,numtrial));
        lowmem = false;
    catch ME
        if strcmp(ME.identifier,'MATLAB:nomem')
            fprintf('\nOut of memory, trying to use single precision... ')
            FF = complex(zeros(numfrq,numsens,numtaper,nbt,numtrial,'single'));
            fprintf('ok.\n')
        else
            rethrow(ME)            
        end
    end 
catch ME
    lowmem = strcmp(ME.identifier,'MATLAB:nomem');
    if lowmem
        fprintf('\nOut of memory, using temporary disk files.\n')
        if ~exist('virtfile','var')     % create random temp file name to avoid overwrites on clusters
            crs='-0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ_abcdefghijklmnopqrstuvwxyz';
            warning('off','MATLAB:RandStream:ActivatingLegacyGenerators')
            c=clock; rand('seed',ceil((10*c(end)-fix(10*c(end)))*100));
            warning('on','MATLAB:RandStream:ActivatingLegacyGenerators')
            virtfile=crs(ceil(rand(1,8)*64)); clear crs c
        end
        if all2all
            for f=1:numfrq
                fidc(f) = fopen(sprintf('%s%03d.tmp',virtfile,f),'w+');
            end
            sDat = numsens*numtaper;        
        else
            fidc = fopen([virtfile '.tmp'],'w+');
            sDat = numfrq*numsens*numtaper;
        end
    else
        rethrow(ME)
    end
end  

index = 1:ns;
tim=zeros(nbt,1);
for t=1:nbt
    for kkk=1:numtrial
        tmp = fft(wdw .* repmat(detrend(meg.data(index,:,kkk)),[1 1 numtaper]), nfft, 1); 
        if lowmem
            if all2all
                for f=1:numfrq
                    fwrite(fidc(f),real(tmp(select(f),:,:)),'double');
                    fwrite(fidc(f),imag(tmp(select(f),:,:)),'double');
                end
            else
                fwrite(fidc,real(tmp(select,:,:)),'double');
                fwrite(fidc,imag(tmp(select,:,:)),'double');
            end
        else
            FF(:,:,:,t,kkk) = tmp(select,:,:);
        end
    end
    tim(t) = mean(lat(index([1 end])));
    index = index + stps;  
end
clear meg
if isext            
    EE = complex(zeros(numfrq,size(ext,2),numtaper,nbt,numtrial,class(FF)));
    index = 1:ns;
    for t=1:nbt
        for kkk=1:numtrial
            tmp = fft(wdwe .* repmat(detrend(ext(index,:,kkk)),[1 1 numtaper]), nfft, 1); 
            EE(:,:,:,t,kkk) = tmp(select,:,:);
        end
        index = index + stps;
    end
end
clear tmp ext index

fprintf('done.\n')

% Source connectivity
% -------------------
if all2all
    % For large numbers of connections, this is much much much much .... faster
    fprintf('Calculating source complex coherency: ')
    
    if lowmem
        for f=1:numfrq
            fseek(fidc(f),0,-1);  % rewind
        end 
    end
    if strcmpi(ftype,'b')
        Coh=zeros(ncomp,nbt,numband);
    else
        Coh=zeros(ncomp,nbt,numfrq);
    end
    perc10=floor((numband*nbt)/10);
    if (perc10>2 || mtaper), numstep=10; step10=round(linspace(perc10,numband*nbt,numstep));
    else step10=numband*nbt; numstep=1;
    end, clear perc10
    perccount=0;
    
    for t=1:nbt
        
        for b=1:numband
        if strcmpi(ftype,'b')
            fXY=zeros(numvirt); fX=zeros(numvirt,1);
        end
        for f=bandsel{b}.'
        
            x = zeros(numtrial,numvirt,numtaper); 
            if lowmem
                for k=1:numtrial
                    FF = reshape(fread(fidc(f),sDat,'double'),[1 numsens numtaper]) + 1i*reshape(fread(fidc(f),sDat,'double'),[1 numsens numtaper]);
                    for ee=1:numtaper 
                        x(k,:,ee) = FF(1,:,ee) * W;
                    end
                end
            else
                for ee=1:numtaper 
                    x(:,:,ee) = permute(FF(f,:,ee,t,:),[5 2 1 3 4]) * W; 
                end
            end
            
            if mtaper
                XY=zeros(numvirt); X=zeros(numvirt,1); 
                for k=1:numtrial
                    %tXY=zeros(numvirt); tX=zeros(numvirt,1); 
                    %for ee=1:numtaper
                    tmp = permute(x(k,:,:),[3 2 1]);
                    tXY = tmp' * tmp;  
                    tmp = abs(tmp);
                    tX  = diag(tmp.'*tmp);
                    %tX = tX + tmp; %sum(abs(x).^2, 2);
                    %end
                    XY=XY+tXY./numtaper; 
                    X =X +tX./numtaper; 
                    clear tXY tX
                end
            else
                XY = x' * x;   % x' is the same as conj(x).'
                x = abs(x);
                X = diag(x.'*x);
            end

            if strcmpi(ftype,'b')
                fXY = fXY+XY;
                fX  = fX+X;
            else
                % We bring the matrix to vector form to save more than half of
                % memory
                XY = XY./sqrt(X*X.');
                for c=1:ncomp
                    Coh(c,t,f) = XY(ch(c,1),ch(c,2)); % complex coherence 
                end      
            end
        end  % f
        clear X XY
        if strcmpi(ftype,'b')
            fXY = fXY./sqrt(fX*fX.');
            for c=1:ncomp
                Coh(c,t,b) = fXY(ch(c,1),ch(c,2)); % complex coherence 
            end                 % complex coherence 
        end
        if any( step10==((t-1)*numband+b) )
            perccount=perccount+100/numstep;
            fprintf('...%d%%',perccount)
        end   
        end  % numband
    end  % nbt
    fprintf('\n')
    
    if lowmem
        for f=1:numfrq
            fclose(fidc(f));
        end
        delete([virtfile '*.tmp']); 
    end
    
else  % if only few connections need to be calculated, it's faster this way
    
    if dosavevirt
        vchn = unique(ch(:));
        nvchn= length(vchn);
        if ncomp/nvchn>2
            fprintf('Saving virtual data as temporary files: ')
            [pTMP,fTMP]=fileparts(tempfile);
            domkdir = ( ~isempty(pTMP) && ~exist(pTMP,'dir') );
            if domkdir, mkdir(pTMP), end

            perc10=floor(nvchn/10);
            if perc10>10, step10=[perc10:perc10:nvchn]; numstep=10;
            else step10=nvchn; numstep=1;
            end, clear perc10
            perccount=0;
            for c=1:nvchn
                x = zeros(numfrq,nbt,numtrial,numtaper); 
                if lowmem
                    fseek(fidc,0,-1);   % rewind
                    for t=1:nbt      
                        for kkk=1:numtrial
                            FF = reshape(fread(fidc,sDat,'double'),[numfrq numsens numtaper]) + 1i*reshape(fread(fidc,sDat,'double'),[numfrq numsens numtaper]);
                            for ee=1:numtaper 
                                x(:,t,kkk,ee) = FF(:,:,ee) * W(:,vchn(c));
                            end
                        end 
                    end         
                else
                    for t=1:nbt      
                        for kkk=1:numtrial
                            for ee=1:numtaper 
                                x(:,t,kkk,ee) = FF(:,:,ee,t,kkk) * W(:,vchn(c));
                            end
                        end 
                    end                
                end
                fidv = fopen([tempfile '_' sprintf('%03d',vchn(c)) '.tmp'],'w');
                fwrite(fidv,real(x),'double');
                fwrite(fidv,imag(x),'double');
                fclose(fidv);

                if any(step10==c)
                    perccount=perccount+100/numstep;
                    fprintf('...%d%%',perccount)
                end                
            end
            vDat=prod([numfrq nbt numtrial numtaper]);
            fprintf('\n')
        else
            dosavevirt=false;
        end
        clear vchn 
    end

    fprintf('Calculating source complex coherency: ')
    if strcmpi(ftype,'b')
        Coh=zeros(ncomp,nbt,numband);
    else
        Coh=zeros(ncomp,nbt,numfrq);
    end
    perc10=floor(ncomp/10);
    if perc10>10, step10=[perc10:perc10:ncomp]; numstep=10;
    else step10=ncomp; numstep=1;
    end, clear perc10
    perccount=0;
    for c=1:ncomp
    
        if ~dosavevirt
            x = zeros(numfrq,nbt,numtrial,numtaper); y=x;
            if lowmem
                fseek(fidc,0,-1);   % rewind
                for t=1:nbt      
                    for kkk=1:numtrial
                        FF = reshape(fread(fidc,sDat,'double'),[numfrq numsens numtaper]);
                        for ee=1:numtaper 
                            x(:,t,kkk,ee) = FF(:,:,ee) * W(:,ch(c,1));
                            if ~isext, y(:,t,kkk,ee) = FF(:,:,ee) * W(:,ch(c,2)); end
                        end
                    end 
                end
            else
                for t=1:nbt      
                    for kkk=1:numtrial
                        for ee=1:numtaper 
                            x(:,t,kkk,ee) = FF(:,:,ee,t,kkk) * W(:,ch(c,1));
                            if ~isext, y(:,t,kkk,ee) = FF(:,:,ee,t,kkk) * W(:,ch(c,2)); end
                        end
                    end 
                end
            end
        else
            fidv = fopen([tempfile '_' sprintf('%03d',ch(c,1)) '.tmp'],'r');
            realx = fread(fidv,vDat,'double');
            imagx = fread(fidv,vDat,'double');
            fclose(fidv);
            x = reshape(complex(realx,imagx),[numfrq nbt numtrial numtaper]); clear realx imagx
            if ~isext
                fidv = fopen([tempfile '_' sprintf('%03d',ch(c,2)) '.tmp'],'r');
                realy = fread(fidv,vDat,'double');
                imagy = fread(fidv,vDat,'double');
                fclose(fidv);
                y = reshape(complex(realy,imagy),[numfrq nbt numtrial numtaper]); clear realy imagy
            end        
        end      
        if isext            
            y = permute(EE(:,ch(c,2),:,:,:),[1 4 5 3 2]);  
        end

        if mtaper
            XY = sum(mean(y.*conj(x), 4), 3);             
            X = sum(mean(abs(x).^2, 4), 3);
            Y = sum(mean(abs(y).^2, 4), 3);    
        else
            XY = sum(y.*conj(x), 3);           
            X = sum(abs(x).^2, 3);
            Y = sum(abs(y).^2, 3);     
        end

        if strcmpi(ftype,'b')
            for b=1:numband
                bX = sum(X(bandsel{b},:),1);
                bY = sum(Y(bandsel{b},:),1);
                bXY = sum(XY(bandsel{b},:),1); 
                Coh(c,:,b) = bXY./sqrt(bX.*bY);             % complex coherence 
            end
        else
            Coh(c,:,:) = (XY./sqrt(X.*Y)).';         % complex coherence 
        end
        if any(step10==c)
            perccount=perccount+100/numstep;
            fprintf('...%d%%',perccount)
        end    
    end % ncomp
    fprintf('\n')
    
    if lowmem, fclose(fidc); delete(virtfile); end
    if dosavevirt
        fprintf('Deleting temporary files: ')
        delete(fullfile(pTMP,[fTMP '*.tmp']))
        if ( exist('domkdir','var') && domkdir )
            d=dir(pTMP); if length(d)<3, rmdir(pTMP); end
        end
        fprintf('done.\n')
    end
    
end % if all2all

% Output
% ------
CC.coh=Coh; clear Coh
CC.method = 'ccohere';
CC.diminfo = 'connection * time * frequency';
CC.comps=ch;
CC.frq=fv;
CC.time=[tim-1000*stp/2 tim+1000*stp/2];
if exist('baseline','var')
    CC.baseline = baseline;
end
CC.N=numtrial;
CC.len=lentrial;
if mtaper
    CC.taperwidth=taper;
    CC.numtaper=numtaper;
end
if useroi
    CC.roidef = roidef;
    CC.roifile = roifile;
end

if ~exist('ccohere','dir'), mkdir ccohere, end
save(fullfile('ccohere',['CC' chfile(2:end)]),'CC')

fprintf('Done (sccohere3).\n')
tempus=toc;
fprintf('Elapsed time is %d minutes and %d seconds.\n',round(tempus/60),round(rem(tempus,60)))
