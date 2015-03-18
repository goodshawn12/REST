
function sglcohere31(filtfile,wfile,chfile,ftype,lofrq,hifrq,nfft,varargin)
%SGLCOHERE31  Calculates source general lagged coherence for non-event-locked multi-trial data  
%             in which each trial is an artifact-free segment of a contiuous signal.
%
%  sglcohere31(filtfile,weightfile,connections,ftype,lofrq,hifrq,nfft,¦options¦)
%                         (variables in ¦¦ are optional)
% 
% Input:
%  filtfile     file containing filtered sensor data (filtdata_*.mat).
%  weightfile   file containing weight matrix from inverse solution (W_*.mat)
%  connections  can be a Nx2 matrix containing the indices of voxel pairs
%               to be analyzed (created with fcm_conndef) or the name of a text 
%               file containing these indices (created with fcm_preparejobs).
%  ftype        's'=frequency [s]pectrum, 'b'=frequency [b]and.
%  lofrq        lower bound of frequencies to be analyzed (in Hz).
%  hifrq        upper bound of frequencies to be analyzed (in Hz).
%  nfft         indicates the number of frequency bins created by FFT.
%               Should be a power of 2.
%  options      parameter/value pairs to specify additional properties:
%               For multi-tapering with Slepian sequences: 
%                   sglcohere31(...,'taper',freqwidth_in_Hz)
%               For calculating complex coherence with extracerebral channels
%               (e.g., EMG for cortico-muscular coherence)
%                   sglcohere31(...,'ext',extdata), with extdata being the non-
%                   event-locked multi-trial extracerebral data (samples X
%                   channels X trials) or a corresponding matlab file containing 
%                   variable "ext".
%               To calculate FC between anatomical ROIs:
%                   sglcohere31(...,'roi',roifile), where roifile contains a structure 
%                   R which can be created with fcm_voxel2roi.
%               To use independent components:
%                   sglcohere31(...,'ica',icadata), where icadata contains
%                   variable IC, the non-artifact independent components and 
%                   variable A, the mixing matrix.
%
% Output:
%  results are saved in a directory "glcohere".
%
% Reference:
%  - Pascual-Marqui. Coherence and phase synchronization: generalization to pairs of multivariate time series,
%    and removal of zero-lag contributions. arXiv:0706.1776v3 [stat.ME] 12 July 2007. (http://arxiv.org/pdf/0706.1776)
%  - Pascual-Marqui: Instantaneous and lagged measurements of linear and nonlinear dependence between groups of
%    multivariate time series: frequency decomposition. arXiv:0711.1455 [stat.ME], 2007-November-09, http://arxiv.org/abs/0711.1455

if nargin<7
    fprintf('\nusage: sglcohere31 filtfile wfile chfile ftype lofrq hifrq nfft [options]\n\n')
    return
end

tic;

if ischar(lofrq), lofrq=str2double(lofrq); end
if ischar(hifrq), hifrq=str2double(hifrq); end
if nargin<7 || isempty(nfft)
    nfft=2^nextpow2(fs);
elseif ischar(nfft)
    nfft=str2double(nfft);
end
if ( ~isscalar(lofrq) || ~isscalar(hifrq) )
    error('Multiple frequency bands are not supported for this datatype.')
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
[lentrial,numsens,numtrial] = size(meg.data);  
[dum,dum,eW]=fileparts(wfile); clear dum
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

% Info Display
fprintf('\nsglcohere31 %s %s %s %s %3.1f %3.1f %d\n\n',filtfile,wfile,chfile,ftype,lofrq,hifrq,nfft)

% Deal with variable inputs.
mtaper=false; isext=false; dosavevirt=false; useroi=false;
if nargin>7
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
    
    ii=find(strcmpi('savevirt',varargin),1);
    dosavevirt=~isempty(ii);
    if dosavevirt
        tempfile=varargin{ii+1};
    end      
    
    ii=find(strcmpi('tmp',varargin),1);
    if ~isempty(ii)
        virtfile = varargin{ii+1};
    end    
end
clear varargin

% Prepare
fs = meg.srate;
tim = [meg.latency(1) meg.latency(end)];

frq=[lofrq min(hifrq,fs/2)];    % check input freq limits
ffft=[0:fs/nfft:fs/2]';
select = find(ffft>=frq(1) & ffft<=frq(2));
numfrq = length(select);
numvirt= size(W,2);
if strcmpi(ftype,'b'), fv = [ffft(select(1)) ffft(select(end))];
else   fv  = ffft(select);
end
fd = (fs/nfft)/2;

if mtaper
    taper = [lentrial/fs taper];
    numtaper = floor( 2*taper(2)*taper(1) - 1);
    NW = taper(1)*taper(2);
    N  = taper(1)*fs;     
    wdw = dpss(N, NW, 'calc');
    wdw = wdw(:,1:numtaper);
    wdw = permute( repmat(wdw,[1 1 numsens]), [1 3 2]);  
    fprintf('Using %d tapers.\n',numtaper)
else
    wdw=repmat(hanning(lentrial),[1 numsens]);
    numtaper=1;
end
if isext, wdwe = wdw(:,size(ext,2),:); end

% Sensor Fourier Transform
%-------------------------
fprintf('Calculating sensor Fourier transform: ')

FF = zeros(numfrq,numsens,numtaper,numtrial);

for k=1:numtrial
    meg.data(:,:,k) = detrend(meg.data(:,:,k));  
    tmp = fft(wdw .* repmat(meg.data(:,:,k),[1 1 numtaper]), nfft, 1); 
    tmp = tmp(select,:,:);
    for ee=1:numtaper   % Normalize to remove amplitudes
        tmp2=tmp(:,:,ee);
        tmp(:,:,ee) = tmp2 * (sqrt(diag(diag(tmp2' * tmp2))));
    end
    FF(:,:,:,k) = tmp;
end
clear meg
if isext            
    EE = zeros(numfrq,size(ext,2),numtaper,numtrial);
    for k=1:numtrial
        ext(:,:,k) = detrend(ext(:,:,k));
        tmp = fft(wdwe .* repmat(ext(:,:,k),[1 1 numtaper]), nfft, 1); 
        tmp = tmp(select,:,:);
        for ee=1:numtaper
            tmp2=tmp(:,:,ee);
            tmp(:,:,ee) = tmp2 * (sqrt(diag(diag(tmp2' * tmp2))));
        end
        EE(:,:,:,k) = tmp;
    end
end
clear tmp tmp2 ext
fprintf('done.\n')

% Source connectivity
% -------------------
if ~isext && (ncomp>numfrq)
    % For large numbers of connections, this is much much much much .... faster

    fprintf('Calculating source general lagged coherence: ')
    if strcmpi(ftype,'b')
        fXY=zeros(numvirt); fX=zeros(numvirt,1); 
    else
        Coh=zeros(ncomp,1,numfrq);
    end
    perc10=floor(numfrq/10);
    if (perc10>=3 || mtaper), numstep=10; step10=round(linspace(perc10,numfrq,numstep));
    else step10=numfrq; numstep=1;
    end, clear perc10
    perccount=0;
    for f=1:numfrq
        x = zeros(numtrial,numvirt,numtaper); 
        for ee=1:numtaper 
            x(:,:,ee) = permute(FF(f,:,ee,:),[4 2 1 3]) * W;
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
            XY = x' * x;         % x' = conj(x).' 
            x = abs(x);
            X = diag(x.'*x);
        end

        if strcmpi(ftype,'b')
            fXY = fXY+XY; 
            fX = fX+X;
        else
            % We bring the matrix to vector form to save more than half of
            % the memory
            %XY = XY./sqrt(X*X.');
            XY = imag(XY).^2 ./ ( X*X.' - (real(XY).^2) );  % eq 25 in Pascual-Marqui 2007
            for c=1:ncomp
                Coh(c,1,f) = XY(ch(c,1),ch(c,2)); % complex coherence 
            end      
        end
        if any( step10==f )
            perccount=perccount+100/numstep;
            fprintf('...%d%%',perccount)
        end   
    end % for f=1:numfrq
    clear X XY
    if strcmpi(ftype,'b')
        %fXY = fXY./sqrt(fX*fX.');
        fXY = imag(fXY).^2 ./ ( fX*fX.' - (real(fXY).^2) );  % eq 25 in Pascual-Marqui 2007
        Coh = zeros(ncomp,1);
        for c=1:ncomp
            Coh(c) = fXY(ch(c,1),ch(c,2)); % complex coherence 
        end               
    end
    fprintf('\n')
    
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
                x = zeros(numfrq,numtrial,numtaper); 
                for k=1:numtrial      
                    for ee=1:numtaper 
                        x(:,k,ee) = FF(:,:,ee,k) * W(:,vchn(c));
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

            vDat=prod([numfrq numtrial numtaper]);
            fprintf('\n')
        else
            dosavevirt=false;
        end
        clear vchn 
    end
    
    fprintf('Calculating source general lagged coherence: ')
    if strcmpi(ftype,'b')
        Coh=zeros(ncomp,1);
    else
        Coh=zeros(ncomp,numfrq);
    end
    perc10=floor(ncomp/10);
    if perc10>10, step10=[perc10:perc10:ncomp]; numstep=10;
    else step10=ncomp; numstep=1;
    end, clear perc10
    perccount=0;
    for c=1:ncomp

        if ~dosavevirt
            x = zeros(numfrq,numtrial,numtaper); y=x;
            for k=1:numtrial      
                for ee=1:numtaper 
                    x(:,k,ee) = FF(:,:,ee,k) * W(:,ch(c,1));
                    if ~isext, y(:,k,ee) = FF(:,:,ee,k) * W(:,ch(c,2)); end
                end 
            end
        else
            fidv = fopen([tempfile '_' sprintf('%03d',ch(c,1)) '.tmp'],'r');
            realx = fread(fidv,vDat,'double');
            imagx = fread(fidv,vDat,'double');
            fclose(fidv);
            x = reshape(complex(realx,imagx),[numfrq numtrial numtaper]); clear realx imagx
            if ~isext
                fidv = fopen([tempfile '_' sprintf('%03d',ch(c,2)) '.tmp'],'r');
                realy = fread(fidv,vDat,'double');
                imagy = fread(fidv,vDat,'double');
                fclose(fidv);
                y = reshape(complex(realy,imagy),[numfrq numtrial numtaper]); clear realy imagy
            end        
        end
        if isext            
            y = permute(EE(:,ch(c,2),:,:),[1 4 3 2]);  
        end

        if mtaper
            XY = sum(mean(y.*conj(x), 3), 2);             
            X = sum(mean(abs(x).^2, 3), 2);
            Y = sum(mean(abs(y).^2, 3), 2);    
        else
            XY = sum(y.*conj(x), 2);           
            X = sum(abs(x).^2, 2);
            Y = sum(abs(y).^2, 2);     
        end

        if strcmpi(ftype,'b')
            X = sum(X);
            Y = sum(Y);
            XY = sum(XY); 
            Coh(c) = imag(XY).^2 ./ ( X.*Y - (real(XY).^2) );  % eq 25 in Pascual-Marqui 2007   
        else
            Coh(c,:) = imag(XY).^2 ./ ( X.*Y - (real(XY).^2) );  % eq 25 in Pascual-Marqui 2007   
        end
        if any(step10==c)
            perccount=perccount+100/numstep;
            fprintf('...%d%%',perccount)
        end    
    end
    if ~strcmpi(ftype,'b')
        Coh=reshape(Coh,[ncomp 1 numfrq]);
    end
    fprintf('\n')

    if dosavevirt
        fprintf('Deleting temporary files: ')
        delete(fullfile(pTMP,[fTMP '*.tmp']))
        if ( exist('domkdir','var') && domkdir )
            d=dir(pTMP); if length(d)<3, rmdir(pTMP); end
        end
        fprintf('done.\n')
    end
end  % if all2all

% Output
CC.coh=Coh;
CC.method = 'glcohere';
CC.diminfo = 'connection * time * frequency';
CC.comps=ch;
CC.frq=[fv(:,1)-fd fv(:,end)+fd];
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

if ~exist('glcohere','dir'), mkdir glcohere, end
save(fullfile('glcohere',['CC' chfile(2:end)]),'CC')

fprintf('Done (sglcohere31).\n')
tempus=toc;
fprintf('Elapsed time is %d minutes and %d seconds.\n',round(tempus/60),round(rem(tempus,60)))
