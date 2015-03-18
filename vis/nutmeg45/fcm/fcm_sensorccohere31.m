function SC = fcm_sensorccohere31(t,sig,ch,ftype,frq,nfft,taper)
%FCM_SENSORCCOHERE31  Complex coherency of 3D sensor data. 
%Each trial is a segment of a continuous signal.
%
%   SC = FCM_SENSORCCOHERE31(FS,SIG,CONNECTIONS,FTYPE,FREQ,¦NFFT¦,¦TAPER¦) 
%   ( parameters in ¦¦ are optional )
%
%   SC          output structure containing the complex valued coherence
%               data.
%   FS          time vector or sampling rate
%   SIG         must be a 3-dimensional signal (time*channel*trial).
%   CONNECTIONS must be a Nx2 matrix indicating the indices of N channel pairs
%               to be analyzed, or the string 'all' for all possible channel pairs.
%   FTYPE       's'=frequency [s]pectrum, 'b'=frequency [b]and.
%   FREQ        If FREQ is a 2-element vector, it indicates the upper and lower limit
%               of the frequency range to be analysed. If FREQ contains a scalar, all 
%               frequencies up to FREQ Hz are analysed. 
%   NFFT        indicates the number of bins used for FFT. Default is the next
%               power of 2 of the trial length.
%   TAPER       For MULTI TAPERING: freq width in Hz.
%
%   (c) 2004-2010  Adrian G. Guggisberg


if nargin <5, help fcm_sensorccohere31, Coh=[]; return, end
[len,numsens,numtrial]=size(sig);
if nargin < 6 || isempty(nfft), nfft=2^nextpow2(len); end
if nargin < 7, taper=[]; end

if isscalar(t)
    Fs=t;
    t=[0 ((len-1)/Fs)*1000];
else
    Fs=round(1/(t(2)-t(1)));
end

if ischar(ch) 
    if strcmpi(ch,'all')
        ncomp=(numsens*(numsens-1))/2;
        ch=zeros(ncomp,2,'uint32');
        currstart=1;
        for cc=1:numsens
            currlength=length(cc+1:numsens);
            currend=currstart+currlength-1;
            ch(currstart:currend,:)=[cc*ones(currlength,1) [cc+1:numsens]'];
            currstart=currend+1;
        end
        clear curr* 
    else
        error('CONNECTIONS must be Nx2 numeric matrix or the string ''all''.')
    end
else
    ncomp=size(ch,1);
    if ( max(ch(:)) > size(sig,2) )
        error('Channel numbers in input parameter CONNECTIONS do not exist.');
    end
end

if length(frq)==2, frq=[max(frq(1),1/(len/Fs)) min(frq(2),Fs/2)];
elseif length(frq)<2, frq=[1/(len/Fs) min(frq,Fs/2)];
else, error('FREQ must be a 1 or 2 element vector.')
end

mtaper = ~isempty(taper);
if mtaper
    taper = [len/Fs taper];
    numtaper = floor( 2*taper(2)*taper(1) - 1);
    NW = taper(1)*taper(2);
    N  = taper(1)*Fs;     
    wdw = dpss(N, NW, 'calc');
    wdw = wdw(:,1:numtaper);
    wdw = permute( repmat(wdw,[1 1 numsens]), [1 3 2]);  
    fprintf('Using %d tapers.\n',numtaper)
else
    wdw=repmat(hanning(len),[1 numsens]);
    numtaper=1;
end

% Main part
ffft=0:Fs/nfft:Fs/2;
select = find(ffft>=frq(1) & ffft<=frq(2));
numfrq = length(select);
%select = [1+floor(frq(1)*nfft/Fs):1+ceil(frq(2)*nfft/Fs)];
%if any(any(imag([x y])~=0)),                            % if x and y are complex
%    select = [select fliplr(1+nfft-select)];
%end

fprintf('Calculating sensor Fourier transform: ')
FF = zeros(numfrq,numsens,numtaper,numtrial);
for k=1:numtrial
    tmp = detrend(sig(:,:,k));  
    tmp = fft(wdw .* repmat(tmp,[1 1 numtaper]), nfft, 1); 
    FF(:,:,:,k) = tmp(select,:,:);
end
clear tmp 
fprintf('done.\n')

fprintf('\nCalculating Coherence: ')

if strcmpi(ftype,'b')
    fXY=zeros(numsens); fX=zeros(numsens,1); 
else
    Coh=zeros(ncomp,length(select));
end

if ncomp>numfrq
        
    perc10=floor(numfrq/10);
    if perc10>=5, numstep=10; step10=round(linspace(perc10,numfrq,numstep));
    else numstep=1; step10=numfrq; 
    end, clear perc10
    perccount=0;
    for f=1:numfrq  
        if mtaper
            x = permute(FF(f,:,:,:),[3 2 4 1]);   
            XY=zeros(numsens); X=zeros(numsens,1); 
            for k=1:numtrial
                tmp = x(:,:,k);
                tXY = tmp' * tmp;  
                tmp = abs(tmp);
                tX  = diag(tmp.'*tmp);
                XY=XY+tXY./numtaper; 
                X =X +tX./numtaper; 
                clear tXY tX 
            end
        else
            x = permute(FF(f,:,:,:),[4 2 3 1]);   
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
            XY = XY./sqrt(X*X.');
            for c=1:ncomp
                Coh(c,f) = XY(ch(c,1),ch(c,2)); % complex coherence 
            end      
        end
        if any( step10==f )
            perccount=perccount+100/numstep;
            fprintf('...%d%%',perccount)
        end   
    end % for f=1:numfrq
    clear X XY
    if strcmpi(ftype,'b')
        fXY = fXY./sqrt(fX*fX.');
        Coh = zeros(ncomp,1);
        for c=1:ncomp
            Coh(c) = fXY(ch(c,1),ch(c,2)); % complex coherence 
        end               
    end
    fprintf('\n')
    
else

    perc10=floor(ncomp/10);
    if perc10>9, step10=[perc10:perc10:ncomp]; numstep=10;
    else step10=ncomp; numstep=1;
    end, clear perc10
    perccount=0;
    for c=1:ncomp
        x = permute(FF(:,ch(c,1),:,:),[1 4 3 2]);
        y = permute(FF(:,ch(c,2),:,:),[1 4 3 2]);

        if mtaper
            Xy = mean(y.*conj(x), 3);
            Xx = mean(abs(x).^2, 3);
            Yy = mean(abs(y).^2, 3);        
        else
            Xy = y.*conj(x);
            Xx = abs(x).^2;
            Yy = abs(y).^2;
        end
        Pxx = sum(Xx(select,:),2);
        Pyy = sum(Yy(select,:),2);
        Pxy = sum(Xy(select,:),2);  
        %end
        if strcmpi(ftype,'b')
            Pxx = sum(Pxx);
            Pyy = sum(Pyy);
            Pxy = sum(Pxy); 
            Coh(c) = Pxy/sqrt(Pxx*Pyy);             % coherence function estimate
        else
            Coh(c,:) = Pxy./sqrt(Pxx.*Pyy);        % coherence function estimate 
        end
        if any(step10==c)
            perccount=perccount+100/numstep;
            fprintf('...%d%%',perccount)
        end       
    end
    fprintf('\n')
end
    
% freq vectors
if strcmpi(ftype,'b'), ff = frq;
else,   ff  = ffft(select);
end

SC.coh=Coh;
SC.method='ccohere';
SC.frq=ff;
SC.comps=ch;
SC.time=[t(1) t(end)];
SC.N = numtrial;
SC.len = len;
if mtaper
    SC.taperwidth=taper(2);
    SC.numtaper = numtaper;
end
    
    
