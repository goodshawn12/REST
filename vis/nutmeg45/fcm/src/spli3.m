function spli3(filtfile,wfile,chfile,n,stp,varargin)
%SPLI3  Calculates source phase lag index for event-locked multi-trial data.
%
%   spli3(filtfile,weightfile,connections,len,step,options)
% 
% Input:
%  filtfile     file containing filtered sensor data (filtdata_*.mat).
%  weightfile   file containing weight matrix from inverse solution (W_*.mat)
%  connections  can be a Nx2 matrix containing the indices of voxel pairs
%               to be analyzed (created with fcm_conndef) or the name of a text 
%               file containing these indices (created with fcm_preparejobs).
%  len          the length of one data segment in sec. 
%  step         the time step between data segments in sec. 
%  options      parameter/value combinations to specify additional properties:
%               For calculating complex coherence with extracerebral channels
%               (e.g., EMG for cortico-muscular coherence)
%                   spli3(...,'ext',extdata), with extdata being the 
%                   event-locked multi-trial extracerebral data (samples X
%                   channels X trials) or a corresponding matlab file containing 
%                   variable "ext".
%               You can specify the method to calculate the PLI:
%                   spli31(...,'bandmethod','hilbert') is the official but very slow way
%                   spli31(...,'bandmethod','fasthilbert') also uses a Hilbert transform,  
%                       but then estimates the phase difference in each data segment rather 
%                       than at each sample (default).
%               To re-bandpass filter the data before running connectivities.
%                   spli3(...,'b',lofrq,hifrq)
%               You can specify a parameter file for the bandpass filter
%                   spli3(...,'filt',filtconfigfile)
%                   filtconfigfile can be created with nut_tfbf_gui.
%                   Default is 'ellipauto.mat'.
%               To obtain additional baseline corrected outputs:
%                   spli3(...,'baseline',basestart,basestop), where
%                   basestart and basestop are seconds relative to the
%                   beginning of each trial.
%               To use independent components:
%                   spli3(...,'ica',icadata), where icadata contains
%                   variable IC, the non-artifact independent components and 
%                   variable A, the mixing matrix.
%               To calculate FC between anatomical ROIs:
%                   spli3(...,'roi',roifile), where roifile contains a structure 
%                   R which can be created with fcm_voxel2roi.
%
% Output:
%  results are saved in a directory "pli".
%
% References:
%   - Stam et al. Hum Brain Mapp 2007; 28: 1178–1193


if nargin<5
    fprintf('\nusage:  spli3 filtfile weightfile connectionfile length_timewin step [options]\n\n')
    return
end

if ischar(n)
    n=str2double(n);
end
if ischar(stp)
    stp=str2double(stp);
end

tic;

% Open all files
if ischar(chfile)
    [chpath,chfile,chext]=fileparts(chfile);
    if isempty(chext), chext='.txt'; end
    ch = load(fullfile(chpath,[chfile chext]));
else    % if comps given as matrix
    ch=chfile;
    chfile='comps';
end

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

% Info Display
fprintf('\nspli3 %s %s %s %d %d\n\n',filtfile,wfile,chfile,n,stp)

% Deal with variable inputs.
isext=false; useroi=false; ds=1; bandmethod='fasthilbert'; %defaults
if nargin>5
    ii=find(strcmpi('ext',varargin));
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
    
    ii=find(strcmpi('filt',varargin));
    if ~isempty(ii)
        load(varargin{ii+1});   % must contain filt structure with filter settings
        if ~exist('filt','var'), error('Invalid filter configuration file.'), end
    end

    if any(strcmp('s',varargin))
        error('Frquency spectrum is not implemented for this configuration. Does it really make sense?')
    end

    ii=find(strcmp('b',varargin),1);
    if ~isempty(ii)
        if length(varargin)>=ii+2
            lofrq=varargin{ii+1}; hifrq=varargin{ii+2};
            if ~isempty(lofrq) && ( lofrq-frq(1)>1 || hifrq-frq(2)<-1 )        
                if ~exist('filt','var')  
                    load ellipauto      % default filter parameters
                end
                fprintf('Filtering sensor data: ')
                meg.data = nut_filter2(meg.data,filt.class,filt.type,filt.order,lofrq,hifrq,meg.srate,1);
                if isext
                    ext = nut_filter2(ext,filt.class,filt.type,filt.order,lofrq,hifrq,meg.srate,1);
                end
                frq = [lofrq hifrq];
                fprintf('done.\n')    
            end
        end
    end
    
    ii=find(strcmpi('bandmethod',varargin),1);
    if ~isempty(ii) 
        bandmethod=lower(varargin{ii+1}); 
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

end

% Prepare
fs = meg.srate;
lat= meg.latency;
ns = floor(n*fs);
stp= floor(stp*fs);
nbt= floor((lentrial-ns)/stp) + 1;    % number of time points
ncomp=size(ch,1);
numvirt=size(W,2);
all2all = (~isext && ncomp>numtrial*nbt);

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


% Main Part

fprintf('Calculating sensor Fourier transform: ')
% Note that projecting the Fourier transform instead of the data itself
% speeds up a lot the algorithm, because only half of the data (the positive
% Fourier coefficients) has to be projected. We win about 30-40% of time!

try
    FF=complex(zeros(numfrq,numsens,nbt,numtrial),zeros(numfrq,numsens,nbt,numtrial));
catch ME
    if strcmp(ME.identifier,'MATLAB:nomem');
        fprintf('\nOut of memory, trying to use single precision... ')
        FF=complex(zeros(numfrq,numsens,nbt,numtrial,'single'),zeros(numfrq,numsens,nbt,numtrial,'single'));
        fprintf('ok.\n')
    else
        rethrow(ME)
    end
end

tim = zeros(nbt,2);
for k=1:numtrial
    meg.data(:,:,k) = detrend(meg.data(:,:,k));  
    index = 1:ns;
    for tt=1:nbt
        tmp = fft(meg.data(index,:,k), ns, 1); 
        FF(:,:,tt,k) = tmp(select,:);          
        if k==1, tim(tt,:) = lat(index([1 end])).'; end
        index = index + stp; index=index(index<=lentrial);
    end
end
clear meg
if isext            
    EE = complex(zeros(numfrq,size(ext,2),numtrial,class(FF)),zeros(numfrq,size(ext,2),numtrial,class(FF)));
    for k=1:numtrial
        ext(:,:,k) = detrend(ext(:,:,k));  
        index = 1:ns;
        for tt=1:nbt
            tmp = fft(ext(index,:,k), ns, 1); 
            EE(:,:,tt,k) = tmp(select,:);          
            index = index + stp; index=index(index<=lentrial);
        end
    end
end
clear tmp ext
fprintf('done.\n')

% Try to speedup some more
% We can remove the filtered frequencies for now and then add them back as 0 later.
D=mean(mean(abs(diff(log10(abs(FF(:,1:3:end,1:3:end))+1e-12))),3),2);     % find bins where spectrogram approaches 0
f=(D>0.02); clear D
if (numfrq-sum(f)>10) && (sum(abs(diff(f)))<3)  % only possible if filtered data present
    select=[1;find(f)+1];  clear f
    fprintf('Speedup on: working only on %d of %d freq bins.\n',length(select),lentrial)    
    FF = FF(select,:,:,:);
    if isext, EE = EE(select,:,:,:); end
end

fprintf('Calculating source PLI: ')
Coh = zeros(ncomp,nbt);
if all2all
    perc10=floor(nbt/10);
    if perc10>=1, numstep=10; step10=round(linspace(perc10,nbt,numstep));
    else step10=nbt; numstep=1;
    end, clear perc10
    perccount=0;    

    for tt=1:nbt
        XY = zeros(numvirt);
        for k=1:numtrial
            x = zeros(numfrq,numvirt); 
            x(select,:) = FF(:,:,tt,k) * W;
            
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
            end
%             XY2 = zeros(numvirt);
%             for t=1:ns
%                 X = repmat(x(t,:),[numvirt 1]);
%                 XY2 = XY2 + sign(sin(X-X.'));
%             end       
%             XY2 = XY2 ./ ns; % mean  
%             XY = XY + XY2;
        end
        
        switch bandmethod
            case 'hilbert'
                XY = XY ./ (numtrial*length(1:ds:ns));
            otherwise
                XY =XY ./ numtrial;
        end
            
        for c=1:ncomp
            Coh(c,tt) = XY(ch(c,1),ch(c,2)); 
        end
        
        if any( step10==tt )
            perccount=perccount+100/numstep;
            fprintf('...%d%%',perccount)
        end
            
    end
              
    clear x XY
    
else

    perc10=floor(ncomp/10);
    if perc10>10, step10=[perc10:perc10:ncomp]; numstep=10;
    else step10=ncomp; numstep=1;
    end, clear perc10
    perccount=0;
    for c=1:ncomp

        x = zeros(numfrq,numtrial); y=x;
        for tt=1:nbt
            for k=1:numtrial      
                x(select,k) = FF(:,:,tt,k) * W(:,ch(c,1));
                if ~isext, y(select,k) = FF(:,:,tt,k) * W(:,ch(c,2)); end
            end
            if isext            
                y(select,:) = permute(EE(:,ch(c,2),tt,:),[1 4 2 3]);  
            end

            % this is the same as X=mean(angle(hilbert(data)))
            x = angle(ifft(cat(1,h.*x,neg0))); 
            y = angle(ifft(cat(1,h.*y,neg0))); 
            Coh(c,tt) = mean(mean(sign(sin(x-y)),1),2);

            if any(step10==c)
                perccount=perccount+100/numstep;
                fprintf('...%d%%',perccount)
            end
        end
    end % for 1:ncomp
end % if all2all

fprintf('\n')

% Output
CC.coh=Coh; 
CC.method='pli';
CC.comps=ch;
CC.frq=frq;
CC.time=tim;
if exist('baseline','var')
    CC.baseline = baseline;
end
CC.N=numtrial;
CC.len=lentrial;
if useroi
    CC.roidef = roidef;
    CC.roifile = roifile;
end

if ~exist('pli','dir'), mkdir pli, end
save(fullfile('pli',['CC' chfile(2:end)]),'CC')

fprintf('Done (spli3).\n')
tempus=toc;
fprintf('Elapsed time is %d minutes and %d seconds.\n',round(tempus/60),round(rem(tempus,60)))
