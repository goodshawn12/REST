function [hec,shec,coh,icoh]=nut_sourceConnectivity(conset);
% function nut_sourceConnectivity(conset)

if ~isfield(conset,'W')
    load(conset.inputfile);
else
    W=conset.W;
    conset=rmfield(conset,'W');
end
if ~exist('W','var')
    error('wrong input file. enter weights file')
end
freq=conset.filterband;

[datapath,weightsfile,ext]=fileparts(conset.inputfile);

sessionfile=[datapath '/' weightsfile(9:(strfind(weightsfile,'firls')-2)) '.mat'];
filtclass='firls';
if ~exist(sessionfile,'file')
    sessionfile=[weightsfile(9:(strfind(weightsfile,'butter')-2)) '.mat'];
    filtclass='butter';
end
if ~exist(sessionfile,'file')
    error('what filter type did you use for creating weights?')
end
coreg=load(sessionfile,'coreg');
try
    size(coreg{1});
    voxels=load(sessionfile,'voxels');
catch
    load(sessionfile,'nuts');
    coreg=nuts.coreg;
    meg=nuts.meg;
    voxels=nuts.voxels;
    clear nuts
end

if conset.filter==1
    megfile=[weightsfile(9:(strfind(weightsfile,'Hz')+1)) '.mat'];
    try
        load([datapath '/' megfile]);
    catch
        conset.filter=0;
    end
end
if conset.filter==0
    if ~exist('meg','var')
        meg=load(sessionfile,'meg');
    end
    filt.order=200;
    if size(meg.data,1)<600
        filt.order=round(.33*size(meg.data,1));
    end
    meg.data = nut_filter2(meg.data,filtclass,'bp',filt.order,freq(1),freq(2),meg.srate,1);
end

if conset.trials==0
    trialuse=1:size(meg.data,3);
else
    trialuse=conset.trials;
end
if ~isnumeric(trialuse)
    error('wrong specification of trials to use.  either vector or `all`')
end
meg.data=meg.data(:,:,trialuse);

beam.timepts=meg.latency;
fs=meg.srate;

if length(meg.latency)<100
    error('probably not enough time points in this dataset');
end

if conset.seed.type==1
    seednode=conset.seed.index;
elseif conset.seed.type==2
    seednode=dsearchn(voxels,nut_mri2meg(nut_mni2mri(conset.seed.index,coreg)));
end


chunk=round(conset.chunktime*fs/1000);
schunk=round(conset.schunktime*fs/1000);
windows=conset.windows;

for ii=1:size(windows,1)
    numchunk(ii)=floor((windows(ii,2)-windows(ii,1))/conset.chunktime);
    numschunk(ii)=conset.chunktime*numchunk(ii)/conset.schunktime;
    timetot(ii)=numchunk(ii)*chunk;
    indwind=dsearchn(beam.timepts,windows(ii,1:2)');
    if timetot(ii)~=diff(indwind)
        indwind(2)=indwind(1)+timetot(ii)-1;
    end
end

%% Hilbert Envelope Correlation
if conset.do.hec || conset.do.shec
    for ii=1:size(windows,1)
        seeds=reshape(W(:,seednode)'*reshape(permute(meg.data(indwind(1):(indwind(2)-1),:,:),[2 1 3]),size(meg.data,2),timetot(ii)*size(meg.data,3)),[timetot(ii) 1 size(meg.data,3)]);
        seedahs=nut_abshilbert(seeds);
        for kk=1:size(meg.data,3)
            s=W'*meg.data(indwind(1):(indwind(2)-1),:,kk)';
            for ll=1:size(W,2)
                ahs=nut_abshilbert(s(ll,:));
                if conset.do.shec
                    for nn=1:numschunk
                        mahs(ll,nn,kk)=mean(ahs(((nn-1)*schunk+1):nn*schunk));
                    end
                end
                if conset.do.hec
                    for nn=1:numchunk
                        hectmp(ll,nn,kk)=corr(ahs(((nn-1)*chunk+1):nn*chunk)',seedahs(((nn-1)*chunk+1):nn*chunk,1,kk));
                    end
                end
            end
            clear ahs
            if conset.do.shec
                shectmp(:,kk)=corr(mahs(seednode,:,kk)',mahs(:,:,kk)');
            end
        end
        clear s mahs
        if conset.do.hec
            hec(:,ii)=mean(reshape(hectmp,[size(hectmp,1) size(hectmp,2)*size(hectmp,3)]),2);
            clear hectmp
        end
        if conset.do.shec
            shec(:,ii)=mean(shectmp,2);
            clear shectmp
        end
    end
end
if ~exist('hec','var')
    hec=[];
end
if ~exist('shec','var')
    shec=[];
end

%% Coherence
if conset.do.coh || conset.do.icoh
    for ii=1:size(windows,1)
        linf=linspace(0,fs,chunk);
        hanning_wdw=hanning(fs*conset.chunktime/1000);
        frequse=[dsearchn(linf',freq(1)):dsearchn(linf',freq(2))];
        cohtmp=zeros(size(W,2),size(numchunk(ii)));
        icohtmp=zeros(size(W,2),size(numchunk(ii)));
        
        seeds=1e15*reshape(W(:,seednode)'*reshape(permute(meg.data(indwind(1):(indwind(2)-1),:,:),[2 1 3]),size(meg.data,2),timetot(ii)*size(meg.data,3)),[timetot(ii) 1 size(meg.data,3)]);
        for ll=1:size(W,2)
            s=1e15*reshape(W(:,ll)'*reshape(permute(meg.data(indwind(1):(indwind(2)-1),:,:),[2 1 3]),size(meg.data,2),timetot(ii)*size(meg.data,3)),[timetot(ii) 1 size(meg.data,3)]);
            for nn=1:numchunk(ii)
                seedsf=fft(repmat(hanning_wdw,1,size(meg.data,3)).*detrend(squeeze(seeds((1+(nn-1)*chunk):(1+nn*chunk-1),:,:))),[],1);
                sf=fft(repmat(hanning_wdw,1,size(meg.data,3)).*detrend(squeeze(s((1+(nn-1)*chunk):(1+nn*chunk-1),:,:))),[],1);
                [tmp1,tmp2]=nut_coherence([permute(seedsf,[3 1 2]); permute(sf,[3 1 2])],frequse);
                cohtmp(ll,nn)=tmp1(1,2);
                icohtmp(ll,nn)=tmp2(1,2);
            end
        end
        coh(:,ii)=mean(cohtmp,2);
        icoh(:,ii)=mean(icohtmp,2);
    end
end

if ~exist('coh','var')
    coh=[];
end
if ~exist('icoh','var')
    icoh=[];
end

%% saving out results
beam=load(['s_beamtf_' weightsfile(9:end)]);
if isfield(beam,'beam')
    beam=beam.beam;
end

beam.timepts=1;
beam.timewindow=[0.5 1.5];
beam.srate=1;
beam.bands=conset.freq;
if conset.do.hec
    for ii=1:size(hec,2)
        beam.s{ii}=hec(:,ii);
    end
    save(['s_beam_hec_' weightsfile(9:end)],'-struct','beam');
end
if conset.do.shec
    for ii=1:size(shec,2)
        beam.s{ii}=shec(:,ii);
    end
    save(['s_beam_shec_' weightsfile(9:end)],'-struct','beam');
end
if conset.do.coh
    for ii=1:size(coh,2)
        beam.s{ii}=coh(:,ii);
    end
    save(['s_beam_coh_' weightsfile(9:end)],'-struct','beam');
end
if conset.do.icoh
    for ii=1:size(icoh,2)
        beam.s{ii}=icoh(:,ii);
    end
    save(['s_beam_icoh_' weightsfile(9:end)],'-struct','beam');
end






