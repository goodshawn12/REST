function [data,sim,L,comp]=nut_sources2sensors(s,t,L,sim);
% function [data,sim,L,comp]=nut_sources2sensors(s,t,L,sim);
%
% INPUTS:
% s = source times courses (can be output from nut_create_sim_sources)
% t = time vector for sources (can be output from nut_create_sim_sources)
% L (input) = lead field precomputed
% sim.nk = number sensors for inputted lead field
% sim.nt = number of time points
% sim.numsens = desired number of sensors to output data onto
% sim.dist = distance contstraint (broken right now)
% sim.iv0 = voxel indices of locations of sources, from a nuts.voxels that the lead field corresponds to.  Can be left empty, and random locations will be chosen.
% sim.numdip = number of source locations
% sim.numint = number of interference locations
% sim.addmean = if want to add mean to poststim y1
% sim.nv = number of voxels in nuts.voxels corresponding to inputted lead field
% sim.nf = number of orientations of lead field per voxel location (normally 1, 2, or 3)
% sim.etacorr = correlation between directions/orientations of source (=1 if fixed dipole, 0 if x,y,z orientations are totally independent)
% sim.alp = correlation of first two sources.  if more than 2 sources, then additional sources not forced to be correlated or not.
% sim.inter = if 'real, then sim.realdata needs to exist (you bring in from elsewhere), or 'inter', 'sin', or 'gauss' will project those time courses onto sensor data from random locations
% sim.snr = in dB, signal to noise ratio, of just sensor noise
% sim.snir = in db, signal to noise-plus-interference ratio, including real brain noise or simulated inteference 
% sim.seed (and all the substructures within) = integer to get the same random number every time
% example:  sim.seed.sourceloc=13;sim.seed.eta=17;sim.seed.sourcegauss=19;sim.seed.noise=23;sim.seed.interts=29;sim.seed.interloc=31;sim.seed.interorient=37;sim.seed.sourcefreqphase=41;sim.seed.interfreqphase=43;sim.seed.sourcewindow=47;
%
% OUTPUTS:
% L (output) = lead field with possibly reduced number of sensors
% data = sensor data as you created it!
% sim = new structures added or modified from what you input to keep a record of what you did
% sim.iv1 = actual voxel locations used for source locations
% sim.sources = final record of source time course, after figuring fixed/rotation and correlated source parameters
% sim.infv = random location of interference sources (if simulated)

if sim.numsens<sim.nk
    origsens=randperm(sim.nk);
    sim.sensorder=origsens(1:sim.numsens);
    L=L(sim.sensorder,:,:);
    sim.realdata=sim.realdata(sim.sensorder,:);
    sim.nk=sim.numsens;
    sim.nm_sefa=ceil(sim.numsens/2);
else
    sim.sensorder=1:sim.nk;
end

if sim.seed.sourceloc
    rand('state',sim.seed.sourceloc);
else
    rand('state',sum(100*clock));
end

if (length(sim.dist)==2) & (length(sim.iv0)==0) & (sim.numdip==2)
    error('this option broken right now.  can"t specify distance constraint')
    sim.iv1(1)=round(rand(1)*(nv-1)+1);
    dip.voxelsMRI=nut_coordtfm(dip.voxels,dip.coreg.meg2mri_tfm);
    distances = nut_rownorm(nut_coord_diff(dip.voxelsMRI,dip.voxelsMRI(sim.iv1(1),:)));
    ind=find(distances<=sim.dist(2) & distances>=sim.dist(1));
    sim.iv1(2)=ind(round(rand*(length(ind)-1)+1));
    % elseif isempty(sim.iv0)
    %     for ii=1:sim.numdip
    %         sim.iv1(ii)=round(rand(1)*(sim.nv-1)+1);
    %     end
elseif length(sim.iv0)<sim.numdip
    sim.iv1=sim.iv0;
    while length(sim.iv1)<sim.numdip
        sim.iv1(length(sim.iv1)+1)=round(rand(1)*(sim.nv-1)+1);
    end
    if length(unique(sim.iv1))<length(sim.iv1)
        error('haven"t fixed this yet. trying to put two sources at the same place')
    end
elseif length(sim.iv0)>sim.numdip
    sim.iv1=sim.iv0(1:sim.numdip)
    disp('not using all specified voxels')
elseif length(sim.iv0)==sim.numdip
    sim.iv1=sim.iv0;
    %     else
    %         sim.iv1=sim.iv0;
    %         while length(sim.iv1)<sim.numdip
    %             sim.iv1(length(sim.iv1)+1)=round(rand(1)*(sim.nv-1)+1);
    %         end
    %     end
end

if sim.seed.eta
    rand('state',sim.seed.eta);
else
    rand('state',sum(100*clock));
end
eta0=2*rand(sim.numdip,sim.nf)-1;
for ii=1:sim.numdip
    sim.eta0(ii,:)=eta0(ii,:)/nut_rownorm(eta0(ii,:));
    %         aeta0(:,ii)=squeeze(L(:,1,sim.iv1(ii)))*sim.eta0(ii,1)+squeeze(L(:,2,sim.iv1(ii)))*sim.eta0(ii,2);
end
t0=dsearchn(t',0);
postnt=t(t0:end);
prent=t(1:t0-1);

% create source data
%     if sim.sken
sim.source=s(1:sim.nf*sim.numdip,:,:);
%     else
%         if sim.sourcegauss_seed
%             randn('state',sim.sourcegauss_seed);
%         else
%             randn('state',sum(100*clock));
%         end
%         x0=[zeros(sim.nf*sim.numdip,length(prent)) randn(sim.nf*sim.numdip,length(postnt))];
%         sim.source=x0;
%     end
if ceil(sim.etacorr)
    etacorrb=sqrt(1-sim.etacorr^2);
    if sim.nf>1
        sim.source(2:sim.nf:sim.nf*sim.numdip,:,:)=sim.etacorr*sim.source(1:sim.nf:sim.nf*sim.numdip,:,:)+etacorrb*sim.source(2:sim.nf:sim.nf*sim.numdip,:,:);
    end
    if sim.nf>2
        sim.source(sim.nf:sim.nf:sim.nf*sim.numdip,:,:)=sim.etacorr*sim.source(1:sim.nf:sim.nf*sim.numdip,:,:)+etacorrb*sim.source(2:sim.nf:sim.nf*sim.numdip,:,:);
    end
end
if sim.numdip>1
    beta=sqrt(1-sim.alp^2); %correlation of sources
    %         sim.source(3,:)=sim.alp*sim.source(1,:)+beta*sim.source(3,:);
    %         sim.source(4,:)=sim.alp*sim.source(2,:)+beta*sim.source(4,:);
    for kk=1:sim.nf
        sim.source(sim.nf+kk,:,:)=sim.alp*sim.source(kk,:,:)+beta*sim.source(sim.nf+kk,:,:);
    end
end
y1=zeros(sim.nk,sim.nt,sim.numtrial);
% fixme.  eta just scaling, not rotating here.
for jj=1:size(sim.source,3)
    for ii=1:sim.numdip
        sim.etas(sim.nf*ii-sim.nf+1:sim.nf*ii,:,jj)=nut_rowmult_jz(sim.eta0(ii,:)',sim.source(sim.nf*ii-sim.nf+1:sim.nf*ii,:,jj));
        y1(:,:,jj)=y1(:,:,jj)+squeeze(L(:,:,sim.iv1(ii)))*sim.etas(sim.nf*ii-sim.nf+1:sim.nf*ii,:,jj);
    end
end
y1(:,t0:end,:)=y1(:,t0:end,:)+sim.addmean;

% create interference and noise
if sim.seed.noise
    randn('state',sim.seed.noise);
else
    randn('state',sum(100*clock));
end
n1=randn(size(y1));
u1=zeros(size(y1));
%     if sim.interts_seed
%         randn('state',sim.interts_seed);
%     else
%         randn('state',sum(100*clock));
%     end
%     u0=randn(sim.numint,nt);
%
if sim.seed.interloc
    rand('state',sim.seed.interloc);
else
    rand('state',sum(100*clock));
end
if sim.seed.interorient
    randn('state',sim.seed.interorient);
else
    randn('state',sum(100*clock));
end

switch sim.inter
    case 'none'  % don't need to do anything further here
        clear u0
        %         case 'gauss'
        %             for ii=1:sim.numint
        %                 sim.infv(ii)=round(rand(1)*(sim.nv-1)+1);
        %                 etai=randn(sim.nf,1);
        %                 %                 aetai=squeeze(L(:,1,sim.infv(ii)))*etai(1)+squeeze(L(:,2,sim.infv(ii)))*etai(2);
        %                 aetai=squeeze(L(:,:,sim.infv(ii)))*etai;
        %                 u1=u1+aetai*u0(ii,:);
        %             end
    case {'sin','gauss','laplacian'}
        clear u0
        for jj=1:sim.numtrial
            for ii=1:sim.numint
                sim.infv(ii)=round(rand(1)*(sim.nv-1)+1);
                etai=randn(sim.nf,1);
                %                 aetai=squeeze(L(:,1,sim.infv(ii)))*etai(1)+squeeze(L(:,2,sim.infv(ii)))*etai(2);
                aetai=squeeze(L(:,:,sim.infv(ii)))*etai;
                u1(:,:,jj)=u1(:,:,jj)+aetai*s(ii+sim.nf*sim.numdip,:,jj);
            end
        end
    case 'real'
        clear n1 u1 u0
        i1=sim.realdata;
    otherwise
        errordlg('you''re not interfering now')
end

% get correct SNR and SNIR
switch sim.inter
    case {'none','gauss','sin'}
        SNR = 10*log10(sum(sum(mean(y1(:,t0:end,:),3).^2,2))/sum(sum(mean(n1(:,t0:end,:),3).^2,2)));
        ngain=sqrt(10^((sim.snr-SNR)/10));
        n1=n1/ngain;

        switch sim.inter
            case 'none'
                i1=n1+u1;
                y=y1+i1;
                SNR = 10*log10(sum(sum(mean(y1(:,t0:end,:),3).^2,2))/sum(sum(mean(n1(:,t0:end,:),3).^2,2)))
            case {'gauss','sin'}
                i1=u1+n1;
                aquad=sum(sum(mean(u1(:,t0:end,:),3).^2,2));
                bquad=2*sum(sum(mean(u1(:,t0:end,:),3).*mean(n1(:,t0:end,:),3),2));
                cquad=sum(sum(mean(n1(:,t0:end,:),3).^2,2))-sum(sum(mean(y1(:,t0:end,:),3).^2,2))/10^(sim.snir/10);
                snir_quadp=(-bquad+sqrt(bquad^2-4*aquad*cquad))/(2*aquad);
                if ~isreal(snir_quadp)
                    errordlg('you chose the SNR too low or the SNIR too high')
                end
                u1=snir_quadp*u1;
                i1=u1+n1;
                y=y1+i1;
                SNR = 10*log10(sum(sum(mean(y1(:,t0:end,:),3).^2,2))/sum(sum(mean(n1(:,t0:end,:),3).^2,2)))
                SIR = 10*log10(sum(sum(mean(y1(:,t0:end,:),3).^2,2))/sum(sum(mean(u1(:,t0:end,:),3).^2,2)))
                SNIR = 10*log10(sum(sum(mean(y1(:,t0:end,:),3).^2,2))/sum(sum(mean(i1(:,t0:end,:),3).^2,2)))
            otherwise
                errordlg('how on earth...');
        end

    case 'real'
        SNIR = 10*log10(sum(sum(mean(y1(:,t0:end,:),3).^2,2))/sum(sum(mean(i1(:,t0:end,:),3).^2,2)));
        snirgain=sqrt(10^((sim.snir-SNIR)/10));
        i1=i1/snirgain;
        SNIR = 10*log10(sum(sum(mean(y1(:,t0:end,:),3).^2,2))/sum(sum(mean(i1(:,t0:end,:),3).^2,2)))
        y=y1+i1;
    otherwise
        errordlg('you''re messing up the SNR')
end
data=y;
comp.y1=y1;
comp.i1=i1;
switch sim.inter
    case {'none','gauss','sin'}
        comp.u1=u1;
        comp.n1=n1;
end

