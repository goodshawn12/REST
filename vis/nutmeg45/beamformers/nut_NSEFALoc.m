function nsefaloc_out=nut_NSEFALoc(Lp,data,flags);
% function [lqv,g]=algo1a(y,f,sig0,lam0,nem,phi,nm_x,plotflag);
% model y=F*G*phi+A*x+v

% global bolts

if flags.qsub
    voxperjob=500; %feel free to change
    pauselength=30; % seconds
    [nk,nf,nv]=size(Lp);
    numjobs=ceil(nv/voxperjob);
    [crap,whereami] = unix('id -gn');
    if(strcmp(whereami(1:3),'bil'))
        pristring = '-p -700 ';
    else
        pristring = '';
    end
    newdir=['tmpnsefaloc' num2str(flags.tb.nl) '_' num2str(flags.tb.nm_x)];    
    mkdir(newdir)
    mkdir('qout');
    cd(newdir)
    meg=data.y;
    latency=data.latency;
    save('nsefalocin.mat','meg','latency','flags');
    Lpfull=Lp;clear Lp
    for ii=1:numjobs
        if ii==numjobs
        Lp=Lpfull(:,:,(ii-1)*voxperjob+1:end);
        else
        Lp=Lpfull(:,:,(ii-1)*voxperjob+1:(ii-1)*voxperjob+voxperjob);
        end
        save(['nsefalocin' num2str(ii) '.mat'],'Lp')
%         nsefalocrun(pwd,ii)
        unix(['qsub ' pristring ' /home/johannaz/matlab/qnsefaloc.csh ' pwd ' ' num2str(ii)]);
        %         unix(['qsub ' pristring '-t 1-' num2str(numjobs) ' /data/research_meg/tfbf/bin/qtfrun.csh ' sessionfile ' ' outname ' ' algo{jj}]);
    end
    done=zeros(1,numjobs);
    while ~prod(done)
        pause(pauselength)
        for ii=1:numjobs
            done(ii)=exist(['nsefalocout' num2str(ii) '.mat']);
        end
    end
    nsefaloc_out.lmap=zeros(nv,2);
    nsefaloc_out.g=zeros(nf,flags.tb.nl,nv,2);
    nsefaloc_out.W1=zeros(nk,nf,nv,2);
    for ii=1:numjobs
        tbout{ii}=load(['nsefalocout' num2str(ii) '.mat']);
        nsefaloc_out.tbflags=tbout{1}.tb1aout.tbflags;
        nsefaloc_out.phi=tbout{1}.tb1aout.phi;
        nsefaloc_out.Ws=tbout{1}.tb1aout.Ws;
        if ii==numjobs
        nsefaloc_out.lmap((ii-1)*voxperjob+1:end,:)=tbout{ii}.tb1aout.lmap;
        nsefaloc_out.g(:,:,(ii-1)*voxperjob+1:end,:)=tbout{ii}.tb1aout.g;
        nsefaloc_out.W1(:,:,(ii-1)*voxperjob+1:end,:)=tbout{ii}.tb1aout.W1;
        else
        nsefaloc_out.lmap((ii-1)*voxperjob+1:(ii-1)*voxperjob+voxperjob,:)=tbout{ii}.tb1aout.lmap;
        nsefaloc_out.g(:,:,(ii-1)*voxperjob+1:(ii-1)*voxperjob+voxperjob,:)=tbout{ii}.tb1aout.g;
        nsefaloc_out.W1(:,:,(ii-1)*voxperjob+1:(ii-1)*voxperjob+voxperjob,:)=tbout{ii}.tb1aout.W1;
        end
    end
    !rm nsefaloc*.mat
    cd ..
    rmdir(newdir);
else
    nsefaloc_out=nut_nsefaloc_sub(Lp,data,flags);
end
