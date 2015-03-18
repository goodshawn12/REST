function sakeout=nut_Saketini(Lp,data,flags);
%
% dipole source localization for the model y = F*s + A*x + B*u + v
% s~N(0,phi) , x~N(0,I), v~N(0,lam)

% global bolts

if flags.qsub
    voxperjob=500;
    pauselength=30; % seconds
    [nk,nf,nv]=size(Lp);
    nt=size(data.y,1);
%     ntpost=size(bolts.meg(dsearchn(bolts.latency,0):end,:),1);
    numjobs=ceil(nv/voxperjob);
    [crap,whereami] = unix('id -gn');
    if(strcmp(whereami(1:3),'bil'))
        pristring = '-p -700 ';
    else
        pristring = '';
    end
    newdir=['tmpsake' num2str(flags.sake.nl) '_' num2str(flags.sake.nm)];    
    mkdir(newdir)
    mkdir('qout');
    cd(newdir)
    meg=data.y;
    latency=data.latency;
    if ~isfield(flags.sake,'eigs')
%         flags.sake.eigs=max(bolts.params.signalspace);
        error('need to specifiy flags.sake.eigs');
    end
    save('sakein.mat','meg','latency','flags');
    Lpfull=Lp;clear Lp
    for ii=1:numjobs
        if ii==numjobs
        Lp=Lpfull(:,:,(ii-1)*voxperjob+1:end);
        else
        Lp=Lpfull(:,:,(ii-1)*voxperjob+1:(ii-1)*voxperjob+voxperjob);
        end
        save(['sakein' num2str(ii) '.mat'],'Lp')
%         sakerun(pwd,ii)
        unix(['qsub ' pristring ' /home/johannaz/matlab/qsake.csh ' pwd ' ' num2str(ii)]);
        %         unix(['qsub ' pristring '-t 1-' num2str(numjobs) ' /data/research_meg/tfbf/bin/qtfrun.csh ' sessionfile ' ' outname ' ' algo{jj}]);
    end
    done=zeros(1,numjobs);
    while ~prod(done)
        pause(pauselength)
        for ii=1:numjobs
            done(ii)=exist(['sakeout' num2str(ii) '.mat']);
        end
    end
    sakeout.pmap=zeros(nf,nf,nv);
    sakeout.pmap0=zeros(nf,nf,nv);
    sakeout.lmap0=zeros(nv,1);
    sakeout.lmap=zeros(nv,1);
    sakeout.lratmap=zeros(nv,1);
%     sakeout.sbar0=zeros(nv,ntpost,nf);
%     sakeout.sbar_out=zeros(nf,ntpost,nv);
    sakeout.pmapt0=zeros(nv,1);
    sakeout.pmapt=zeros(nv,1);
    sakeout.W1=zeros(nk,nf,nv);
    for ii=1:numjobs
        saout{ii}=load(['sakeout' num2str(ii) '.mat']);
        sakeout.sakeflags=saout{1}.sakeout.sakeflags;
        if ii==numjobs
            sakeout.lmap((ii-1)*voxperjob+1:end)=saout{ii}.sakeout.lmap;
            sakeout.lmap0((ii-1)*voxperjob+1:end)=saout{ii}.sakeout.lmap0;
            sakeout.lratmap((ii-1)*voxperjob+1:end)=saout{ii}.sakeout.lratmap;
            sakeout.pmap(:,:,(ii-1)*voxperjob+1:end,:)=saout{ii}.sakeout.pmap;
            sakeout.pmap0(:,:,(ii-1)*voxperjob+1:end,:)=saout{ii}.sakeout.pmap0;
%             sakeout.sbar0((ii-1)*voxperjob+1:end,:,:)=saout{ii}.sakeout.sbar0;
%             sakeout.sbar_out(:,:,(ii-1)*voxperjob+1:end)=saout{ii}.sakeout.sbar_out;
            sakeout.pmapt((ii-1)*voxperjob+1:end)=saout{ii}.sakeout.pmapt;
            sakeout.pmapt0((ii-1)*voxperjob+1:end)=saout{ii}.sakeout.pmapt0;
            sakeout.W1(:,:,(ii-1)*voxperjob+1:end)=saout{ii}.sakeout.W1;
        else
            sakeout.lmap((ii-1)*voxperjob+1:(ii-1)*voxperjob+voxperjob)=saout{ii}.sakeout.lmap;
            sakeout.lmap0((ii-1)*voxperjob+1:(ii-1)*voxperjob+voxperjob)=saout{ii}.sakeout.lmap0;
            sakeout.lratmap((ii-1)*voxperjob+1:(ii-1)*voxperjob+voxperjob)=saout{ii}.sakeout.lratmap;
            sakeout.pmap(:,:,(ii-1)*voxperjob+1:(ii-1)*voxperjob+voxperjob,:)=saout{ii}.sakeout.pmap;
            sakeout.pmap0(:,:,(ii-1)*voxperjob+1:(ii-1)*voxperjob+voxperjob,:)=saout{ii}.sakeout.pmap0;
%             sakeout.sbar0((ii-1)*voxperjob+1:(ii-1)*voxperjob+voxperjob,:,:)=saout{ii}.sakeout.sbar0;
%             sakeout.sbar_out(:,:,(ii-1)*voxperjob+1:(ii-1)*voxperjob+voxperjob)=saout{ii}.sakeout.sbar_out;
            sakeout.pmapt((ii-1)*voxperjob+1:(ii-1)*voxperjob+voxperjob)=saout{ii}.sakeout.pmapt;
            sakeout.pmapt0((ii-1)*voxperjob+1:(ii-1)*voxperjob+voxperjob)=saout{ii}.sakeout.pmapt0;
            sakeout.W1(:,:,(ii-1)*voxperjob+1:(ii-1)*voxperjob+voxperjob)=saout{ii}.sakeout.W1;
        end
    end
    !rm sake*.mat
    cd ..
    rmdir(newdir);
else
    sakeout=nut_saketini_sub(Lp,data,flags);
end
