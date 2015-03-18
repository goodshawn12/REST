function [hec,shec]=nut_sourceConnectivity_HilbertEnvelope(inputfile,seednode,chunktime,schunktime,windows)
% function [hec,shec]=nut_sourceConnectivity_HilbertEnvelope(inputfile,seednode,chunktime,schunktime,windows)
%
% inputfile either weights or s_beam_* file
%           you'd use s_beam_* file if number of voxels*timepts is small
%           otherwise use weights file (for memory reasons)
% seednode: voxel index for this subject (call this function in a loop over
%           seed voxels)
% suggested: chunktime=3000; in ms
% suggested: schunktime=500; in ms
% windows is size conditions x 3; for each condition (active/control, or
%           cond1, cond2, etc,:  (1) start and (2)end points in
%           milliseconds of time window, plus (3) which s{?} structure to use
%           (set to 1 if inputfile is weights)

beam=load(inputfile);
if isfield(beam,'coreg')
    beam=nut_beam_legacy_compatibility(beam);
    if ~isfield(beam,'s')
        error('Not a valid s_beam file.')
    end
    if size(beam.s{1},2)<100
        error('probably not enough time points in this s_beam_* file');
    end
    inputtype=1;
elseif isfield(beam,'W')
    W=beam.W;
    clear beam
    megfile=[inputfile(9:(strfind(inputfile,'Hz')+1)) '.mat'];
    load(megfile);
    inputtype=2;
    beam.timepts=meg.latency;
else
    error('wrong input file. enter either weights or s_beam_* file')
end

fs=round(1000/(beam.timepts(2)-beam.timepts(1)));
chunk=chunktime*fs/1000;
schunk=schunktime*fs/1000;


for ii=1:size(windows,1)
    numchunk(ii)=floor((windows(ii,2)-windows(ii,1))/chunktime);
    numschunk(ii)=chunktime*numchunk(ii)/schunktime;
    timetot(ii)=numchunk(ii)*chunk;
    indwind=dsearchn(beam.timepts,windows(ii,1:2)');
    if timetot(ii)~=diff(indwind)
        indwind(2)=indwind(1)+timetot(ii)-1;
    end
    %
    if inputtype==1
        ahs=permute(nut_abshilbert(permute(beam.s{windows(ii,3)}(:,indwind(1):(indwind(2)-1),:),[2 1 3])),[2 1 3]);
        hectmp=zeros(size(ahs,1),numchunk,size(ahs,3));
        mahs=zeros(size(ahs,1),numschunk,size(ahs,3));
        shectmp=zeros(size(ahs,1),size(ahs,3));
        hec=zeros(size(ahs,1),size(windows,1));
        shec=zeros(size(ahs,1),size(windows,1));
        for kk=1:size(ahs,3)
            for nn=1:numchunk
                hectmp(:,nn,kk)=corr(ahs(seednode,((nn-1)*chunk+1):nn*chunk,kk),ahs(:,((nn-1)*chunk+1):nn*chunk,kk));
            end
            for nn=1:numschunk
                mahs(:,nn,kk)=mean(ahs(:,((nn-1)*schunk+1):(nn*schunk),kk),2);
            end
            shectmp(:,kk)=corr(mahs(seednode,:,kk)',mahs(:,:,kk)');
        end
        hec(:,ii)=mean(reshape(hectmp,[size(hectmp,1) size(hectmp,2)*size(hectmp,3)]),2);
        clear hectmp
        shec(:,ii)=mean(shectmp,2);
        clear shectmp
        %
    elseif inputtype==2
        seeds=reshape(W(:,seednode)'*reshape(permute(meg.data(indwind(1):(indwind(2)-1),:,:),[2 1 3]),size(meg.data,2),timetot(ii)*size(meg.data,3)),[timetot(ii) 1 size(meg.data,3)]);
        seedahs=nut_abshilbert(seeds);
        for kk=1:size(meg.data,3)
            s=W'*meg.data(indwind(1):(indwind(2)-1),:,kk)';
            for ll=1:size(W,2)
                ahs=nut_abshilbert(s(ll,:));
                for nn=1:numschunk
                    mahs(ll,nn,kk)=mean(ahs(((nn-1)*schunk+1):nn*schunk));
                end
                for nn=1:numchunk
                    hectmp(ll,nn,kk)=corr(ahs(((nn-1)*chunk+1):nn*chunk)',seedahs(((nn-1)*chunk+1):nn*chunk,1,kk));
                end
            end
            clear ahs
            shectmp(:,kk)=corr(mahs(seednode,:,kk)',mahs(:,:,kk)');
        end
        clear s mahs
        hec(:,ii)=mean(reshape(hectmp,[size(hectmp,1) size(hectmp,2)*size(hectmp,3)]),2);
        clear hectmp
        shec(:,ii)=mean(shectmp,2);
        clear shectmp
    end
end




