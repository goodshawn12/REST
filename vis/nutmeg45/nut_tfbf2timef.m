function nut_tfbf2timef(sessionfile,algo,controlperiod)
% NUT_TFBF2TIMEF aggregates individual TFBF outputs into s_beamtf composite
%
%   nut_tfbf2timef(sessionfile_filt_beamparam,algo)
%
% sessionfile_filt_beamparam      string containing the filenames of the
%               session file, the filter parameter file, and the beamformer 
%               parameter file
% algo          string indicating the beamformer algorithm (optional,
%               default is 'SAM')

if nargin<2, algo='SAM'; end

[sessionpath,sessionname,ext]=fileparts(sessionfile);


% matfiles = dir(['s_beamtf_' sessionname '_*to*_*to*ms_' algo '.mat']);
if exist('controlperiod')
    matfiles = dir(['s_beamtf_' sessionname '_*to*Hz_' controlperiod '_*to*ms_' algo '.mat']);
else
    matfiles = dir(['s_beamtf_' sessionname '_*to*Hz_*to*ms_' algo '.mat']);
end% [SAMimage,coords,SAMparams] = nut_read_svl(svlfiles(1).name);
% temp=load(matfiles(1).name);

for ii=1:size(matfiles,1)
    if exist('controlperiod')
        %     TFinfo(ii,:)=sscanf(matfiles(ii).name,['s_beamtf_' sessionname '_%fto%f_' controlperiod '_%dto%dms_' algo '.mat'])';
        TFinfo(ii,:)=sscanf(matfiles(ii).name,['s_beamtf_' sessionname '_%fto%fHz_' controlperiod '_%dto%dms_' algo '.mat'])';
    else
        %     TFinfo(ii,:)=sscanf(matfiles(ii).name,['s_beamtf_' sessionname '_%fto%f_%dto%dms_' algo '.mat'])';
        TFinfo(ii,:)=sscanf(matfiles(ii).name,['s_beamtf_' sessionname '_%fto%fHz_%dto%dms_' algo '.mat'])';
    end
end

timewins = [unique(TFinfo(:,3:4),'rows')];
bands = [unique(TFinfo(:,1:2),'rows')];

if(size(timewins,1)>1)
    freqbin=1; %this is just using the first frequency band to find timestep
    if exist('controlperiod')
        tmpfiles = dir(['s_beamtf_' sessionname '_' num2str(bands(freqbin,1)) 'to' num2str(bands(freqbin,2)) 'Hz_' controlperiod '_*to*ms_' algo '.mat']);
    else
        tmpfiles = dir(['s_beamtf_' sessionname '_' num2str(bands(freqbin,1)) 'to' num2str(bands(freqbin,2)) 'Hz*to*ms_' algo '.mat']);
    end
    clear TFinfo
    
    for ii=1:size(tmpfiles,1)
        if exist('controlperiod')
            %             TFinfo(ii,:)=sscanf(matfiles(ii).name,['s_beamtf_' sessionname '_%fto%f_' controlperiod '_%dto%dms_' algo '.mat'])';
            TFinfo(ii,:)=sscanf(tmpfiles(ii).name,['s_beamtf_' sessionname '_%fto%fHz_' controlperiod '_%dto%dms_' algo '.mat'])';
        else
            %             TFinfo(ii,:)=sscanf(matfiles(ii).name,['s_beamtf_' sessionname '_%fto%f_%dto%dms_' algo '.mat'])';
            TFinfo(ii,:)=sscanf(tmpfiles(ii).name,['s_beamtf_' sessionname '_%fto%fHz_%dto%dms_' algo '.mat'])';
        end
    end

    timewins = [unique(TFinfo(:,3:4),'rows')];
    timestep = timewins(2,1)-timewins(1,1);
    
    if(timewins(1,2) == timewins(2,1)) % no overlap
        beam.timewindow = timewins;
    else
        % the following is a check for consistency of 2 different ways of computing beam.timewindow
        tw1=[(timewins(1,2):timestep:timewins(end-1,1))' (timewins(1,2):timestep:timewins(end-1,1))'+timestep];
        tw2=[((timewins(1,1)+timestep):timestep:(timewins(end,2)-2*timestep))' ((timewins(1,1)+2*timestep):timestep:(timewins(end,2)-timestep))'];
        if size(tw1)==size(tw2)
            if tw1==tw2
                % this will be true if you arrange time windows across frequency bands the way that sarang and adrian originally did
            end
        else
            disp('warning: check if the final timesteps plotted in results_viewer are what you expected? (jz ignore)');
        end
        beam.timewindow=tw2; %we think this should work for everyone, however may or may be equal to tw1 depending what you've done
    end
else
    beam.timewindow = timewins;
end

timepts = mean(beam.timewindow,2);
beam.timepts = timepts;
if(size(timepts,1)==1)
    beam.srate = 1;
else
%     beam.timewindow = timewins;
    beam.srate = 1000/timestep;
end
beam.bands = bands;


% beam.tf = zeros(length(SAMimage),size(timepts,1),size(bands,1));

% load in first matfile just to set dimensions on Sact, Scon, beam.tf
% temp=load(matfiles(1).name,'beam');
temp.beam=load(matfiles(1).name);
if isfield(temp.beam,'beam')
    temp.beam=temp.beam.beam;
end
% temp=load(matfiles(1).name,'beam','Sact','Scon');
% beam.tf = zeros(size(temp.beam.tf,1),size(timepts,1),size(bands,1));
% beam.Sact = beam.tf;
% beam.Scon = beam.tf;
% beam.s{1} = zeros(size(temp.beam.tf,1),size(timepts,1),size(bands,1));
% beam.s{2} = beam.s{1};
beam.s{1} = zeros(size(temp.beam.s{1}));
beam.s{2} = zeros(size(temp.beam.s{2}));
if(length(temp.beam.s) == 3)
    noiseflag = true;
    beam.s{3} = zeros(size(temp.beam.s{3}));
else
    noiseflag = false;
end
if isfield(temp.beam,'z')
    zflag=1;
    if size(temp.beam.z,1)==1
        beam.z{1}=zeros(size(temp.beam.z{1}));
    else
        beam.z{1}=zeros(size(temp.beam.z));
    end
else
    zflag=0;
end

if(isfield(temp.beam,'p'))
    pflag = 1;
    if size(temp.beam.z,1)==1
        beam.p{1}=zeros(size(temp.beam.p{1}));
    else
        beam.p{1}=zeros(size(temp.beam.p));
    end
else
    pflag = 0;
end


if isfield(temp.beam,'sasd')
    varflag=1;
    beam.sasd{1}=zeros(size(temp.beam.sasd{1}));
    beam.scsd{1}=zeros(size(temp.beam.scsd{1}));
else
    varflag=0;
end

for freqbin=1:size(bands,1)
    freqbin

    if exist('controlperiod')
%     matfiles = dir(['s_beamtf_' sessionname '_' num2str(bands(freqbin,1)) 'to' num2str(bands(freqbin,2)) '_' controlperiod '*to*ms_' algo '.mat']);
    matfiles = dir(['s_beamtf_' sessionname '_' num2str(bands(freqbin,1)) 'to' num2str(bands(freqbin,2)) 'Hz_' controlperiod '*to*ms_' algo '.mat']);
    else
%     matfiles = dir(['s_beamtf_' sessionname '_' num2str(bands(freqbin,1)) 'to' num2str(bands(freqbin,2)) '*to*ms_' algo '.mat']);
    matfiles = dir(['s_beamtf_' sessionname '_' num2str(bands(freqbin,1)) 'to' num2str(bands(freqbin,2)) 'Hz*to*ms_' algo '.mat']);
    end
    clear TFinfo
    
    for ii=1:size(matfiles,1)
        if exist('controlperiod')
%             TFinfo(ii,:)=sscanf(matfiles(ii).name,['s_beamtf_' sessionname '_%fto%f_' controlperiod '_%dto%dms_' algo '.mat'])';
            TFinfo(ii,:)=sscanf(matfiles(ii).name,['s_beamtf_' sessionname '_%fto%fHz_' controlperiod '_%dto%dms_' algo '.mat'])';
        else
%             TFinfo(ii,:)=sscanf(matfiles(ii).name,['s_beamtf_' sessionname '_%fto%f_%dto%dms_' algo '.mat'])';
            TFinfo(ii,:)=sscanf(matfiles(ii).name,['s_beamtf_' sessionname '_%fto%fHz_%dto%dms_' algo '.mat'])';
        end
    end
    
    
    [timewins,timeidx] = sortrows([TFinfo(:,3) TFinfo(:,4)]);
%     tf = zeros(length(SAMimage),size(timepts,1));

%     Sact = zeros(size(beam.Sact,1),size(timewins,1));
    Sact = zeros(size(beam.s{1},1),size(timewins,1));
    Scon = Sact;
    if(noiseflag)
        noise = Sact;
    end
    if(zflag)
        zz = Sact;
    end
    if(pflag)
        pp = Sact;
    end
    if(varflag)
        sasd = Sact;
        scsd = Scon;
    end
    
    for timebin=1:size(timewins,1)
        filename = matfiles(timeidx(timebin)).name;
%         temp=load(filename,'beam','Sact','Scon');
%         temp=load(filename,'beam');
        temp.beam=load(filename);
        if isfield(temp.beam,'beam')
            temp.beam=temp.beam.beam;
        end
%         if ~isfield(temp,'beam')
%             temp.beam=load(filename);
%         end
        Sact(:,timebin) = temp.beam.s{1};
        Scon(:,timebin) = temp.beam.s{2};
        if(noiseflag)
            noise(:,timebin) = temp.beam.s{3};
        end
        if(zflag)
            zz(:,timebin) = temp.beam.z{1};
        end
        if(pflag)
            pp(:,timebin) = temp.beam.p{1};
        end
        if(varflag)
            sasd(:,timebin) = temp.beam.sasd{1};
            scsd(:,timebin) = temp.beam.scsd{1};
        end
    end

    for timebin=1:size(timepts,1)
        timeselect = find((timepts(timebin) > timewins(:,1)) & (timepts(timebin) < timewins(:,2)));
%         beam.Sact(:,timebin,freqbin) = mean(Sact(:,timeselect),2);
%         beam.Scon(:,timebin,freqbin) = mean(Scon(:,timeselect),2);
        beam.s{1}(:,timebin,freqbin) = mean(Sact(:,timeselect),2);
        beam.s{2}(:,timebin,freqbin) = mean(Scon(:,timeselect),2);
        if(noiseflag)
            beam.s{3}(:,timebin,freqbin) = mean(noise(:,timeselect),2);
        end
        if(zflag)
            beam.z{1}(:,timebin,freqbin) = mean(zz(:,timeselect),2);
            %is this legit to do?
        end
        if(pflag)
            beam.p{1}(:,timebin,freqbin) = mean(pp(:,timeselect),2);
            %is this legit to do?
        end
        if(varflag)
            beam.sasd{1}(:,timebin,freqbin) = mean(sasd(:,timeselect),2);
            %is this legit to do?
            beam.scsd{1}(:,timebin,freqbin) = mean(scsd(:,timeselect),2);
            %is this legit to do?
        end
    end
    
   
%     for timebin=1:size(timepts,1)
%         timeselect = find((timepts(timebin) > timewins(:,1)) & (timepts(timebin) < timewins(:,2)));
%         beam.tf(:,timebin,freqbin) = mean(tf(:,timeselect),2);
%     end
end

if(0) % let's just do this in nut_timef_viewer
    switch(units)
        case 'dB'
            beam.tf = 10*log10(abs(beam.Sact./beam.Scon));
        case 'F'
            beam.tf = beam.Sact./beam.Scon;
        case 't'
            beam.tf = beam.Sact - beam.Scon;
        case 'PW'
            % Scon -> ERD, Sact -> ERS
            beam.tf = beam.Sact;
            select = find(beam.Scon>beam.Sact);
            beam.tf(select)=-beam.Scon(select);
        case 'raw'
            beam.tf = beam.Sact;
    end
end


% beam.s_timef(:,1,1) = SAMimage;

%beam.voxels = coords*1000; % convert from m to mm
%beam.voxelsize = repmat(SAMparams.stepsize*1000,[1 3]);
beam.voxels = temp.beam.voxels;
beam.voxelsize = temp.beam.voxelsize;
beam.params = temp.beam.params;
% beam.params.units = units;
beam.coreg = temp.beam.coreg;

if exist('controlperiod')
    outputfile = fullfile(sessionpath,['s_beamtf_' sessionname '_' algo '_' controlperiod '_all.mat']);
else
    outputfile = fullfile(sessionpath,['s_beamtf_' sessionname '_' algo '_all.mat']);
end
save(outputfile,'-struct','beam','-V7.3');


function F=weirdFtostandardF(Fweird)
% converts CTF's weird F value to a standard F ratio
F = zeros(size(Fweird));
F(Fweird>=0) = Fweird(Fweird>=0)+1;
F(Fweird<0) = 1./(1-Fweird(Fweird<0));
