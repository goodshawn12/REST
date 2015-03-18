function nut_svl2timef(SAMdir)
% nut_svl2timef([time1start time1end; time2start time2end;...],
%               [band1start band1end; band2start band2end;...],
%               {{svl11,svl12,...},{svl21,svl22,...},...})

% converts SAM svl volume into s_beam

cd(SAMdir);

svlfiles = dir('*.svl');
[SAMimage,coords,SAMparams] = nut_read_svl(svlfiles(1).name);

for ii=1:size(svlfiles,1)
    SAMinfo(ii,:)=sscanf(svlfiles(ii).name,'%d,%dto%dms,%d-%dHz,D3.svl')';
end

bands = [unique(SAMinfo(:,4:5),'rows')];
timewins = [unique(SAMinfo(:,2:3),'rows')];

if(size(timewins,1)>1)
    timestep = timewins(2,1)-timewins(1,1);
    beam.timewindow = [(timewins(1,2):timestep:timewins(end-1,1))' (timewins(1,2):timestep:timewins(end-1,1))'+timestep];
else
    timestep = 100; % otherwise, arbitrary
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


beam.tf = zeros(length(SAMimage),size(timepts,1),size(bands,1));

for freqbin=1:size(bands,1)
    svlfiles = dir(['*' num2str(bands(freqbin,1)) '-' num2str(bands(freqbin,2)) 'Hz*.svl']);
    clear SAMinfo
    for ii=1:size(svlfiles,1)
        SAMinfo(ii,:)=sscanf(svlfiles(ii).name,'%d,%dto%dms,%d-%dHz,D3.svl')';
    end
    timewins = [SAMinfo(:,2) SAMinfo(:,3)];

    tf = zeros(length(SAMimage),size(timepts,1));

    for timebin=1:size(timewins,1)
        filename = svlfiles(timebin).name;
        [SAMimage,coords,SAMparams] = nut_read_svl(filename);
        bw = [SAMparams.hpfreq SAMparams.lpfreq];
        if(any(bw ~= bands(freqbin,:)))
            error('bandwidth mismatch.');
        end
        
        tf(:,timebin) = 10*log10(weirdFtostandardF(SAMimage));
%         tf(:,timebin,:)=SAMimage;
    end
    
    for timebin=1:size(timepts,1)
        timeselect = find((timepts(timebin) > timewins(:,1)) & (timepts(timebin) < timewins(:,2)));
        beam.tf(:,timebin,freqbin) = mean(tf(:,timeselect,:),2);

        % convert to dB
    end
end

% beam.s_timef(:,1,1) = SAMimage;

beam.voxels = coords*1000; % convert from m to mm
beam.voxelsize = repmat(SAMparams.stepsize*1000,[1 3]);
beam.params.beamformertype='CTF-SAM';



% dummy info
beam.coreg.mripath = '/data/advaita/nutmeg/blank.img';
beam.coreg.meg2mri_tfm = eye(4);

[pathname,basename,ext]=fileparts(filename);
outputfile = fullfile(pathname,['s_beamtf_' basename '.mat']);
save(outputfile,'beam');


function F=weirdFtostandardF(Fweird)
% converts CTF's weird F value to a standard F ratio
F = zeros(size(Fweird));
F(Fweird>=0) = Fweird(Fweird>=0)+1;
F(Fweird<0) = 1./(1-Fweird(Fweird<0));
