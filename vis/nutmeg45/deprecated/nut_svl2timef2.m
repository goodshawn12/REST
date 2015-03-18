function nut_svl2timef2(SAMdir,coreg,outbase,inbase,isavg)
% nut_svl2timef2(SAMdir,{coreg},{output_filename},{svlstart})
%
% converts SAM svl volume into s_beam.
% Configuration: Timewindows must have the same "middle frequency" and the same
%                time step size across frequency bands, and all frequency bands 
%                must have the same number of timewindows. Eg:
%                                    "Middle frequency":   150ms     250ms
%                                    Time step size:            100ms
%                Low frequency band (time width 300ms):   0-300ms  100-400ms 
%                Middle freq band (time width 200ms):    50-250ms  150-350ms
%                High freq band (time width 100ms):     100-200ms  200-300ms
%
% INPUT [variables in {} are optional; use can use [] to set to default]
% coreg:           Optional input containing a structure with MRI coregistration.
% output_filename: Optional input containing filename of new file to be created.
% svlstart:        Optional string to differentiate svl files from different settings
%                  in the same SAM folder. The filenames are expected to be named, e.g.:
%                  svlstart,01,-100to+100ms,15-30Hz,D3.svl

if nargin<2; coreg=[]; end
if nargin<3, outbase=[]; end
if nargin<4, 
    inbase=[]; 
elseif ~isempty(inbase)
    inbase=[inbase ','];
end
if nargin<5, isavg=false; end
dum=pwd;
cd(SAMdir);

svlfiles = dir([inbase '*.svl']);
[SAMimage,coords,SAMparams] = nut_read_svl(svlfiles(1).name);

for ii=1:size(svlfiles,1)
    SAMinfo(ii,:)=sscanf(svlfiles(ii).name,[inbase '%d,%dto%dms,%d-%dHz,D3.svl'])';
end

bands = [unique(SAMinfo(:,4:5),'rows')];
timewins = [unique(SAMinfo(:,2:3),'rows')];

timepts = unique(mean(timewins,2));
if length(unique(diff(timepts)))>1, error('Middle frequencies and time steps of timewindows must not differ across frequency bands.'), end
beam.timepts = timepts;

if (size(timewins,1)>1)
    timestep = timewins(2,1)-timewins(1,1);
    beam.timewindow = [timepts-timestep/2 timepts+timestep/2];
else
    timestep = 100; % otherwise, arbitrary
    beam.timewindow = timewins;
end

if(size(timepts,1)==1)
    beam.srate = 1;
else
%     beam.timewindow = timewins;
    beam.srate = 1000/timestep;
end
beam.bands = bands;


beam.s{1} = zeros(length(SAMimage),size(timepts,1),size(bands,1));

for freqbin=1:size(bands,1)
    svlfiles = dir([inbase '*' num2str(bands(freqbin,1)) '-' num2str(bands(freqbin,2)) 'Hz*.svl']);
    clear SAMinfo
    for ii=1:size(svlfiles,1)
        SAMinfo(ii,:)=sscanf(svlfiles(ii).name,[inbase '%d,%dto%dms,%d-%dHz,D3.svl'])';
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
        
        tf(:,timebin) = 10*log10(weirdFtostandardF(SAMimage));   % convert to dB
%         tf(:,timebin,:)=SAMimage;
    end
    
    if isavg    % normally false, but this is done by other version of nut_svl2timef
        for timebin=1:size(timepts,1)
            timeselect = find((timewins(:,1) < timepts(timebin)) & (timewins(:,2) > timepts(timebin)));
            beam.s{1}(:,timebin,freqbin) = mean(tf(:,timeselect),2);
        end
    else
        beam.s{1}(:,:,freqbin)=tf;
    end
end

% beam.s_timef(:,1,1) = SAMimage;

beam.voxels = coords*1000; % convert from m to mm
beam.voxelsize = repmat(SAMparams.stepsize*1000,[1 3]);
beam.params.beamformertype='CTF-SAM';


% dummy info
if isempty(coreg)
    beam.coreg.mripath = which('blank.img');
    beam.coreg.meg2mri_tfm = eye(4);
else
    beam.coreg=coreg;
end

if isempty(outbase)
    [pathname,outbase,ext]=fileparts(filename);
end
outputfile = ['s_beamtf_' outbase '.mat'];
save(outputfile,'beam');
cd(dum)

% ------------------------------------
function F=weirdFtostandardF(Fweird)
% converts CTF's weird F value to a standard F ratio
F = zeros(size(Fweird));
F(Fweird>=0) = Fweird(Fweird>=0)+1;
F(Fweird<0) = 1./(1-Fweird(Fweird<0));
