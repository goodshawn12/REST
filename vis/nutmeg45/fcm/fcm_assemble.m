function CC = fcm_assemble(confilename,killfolder)
% FCM_ASSEMBLE  loads and assembles source connectivity data calculated on a Linux cluster.
%  
%  fcm_assemble(confilename,killfolder)
%
% Data to be read must be in folder 'ccohere'/'pli'/'ampcorr' in current working directory.
% confilename     name of assembled file to be created.
% killfolder      (optional) if true, the directory 'ccohere'/'pli'/'ampcorr' is deleted after
%                 reading. Default is false, except if the directory contains only one file, in
%                 which case the directory is always deleted.

if nargin<2, killfolder=false; end

global fuse

ICI.coh=[];
ICI.comps=uint32([]);

fona = fuse.funconn;
if strcmp(fuse.funconn,'ccohere') && exist('comcoh','dir')      % legacy compatibility
    if ~exist('ccohere','dir')
        fona='comcoh';
    else
        warning('Reading folder "ccohere" and NOT "comcoh"!')
    end
end

cd(fona)

d=dir('CC*.mat');
numjobs=length(d);
if numjobs==0, error('No files.'), end
disp(int2str(numjobs))

if numjobs==1
    if (length(confilename)<4 || ~strcmpi(confilename(end-3:end),'.mat')), confilename=[confilename '.mat']; end
    if ispc
        dos(['move ' d(1).name ' ..\' confilename]);
    else
        unix(['mv ' d(1).name ' ../' confilename]);    
    end
    cd ..
    rmdir(fona,'s')
    return
end

for k=1:numjobs
    if isempty(strfind(d(k).name,sprintf('%02d',k)))
        error('Missing file nr %d',k)
    end
    load(d(k).name);
    if ~isfield(CC,'coh') || isempty(CC.coh)
        error('Invalid data file.')
    end
    if isfield(CC,'method') && ~strcmpi(fuse.funconn,CC.method)
        error('Mismatch between the functional connectivity method set in your current configuration and the method used for the file.')
    end
    ICI.coh=cat(1,ICI.coh,CC.coh);
    ICI.comps=cat(1,ICI.comps,CC.comps);
end

ICI.frq=CC.frq;
ICI.time=CC.time;
ICI.N = CC.N;
if isfield(CC,'len'),
    ICI.len = CC.len;
end
ICI.method = fuse.funconn;
if isfield(CC,'taperwidth')
    ICI.taperwidth=CC.taperwidth;
    ICI.numtaper=CC.numtaper;
end
CC=ICI; clear ICI

cd ..
fprintf('Saving as %s.\n',confilename)
save(confilename,'CC');


if killfolder
    rmdir(fona,'s')
end

