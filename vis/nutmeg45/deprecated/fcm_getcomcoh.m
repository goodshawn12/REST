function ICI=fcm_getcomcoh(killcomcohfolder)
% FCM_GETCOMCOH  loads and assembles complex coherence data. 
%
%  CC = fcm_getcomcoh(killcomcohfolder)
%
% Data to be read must be in folder 'comcoh' in current working directory.
% killcomcohfolder      (optional) if true, the directory "comcoh" is deleted after
%                       reading. Default is false.

warning('FCM_GETCOMCOH has been replaced by FCM_ASSEMBLE. It will be removed in the near future.')

if nargin<1, killcomcohfolder=false; end

ICI.coh=[];
ICI.comps=[];

cd comcoh

d=dir('CC*.mat');
numjobs=length(d);
if numjobs==0, error('No files.'), end
disp(int2str(numjobs))

for k=1:numjobs
    if isempty(strfind(d(k).name,sprintf('%02d',k)))
        error('Missing file nr %d',k)
    end
    load(d(k).name);
    if ~isfield(CC,'coh') || isempty(CC.coh)
        error('Invalid data file.')
    end
    if size(CC.coh,1)~=size(CC.comps,1)     % Legacy compatibility
        CC.coh=permute(CC.coh,[2 3 1]); 
    end
    ICI.coh=cat(1,ICI.coh,CC.coh);
    ICI.comps=cat(1,ICI.comps,CC.comps);
end

ICI.frq=CC.frq;
ICI.time=CC.time;
ICI.N = CC.N;
if isfield(CC,'taperwidth')
    ICI.taperwidth=CC.taperwidth;
    ICI.numtaper=CC.numtaper;
end

cd ..

if killcomcohfolder
    rmdir('comcoh','s')
end

