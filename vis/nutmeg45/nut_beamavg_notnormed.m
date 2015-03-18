function nut_beamavg_notnormed(varargin)
% nut_beamavg_notnormed(varargin)
% varargin should be names of as many beam files as desired to average
% saves output as name of first file with _avg5 appened (where 5 is the
% number of files averaged, for example
% note that all files must contain the same voxels

ignoreorientation = true; % i can't imagine our sketchy orientations are useful for this; % depends on the method!
avgovertime = true;

for ii=1:nargin
    load(varargin{ii});
    if avgovertime
        t0=dsearchn(beam.timewindow,0);
        s_th{ii}=mean(beam.s_th(:,t0:end),2);
        s_ph{ii}=mean(beam.s_ph(:,t0:end),2);
        s_z{ii}=mean(beam.s_z(:,t0:end),2);
    else
        s_th{ii}=beam.s_th;
        s_ph{ii}=beam.s_ph;
        s_z{ii}=beam.s_z;
    end
    s_beam{ii} = sqrt(s_th{ii}.^2 + s_ph{ii}.^2 + s_z{ii}.^2);
    if ii==1
        beamfirst=beam;
    end
    if ignoreorientation
        maxblob{ii}=max(s_beam{ii});
        minblob{ii}=min(s_beam{ii});
        scalefactor{ii}=1000/max(abs(maxblob{ii}),abs(minblob{ii}))
    else
        maxblob{ii}=max([max(s_th{ii}),max(s_ph{ii}),max(s_z{ii})]);
        minblob{ii}=min([min(s_th{ii}),min(s_ph{ii}),min(s_z{ii})]);
        scalefactor{ii}=1000/max(abs(maxblob{ii}),abs(minblob{ii}))
    end
end
clear beam
beam=beamfirst;
beam.s_th=zeros(size(s_th{ii}));
beam.s_ph=zeros(size(s_ph{ii}));
beam.s_z=zeros(size(s_z{ii}));
if avgovertime
    beam.timewindow=1;
end
if(ignoreorientation)
    for ii=1:nargin
        beam.s_th=beam.s_th+scalefactor{ii}*s_beam{ii};
        beam.s_ph=zeros(size(beam.s_th));
        beam.s_z=zeros(size(beam.s_th));
    end
else
    for ii=1:nargin
        beam.s_th=beam.s_th+s_th{ii}*scalefactor{ii};
        beam.s_ph=beam.s_ph+s_ph{ii}*scalefactor{ii};
        beam.s_z=beam.s_z+s_z{ii}*scalefactor{ii};
    end
end
beam.s_th=beam.s_th/nargin;
beam.s_ph=beam.s_ph/nargin;
beam.s_z=beam.s_z/nargin;

[pathname,filename,fileext] = fileparts(varargin{1});
save([filename '-avg' num2str(nargin) '.mat'],'beam');