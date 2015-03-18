function nut_constructSAMparams(startpt,endpt,windowsize,stepsize,controlstart,controlend,activemarker,controlmarker,condname)
% Construct time window parameters for CTF SAM space-time-frequency analysis.
% nut_constructSAMparams(startpt,endpt,windowsize,stepsize,baselinestart,baselineend, ...
%            {activemarker},{controlmarker},{condname})
%
% Arguments in {} are optional.
% Note that active and baseline time windows must have the same length; therefore
% baselineend-baselinestart must be equal to windowsize.
% activemarker and controlmarker can contain the name of a CTF marker. If empty, the
%           analysis is relative to time 0 in all trials.
% condname  is an arbitrary condition identifier, that corresponds to the svlstart
%           argument in nut_svl2timef2


if (nargin<7 || isempty(activemarker))
    activemarker='_time_';
end
if (nargin<8 || isempty(controlmarker))
    controlmarker='_time_';
end
if nargin<9
    condname='';
else
    condname=[condname ','];
end

t1 = (startpt:stepsize:(endpt-windowsize));
t2 = ((startpt+windowsize):stepsize:endpt);

for ii=1:length(t1)
    filestring = sprintf('%s%02d,%+dto%+dms',condname,ii,t1(ii),t2(ii));
    fid=fopen(filestring,'w');
    
    
    fprintf(fid,'1\n%s   %f   %f\n',activemarker,t1(ii)/1000,t2(ii)/1000);
    fprintf(fid,'1\n%s   %f   %f\n',controlmarker,controlstart/1000,controlend/1000);
   
    fclose(fid);
end