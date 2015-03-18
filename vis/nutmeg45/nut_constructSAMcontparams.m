function nut_constructSAMcontparams(startpt,endpt,windowsize,stepsize,activemarker,controlmarker)
% Construct time window parameters for CTF SAM space-time-frequency analysis with SAM
% constasts between 2 conditions.
% nut_constructSAMcontparams(startpt,endpt,windowsize,stepsize,{activemarker},{controlmarker})
%
% Arguments in {} are optional.

if (nargin<5 || isempty(activemarker))
    activemarker='_time_';
end
if nargin<6
    controlmarker='_time_';
end

t1 = (startpt:stepsize:(endpt-windowsize));
t2 = ((startpt+windowsize):stepsize:endpt);

for ii=1:length(t1)
    filestring = sprintf('%02d,%+dto%+dms',ii,t1(ii),t2(ii));
    fid=fopen(filestring,'w');
    
    
    fprintf(fid,'1\n%s   %f   %f\n',activemarker,t1(ii)/1000,t2(ii)/1000);
    fprintf(fid,'1\n%s   %f   %f\n',controlmarker,t1(ii)/1000,t2(ii)/1000);
   
    fclose(fid);
end