function fcm_preparejobs(comps,numjobs)
% FCM_PREPAREJOBS creates text files for parallel job distribution on Linux cluster.
%  fcm_preparejobs(comps,numjobs)
% COMPS     connection definition obtained from fcm_conndef.
% NUMJOBS   number of computers that you want to work on your data.

ncomps=size(comps,1);

step=ceil(ncomps/numjobs);
index=1:step;

if ~exist('comps','dir')
    mkdir comps
else
    error('A directory named "comps" already exists.')
end
cd comps
for k=1:numjobs
    CI=comps(index,:);
%    fid=fopen(sprintf('CI%02d.txt',k),'wt');
%     for lin=1:length(index)
%         fprintf(fid,'%d %d\n',CI(lin,:));
%     end
%    fclose(fid);
    dlmwrite(sprintf('CI%02d.txt',k),CI,' ');
    index=index+step;
    index=index(index<=ncomps);
end
