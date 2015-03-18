function [pnt, tri]=read_om_tri(trifile)

%%
fid = fopen(trifile,'r');
vnum=abs(str2num(fgetl(fid)));


%%

pnt=cell2mat(textscan(fid,'%f %f %f %f %f %f'));
pnt=pnt(:,1:3);
tnum=abs(str2num(fgetl(fid)));
tri=cell2mat(textscan(fid,'%f %f %f')) + 1; % nerds count from 0, matlab counts from 1

%%
fclose(fid);


%%


