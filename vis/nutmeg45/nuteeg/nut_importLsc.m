%%
% Imports localspheres exported from NUTEEG.
% @author Daniel D.E. Wong
%%
function nut_importLsc

global nuts
[file,path] = uigetfile('*.mat','Open Localsphere File');
if isequal(file,0); return; end;
load(strcat(path,file));

nuts.meg.lsc = lsc;
nuts.meg.lsc_sensor_labels = lsc_sensor_labels;