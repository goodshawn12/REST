%%
% Exports localsphere centers.
% @author Daniel D.E. Wong
%%
function nut_exportLsc

global nuts
lsc = nuts.meg.lsc;
lsc_sensor_labels = nuts.meg.lsc_sensor_labels;

[file,path] = uiputfile('*.mat','Save BEM matrices');  % BEM points are saved in MRI coordinates
if ~isequal(file,0) & ~isequal(path,0)
    save(strcat(path,file),'lsc','lsc_sensor_labels');
end