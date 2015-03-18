function nut_exportmultimg(name)
% NUT_EXPORTMULTIMG  export an analyze image for each time window.
% 
%  nut_exportmultimg('name_of_exported_files')
%
% nut_timef_viewer must be open, showing the data with the threshold
% settings you want to use for export. The frequency band you want to export must 
% be selected.

global beam

for k=1:length(beam.timepts)
    nut_mat2img(k,sprintf('%s%02d.img',name,k));
end