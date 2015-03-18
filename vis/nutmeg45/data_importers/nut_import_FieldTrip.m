function nut_import_FieldTrip

global nuts ndefaults
[filename, pathname]=uigetfile('*.mat','Open FieldTrip "data" file...');
if isequal(filename,0)|isequal(pathname,0)
  return;
end
[crap,datafile,morecrap] = fileparts(filename);
nuts.meg.filename = datafile;

load(fullfile(pathname,filename));

[filename, pathname]=uigetfile('*.mat','Open FieldTrip "grid" file. Select "Cancel" if not needed.');
if isequal(filename,0)|isequal(pathname,0)
  % no grid file selected
  nuts = nut_ft2nuts(data);
else
  grid=load(fullfile(pathname,filename));
  gridname=fieldnames(grid);
  grid=grid.(gridname{1});  
  nuts = nut_ft2nuts(data,grid);
end
