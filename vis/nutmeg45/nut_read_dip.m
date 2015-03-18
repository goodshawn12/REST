function dip = nut_read_dip(filename);
% NUT_READ_DIP   reads CTF dipole files (*.dip).
%
%   dip = nut_read_dip(filename);
%
%   dip         structure containing dipoles, model, and underlying structural MRI
%   filename    name of .dip file

fid = fopen(filename, 'rt');
if fid==-1
  error(sprintf('NUT_READ_DIP: Could not open file %s', filename));
end

currlin = '';
while ~strcmp(currlin,'Model')                          % Model
    currlin=fgetl(fid);
end
currlin=fgetl(fid);                                     % {
currlin=fgetl(fid);                                     % MODEL_TYPE:   ........
[dum, value] = strtok(currlin, ':');
dip.model = value(1+find(~isspace(value(2:end))));      

currlin=fgetl(fid);                                     % MRI_FILE:     ......
[dum, value] = strtok(currlin, ':');
dip.mri = value(1+find(~isspace(value(2:end))));

currlin = '';
while ~strcmp(currlin,'Dipoles')                        % Dipoles
    currlin=fgetl(fid);
end
for k=1:3                                               % { , // , //
    currlin=fgetl(fid);
end

currlin=fgetl(fid);                                     % 1: .......
[d,value]=strtok(currlin,':');
d=str2num(d);
while ~isempty(d)
    [x_cm,y_cm,z_cm,xo,yo,zo,mom_nA,lab]=strread(value,':%f%f%f%f%f%f%f%s');
    dip.dipole(d)=struct('x_cm',x_cm,'y_cm',y_cm,'z_cm',z_cm,'xo',xo,'yo',yo,'zo',zo,'mom_nA',mom_nA,'label',char(lab));
    %[dip.dipole(d).x_cm,dip.dipole(d).y_cm,dip.dipole(d).z_cm,dip.dipole(d).xo,dip.dipole(d).yo,dip.dipole(d).zo,dip.dipole(d).mom_nA,dip.dipole(d).label]=strread(value,':%f%f%f%f%f%f%f%s'); 
    %dip.dipole(d).label=char(dip.dipole(d).label);
    currlin=fgetl(fid);                                 % 2:end: .......
    [d,value]=strtok(currlin,':');
    d=str2num(d);
end

fclose(fid);
