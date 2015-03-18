function Lp = nut_import_smacleadfield(fn)
% Lp = nut_import_smacleadfield('leadfieldfile.lf');
% Adapted from SMAC's ioReadLeadfield


leadfield = double([]);
machineFormat = 'native';

% Open file for reading only
[fid, message] = fopen(fn, 'r', machineFormat);
if fid == -1
    disp(sprintf('nut_import_smacleadfield: %s', message));
    return;
end

% Header : nbRows(Ne) and nbCols(3*M)
leadfieldSize = fread(fid, 2, 'uint32');
% Lead field : nbRows*nbCols of 64bit floating point reals, 
%             stored ROW BY ROW (while natural Matlab storage
%             is column by column) => TRANSPOSITION
[leadfield, count] = fread(fid, prod(leadfieldSize), 'float64');
if (count ~= prod(leadfieldSize))
    error('nut_import_smacleadfield:OpenFile', sprintf('Corrupted leadfield file'));
end
leadfield = reshape(leadfield, leadfieldSize(2), leadfieldSize(1))';

% Close file
fclose(fid);

% bring to nutmeg convention
Lp = zeros(size(leadfield,1),3,size(leadfield,2)/3);
for k=1:3
    Lp(:,k,:)=leadfield(:,k:3:end);
end