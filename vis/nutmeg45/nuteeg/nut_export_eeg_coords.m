function nut_export_eeg_coords

global nuts

if ~isfield(nuts,'meg') || ~isfield(nuts.meg,'sensorCoord')
    errordlg('No coordinates to export!'), return
end

[fi, pa] = uiputfile( ...
   {'*.sfp;*.xyz' ,'Supported Formats'; ...
    '*.sfp'       ,'BESA'; ...
    '*.xyz'       ,'Cartool electrode coordinates'}, ...
    'Export Sensor Coordinates...');
if isequal(fi,0) || isequal(pa,0), return, end
outfile = fullfile(pa,fi);
[dum,dum,ext]=fileparts(outfile);

numsens=length(nuts.meg.sensorCoord);
sensCoord = nut_meg2mri(nuts.meg.sensorCoord);

fid=fopen(outfile,'wt');
switch ext
    case '.xyz'
        fprintf(fid,'%d\t1<n',numsens);
        for k=1:numsens
            fprintf(fid,'%1.4f\t%1.4f\t%1.4f\t%s\n',sensCoord(k,:),nuts.meg.sensor_labels{k});
        end
    case '.sfp'
        %TODO: allow exporting fiducials in same file if available
        for k=1:numsens
            fprintf(fid,'%s\t%1.4f\t%1.4f\t%1.4f\n',nuts.meg.sensor_labels{k},sensCoord(k,:));
        end
    otherwise
        fclose(fid); delete(outfile);
        errordlg('Only Cartool''s .xyz and BESA''s .sfp formats are currently supported for export.'), return
end
fclose(fid);

msgbox(sprintf('Electrode coordinates were exported as %s.',outfile))
        