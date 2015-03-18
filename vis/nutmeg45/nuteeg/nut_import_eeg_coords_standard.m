function nut_import_eeg_coords_standard
% Imports coordinates of standard EEG coordinates, if individual headpoints or
% electrode positions are not available. The standard coordinates are then
% adjusted to the individual head shape (just a first approximation, it is
% recommended to project the sensors to the scalp in a second step, using
% nut_eegmralign).
% Electrode coordinate files must indicate Cartesian coordinates, organized
% as follows:
%   xcoord  ycoord  zcoord  eleclabel
% with the nasion being +X, the left hemisphere +Y, the right hemisphere -Y, 
% and the vertex [0 0 +Z]. 
% The individualized coordinates are stored in nuts.meg.sensorCoord.


global nuts st

if ~isfield(nuts,'coreg'), errordlg('Load MRI first.'), return, end
if ~isfield(nuts.coreg,'fiducials_mri_mm'), errordlg('You need to define the fiducials.'), return, end

% Import sensor locations and headshape
nutroot = fileparts(which('nutmeg.m'));
cwd = pwd;
cd([nutroot filesep 'templates'])
[efile, epath] = uigetfile({'*.txt;*.els;*.xyz;*.csv', 'Electrode location files (*.txt, *.els, *.xyz, *.csv)'; ...
                                        '*.xyz;*.els' , 'Cartool electrode coordinates (*.xyz, *.els)'; ...
                                        '*.txt' , 'Text file'; ...
                                        '*.csv' , 'Comma separated file'},'Open electrode location file');
cd(cwd);
if isequal(efile,0)
    return;
end

% Import
[dum,dum,ext] = fileparts(efile);
switch lower(ext)
    case '.xyz'
        [x y z lab] = textread(fullfile(epath,efile),'%f%f%f%s','headerlines',1);
        
    case {'.txt','.csv'}
        % Data must be saved as: x y z eleclabel, without header lines
        [x y z lab] = textread(fullfile(epath,efile),'%f%f%f%s');
        
    case '.els'
        [x y z lab] = textread(fullfile(epath,efile),'%f%f%f%s','headerlines',6);
        % TODO: this does not work if several clusters of electrodes are
        % defined.
    otherwise
        errordlg('This data format is not supported yet, but it may be easy to code...')
        return
end
       

loc = [x y z]; 
if isfield(nuts.meg,'sensor_labels') 
    if ~isequal(length(nuts.meg.sensor_labels),size(loc,1))
        errordlg('The number of electrodes in this location file do not match the number of data channels.')
        return
    end
else
    nuts.meg.sensor_labels=lab;
end
clear x y z lab

% Convert units
% unit = questdlg('Select units of the file', 'Import xyz', 'normalized', 'mm', 'normalized');
% switch unit
%     case 'normalized'
% Blow up unitary shere to arbitrary 100 mm 
%loc = loc.* 100;
%     otherwise
% end

% Get mri fiducials in MEG space
mrifidmeg = nut_mri2meg(nuts.coreg.fiducials_mri_mm);

% Obtain indices of landmark electrodes
czidx=[]; fpzidx=[]; t7idx=[]; t8idx=[]; ozidx=[];
if ( strcmpi(nuts.meg.system,'Biosemi') && length(nuts.meg.sensor_labels)>=128 && length(nuts.meg.sensor_labels)<160)  % Biosemi 128 channel system
    czidx  = strmatch('A1',nuts.meg.sensor_labels,'exact');
    fpzidx = strmatch('C17',nuts.meg.sensor_labels,'exact');
    t7idx = strmatch('D23',nuts.meg.sensor_labels,'exact');
    t8idx = strmatch('B26',nuts.meg.sensor_labels,'exact');
    ozidx = strmatch('A23',nuts.meg.sensor_labels,'exact');
elseif ( strcmpi(nuts.meg.system,'EGI') && length(nuts.meg.sensor_labels)>=204 )    % EGI 257 channel system
    czidx = strmatch('REF',nuts.meg.sensor_labels,'exact');
    fpzidx = strmatch('30',nuts.meg.sensor_labels,'exact');
    t7idx = strmatch('69',nuts.meg.sensor_labels,'exact');
    t8idx = strmatch('194',nuts.meg.sensor_labels,'exact');
    ozidx = strmatch('127',nuts.meg.sensor_labels,'exact');
% elseif
%   ADD your system's landmark elec labels here
else % check if the elec names are according to the 1020 system
    czidx  = strmatch('Cz',nuts.meg.sensor_labels);
    fpzidx = strmatch('Fpz',nuts.meg.sensor_labels);
    t7idx = strmatch('T7',nuts.meg.sensor_labels);
    t8idx = strmatch('T8',nuts.meg.sensor_labels);
    ozidx = strmatch('Oz',nuts.meg.sensor_labels);
end
    
if ~isscalar(czidx)
    resp = inputdlg({'Indicate name of vertex electrode (Cz):'},'Import standard elec positions',1);
    czidx = strmatch(resp{1},nuts.meg.sensor_labels,'exact');
end
if ~isscalar(fpzidx)
    resp = inputdlg({'Indicate name of electrode closest to nasion (Fpz):','Distance from nasion (in mm):'},'Import standard elec positions',1,{'','35'});
    fpzidx = strmatch(resp{1},nuts.meg.sensor_labels,'exact');
    fpzdist = str2num(resp{2});
else
    resp = inputdlg(sprintf('Distance of electrode %s from nasion (in mm):',nuts.meg.sensor_labels{fpzidx}),'Import standard elec positions',1,{'35'});
    fpzdist = str2num(resp{1});
end
if ~isscalar(t7idx)
    resp = inputdlg({'Indicate name of electrode closest to left ear (T7):'},'Import standard elec positions',1);
    t7idx = strmatch(resp{1},nuts.meg.sensor_labels,'exact');
end
if ~isscalar(t8idx)
    resp = inputdlg({'Indicate name of electrode closest to right ear (T8):'},'Import standard elec positions',1);
    t8idx = strmatch(resp{1},nuts.meg.sensor_labels,'exact');
end
if ~isscalar(ozidx)
    resp = inputdlg({'Indicate name of electrode closest to inion (Oz):'},'Import standard elec positions',1);
    ozidx = strmatch(resp{1},nuts.meg.sensor_labels,'exact');
end

% Obtain scalp surface
resp = questdlg('I need the scalp surface to adjust standard coordinates to individual head shape. How should I obtain it?', ...
    'Adjust coordinates.', 'Load surface file', 'Calculate from MRI', 'Load surface file');
switch resp
    case 'Load surface file'
        if ~isfield(nuts.coreg,'volfile')
            [sfile, spath] = uigetfile({'*_vol.mat', 'NUTEEG head volume file (*_vol.mat)'},'Open surface file');
            if isequal(sfile,0), return, end
            nuts.coreg.volfile = fullfile(spath,sfile);
        end
        load(nuts.coreg.volfile);
        position = nut_mri2meg(vol.bnd(1).vertices);
        clear vol        
    case 'Calculate from MRI'
        thresh = nut_find_thresh(spm_read_vols(st.vols{1}));
        position=nut_vol2cloud6(thresh,spm_read_vols(st.vols{1}));
        position = nut_mri2meg(nut_voxels2mm(position));
    otherwise
        return
end

% Match with scalp surface
% scalp size
t0 = max(position(:,3));    % vertex
tb = t0 - fpzdist;
l0 = mrifidmeg(1,2);
r0 = mrifidmeg(2,2);
a0 = mrifidmeg(3,1);
p0 = min(position(:,1));

% stretch eeg coordinates
%loc(:,2) = loc(:,2) - loc(czidx,2);
loc(:,3) = loc(:,3)+100;  % arbitary shift to all positive values
zgain = 1/ ((loc(czidx,3)-loc(fpzidx,3)) / tb);
loc(:,3) = loc(:,3).*zgain;
loc(:,3) = loc(:,3) - mean([ loc(fpzidx,3)-fpzdist   loc(czidx,3)-t0 ]); % shift back to correct positon.

lgain = 1/ ( loc(t7idx,2) / l0 );
eidx = (loc(:,2)>0);
loc(eidx,2)=loc(eidx,2).*lgain;

rgain = 1/ ( loc(t8idx,2) / r0 );
eidx = (loc(:,2)<0);
loc(eidx,2)=loc(eidx,2).*rgain;

again = 1/ ( loc(fpzidx,1) / a0 );
eidx = loc(:,1)>0;
loc(eidx,1)=loc(eidx,1).*again;

pgain = 1/ ( loc(ozidx,1) / p0 );
eidx = loc(:,1)<0;
loc(eidx,1)=loc(eidx,1).*pgain;

% Add to nuts structure
nuts.meg.sensorCoord = loc;

resp = questdlg('Would you like to project the electrodes to the scalp (recommended but needs the Helsinki BEM library)?', ...
    'Adjust coordinates', 'Yes', 'No', 'Yes');
if strcmp(resp, 'Yes')
    nut_eegmralign;
end
