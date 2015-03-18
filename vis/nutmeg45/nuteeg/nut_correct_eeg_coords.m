function nut_correct_eeg_coords(label2corr,templatefile)
% can interpolate single electrodes whose position was badly digitized.
% nut_correct_eeg_coords(label2corr,[templatefile])

global nuts

if isempty(nuts), error('Must load NUTMEG session first.'), end

bad = strmatch(label2corr,nuts.meg.sensor_labels,'exact');
if isempty(bad), error('No electrode %s found.',label2corr), end

if nargin<2
    if strncmp(spm('ver'),'SPM8',4), cwd=pwd; cd([fileparts(which('spm.m')) filesep 'EEGtemplates']); end
    [template_filename, template_path]=uigetfile( ...
       {'*.sfp;*.elc'       ,'Supported Formats'; ...
        '*.sfp'             ,'BESA'; ...
        '*.elc'             ,'ASA'}, ...
        'Load MNI template coordinates of your system...');
    if isequal(template_filename,0) || isequal(template_path,0)
        return
    end
    templatefile = fullfile(template_path,template_filename);
    if strncmp(spm('ver'),'SPM8',4), cd(cwd), clear cwd, end
end

elec = ft_read_sens(templatefile);
tind = ismember(elec.label,nuts.meg.sensor_labels);
elec.pnt = elec.pnt(tind,:); elec.label = elec.label(tind); clear tind
badintemp = strmatch(label2corr,elec.label,'exact');
if isempty(badintemp), error('No electrode %s found in template.',label2corr), end

elec.pnt = nut_mri2meg(nut_mni2mri(elec.pnt));
[dum,surr] = sort(nut_rownorm(nut_coord_diff(elec.pnt,elec.pnt(badintemp,:))));

surr = surr(2:5); % index of 4 closest electrodes

ratio = repmat(elec.pnt(badintemp,:),[4 1]) - elec.pnt(surr,:);

surrind = ismember(nuts.meg.sensor_labels,elec.label(surr));
nuts.meg.sensorCoord(bad,:) = mean(nuts.meg.sensorCoord(surrind,:) + ratio,1);







