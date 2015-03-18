function nut_batch(sessionfile,dsfile,lscfile,voxelsfile,voxelsize,sbeamname,preprocessing,whichsensors,nwd);
%function nut_batch(sessionfile,dsfile,lscfile,voxelsfile,voxelsize,sbeamname,preprocessing,whichsensors,nwd)
%
% sessionfile should have MRI and fiducials loaded already, and optional VOI
%
% example of inputs:
% nwd='/data/mimamsa/johannaz/dti_fmri_meg_data/meg/jz/sommot12dec03/sef_dig1_lvox_asens_2mm/nuts/';
% sessionfile=[nwd 'allinfo.mat'];
% dsfile='/data/mimamsa/johannaz/dti_fmri_meg_data/meg/jz/sommot12dec03/ZumerJohanna_SEF_20031212_01-blav.ds/ZumerJohanna_SEF_20031212_01-blav.meg4';
% lscfile='/data/mimamsa/johannaz/dti_fmri_meg_data/meg/jz/default.hdm';
% basecorrect=2;
% voxelsfile=[nwd 'voxels.mat'];
% voxelsize=5;  % in mm
% numeigs=2;
% timestart=0;     % in ms
% timestop=240;     % in ms
% algorithm='Default Beamformer';
% sbeamname=[nwd  'ZumerJohanna_SEF_20031212_01-blav_0to240_2_batch.mat'];
% %end example

global nuts bolts
cwd=pwd;
if(exist('nwd','var'))
    cd(nwd);
end
nut_defaults;
% nut_opensession(sessionfile)
load(sessionfile);
nut_importmeg(dsfile,lscfile);

nut_preprocessing_defaults;
if exist('whichsensors','var') & ~isempty(whichsensors)
    nuts.meg.goodchannels=whichsensors;
end
if exist('preprocessing','var') & ~isempty(preprocessing)
    preprocfields = fieldnames(preprocessing);
    
    % merge supplied preprocessing structure with defaults
    for i=1:length(preprocfields)
        nuts.preprocessing.(preprocfields{i})=preprocessing.(preprocfields{i});
    end
end
% nuts.preprocessing.bf_timeinterval=[timestart timestop];
% nuts.preprocessing.signalspace=numeigs;
if exist('voxelsfile','var')   %voxelsfile is optional;  if not used, then VOIvoxels should be in sessionfile already.
    load(voxelsfile);
    [nuts.Lp,nuts.voxels]=nut_compute_lead_field(voxelsize,VOIvoxels);  
else
    [nuts.Lp,nuts.voxels]=nut_compute_lead_field(voxelsize);
end
bfguihandle = nut_beamforming_gui('nut_baseline_Callback');
delete(bfguihandle);

%data = nut_filter(data,beam.timewindow,nuts.meg.srate,nuts.preprocessing.baseline,nuts.preprocessing.notch,nuts.preprocessing.bpf,nuts.preprocessing.bpf_low_cutoff,nuts.preprocessing.bpf_high_cutoff);

% numeigs = nuts.preprocessing.signalspace;
time1 = nuts.preprocessing.bf_timeinterval(1); % assumes time intervals of interest are those of the beamformer tool
time2 = nuts.preprocessing.bf_timeinterval(2);
% algorithm = nuts.preprocessing.beamformertype;
timept1 = dsearchn(nuts.meg.latency,time1);
timept2 = dsearchn(nuts.meg.latency,time2);
handles=[];
denoise_str = strrep(['nut_' nuts.preprocessing.denoisertype],' ','_');
feval(denoise_str,handles,nuts.preprocessing.signalspace,timept1,timept2); % update bolts.meg, bolts.params.Rzz1 (bolts.meg should be updated from nuts.meg.data every time Beamformer Tool is opened otherwise the denoising will be applied multiple times)

% nut_cov_eigen(numeigs,timestart,timestop)
% nut_cov_eigen(nuts.preprocessing.invtype);
bolts.params.InvRzz1=nut_inv(double(bolts.params.Rzz1),nuts.preprocessing.invtype);

nut_generate_beamforming_activations(sbeamname,timept1,timept2,nuts.preprocessing.beamformertype);
cd(cwd);
%close all hidden
clear global nuts beam bolts
% clear all global

