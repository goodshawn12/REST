function nut_defaults
% defaults
% please change them depending on your preferences.
% make an individual copy for each user, putting it ahead of the other
% nutmeg code in your path.  

disp(['Using NUTMEG settings located here: ' which('nut_defaults')]);

global ndefaults

% helps maintain similar GUI appearance/spacing across platforms
set(0,'defaultuicontrolfontname','Sans Serif');

% MNI to Talairach transform method
% ='icbm' for MNI-Talairach conversion based on ICBM (Lancaster) method;
% ='brett' for original Matthew Brett method
ndefaults.mni2tal='brett'; 

% VOI selection
ndefaults.voi.selectmethod=1;  % =1 for fancy roipoly select, =2 for simpler drawrect, =3 for non-graphical selection
ndefaults.voi.reposition=1;  % =1 to reposition to center of blob, =0 to not reposition,

ndefaults.meg.btitype=0; % =0 for 37ch/74ch Bti system, =1 for ~150+ channel BTi system
ndefaults.meg.single=1; %=1 for keep data single; =0 want double

ndefaults.ctf.gradcorr = 0;

ndefaults.lf.voxelMRIalign=2; % =0, align to MEG head coord grid, =1, align to MRI coord grid, =2, to align to MNI coord grid
ndefaults.lf.numcomp=3; % =3 for 3-component (x,y,z) lead fields; =2 for 2-component (theta,phi) lead fields -- valid for single sphere only

ndefaults.postbc=0; % 1 for post-filtering baseline correction
ndefaults.rzzcov=1; %0 means y*y', =1 means cov(y) (by default)

ndefaults.bf.signalspace=1; % =1 for 'number' to mean 1:number unless specificy num:num ('3' mean 1:3), =0 to exactly specify eig signalspace ('3' means just 3)
ndefaults.bf.viewwaitbar=1;
ndefaults.bf.viewinv=1; %=1 default

% set to approx size of RAM nutmeg is allowed to use (in GB)
% this will be used to guide loading of large multiple-trial datasets
ndefaults.RAM = 1;  % GB

ndefaults.vwr.plotweights=0;

ndefaults.sliders=0;

ndefaults.tfbf.border = 0;  % =1 to draw border around tfbf patches
ndefaults.tfbf.normintens = 0;  % normalize intensities in time-frequency mode of nut_results_viewer, 1=yes, 0=no (default)
ndefaults.tsbf.normintens = 0;  % normalize intensities in time-series mode of nut_results_viewer

ndefaults.FDR = 0.05;      % default false discovery rate (FDR)