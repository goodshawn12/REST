% script NUTMEG - TEST

ft_defaults;

addpath('/i2bm/platform/spm2');
addpath('/neurospin/meg_tmp/ResonanceMeg_Baptiste_2009/SCRIPTS/NUTMEG/scripts');
addpath(genpath('/neurospin/meg_tmp/ResonanceMeg_Baptiste_2009/SCRIPTS/NUTMEG/nutmeg'));

load('/neurospin/meg_tmp/ResonanceMeg_Baptiste_2009/DATA/20100205_vvw/FIELDTRIP/MATLAB2/cond100Run1AllChansNoPCA.mat');
load('/neurospin/meg_tmp/ResonanceMeg_Baptiste_2009/DATA/20100205_vvw/SEG_BRAIN/all.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPARE_SINGLESHELL creates a simple and fast method for the MEG forward
% calculation for one shell of arbitrary shape
cfg               = [];
cfg.smooth        = 5;
cfg.mriunits      = 'mm'; 
cfg.sourceunits   = 'cm'; 
cfg.threshold     = 0.5;
[vol, cfg]        = prepare_singleshell(cfg, segmentedmri);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPARE_LEADFIELD computes the forward model for many dipole locations
% on a regular 2D or 3D grid and stores it for efficient inverse modelling
cfg                 = [];
cfg.channel         = {'all'};
cfg.grad            = cond50Run1Mags.grad ;
cfg.grid.xgrid      ='auto';
cfg.grid.ygrid      ='auto';
cfg.grid.zgrid      ='auto'; 
cfg.grid.resolution = 1;
cfg.vol             = vol;
cfg.reducerank      = 2;
cfg.normalize       = 'no';
cfg.normalizeparam  = 0.5;
[testgrid] = prepare_leadfield(cfg,cond50Run1Mags);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nuts = nut_ft2nuts(cond50Run1Mags,testgrid);

PrevFolder = cd('/neurospin/meg_tmp/ResonanceMeg_Baptiste_2009/DATA/20100205_vvw/NUTMEG');
save mynutmegsession nuts
cd(PrevFolder)
