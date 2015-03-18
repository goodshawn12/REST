% FCM  the Functional Connectivity Map toolbox of NUTMEG
% Overview of the different functions and graphical user interfaces (GUIs):
%
% Command                         Comments
% -------                         --------
%
% GUI
% ---
% fcm_gui                         The GUI starting point to all functions.
% fcm_across_subject_gui          GUI for across subject stats (P-image, correlations)
% fcm_roiselect_gui               GUI for selection of anatomical ROIs
%
% Command line usage
% ------------------
% Before starting the FCM toolbox:
% 1. nutmeg                       Start NUTMEG and read in coregistration and data, calculate 
%                                 lead field.
% 2. nut_tfbf_gui  or  tfbf       Run a single-state time-frequency beamformer to obtain the 
%                                 source dipole weights. 
%                                 Make sure to save the weight and filtered data outputs.  
%
% FCM workflow:
%  [Set config]
% 1. fcm_start                    Sets the FCM configuration: seed voxels, connections, output, etc.        
% 2. fcm_anatroi                  (optional) Select seed regions of interest
%    fcm_paintroi
%    fcm_functroi
%    fcm_roiidx
% 3. fcm_conndef                  Defines the connections between voxels to be analyzed.
% 
%  [Calculate functional connectivity on a single computer]
% 4. fcm_sourceconnectivity       Calculates functional connectivity metrics for voxels
% 
%  [Calculate functional connectivity on a linux cluster]
% 4a. fcm_preparejobs             Creates text files with job descriptions for Linux cluster. 
%                                 Output is stored in a new directory named "comps".
% 4b. sccohere1/3/31, spli1/3/31, sampcorr1/3/31  
%                                 Calculate source connectivity on linux cluster (compiled binary)
% 4c. fcm_assemble                Assembles the cluster data 
%
%  [Calculate and visualize maps]
% 5. fcm_comcoh2beam              Creates maps from complex coherence results.
% 6. nut_results_viewer           View maps.
%
%  [Across subject analyses]
% - fcm_beam2Pimage               Calculates "P-images".
% - fcm_corr                      Correlates maps with behavioral/clinical variables.
%
% See also  TFBF, FCM_START, FCM_ANATROI, FCM_PAINTROI, FCM_FUNCTROI, FCM_ROIIDX
% FCM_CONNDEF, FCM_PREPAREJOBS, FCM_SOURCECONNECTIVITY, FCM_ASSEMBLE, FCM_COMCOH2BEAM,
% NUT_RESULTS_VIEWER, FCM_BEAM2PIMAGE, FCM_CORR

help fcm