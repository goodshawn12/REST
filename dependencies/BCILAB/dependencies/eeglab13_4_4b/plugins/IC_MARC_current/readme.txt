Readme file for IC_MARC
------------------------------------------------------------
This file contains information on the feature extraction functions found in IC_MARC (classification of Independent Components of EEG into Multiple ARtifact Classes). The functions are most easily used through the main function ic_feature_extraction.m, which performs preprocessing and then calls the functions that calculate features, or by placing the folder IC_MARC in the dependencies folder in EEGLab. If placed in the dependencies folder in EEGLab, IC_MARC can be run as an EEGLab plugin via a menu item in the tools menu.

Please contact laura dot frolich at gmail dot com if errors occur in the code,  if you have suggestions for improvements, if you have problems using the code, or if you have other comments.

****************************
Special fields in common spatial and temporal space
****************************
The feature extraction functions require some special fields to be present in the EEGLab data structure given as input. These fields are calculated and added by the main function ic_feature_extraction. The fields contain representations of the filtered and downsampled time series and downsampled scalp maps. This ensures that ICs from different studies can be directly compared since they lie in a common space, and that a classifier will generalize across studies more easily. The fields concerned are the following:

virtual_srate: Sampling rate that IC time series are downsampled to 200 Hz.

virtual_pnts: Number of points in the downsampled IC time series.

icaact_unfiltered: IC time series downsampled to virtual_srate (200Hz) (and filtered by the function resample before downsampling to avoid aliasing)

icaact_filtered_resampled: IC time series that have both been downsampled to 200Hz and bandpass filtered between 3Hz and 90Hz.

virtual_topography: Matrix containing downsampled scalp maps of ICs in rows. See the function virtual_topography for details.

virtual_chanlocs: The channel locations corresponding to the downsampled scalp maps, the standard 64 electrode 10-20 system placements.

virtual_nbchan: Number of channels in virtual_chanlocs (64).

****************************
Dependencies
****************************
Almost all functions rely on functionality in EEGLab (A. Delorme and S. Makeig).

The function hurst_exponents requires the Matlab Wavelet toolbox.

The functions spatially_normalize_icdecomp, temporally_normalize_icdecomp, and virtual_topography require the function zscore in the Matlab Statistics toolbox.
****************************
Data files in this folder
****************************
Multinomial regression models using the three feature sets discussed in the accompanying paper were trained on emotion data. The models are contained in the following .mat files. The mean and standard deviation of the data used to train the models are also included. New data should have this mean subtracted and be divided by the standard deviation before passing feature values to the classifier.

established_features.mat: The variable mod contains the multinomial regression model trained using features in the feature set established in the accompanying paper. Mean and standard deviations of these features are included as mu and sigma, respectively.

spatial_established_features.mat: The variable mod contains the multinomial regression model trained using spatial features selected often in the leave-one-subject-out cross-validation described in the accompanying paper. Mean and standard deviations of these features are included as mu and sigma, respectively. This feature set consists of the following feaures:
GD, SED, var_front, var_back, vertical_polarity, central, xcoord, zcoord, dipole_residual_variance, central_activation, spat_dist_extrema, scalp_entropy.

spatial2.mat: The variable mod contains the multinomial regression model trained using the features calculated with the argument value {'spatial2'} for the argument 'features'. These features do not contain dipole fit features, which are computationally demanding due to the dipole fit. This feature set was optimized heuristically by hand. Mean and standard deviations of these features are included as mu and sigma, respectively. This feature set consists of the following features:
GD, SED, frontal, posterior, left, right, abs_med_topog, cdn, 2ddft, log_range_spatial, spat_dist_extrema, scalp_entropy'

Files in ic_feature_extraction_functions needed by the feature extraction functions:
dipolfit_matrix.mat: The function current_density_norm requires data contained in dipolfit_matrix.mat, published with (I. Winkler et al.), to be present in the same directory as current_density_norm.m.

standard_BESA.mat: The file standard_BESA.mat from EEGLab (A. Delorme and S. Makeig) is used by the functions dipfit_features and dipole_residual_variance.

standard-10-5-cap385.elp: The file standard-10-5-cap385.elp from EEGLab (A. Delorme and S. Makeig) is used by the functions calc_central_and_border_activation, calc_central_and_border_activation, scalpmap_features, dipole_residual_variance, and virtual_topography.

****************************
Features not included during model construction
****************************
Some of the features that can be calculated by functions in this folder were not included in the initial pool of features on which analyses in the paper "Classification of independent components of EEG into multiple artifact classes" were based. The concerned features are the following:

mean activations around eyes from, from which SED is calculated: These two features were not included since they are not linearly independent when SED is also included as a feature.

mean activations in the frontal and posterior scalp regions, from which SAD is calculated: These two features were not included since they are not linearly independent when SAD is also included as a feature.

variance and mean of skewness over 15 second intervals: Since the recordings in one of our data sets were only four seconds long, these two features could not be calculated. It was one of the six final features in the feature set found in (I. Winkler et al.).

variance and mean of variance over 15 second intervals: Since the recordings in one of our data sets were only four seconds long, these two features could not be calculated.

****************************
References
****************************
I. Winkler et al.:  Irene Winkler and Stefan Haufe and Michael Tangermann (2011), "Automatic Classification of Artifactual ICA-Components for Artifact Removal in EEG Signals", Behavioral and Brain Functions, 7:30 doi:10.1186/1744-9081-7-30, http://www.behavioralandbrainfunctions.com/content/7/1/30/

A. Delorme and S. Makeig: Arnaud Delorme and Scott Makeig (2004) "EEGLAB: an open source toolbox for analysis of single-trial EEG dynamics", Journal of Neuroscience Methods 134:9-21
