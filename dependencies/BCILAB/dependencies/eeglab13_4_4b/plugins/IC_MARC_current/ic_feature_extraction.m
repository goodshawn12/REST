function [features, outeeg, sp_icas, feature_names, function_timing] = ic_feature_extraction(varargin)
% [features, outeeg, sp_icas, feature_names, function_timing] = ic_feature_extraction(varargin)
% Calculate features of ICs.
%
% Input:
% eeg: EEGLab data structure. Mandatory input.
%
% Optional input:
% features: The features to calculate (default: {'all'}). Possible values:
%   Vector of integers between 1 and 73, denoting the 73 features that can
%   be calculated by this function.
%   Cell array of vectors of integers between 1 and 73 and/or one or more
%   of the following strings:
%   'all', 'All': All features are calculated
%   'established_feature_set', 'established_features': Calculate the features in the established
%       feature set consisting of 14 features as described in the
%       accompanying paper.
%   'established_spatial_features': Calculate the features in the
%       established spatial feature set (described in the
%       accompanying paper.)
%   'spatial2': Calculate a set of spatial features without the
%   computationally demanding dipole features. This feature set was not
%   optimized systematically, but by some fiddling by hand.
%
% force_preprocessing: Force preprocessing steps to be performed even if
% unnecessary. If no temporal or spectral features are requested, it is not
% necessary to downsample and filter time series. However, if the outputs
% produced during preprocessing are desired, force_preprocessing can be set
% to true. (default: false)
%
% varargin: Cell array of named parameter/value pairs to pass to functions
% called by ic_feature_extraction. All parameters are passed on to all
% called functions. Each function uses the parameters with names that it
% takes as input and ignores all other parameters.
%
% Output:
% features: Calculated feature values.
%
% outeeg: The EEGlab data structure given as input, with the dipole fit found
% in this function and virtual fields.
%
% sp_icas: A cell array which contains, for each IC, a vector, f, giving the
% midpoints of frequency bands and a matrix, pxx, giving   the band power
% in one second intervals in the theta (4-7Hz), alpha
% (8-13Hz), beta (14-20Hz), low gamma (21-30Hz), medium gamma (30-45Hz),
% power grid voltage (46-65Hz), and high gamma (66-80Hz) bands. Time moves
% over columns and each row corresponds to a frequency band. The first row
% is the theta band and the last row the high gamma band.
%
% feature_names: Names of the features calculated.
%
% Examples of calls of ic_feature_extraction:
%
% Calculate only the first feature (GD)
% [ic_feats eegout speccomp feature_names] = ic_feature_extraction(eeg, 1);
%
% Calculate only the first feature (GD), but force temporal and spectral
% preprocessing
% [ic_feats eegout speccomp feature_names] = ic_feature_extraction(eeg, 1, 'force_preprocessing', true);
%
% Calculate the feature set established in the accompanying paper
% [ic_feats eegout speccomp feature_names] = ic_feature_extraction(eeg, {'established_feature_set'});
%
% Calculate the feature set established in the accompanying paper as well
% as features 68 and 69, and pass on the path to the file standard_BESA.mat
% to the feature extraction functions.
% [ic_feats eegout speccomp feature_names] = ic_feature_extraction(eeg, {'established_feature_set', 68:69},...
% 'varargs', {'besapath', '~/Documents/bcilab-1.0/dependencies/eeglab_10_0_1_0x/plugins/dipfit2.2/standard_BESA/'});

% Copyright (C) 2013  Laura Froelich (laura dot frolich at gmail dot com)
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.




if ~usejava('jvm')
    error([mfilename ' requires Java to run.']);
end

mfile_location=mfilename('fullpath');
ICMARC_directory = mfile_location(1:(strfind(mfile_location, 'ic_feature_extraction')-1));
frombcilab_location = [ICMARC_directory 'private' filesep 'from_bcilab'];
addpath(frombcilab_location) % the folder private/from_bcilab contains functions from BCILab that handle passing
% of arguments between functions.


arg_define(2,varargin, ...
    arg({'eeg'},[],[],'EEGLab struct containing EEG data.', 'cat', 'Core Parameters'), ...
    arg({'features'}, {'all'}, [],...
    'Features to calculate.','cat','Core Parameters'), ...
    arg({'force_preprocessing'},false,[true, false],'Force pre-processing steps even if unnecessary.'),...
    arg('varargs', {},[], 'Arguments to functions called by ic_feature_extraction.'));
% save the data set that was given as input so that it can be returned
% unaltered, except for added virtual fields
outeeg = eeg;

talairachjar_location = [ICMARC_directory 'talairach.jar'];

% save globals to local space since these can be cleared by the javaaddpath
% call

baseVarsNames=evalin('base', 'who()');
globalVarsNames = who('global');
globalVarsArray = cell(1, numel(globalVarsNames));
for iVar = 1:numel(globalVarsNames)
    %eval(sprintf('global %s', globalVars{iVar}));
    eval(sprintf(['global ' globalVarsNames{iVar} ';']));
  eval(sprintf(['globalVarsArray{' num2str(iVar) '}=' globalVarsNames{iVar} ';']));
  
end

javaaddpath(talairachjar_location ) % for add_anatomical_dipoles

% resave your variables to global space
for iVar = 1:numel(globalVarsNames)
    eval(sprintf(['global ' globalVarsNames{iVar}]));
  eval(sprintf([globalVarsNames{iVar} '= globalVarsArray{' num2str(iVar) '};' ]));
  if sum(strcmp(globalVarsNames{iVar}, baseVarsNames))
      evalin('base', sprintf(['global ' globalVarsNames{iVar}]))
  end
end

feature_names = {
    % Start of spatial features
    'GD',... % 1
    'SED',...
    'mean_left',...
    'mean_right',...
    'SAD',... % 5
    'var_front',...
    'var_back',...
    'mean_front',...
    'mean_back',...
    'lateral_eyes',... % 10
    'vertical_polarity',...
    'lefteye',...
    'righteye',...
    'frontal',...
    'central',... % 15
    'posterior',...
    'left',...
    'right',...
    'abs_med_topog',...
    'cdn',... % 20
    'xcoord',...
    'ycoord',...
    'zcoord',...
    'ndipole_labels',...
    'dipole_residual_variance',... % 25
    '2ddft',...
    'central_activation',...
    'border_activation',...
    'log_range_spatial',...
    'spat_dist_extrema',... % 30
    'scalp_entropy',...
    ... % start of spectral features
    'theta',...
    'alpha',...
    'beta',...
    'gamma',... % 35
    'gammamed',...
    'gammaelec',...
    'gammah',...
    'vartheta',...
    'varalpha',... % 40
    'varbeta',...
    'vargamma',...
    'vargammamed',...
    'vargammaelec',...
    'vargammah',... % 45
    'spectral_entropy_mean',...
    'spectral_entropy_var',...
    'low_frequent_power_mean',...
    'low_frequent_power_var',...
    ... % start of temporal features
    'avgskew1s_mean',... % 50
    'avgskew15s_mean',...
    'avgskew1s_var',...
    'avgskew15s_var', ...
    'log_range_temporal_mean',...
    'log_range_temporal_var',... % 55
    'kurtosis_mean',...
    'kurtosis_var', ...
    'hurst1_mean',...
    'hurst2_mean',...
    'hurst3_mean',... % 60
    'hurst1_var',...
    'hurst2_var',...
    'hurst3_var', ...
    'var1s_mean',...
    'var1s_var',... % 65
    'var15s_mean',...
    'var15_var',...
    'max_first_deriv_mean',...
    'max_first_deriv_var',...
    'max_ampl_mean',... % 70
    'max_ampl_var',...
    'time_entropy_mean',...
    'time_entropy_var'};

functions = {...
    'localized_discontinuity_measure_wrapper',... %1
    'computeSED_NOnorm_wrapper',...
    'computeSED_NOnorm_wrapper',...
    'computeSED_NOnorm_wrapper',...
    'computeSAD_wrapper',... %5
    'computeSAD_wrapper',...
    'computeSAD_wrapper',...
    'computeSAD_wrapper',...
    'computeSAD_wrapper',...
    'scalpmap_features',... %10
    'scalpmap_features',...
    'scalpmap_features',...
    'scalpmap_features',...
    'scalpmap_features',...
    'scalpmap_features',... %15
    'scalpmap_features',...
    'scalpmap_features',...
    'scalpmap_features',...
    'scalpmap_features',...
    'current_density_norm',... % 20
    'dipfit_features',...
    'dipfit_features',...
    'dipfit_features',...
    'dipfit_features',...
    'dipole_residual_variance',... %25
    'calc_2ddft',...
    'calc_central_and_border_activation_wrapper',...
    'calc_central_and_border_activation_wrapper',...
    'log_range_spatial',...
    'spatial_distance_extrema',... %30
    'scalp_entropy',...
    'mean_log_band_power',...
    'mean_log_band_power',...
    'mean_log_band_power',...
    'mean_log_band_power',... %35
    'mean_log_band_power',... % 36 gammamed
    'mean_log_band_power',... % 37 gamma_elec
    'mean_log_band_power',...
    'mean_log_band_power',...
    'mean_log_band_power',... %40
    'mean_log_band_power',...
    'mean_log_band_power',...
    'mean_log_band_power',...
    'mean_log_band_power',... % 44 var gammamed
    'mean_log_band_power',... % 45 var gamma_elec
    'spectral_entropy',...
    'spectral_entropy',...
    'low_frequent_power',...
    'low_frequent_power',...
    'avg_skew',... %50
    'avg_skew',...
    'avg_skew',... % 15 second intervals
    'avg_skew',...
    'log_range_temporal',... % 1 second intervals
    'log_range_temporal',... %55
    'kurtosis_temporal',...
    'kurtosis_temporal',...
    'hurst_exponents',...
    'hurst_exponents',...
    'hurst_exponents',... %60
    'hurst_exponents',...
    'hurst_exponents',...
    'hurst_exponents',...
    'local_variance_1s',...
    'local_variance_1s',... % 65
    'local_variance_15s'...
    'local_variance_15s',...
    'max_first_derivative',...
    'max_first_derivative',...
    'max_amplitude',... % 70
    'max_amplitude',...
    'time_entropy',...
    'time_entropy',...
    };



run_functions = containers.Map(functions, repmat({0},1,length(functions)));
feature_calculators = containers.Map(...
    arrayfun(@(x){x}, 1:length(functions)),...
    functions);

if ( (strcmp(class(features), 'double') && all(round(features)==features)) || iscell(features))
    features_to_calculate=parse_feature_specification(features, uint32(feature_calculators.Count));
else
    error(['ic_feature_extraction: Features to be calculated must be passed as a vector of integers, ',...
        'or a cell array containing integers or names of features'])
end

if max(features_to_calculate) > feature_calculators.Count
    error(['ic_feature_extraction: The following numbers do not have corresponding features defined: '...
        num2str(features_to_calculate(features_to_calculate > feature_calculators.Count))])
end

% we only need the channels that were used in the ICA decomposition
eeg = pop_select(eeg, 'channel', eeg.icachansind); 
eeg= check_ica_presence(eeg, 'do_icaact', false);
for i=1:length(eeg.chanlocs)
    % set head radius to 9 cm
    eeg=pop_chanedit(eeg, 'changefield',{i 'sph_radius' '9'},'convert',{'sph2all'});
end

celleegs=[];
sp_icas = [];
eeg.virtual_srate = 200;
eeg.virtual_pnts = [];
eeg.icaact_unfiltered = [];
eeg.icaact_filtered_resampled = [];


if max(features_to_calculate)>31 || force_preprocessing % only perform temporal and spectral
    % preprocessing if temporal or spectral features are requested, or if
    % forced to
    if eeg.trials ==1 && isfield(eeg.event, 'type')% if the input data is continuous and has events, check for boundary events
        boundaries = cellfun(@(x) strcmp(x, 'boundary'), {eeg.event.type});
        if any(boundaries)
            boundary_latencies = ceil(cell2mat({eeg.event(boundaries).latency})); % boundary latencies are between
            % data points, but we need locations of actual data points when
            % splitting the data. The first sample of new data sections
            % are contained in boundary_latencies
            if boundary_latencies(1)~=1
                edgelats = [1 boundary_latencies];
            else
                edgelats = boundary_latencies;
            end
            if boundary_latencies(end)~=eeg.pnts
                edgelats = [edgelats eeg.pnts];
            end
            eeg.event = []; % remove events from temporary EEG structure
            eegcounter = 0;
            usedpnts = 0;
            unusedpnts = 0;
            for ilat = 1:length(edgelats)-1
                interval_pnts = diff([edgelats(ilat) (edgelats(ilat+1)-1)]);
                if interval_pnts > eeg.srate % the part of data 
                    %to be cut out must be longer than one second
                    eegcounter = eegcounter +1;
                    cureeg = pop_select(eeg, 'point', [edgelats(ilat) (edgelats(ilat+1)-1)]);
                    celleegs{eegcounter}= spectral_temporal_preprocessing(cureeg);
                    usedpnts = usedpnts + interval_pnts;
                else
                    unusedpnts = unusedpnts + interval_pnts;
                end
            end
            if isempty(celleegs)
                error(['ic_feature_extraction.m: Temporal and spectral features, '...
                    'which were requested, cannot be calculated since ' ...
                    'no continuous data parts between two boundary events were '...
                    'longer than one second in the data provided.'])
            end
            if length(edgelats)-1 > length(celleegs)
                warning(['ic_feature_extraction.m: Some parts of the continuous data '...
                    'between two boundary events were less '...
                    'than one second (' num2str(100*unusedpnts/(unusedpnts + usedpnts), '%0.3f')...
                    '% of data). Such data will not be used for calculating '...
                    'temporal and spectral features.'])
            end
        else
            eeg= spectral_temporal_preprocessing(eeg);
        end
    else
        eeg= spectral_temporal_preprocessing(eeg);
    end    
    sp_icas = calculate_spectra_wrapper(eeg, 'celleegs', celleegs);
end

spatially_normalized_eeg = spatially_normalize_icdecomp(eeg);

% Safe to remove channels since IC activations have already been
% calculated. The channels are NOT removed in the dataset that is returned,
% but only from the temporary copy used by the feature extraction
% functions.

eeg=remove_channels_wo_locations(eeg); % when using old pop_select code, 
% this call does not modify the values in eeg.icawinv, only removes those
% for the channel(s) removed. With the new pop_select, the values are
% modified
eeg=virtual_topography(eeg);

if ~isempty(celleegs)
    eeg.data = []; % save memory since all data is now stored in celldata
    eeg.icaact = [];
end

if any(strcmp(varargs, 'besapath'))
    besaargind = find(strcmp(varargs, 'besapath'))+1;
    besaarg = varargs{besaargind};
    if strcmp(besaarg(1), '~')
        if ispc; userdir= getenv('USERPROFILE');
        else
            userdir= getenv('HOME');
        end
        besaarg(1)=[];
        besaarg = [userdir besaarg];
        varargs{besaargind} = besaarg;
    end
end


varargs = {varargs{:}, 'topography', spatially_normalized_eeg.icawinv',...
    'sp_icas', sp_icas, 'celleegs', celleegs};

feature_names = feature_names(features_to_calculate);
features = NaN(size(eeg.icawinv,2),length(features_to_calculate));
timing = NaN(1, length(functions)); % added 5/2/2014
nfunc=0; % added 5/2/2014
for ifeat = 1:length(features_to_calculate)
    func = feature_calculators(features_to_calculate(ifeat));
    if ~run_functions(func); % only run the function if the function
        % that calculates the current features has not already been run
        nfunc=nfunc+1; % added 5/2/2014
        
        fhandle = eval(['@' func]);
        noutputs = nargout(fhandle); % get the number of outputs from the
        % current function
        
        % run the current function and save all its outputs
        tic % added 5/2/2014
        [mytemp(1:noutputs).output] = feval(fhandle, eeg, varargs{:});
        timing(nfunc)=toc; % added 5/2/2014
        
        % rename the saved outputs to "(function name)_result"
        eval([func '_result' '= mytemp;']);
        
        % if a function performing dipole fits was called, overwrite the
        % eeg data structure with the eeg data structure containing the
        % dipole fit from the called function
        if strcmp(func, 'dipfit_features') || strcmp(func, 'dipole_residual_variance')
            eeg = eval([func '_result(noutputs).output']);
        end
        
        run_functions(func) = true;
    end
    
    nthoutput = sum(strcmp(functions(1:int16(features_to_calculate(ifeat))), func));
    features(:,ifeat) = eval([func '_result(' num2str(nthoutput) ').output']);
end

outeeg.virtual_srate = eeg.virtual_srate;
outeeg.virtual_pnts = eeg.virtual_pnts;
outeeg.icaact_unfiltered = eeg.icaact_unfiltered;
outeeg.icaact_filtered_resampled = eeg.icaact_filtered_resampled;
outeeg.virtual_topography = eeg.virtual_topography;
outeeg.virtual_chanlocs = eeg.virtual_chanlocs;
outeeg.virtual_nbchan = eeg.virtual_nbchan;
outeeg.dipfit = eeg.dipfit;
function_timing.timing = timing; % added 5/2/2014
function_timing.functions = functions; % added 5/2/2014

end

function feature_vector = parse_feature_specification(features, nfeatures)
% nfeatures is the number of features that can be calculated

feature_vector = [];

establishedfeatureset = [2 5 15 20 23 25 27 29 30 32 48 56 64 72];
feature_map = containers.Map(...
    {'all', 'All','established_features',  'established_feature_set', 'established_spatial_features', 'spatial2'},...
    {1:nfeatures,1:nfeatures, establishedfeatureset, ...
    establishedfeatureset,...
    [1:2, 6:7, 11, 15, 21, 23, 25, 27, 30:31], [1:2 14 16:20 26 29:31]});

if strcmp(class(features), 'double')
    features = {features};
end

for ifeat = 1:length(features)
    current_key = features{ifeat};
    if strcmp(class(current_key), 'double')
        feature_vector = [feature_vector, current_key];
    else
        feature_vector = [feature_vector, feature_map(current_key)];
    end
end

feature_vector = sort(unique(feature_vector));
end
