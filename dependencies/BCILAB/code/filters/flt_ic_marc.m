function [signal, state] = flt_ic_marc(varargin)


if ~exp_beginfun('filter') return; end

% requires epoched data, works best on spatially filtered data
declare_properties('name','IC_MARC', ...
    'depends', 'set_makepos', ...
    'follows', 'flt_orica', ...
    'precedes', 'flt_ica_reproject', ...
    'independent_channels', true, ...
    'independent_trials', false);

icmarcpath = fullfile(fileparts(fileparts(which('bcilab'))), 'IC_MARC');

arg_define(varargin,...
    arg_norep({'signal','Signal'}), ...
    arg({'rejectcls','RejectClasses'}, {'blink', 'lateral', 'muscle', 'heart'}, {'blink', 'neural', 'heart', 'lateral', 'muscle', 'mixed'}, 'Which IC classes to reject.'), ...
    arg({'freq', 'updatefreq','UpdateFreq'}, 6, [], 'How many ICs to classify per second.'), ...
    arg({'icmarcpath', 'IC_MARCPath'}, icmarcpath, [], 'Path to the IC_MARC folder.'), ...
    arg_nogui({'state','State'}));

% state.history: past classifications?
if isempty(state)
    % setup virtual chanlocs
    labels = {'Fp1', 'Fp2', 'F3', 'F4', 'C3', 'C4', 'P3', 'P4', 'O1', 'O2', 'F7', 'F8', 'T7', 'T8', 'P7', 'P8', 'Fz', 'Pz', 'FC1', 'FC2', 'CP1',...
        'CP2', 'FC5', 'FC6', 'CP5', 'CP6', 'F9', 'F10', 'TP9', 'TP10', 'Cz', 'Oz', 'F1', 'F2', 'C1', 'C2', 'P1', 'P2', 'AF3', 'AF4',...
        'FC3', 'FC4', 'CP3', 'CP4', 'PO3', 'PO4', 'F5', 'F6', 'C5', 'C6', 'P5', 'P6', 'AF7', 'AF8', 'FT7', 'FT8', 'TP7', 'TP8', 'PO7', ...
        'PO8', 'AFz', 'FCz', 'CPz', 'POz'};
    virtual_chanlocs_struct = struct('labels', lower(labels));
    virtual_chanlocs = pop_chanedit(virtual_chanlocs_struct,'lookup', 'standard-10-5-cap385.elp');
    for i=1:length(virtual_chanlocs)
        % set head radius to 9 cm
        virtual_chanlocs = pop_chanedit(virtual_chanlocs, ...
            'changefield',{i 'sph_radius' '9'},'convert',{'sph2all'});
    end
    
    % intialize state
    state.virtual_chanlocs = virtual_chanlocs;
    state.cdn_matrix = load(fullfile(icmarcpath, 'spatial2', 'dipolfit_matrix.mat')); 
    state.model = load(fullfile(icmarcpath, 'spatial2.mat'));
    state.class_labels = {'blink', 'neural', 'heart', 'lateral', 'muscle', 'mixed'};
    state.reject = find(ismember(state.class_labels, rejectcls));
    state.cls = zeros(signal.nbchan, 1);
    state.s = 0;
    state.next_ic = 1;
end

% parameters
update_s = max(round(signal.srate / freq), 1);
Winv = [];
if isfield(signal, 'smax')
    smax = signal.smax;
else
    smax = signal.pnts;
end

% run eyeCatch
while smax - state.s >= update_s
    
    % update classifications
    if isempty(Winv)
        Winv = pinv(signal.icaweights * signal.icasphere);
    end
    [predclass, predprob] = runIC_MARC(state, signal, Winv);
    % TODO: figure out how to save posterior probs to a buffer
    % ideas: save to buffer (bad because buffer is not BCILAB compliant...)
    %        save to state (possible, but need to predetermine buffer size)
    
    % update state
    state.cls(state.next_ic) = predclass;
    state.s = state.s + update_s;
    state.next_ic = mod(state.next_ic, signal.nbchan) + 1;
    
end

% update signal
if ~isempty(Winv)
    signal.icawinv = Winv;
end
signal.reject = ismember(state.cls, state.reject);

exp_endfun;

end


% run IC_MARC classifier
function [predclass, predprob] = runIC_MARC(state, signal, icawinv)
% should pass the whole eeg signal for location and in case there are some
% channels being removed.

% standardize icawinv
icawinv = zscore(icawinv(:, state.next_ic));
topog = icawinv';

% interpolate / extrapolate channel activation for 64-ch setting
thetas = cell2mat({state.virtual_chanlocs.sph_theta});
phis =  cell2mat({state.virtual_chanlocs.sph_phi});
sigmas = repmat(0.5, 64,1);
head_radius=9;

n_ics = size(icawinv,2);
virtual_topog = NaN(n_ics,64);

for ic=1:n_ics
    activations = virtual_electrode_activation(phis, thetas, sigmas, signal.chanlocs, topog(ic, :), head_radius);
    virtual_topog(ic,:) = activations;    
end

virtual_topography = zscore(virtual_topog,[],2);

% features extraction using spatial2
features = NaN(size(icawinv,2),12);

% 1st feature: spatial discontinuity
features(:,1) = localized_discontinuity_measure(virtual_topography, state.virtual_chanlocs, n_ics); 

% 2nd feature: SED (spatial eye difference) - lateral difference 
features(:,2) = computeSED_NOnorm_variable_ics_light(virtual_topography,state.virtual_chanlocs,size(virtual_topography, 2), size(virtual_topography, 1));

% 3rd to 6th features: degree of activations in different areas: [front, post, leftarea, rightarea]
[features(:,3), features(:,4), features(:,5), features(:,6)] = scalpmap_features_light(topog, signal.chanlocs, state.virtual_chanlocs);

% 7th feature: median value of icawinv'
features(:,7) = abs(median(topog,2)); % abs_med

% 8th feature: current density norm
features(:,8) = current_density_norm_light(virtual_topography, state.virtual_chanlocs, state.cdn_matrix);

% 9th feature: 2D spatial FFT
features(:,9) = calc_2ddft_light(virtual_topography, state.virtual_chanlocs);

% 10th feature: log of IC topographies
features(:,10) = log_range_spatial_light(virtual_topography);

% 11th feature: log of the 2-norm of the distance between the two 
% electrodes with the highest and lowest activations
features(:,11) = spatial_distance_extrema_light(virtual_topography, state.virtual_chanlocs);

% 12th feature: spatial entropy of scalpmap
features(:,12) = scalp_entropy_light(virtual_topography);

% standardize features 
featsnorm = features - repmat(state.model.mu, size(features,1),1);
featsstand = featsnorm./repmat(state.model.sigma, size(features,1),1);

% classify the IC into 6 different classes: 'blink', 'neural', 'heart', 'lat', 'muscle', 'mixed'
predprob = mnrval(state.model.mod, featsstand);
[~, predclass] = max(predprob, [], 2); 

end
