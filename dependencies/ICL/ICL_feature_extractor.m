% extract features for the ICLabel Classifier

function features = ICL_feature_extractor(EEG, psd_percent)
%% check inputs
if ~exist('psd_percent', 'var') || isempty(psd_percent)
    psd_percent = 100;
end
ncomp = size(EEG.icawinv, 2);

% check for ica
assert(isfield(EEG, 'icawinv'), 'You must have an ICA decomposition to use ICLabel')

% calculate ica activations if missing
if isempty(EEG.icaact)
    EEG.icaact = eeg_getdatact(EEG,'component',1:ncomp);
end

% check ica is real
assert(isreal(EEG.icaact), 'Your ICA decomposition must be real to use ICLabel')

% assuming chanlocs are correct
% assuming reference as desired (* maybe enfore CAR? could cause issue with ica)
% assume dipfit already present

%% calc topo
topo = zeros(32, 32, 1, ncomp);
for it = 1:ncomp
    [~, temp_topo, plotrad] = ...
        topoplotFast(EEG.icawinv(:, it), EEG.chanlocs, 'noplot', 'on');
    temp_topo(isnan(temp_topo)) = 0;
    topo(:, :, 1, it) = temp_topo / max(abs(temp_topo(:)));
end

% cast
topo = single(topo);
    
%% calc psd
psd = eeg_rpsd(EEG, 100);

% extrapolate or prune as needed
nfreq = size(psd, 2);
if nfreq < 100
    psd = [psd, repmat(psd(:, end), 1, 100 - nfreq)];
end

% undo notch filter
for linenoise_ind = [50, 60]
    linenoise_around = [linenoise_ind - 1, linenoise_ind + 1];
    difference = bsxfun(@minus, psd(:, linenoise_around), ...
        psd(:, linenoise_ind));
    notch_ind = all(difference > 5, 2);
    if any(notch_ind)
        psd(notch_ind, linenoise_ind) = mean(psd(notch_ind, linenoise_around), 2);
    end
end

% normalize
psd = bsxfun(@rdivide, psd, max(abs(psd), [], 2));

% reshape and cast
psd = single(permute(psd, [3 2 4 1]));

%% format outputs
features = {0.99 * topo, 0.99 * psd};
