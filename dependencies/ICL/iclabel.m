function EEG = iclabel(EEG)

% extract features
disp 'ICLabel: extracting features...'
features = ICL_feature_extractor(EEG, 15);

% run ICL
disp 'ICLabel: calculating labels...'
labels = run_ICL(features{:});

% save into EEG
disp 'ICLabel: saving results...'
EEG.etc.ic_classification.ICLabel.classes = ...
    {'Brain', 'Muscle', 'Eye', 'Heart', 'Line Noise', 'Channel Noise', 'Other'};
EEG.etc.ic_classification.ICLabel.classifications = labels;
