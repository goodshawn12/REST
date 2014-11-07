% testing script for artifact reduction
% 2014.9.23

% load eegdata and assign ica results
EEGtemp = pop_loadset('filename','Flanker_SH_cleaned_test_icainfo.set','filepath','D:\\Matlab Coding\\OnlineICA\\BCILAB\\userdata\\');
useRUNICA = 1;
useORICA = 1;

if useRUNICA
    EEGrun = EEGtemp;
    EEGrun.icawinv = EEGrun.icawinv_true;
    EEGrun.icaweights = inv(EEGrun.icawinv);
    EEGrun.icasphere = eye(EEGrun.nbchan);
    EEGrun = eeg_checkset(EEGrun);
end
if useORICA
    EEGor = EEGtemp;
    load flanker_b1w8.mat;
    % averaging over last minute
    window = 1; % min
    avg_start = find(orica.time>(20-window)*60*300,1);
    [nChs,~,nPts] = size(orica.icaweights);
    icaweights_all = zeros(nChs);
    for it = avg_start:nPts
        icaweights_all = icaweights_all + orica.icaweights(:,:,it) * orica.icasphere(:,:,it);
    end
    EEGor.icaweights = icaweights_all ./ (nPts - avg_start + 1);
    EEGor.icasphere = eye(EEGor.nbchan);
    EEGor.icawinv = inv(EEGor.icaweights);
    EEGor = eeg_checkset(EEGor);
    clear orica
end


% select and reject components
if useRUNICA
    rejcomp = [18 26 40];
    EEGrun = pop_subcomp(EEGrun, rejcomp, 0);
    EEGrun = eeg_checkset(EEGrun);
end

if useORICA
%     pop_selectcomps(EEGor,rejcomp);
    rejcomp = [4 7 33];
    EEGor = pop_subcomp(EEGor, rejcomp, 0);
    EEGor = eeg_checkset(EEGor);
end


eegplotoptions = { 'winlength', 5, 'events', EEGtemp.event };

eegplot( EEGtemp.data, 'srate', EEGtemp.srate, 'title', 'Scroll channel activities -- eegplot()', ...
    'limits', [EEGtemp.xmin EEGtemp.xmax]*1000 , eegplotoptions{:}, 'data2' , EEGrun.data);

eegplot( EEGtemp.data, 'srate', EEGtemp.srate, 'title', 'Scroll channel activities -- eegplot()', ...
    'limits', [EEGtemp.xmin EEGtemp.xmax]*1000 , eegplotoptions{:}, 'data2' , EEGor.data);

% compare runica and orica results
% it seems that ORICA reconstructed signals are even cleaner
eegplot( EEGrun.data(1:10,:), 'srate', EEGrun.srate, 'title', 'Scroll channel activities -- eegplot()', ...
    'limits', [EEGrun.xmin EEGrun.xmax]*1000 , eegplotoptions{:}, 'data2' , EEGor.data(1:10,:));


% % back projection
% component_keep = setdiff(1:size(EEG.icaweights,1), components);
% compproj = EEG.icawinv(:, component_keep)*eeg_getdatact(EEG, 'component', component_keep, 'reshape', '2d');
% compproj = reshape(compproj, size(compproj,1), EEG.pnts, EEG.trials);
