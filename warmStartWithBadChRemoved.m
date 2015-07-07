
function signal = warmStartWithBadChRemoved(signal)

%% Bad channel detection
% define parameters
min_corr = 0.5;
ignored_quantile = 0.1;
window_len = 2;
max_broken_time = 0.4;
linenoise_aware = true;
rereferenced = false;
protect_channels = [];

% flag channels
if ~exist('removed_channel_mask','var')
    if max_broken_time > 0 && max_broken_time < 1  %#ok<*NODEF>
        max_broken_time = size(signal.data,2)*max_broken_time;
    else
        max_broken_time = signal.srate*max_broken_time;
    end
    
    [nChs,nPts] = size(signal.data);
    window_len = window_len*signal.srate;
    wnd = 0:window_len-1;
    offsets = round(1:window_len:nPts-window_len);
    W = length(offsets);    
    retained = 1:(nChs-ceil(nChs*ignored_quantile));
        
    % optionally ignore both 50 and 60 Hz spectral components...
    if linenoise_aware && signal.srate > 110
        if signal.srate <= 130
            B = design_fir(500,[2*[0 45 50 55]/signal.srate 1],[1 1 0 1 1]);
        else
            B = design_fir(500,[2*[0 45 50 55 60 65]/signal.srate 1],[1 1 0 1 0 1 1]);
        end
        for c=signal.nbchan:-1:1
            X(:,c) = filtfilt_fast(B,1,signal.data(c,:)'); end
    else
        X = signal.data';
    end

    % optionally subtract common reference from data
    if rereferenced
        X = bsxfun(@minus,X,mean(X,2)); end
    
    % for each window, flag channels with too low correlation to any other channel (outside the
    % ignored quantile)
    flagged = zeros(nChs,W);
    for o=1:W
        sortcc = sort(abs(corrcoef(X(offsets(o)+wnd,:))));
        flagged(:,o) = all(sortcc(retained,:) < min_corr);
    end
    
    % mark all channels for removal which have more flagged samples than the maximum number of
    % ignored samples
    removed_channel_mask = sum(flagged,2)*window_len > max_broken_time;
    fprintf('Removing %i channels...',nnz(removed_channel_mask));
    disp({signal.chanlocs(removed_channel_mask==1).labels});
    
    % remove the channels in the protect list
    if ~isempty(protect_channels)
        removed_channel_mask(set_chanid(signal,protect_channels)) = true; end    
end

% annotate the data with what was removed (for visualization)
if ~isfield(signal.etc,'clean_channel_mask')
    signal.etc.clean_channel_mask = true(1,signal.nbchan); end
signal.etc.clean_channel_mask(signal.etc.clean_channel_mask) = ~removed_channel_mask;
signal.etc.badChIndex = find(removed_channel_mask==1);
signal.etc.badChLabels = signal.chanlocs(signal.etc.badChIndex).labels;

%% Whitening
rowmeans = mean(signal.data,2);
data = bsxfun(@minus,signal.data,rowmeans);

[E,D] = eig ( data * data' / nPts );
signal.icasphere = (sqrtm(D)\eye(nChs)) * E';

%% 
