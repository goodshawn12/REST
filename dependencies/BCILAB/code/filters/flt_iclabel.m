function [signal, state] = flt_iclabel(varargin)


if ~exp_beginfun('filter') return; end

% requires epoched data, works best on spatially filtered data
declare_properties('name','ICLabel', ...
    'depends', 'set_makepos', ...
    'follows', 'flt_orica', ...
    'precedes', 'flt_ica_reproject', ...
    'independent_channels', true, ...
    'independent_trials', false);

iclpath = fullfile(fileparts(fileparts(which('bcilab'))), 'ICLabel');

arg_define(varargin,...
    arg_norep({'signal','Signal'}), ...
    arg({'rejectcls','RejectClasses'}, {'Muscle', 'Eye', 'Heart', 'Line Noise', 'Chan Noise',}, {'Brain', 'Muscle', 'Eye', 'Heart', 'Line Noise', 'Chan Noise', 'Other'}, 'Which IC classes to reject.'), ...
    arg({'freq', 'updatefreq','UpdateFreq'}, 6, [], 'How many ICs to classify per second.'), ...
    arg({'iclpath', 'ICLabelPath'}, iclpath, [], 'Path to the ICLabel folder.'), ...
    arg({'psd_sec2samp', 'PSDSecondsToSample'}, 10, [], 'Time in seconds to consider for PSD estimate.'), ...
    arg({'psd_sec4fft', 'PSDSecondsForFFT'}, 1, [], 'Time in seconds to use when computing fft.'), ...
    arg({'psd_poverlap', 'PSDPercentOverlap'}, 0.5, [], 'Percent window overlap for PSD estimate.'), ...
    arg_nogui({'state','State'}));
%     arg({'iclversion', 'ICLabelVersion'}, 'lite', {'default', 'lite', 'beta'}, 'Path to the ICLabel folder.'), ...

% state.history: past classifications?
if isempty(state)
    
    % initialize matconvnet
    w = warning;
    warning('off', 'all');
    vl_setupnn
    warning(w);
    
    % intialize state
    netStruct = load(fullfile(iclpath, 'netICL_lite.mat'));
    state.net = dagnn.DagNN.loadobj(netStruct);
    state.class_labels = {'Brain', 'Muscle', 'Eye', 'Heart', 'Line Noise', 'Chan Noise', 'Other'};
    state.reject = find(ismember(state.class_labels, rejectcls));
    state.cls = zeros(signal.nbchan, 1);
    state.s = 0;
    state.next_ic = 1;
    
    % psd settings
    nwindows = floor((psd_sec2samp / psd_sec4fft - 1) / (1 - psd_poverlap)) + 1;
    winlen = round(signal.srate * psd_sec4fft);
    n_ic = size(signal.icaweights, 1);
    state.psd = struct('nwindows', nwindows, 'sec2samp', psd_sec2samp, ...
        'sec4fft', psd_sec4fft, 'overlap', psd_poverlap, ...
        'unoverlapped', round(winlen * (1 - psd_poverlap)), ...
        'winlen', winlen, 'window', windows('hamming', winlen, 0.54)', ...
        'ffts', nan(floor(winlen/2), nwindows, n_ic), 'winds', nan(nwindows, n_ic), ...
        'buffer', zeros(n_ic, round(signal.srate * psd_sec2samp)), ...
        'bufferind', 1, 'bufferlen', round(signal.srate * psd_sec2samp));
end

% parameters
update_s = max(round(signal.srate / freq), 1);
Winv = [];
if isfield(signal, 'smax')
    smax = signal.smax;
else
    smax = signal.pnts;
end

% run ICLabel
while smax - state.s >= update_s
    
    % format network inputs
    [map, Winv, state] = interp_map(state, signal, Winv);
    [psds, state, signal] = updatePSD(state, signal);
    images = cat(4, map, -map, map(:, end:-1:1, :, :), -map(:, end:-1:1, :, :));
    psds = repmat(psds, [1 1 1 4]);
    input = {
        'in_image', single(images), ...
        'in_psdmed', single(psds)
    };

    % run with mex-files
    state.net.eval(input);
    
    % extract result
    labels = squeeze(state.net.getVar(state.net.getOutputs()).value)';
    labels = reshape(mean(reshape(labels', [], 4), 2), 7, [])';
    
    % TODO: figure out how to save similarity to a buffer
    % ideas: save to buffer (bad because buffer is not BCILAB compliant...)
    %        save to state (possible, but need to predetermine buffer size)
    
    % update state
    [~, state.cls(state.next_ic)] = max(labels);
    state.s = state.s + update_s;
    state.next_ic = mod(state.next_ic, signal.nbchan) + 1;
    
end

% update signal
if ~isempty(Winv)
    signal.icawinv = Winv;
    signal.reject = ismember(state.cls, state.reject);
end

exp_endfun;

end


% linearly interpolate scalp topographies
function [map, Winv, state] = interp_map(state, signal, Winv)

% generate scalp map interpolation matrix (jerry rigged)
if ~isfield(state, 'interp')
    nChan = length(signal.chanlocs);
    in = eye(nChan);
    out = zeros(32^2,nChan);
    [~,Zi,~,Xi,Yi,intx,inty] = topoplotFast_LRBF(zeros(size(signal.chanlocs)), ...
        signal.chanlocs, 'noplot', 'on');
    for it = 1:nChan
        op = rbfcreate(double([inty;intx]),in(:,it)','RBFFunction', 'linear');
        out(:,it) = rbfinterp(double([Xi(:),Yi(:)]'), op);
    end
    state.interp.topoMat = out/in;
    state.interp.topoNaNMask = isnan(Zi);
    state.interp.topoNPixel = size(out,1);
    state.interp.topoMat(state.interp.topoNaNMask,:) = [];
end

% compute mixing matrix
if isempty(Winv)
    Winv = pinv(signal.icaweights * signal.icasphere);
end

% interpolate
map = zeros(state.interp.topoNPixel, 1);
map(~state.interp.topoNaNMask) = state.interp.topoMat * Winv(:, state.next_ic);
map = reshape(map, sqrt(state.interp.topoNPixel),[]);
map = single(0.99 * map / max(abs(vec(map))));

end


function [psd, state, signal] = updatePSD(state, signal)

% buffer data
% calc icaact
if ~isempty(signal.icaact)
    signal.icaact = (signal.icaweights * signal.icasphere) * signal.data;
end

% copy icaact to buffer
state.psd.buffer(:, mod1(state.psd.bufferind + (0:signal.pnts - 1), state.psd.bufferlen)) ...
    = signal.icaact;

% move buffer ind forward
state.psd.bufferind = state.psd.bufferind + signal.pnts;

% determine data available
if state.psd.bufferind == signal.pnts
    ind_newest = nanmax(state.psd.winds);
else
    ind_newest = 0;
end
pnts = state.psd.bufferind - ind_newest;
n_windows = min(ceil((pnts - state.psd.winlen + 1) / state.psd.unoverlapped), ...
                state.psd.nwindows);
ind_new = ind_newest + ((1:n_windows) * state.psd.unoverlapped);

% remove old window if possible
cutoff = state.psd.bufferind - (state.psd.sec2samp * signal.srate);
ind_old = state.psd.winds < cutoff;
state.psd.ffts(:, ind_old, state.next_ic) = nan;
state.psd.winds(ind_old) = nan;

% calc new ffts if enough data available
if length(ind_new)
    % extract new window if possible
    ind_mat = bsxfun(@plus, (1:state.psd.winlen)', ind_new - 1);
    ind_mat = mod1(ind_mat, state.psd.bufferlen);
    windows = reshape(state.psd.buffer(state.next_ic, ind_mat), size(ind_mat));

    % apply fft
    fftest = fft(windows, signal.srate);
    bad = find(isnan(state.psd.winds));
    state.psd.ffts(:, bad(1:length(ind_new)), state.next_ic) = ...
        abs(fftest(2:signal.srate / 2 + 1, :));
end

% compute estimate and set to 1--100Hz
psd = nanmean(state.psd.ffts(:, :, state.next_ic), 2);
if size(psd, 1) < 100
    psd = [psd; repmat(psd(end), 100 - size(psd, 1), 1)];
else
    psd(101:end) = [];
end
psd = psd / max(abs(psd(:)));

end


function out = mod1(x, y)
out = mod(x - 1, y) + 1;
end