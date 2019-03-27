function [signal, state] = flt_eyecatch(varargin)


if ~exp_beginfun('filter') return; end

% requires epoched data, works best on spatially filtered data
declare_properties('name','EyeCatch', ...
    'depends', 'set_makepos', ...
    'follows', 'flt_orica', ...
    'precedes', 'flt_ica_reproject', ...
    'independent_channels', true, ...
    'independent_trials', false);

eyecatch_path = fullfile(fileparts(fileparts(which('bcilab'))), 'eyeCatch', 'libEyeCatch.mat');

arg_define(varargin,...
    arg_norep({'signal','Signal'}), ...
    arg({'cutoff','Cutoff'}, 0.89, [0 1], 'Correlation cutoff for eye-component detection.'), ...
    arg({'freq', 'updatefreq','UpdateFreq'}, 6, [], 'How many ICs to classify per second.'), ...
    arg({'libpath', 'eyecatchlib', 'eyeCatchLib'}, eyecatch_path, [], 'Path to eyeCatch library mat-file.'), ...
    arg_nogui({'state','State'}));

% state.history: past classifications?
if isempty(state)
    state.lib = load(libpath);
    state.cls = false(signal.nbchan, 1);
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
    [map, Winv, state] = interp_map(state, signal, Winv);
    [isEyeIC, similarity] = runEyeCatch(state.lib, map, cutoff);
    % TODO: figure out how to save similarity to a buffer
    % ideas: save to buffer (bad because buffer is not BCILAB compliant...)
    %        save to state (possible, but need to predetermine buffer size)
    
    % update state
    state.cls(state.next_ic) = isEyeIC;
    state.s = state.s + update_s;
    state.next_ic = mod(state.next_ic, signal.nbchan) + 1;
    
end

% update signal
if ~isempty(Winv)
    signal.icawinv = Winv;
end
signal.reject = state.cls;

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

end


% run eyeCatch classifier
function [isEyeIC, similarity] = runEyeCatch(libEyeCatch, map, threshold)

    % normalize map
    normMap = bsxfun(@minus, map,  mean(map));
    normMap = bsxfun(@rdivide, normMap,  std(normMap));
    
    % import library
    lib = libEyeCatch.new_map;   
    
    similarity  = max(abs(lib * normMap)) / length(normMap);
    isEyeIC = similarity > threshold;
end

