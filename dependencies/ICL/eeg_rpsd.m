function [psdmed, psdvar, psdkurt] = eeg_rpsd(EEG, nfreqs)

% clean input cutoff freq
nyquist = floor(EEG.srate / 2);
if ~exist('nfreqs', 'var') || isempty(nfreqs)
    nfreqs = nyquist;
elseif nfreqs > nyquist
    nfreqs = nyquist;
end

% setup constants
ncomp = size(EEG.icaweights, 1);
n_points = min(EEG.pnts, EEG.srate);
window = hamming(n_points)';
cutoff = floor(EEG.pnts / n_points) * n_points;
index = bsxfun(@plus, ceil(0:n_points / 2:cutoff - n_points), (1:n_points)');


psdmed = zeros(ncomp, nfreqs);
if nargout >= 2
    psdvar = zeros(ncomp, nfreqs);
    if nargout >= 3
        psdkurt = zeros(ncomp, nfreqs);
    end
end

% lower memory but slightly slower
for it = 1:ncomp
    % calculate windowed spectrums
    temp = reshape(EEG.icaact(it, index, :), [1 size(index) .* [1 EEG.trials]]);
    temp = bsxfun(@times, temp, window);
    temp = fft(temp, n_points, 2);
    temp = temp .* conj(temp);
    temp = temp(:, 2:nfreqs + 1, :) * 2 / (EEG.srate*sum(window.^2));
    if nfreqs == nyquist
        temp(:, end, :) = temp(:, end, :) / 2; end

    % calculate outputs
    psdmed(it, :) = db(median(temp, 3));
    if nargout >= 2
        psdvar(it, :) = db(var(temp, [], 3), 'power');
        if nargout >= 3
            psdkurt(it, :) = db(kurtosis(temp, [], 3)/4);
        end
    end
end

% % slightly faster but high memory use
% % calculate windowed spectrums
% temp = reshape(EEG.icaact(:, index, :), [ncomp size(index) .* [1 EEG.trials]]);
% temp = bsxfun(@times, temp, window);
% temp = fft(temp, n_points, 2);
% temp = temp .* conj(temp);
% temp = temp(:, 2:nfreqs + 1, :) * 2 / (EEG.srate*sum(window.^2));
% if nfreqs == nyquist
%     temp(:, end, :) = temp(:, end, :) / 2; end
% 
% % calculate outputs
% psdmed = db(median(temp, 3));
% if nargout >= 2
%     psdvar = db(var(temp, [], 3), 'power');
%     if nargout >= 3
%         psdkurt = db(kurtosis(temp, [], 3)/4);
%     end
% end

end