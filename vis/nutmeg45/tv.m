function tv(varargin)
% tv [s_beamtf_results.mat]
% tv is a shortcut to open nut_timef_viewer and nutmeg

global nuts

% disp('calling nut_results_viewer, not nut_timef_viewer!')
if ~isfield(nuts,'fig')
    nutmeg(nuts)
end
% nut_timef_viewer(varargin{:});
nut_results_viewer(varargin{:});
