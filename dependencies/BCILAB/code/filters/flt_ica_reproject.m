function signal = flt_ica_reproject(varargin)
% Spatially project the given icaactivations of an IC decomposition to
% channel data
% Signal = flt_project(Signal, ProjectionMatrix, Componentsubset)
%
% In:
%   Signal           : epoched or continuous EEGLAB data set
%
%   ProjectionMatrix : projection matrix to apply; can be one of the following:
%                      * '.icawinv' : apply ICA back-projection
%
%   ComponentSubset  : subset of components to reproject;
%                      * '.reject' : use selected rejection
%
% Out:
%   Signal : projected EEGLAB data set
% 
% Examples:
%   % project onto a 10-dimensional random subspace, assuming that there are 32 channels in the data
%   eeg = flt_project(eeg,'.icawinv','.reject')

%                                Luca
%                                2018-01-08


if ~exp_beginfun('filter') return; end;

% would be reverted by an ICA
declare_properties('name','Reprojection', 'follows',{'flt_ica', 'flt_orica', 'flt_eyecatch'}, 'independent_channels',false, 'independent_trials',true);

arg_define(varargin,...
    arg_norep({'signal','Signal'}), ...
    arg({'projmat','ProjectionMatrix'}, '.icawinv', {'.icawinv'}, 'Projection matrix. The ICA components are multiplied by this matrix.'),...
    arg({'subcomps','ComponentSubset'}, '.reject', [], 'Component subset. List of component indices which will be used in the reprojection.','shape','row'));

% determine projection matrices
if strcmp(projmat,'.icawinv')
    if isfield(signal,'icawinv') && ~isempty(signal.icawinv)
        subchans = signal.icachansind;
        projmat = signal.icawinv;
    elseif isfield(signal,'icaweights') && ~isempty(signal.icaweights)
        signal.icawinv = pinv(signal.icaweights * signal.icasphere);
        subchans = signal.icachansind;
        projmat = signal.icawinv;
    else
        error('No icaweights or icawinv fields are present.')
    end
end

% ammend channels if subset was used for ica
if ~isempty(subchans)
    subset = set_chanid(signal,subchans);
    if ~isequal(subset,1:signal.nbchan)        
        signal.data = signal.data(subset,:,:,:,:,:,:,:);
        signal.chanlocs = signal.chanlocs(subset);
        signal.nbchan = size(signal.data,1);
    end
end

% determine component rejections
if strcmp(subcomps, '.reject')
    subcomps = ~signal.reject;
elseif isempty(subcomps)
    subcomps = true(1, size(projmat, 1));
end

% project data
[C,S,T] = size(signal.data); %#ok<*NODEF>
if size(projmat,2) ~= C
    error('The given projection matrix needs to have the same number of rows as the data has channels (%i), but had: %i',C,size(projmat,2)); end
if isempty(signal.icaact)
    temp = reshape(signal.data,C,[]);
    signal.icaact = reshape(signal.icaweights * temp,[],S,T);
end
temp = reshape(signal.icaact,C,[]);
signal.data = reshape(projmat(:, subcomps) * temp(subcomps, :),[],S,T);

% rewrite chanlocs
if isempty(subchans) || isequal(subchans, 1:signal.nbchan)
    signal.nbchan = size(signal.data,1);
    if isempty(signal.urchanlocs)
        signal.urchanlocs = signal.chanlocs; end
    signal.chanlocs = signal.chanlocs(subchans);
end

exp_endfun;
