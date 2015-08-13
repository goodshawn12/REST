% Minimum variance dipole fitting
%

function [dipoles, L, moments, rv, state] = dipole_fit(data,K,vertices,state)
nDipoles = 1;

% check inputs
nChan = size(data,1);
nTimes = size(data,2);
if ~any(size(vertices)==3)
    error('Vertices must be a 3 x nVert matrix')
elseif size(vertices,2)==3
    vertices = vertices';
end
nVert = size(vertices,2);
nK = size(K,2);
if exist('state','var') && ~isempty(state)
    flag_state = true;
else
    flag_state = false;
end
if nVert==length(K)
    flag_fixedDir = true;
else
    flag_fixedDir = false;
end
K = double(K); % add check function

% get state
if flag_state
    Kpi_individual = state.Kpi_individual;
else
    % precalculate pseudo inverse of Lead Field Matrix transposed
    Kpi_individual = zeros(size(K'));
    for it = 1:length(Kpi_individual)/3
        Kpi_individual(tri_index(it,nVert,flag_fixedDir),:) = pinv(K(:,tri_index(it,nVert,flag_fixedDir))); end
    state.Kpi_individual = Kpi_individual;
    state.nVert = length(vertices);
end

for it = 1:nTimes
    % initialize possible dipoles
    L = 1:nVert;
    moments = zeros(nK,1);
    
    % calculate moments vector
    for itM = 1:nVert
        moments(tri_index(itM,nVert,flag_fixedDir)) = Kpi_individual(tri_index(L(itM),nVert,flag_fixedDir),:)*data;
    end
    
    % calculate residual variances
    if flag_fixedDir
        y_hat = bsxfun(@times,moments',K(:,tri_index(L,nVert,flag_fixedDir)));
    else
        y_hat = sum(reshape(bsxfun(@times,moments',K(:,tri_index(L,nVert,flag_fixedDir))),[nChan,nVert,3]),3);
    end
    rv = var(bsxfun(@minus,data,y_hat))/var(data);
end

% find lowest residual variance dipole
dipoles = struct('location',[],'moment',[]);
[~,ind] = sort(rv);
for it = 1:nDipoles % add option for multiple dipoles in the future
    dipoles(it).location = vertices(:,L(ind(it)));
    dipoles(it).L = tri_index(ind(it),nVert,flag_fixedDir);
    dipoles(it).moment = moments(dipoles(it).L);
end
rv = rv(ind(1));


function index = tri_index(index,nVert,flag_fixedDir)
if flag_fixedDir
    return
else
    if iscolumn(index)
        index = index'; end
    index = [index,index+nVert,index+2*nVert];
end