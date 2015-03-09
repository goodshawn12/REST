function [dipoles, L, m, w] = bfpf(data,K,vertices,nDipoles,nParticles,Q_location,Q_data,state)
% check inputs
nChan = size(data,1);
nTimes = size(data,2);
if ~any(size(vertices)==3)
    error('Vertices must be a 3 x nVert matrix')
elseif size(vertices,2)==3
    vertices = vertices';
end
nVert = size(vertices,2);
if ~exist('Q_location','var') || isempty(Q_location)
    Q_location = .5*eye(3); end
if ~exist('Q_data','var') || isempty(Q_data)
    Q_data = 1e-5*eye(nChan); end
rng(0);

% precalculate pseudo inverse of Lead Field Matrix transposed
Ktpi = pinv(K');
cQ_location = chol(Q_location);

% initialize weights
weights = ones(1,nParticles)/nParticles;
% initialize dipole locations
L = randi(nVert,1,nParticles);
% iterate
for it = 1:nTimes
    Lold = L;
    % generate random numbers
    w = cQ_location*randn(3,nParticles);
    % displace current locations
    vtemp = vertices(:,Lold) + w;
    % lock new location to existing vertex
    for itPart = 1:nParticles
        [~,L(itPart)] = min(row_pnorm(bsxfun(@minus,vertices,vtemp(:,itPart)),2,'column'));
    end
    % skipping Wo !!!!! this might mean to only take the dipole locations
    W = eye(nChan);
    % calculate moments vector
    m = Ktpi(:,tri_index(Lold,nVert))'*W*data(:,it);
    % update new weights
    x_hat = sum(reshape(bsxfun(@times,m',K(:,tri_index(Lold,nVert))),[nChan,nParticles,3]),3);
    p_ylx = mvnpdf(bsxfun(@minus,data(:,it),x_hat)',zeros(nChan,1)',Q_data);
    weights = weights.*mvnpdf(bsxfun(@minus,data(:,it),x_hat)',zeros(nChan,1)',Q_data)';
    % normalize weights
    weights = weights/sum(weights);
    % resample from L (with replacement) according to importance weights
    if it<nTimes
        L = datasample(L,nParticles,'weights',weights); end
end

dipoles = struct('location',[],'moment',[]);
[~,ind] = sort(weights,2,'ascend');
for it = 1:nDipoles
    dipoles(it).location = vertices(:,L(ind(it)));
    dipoles(it).moment = m(tri_index(ind(it),nParticles));
end


function index = tri_index(index,nVert)
if iscolumn(index)
    index = index'; end
index = [index,index+nVert,index+2*nVert];