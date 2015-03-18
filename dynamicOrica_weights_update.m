alpha = .05;
beta = 1;
delta = .5;
omdelta = 1-delta;

% tbd
A = pinv(weights);
lambda_old = 0.95;
r_old = rand(nChs)/nChs;
mult = zeros([nPts nChs nChs]);
for it = 1:nPts
    mult(it,:,:) = y(:,it)*f(:,it)';
end

r_t1 = bsxfun(@times,omdelta.^dataRange',repmat(permute(r_old,[3 1 2]),[nPts 1 1]));
r_t2 = cumsum(bsxfun(@times,delta*omdelta.^(dataRange'-1),mult),1);
r_t3 = cumsum(bsxfun(@times,delta*omdelta.^(dataRange'-1),repmat(permute(A,[3 1 2]),[nPts 1 1])),1);
rn_k = zeros(1,nPts);
lamda_k = zeros(1,nPts+1);
lambda_k(1) = lambda_old;
for it = 1:nPts
    rn_k(it) =  norm(squeeze(r_t1(it,:,:) + r_t2(it,:,:) - r_t3(it,:,:)),'fro')/nPts;
    lambda_k(it+1) = lambda_k(it)*(1 + alpha*(beta*rn_k(it)-lambda_k(it)));
end

%%

alpha = .05;
beta = 1;
delta = .5;
omdelta = 1-delta;

% tbd
A = pinv(weights);
lambda_old = 0.95;
r_old = rand(nChs)/nChs;
mult = zeros([nChs nChs nPts]);
for it = 1:nPts
    mult(:,:,it) = y(:,it)*f(:,it)';
end
mult = mean(bsxfun(@minus,mult,A),3);
rn = norm(omdelta*r_old + delta*mult,'fro')/nChs;
lambda = lambda_old*(1 + alpha*(beta*rn - lambda_old));
lambda = min(max(0,lambda),1);


r_t1 = bsxfun(@times,omdelta.^dataRange',repmat(permute(r_old,[3 1 2]),[nPts 1 1]));
r_t2 = cumsum(bsxfun(@times,delta*omdelta.^(dataRange'-1),mult),1);
r_t3 = cumsum(bsxfun(@times,delta*omdelta.^(dataRange'-1),repmat(permute(A,[3 1 2]),[nPts 1 1])),1);
rn_k = zeros(1,nPts);
lamda_k = zeros(1,nPts+1);
lambda_k(1) = lambda_old;
for it = 1:nPts
    rn_k(it) =  norm(squeeze(r_t1(it,:,:) + r_t2(it,:,:) - r_t3(it,:,:)),'fro')/nPts;
    lambda_k(it+1) = lambda_k(it)*(1 + alpha*(beta*rn_k(it)-lambda_k(it)));
end

