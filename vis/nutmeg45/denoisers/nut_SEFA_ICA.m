function nut_SEFA_ICA(handles,signalspace,timept1,timept2)
% update bolts.meg, compute Rzz1, InvES1 (InvRzz1 computed in nut_cov_eigen)

global bolts nuts

num_factors = 2;
distribution = 'peaky'; % use 'peaky' or 'sin'
ndxt = 1:num_factors; % use 1:num_factors to use all factors, else specify a single factor

pre = find(nuts.meg.latency < 0);
post = find(nuts.meg.latency >= 0);

if length(pre) < 10 | length(post) < 10
   error('Not enough data in either the pre-stimulus or post-stimulus period for this method')
end

[yproj,xbar,wifa,cy,cysep] = nut_sefaica(bolts.meg(pre,:)',bolts.meg(post,:)',[],0,num_factors,distribution);

if length(ndxt) == 1, Rzz1 = cysep(:,:,ndxt); else, Rzz1 = cy; end
condRzz = cond(Rzz1);
disp(['I just dropped in to see what condition my condition was in: ' num2str((condRzz),'%0.5g')]);
InvES1 = inv(Rzz1);

% update everything EXCEPT InvRzz1 (which is computed in nut_cov_eigen.m)
bolts.params.Rzz1 = Rzz1;
bolts.params.InvES1 = InvES1;
bolts.params.cond = condRzz;

bolts.params.yc_full = zeros(size(bolts.meg'));
for i = ndxt
   bolts.params.yc_full(:,length(pre)+1:end) = bolts.params.yc_full(:,length(pre)+1:end) + yproj(:,:,ndxt(i));
end

return
