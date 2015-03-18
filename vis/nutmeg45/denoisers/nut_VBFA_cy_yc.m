function nut_VBFA_cy_yc(handles,signalspace,timept1,timept2,extra,crap)
% update bolts.meg, compute Rzz1 (InvRzz1 computed in nut_cov_eigen), InvES1

global bolts

num_factors = max(signalspace);
num_iterations = 50;

data=mean(bolts.meg(timept1:timept2,:,:),3)';

[b,lam,sig,yc,cy,bet,weight] = nut_reg_vbfa(data,num_factors,num_iterations,1);

Rzz1=nut_cov(bolts.meg(timept1:timept2,:,:),bolts.flags.avecov);
% Rzz1 = cy;

condRzz = cond(Rzz1);
disp(['I just dropped in to see what condition my condition was in: ' num2str((condRzz),'%0.5g')]);
InvES1 = inv(cy);

% update everything EXCEPT InvRzz1 (which is computed in nut_cov_eigen.m)
bolts.params.Rzz1 = Rzz1;
bolts.params.InvES1 = InvES1;
bolts.params.cond = condRzz;

bolts.params.yc = yc;
bolts.params.yc_full = weight*mean(bolts.meg,3)';
bolts.params.bet = bet;

return
