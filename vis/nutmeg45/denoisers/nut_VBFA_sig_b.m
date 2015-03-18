function nut_VBFA_sig_b(handles,signalspace,timept1,timept2)
% update bolts.meg, compute Rzz1 (InvRzz1 computed in nut_cov_eigen), InvES1

global bolts

num_factors = 5;
num_iterations = 50;

[b,lam,sig,yc,cy,bet,weight] = nut_reg_vbfa(bolts.meg(timept1:timept2,:)',num_factors,num_iterations,1);

Rzz1 = sig;
condRzz = cond(Rzz1);
disp(['I just dropped in to see what condition my condition was in: ' num2str((condRzz),'%0.5g')]);
InvES1 = inv(sig);

% update everything EXCEPT InvRzz1 (which is computed in nut_cov_eigen.m)
bolts.params.Rzz1 = Rzz1;
bolts.params.InvES1 = InvES1;
bolts.params.cond = condRzz;

bolts.params.bet = bet;

return
