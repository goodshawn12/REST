function nut_SEFA_sig_b(handles,signalspace,timept1,timept2)
% update bolts.meg, compute Rzz1 (InvRzz1 computed in nut_cov_eigen), InvES1

global bolts

nl=max(signalspace);
nm=round(size(bolts.meg,2)/3);

[y_clean,cleaning,www,cy,sig] = nut_sefa_wrap(nl,nm);

Rzz1 = sig;
condRzz = cond(Rzz1);
disp(['I just dropped in to see what condition my condition was in: ' num2str((condRzz),'%0.5g')]);
InvES1 = inv(sig);

% update everything EXCEPT InvRzz1 (which is computed in nut_cov_eigen.m)
bolts.params.Rzz1 = Rzz1;
bolts.params.InvES1 = InvES1;
bolts.params.cond = condRzz;

% bolts.params.bet = bet;

return
