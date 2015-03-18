function nut_SEFA_cy_yc(handles,signalspace,timept1,timept2,ptimept1,ptimept2)
% update bolts.meg, compute Rzz1 (InvRzz1 computed in nut_cov_eigen), InvES1

global bolts

nl=max(signalspace);
nm=round(size(bolts.meg,2)/3);

[y_clean,cleaning,www,cy,sig,cyall] = nut_sefa_wrap(bolts.meg,nl,nm,[timept1 timept2],[ptimept1 ptimept2],bolts.flags.avecov);

if bolts.flags.avecov
    Rzz1 = cy;
else
    Rzz1 = cyall;
end
condRzz = cond(Rzz1);
disp(['I just dropped in to see what condition my condition was in: ' num2str((condRzz),'%0.5g')]);
InvES1 = inv(cy);

% update everything EXCEPT InvRzz1 (which is computed in nut_cov_eigen.m)
bolts.params.Rzz1 = Rzz1;
bolts.params.InvES1 = InvES1;
bolts.params.cond = condRzz;

bolts.params.yc = y_clean;
bolts.params.yc_full = cleaning*mean(bolts.meg,3)';
% bolts.params.bet = bet;

return
