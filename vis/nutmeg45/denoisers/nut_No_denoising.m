function nut_No_denoising(handles,signalspace,timept1,timept2,ptimept1,ptimept2)
% update bolts.meg, compute Rzz1 (InvRzz1 computed in nut_cov_eigen), InvES1

global bolts ndefaults

% No need to update bolts.meg for "No denoising"

if(1)  %changed to if(1) on 21/7/08; should be identical to if(0) anyway
%     if ndefaults.rzzcov
%         Rzz1 = cov(double(bolts.meg(timept1:timept2,:)));
%     else
%         Rzz1 = double(bolts.meg(timept1:timept2,:))'*double(bolts.meg(timept1:timept2,:));
%     end
Rzz1=nut_cov(bolts.meg(timept1:timept2,:,:),bolts.flags.avecov);
else
    Rzz1 = bolts.params.Rzz1;
end
condRzz = cond(Rzz1);
disp(['I just dropped in to see what condition my condition was in: ' num2str((condRzz),'%0.5g')]);

if ndefaults.bf.signalspace
   if(length(signalspace)==1)
      signalspace = 1:signalspace;
   end
end

[InvES1,q,u,v]=nut_eiginv(Rzz1,signalspace);

% update everything EXCEPT InvRzz1 (which is computed in nut_cov_eigen.m)
bolts.params.Rzz1 = Rzz1;
bolts.params.InvES1 = InvES1;
bolts.params.cond = condRzz;

bolts.params.signalspace = signalspace;
bolts.params.eigs = q;   % keep eigenvalues around for future reference
bolts.params.u = u;
bolts.params.v = v;

return
