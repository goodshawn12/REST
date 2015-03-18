function nut_cov_eigen(invtype)
% Computes bolts.params.InvRzz only (Rzz1, InvES, et al. computed by denoising method)

global nuts bolts

condRzz = bolts.params.cond;
Rzz1 = double(bolts.params.Rzz1);

condthresh = 1e12; %valid to have thresh same for all types of inversion?

switch(invtype)
   case 'tikhonov'
      if(condRzz > condthresh)    %%%% NOTE: look @ optimum value here
         fprintf(['FYI: using ' invtype ' to estimate inverse of Rzz1; cond = %g\nMAKE SURE YOU KNOW WHAT THIS MEANS!!!\n'],condRzz);
         if isfield(bolts.params,'eigs')
             normRzz = bolts.params.eigs(1);
         else
             normRzz = norm(Rzz1);  % equivalent to first singular value of Rzz1
         end
         gamma = 1e-6*normRzz;
         InvRzz1 = inv(Rzz1+gamma*eye(size(Rzz1)));
      else
         InvRzz1=inv(Rzz1);
      end
   case 'pseudoinverse'
      if(condRzz > condthresh)    %%%% NOTE: look @ optimum value here
         fprintf(['FYI: using ' invtype ' to estimate inverse of Rzz1; cond = %g\nMAKE SURE YOU KNOW WHAT THIS MEANS!!!\n'],condRzz);
         InvRzz1=pinv(Rzz1);
      else
         InvRzz1=inv(Rzz1);
      end
   case 'bayesian'
      if(condRzz > condthresh)    %%%% NOTE: look @ optimum value here
         fprintf(['FYI: using ' invtype ' to estimate inverse of Rzz1; cond = %g\nMAKE SURE YOU KNOW WHAT THIS MEANS!!!\n'],condRzz);
         gamma = var(bolts.meg,1);
         InvRzz1 = inv(Rzz1+(diag(gamma)/(size(nuts.meg.data,1)+1)));
      else
         InvRzz1=inv(Rzz1);
      end
   case 'vanilla'
      if(condRzz > condthresh)    %%%% NOTE: look @ optimum value here
         fprintf(['FYI: using ' invtype ' to estimate inverse of Rzz1; cond = %g\nMAKE SURE YOU KNOW WHAT THIS MEANS!!!\n'],condRzz);
         InvRzz1 = inv(Rzz1);
      else
         InvRzz1=inv(Rzz1);
      end
end           

bolts.params.InvRzz1 = InvRzz1;

return
