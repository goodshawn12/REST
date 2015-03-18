function nutsnew = nut_megrealign(nutsold,nutsnew);
% nuts = nut_megrealign(nutsold,nutsnew)
% uses new lead field from nutsnew to realign data from nutsold.
% requires:
% nutsold.meg.data
% nutsold.Lp
% nutsnew.Lp



lfold = [squeeze(nutsold.Lp(:,1,:)) squeeze(nutsold.Lp(:,2,:)) squeeze(nutsold.Lp(:,3,:))];
lfnew = [squeeze(nutsnew.Lp(:,1,:)) squeeze(nutsnew.Lp(:,2,:)) squeeze(nutsnew.Lp(:,3,:))];

tmp = pinv(lfold);

fprintf('computing interpolation matrix #1\n');
realign = lfnew * tmp;
if 1
  fprintf('computing interpolation matrix #2\n');
  noalign = lfold * tmp;
  fprintf('computing interpolation matrix #3\n');
  bkalign = lfold * pinv(lfnew) * realign;
end

% interpolate the data towards the template gradiometers
for i=1:size(nutsold.meg.data,3)
  fprintf('realigning trial %d\n', i);
  nutsnew.meg.data(:,:,i) = nutsold.meg.data(:,:,i) * realign';
  % also compute the residual variance when interpolating
  if(1)
      rvrealign = rv(nutsold.meg.data(1:100000,:,i), nutsnew.meg.data(1:100000,:,i));
      fprintf('original -> template             RV %.2f %%\n', 100 * mean(rvrealign));
  else
      rvrealign = rv(nutsold.meg.data(:,:,i), nutsnew.meg.data(:,:,i));
      fprintf('original -> template             RV %.2f %%\n', 100 * mean(rvrealign));
      if 1
          datnoalign = nutsold.meg.data(:,:,i) * noalign';
          datbkalign = nutsold.meg.data(:,:,i) * bkalign';
          rvnoalign = rv(nutsold.meg.data(:,:,i), datnoalign);
          rvbkalign = rv(nutsold.meg.data(:,:,i), datbkalign);
          fprintf('original             -> original RV %.2f %%\n', 100 * mean(rvnoalign));
          fprintf('original -> template -> original RV %.2f %%\n', 100 * mean(rvbkalign));
      end
  end
end


function invx = prunedinv(x,r)
[u,q,v]=svd(x,'econ');
q=diag(q);
signalspace = find(q/sum(q) > r);

% compute qinv by truncating and taking reciprocal of q;
% we want to keep original q for future reference and troubleshooting
qinv=zeros(size(q));
qinv(signalspace)=1./q(signalspace);

invx=v*diag(qinv)*u';
