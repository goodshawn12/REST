function [output]=nut_Correlate_Columns(Lp,data, flags) %---------------------------------------------------------
% Lp : lead field ( channel X 3 )
% InvRzz : Inverse data covariance matrix
% InvES : Inverse signal-component-data covariance matrix

% global bolts

if size(Lp,2) ~= 2
   disp(' ');
   error('   This code is (currently) designed only for homogeneous, spherical conductor head models')
end

if 1
   disp('Ignores selected interval and uses all time points instead.')
   t = 1:length(data.latency); % indices
else
   disp('Ignores selected interval and uses post stimulus period instead.')
   t = find(data.latency >= 0);
end

[N,NS] = size(data.y(t,:)); % number of time points, number of sensors

% remove mean
for i = 1:NS
   data.y(:,i) = data.y(:,i) - mean(data.y(t,i));
end

% dimension reduction
% M1_ndx = bolts.params.signalspace; % indices of signal space
M1_ndx = flags.tb.nl; % indices of signal space
M1 = length(M1_ndx); % number of source estimates to calculate
[U,S,V] = svd(data.y(t,:)');
[svalues,pos] = sort(diag(S),'descend');
U = U(:,pos); S = diag(svalues); V = V(:,pos);
w1 = diag(1./svalues(M1_ndx))*U(:,pos(M1_ndx))'*sqrt(N);
b = V(:,pos(M1_ndx))'*sqrt(N); % b = w1*data.y(t,:)'; b is sphered

% blind source separation (removes cross correlation at two lags)
tau = 10; % the second lag (in addition to lag 0)
Rx2 = b(:,1:N-tau)*b(:,tau+1:N)'/(N-tau); Rx2 = 0.5*(Rx2+Rx2');
[w2inv,jnk] = svd(Rx2);
w2 = inv(w2inv);
y = w2*b; % source estimates

w = w2*w1; % estmimated demixing matrix (including dimension reduction)
winv = w1'*w2inv; % estimated mixing matrix

if 1
   disp(' ');
   disp(['   Hit the space bar to view each source estimate.'])
   disp(['   Determine which rows of ''y'' correspond to source (estimates) of interest.'])
   disp(' ');

   figure
   for i = 1:size(y,1)
      plot(data.latency(t),y(i,:));
      disp(['   ndx = ' int2str(i)]);
      pause
   end

   disp(' ');
   disp(['   Set ndx to be a row vector containing the indices of the interesting rows of ''y''. (e.g., ndx = [1 5 7])'])
   disp(['   Then type ''return'' to continue ...'])
   disp(' ');

   keyboard % set ndx, e.g., ndx = [1 5 7];
else
   ndx = 2
end

clear b

M = length(ndx);
winv = winv(:,ndx); % redefines winv

% remove mean
for i = 1:M
   winv(:,i) = winv(:,i) - mean(winv(:,i));
end
winv_std = sqrt(sum(winv.^2,1));

L = zeros(NS,size(Lp,3));

% code comes from localize.m
for i = 1:size(Lp,3)
   Ltmp = Lp(:,:,i);
   for j = 1:size(Ltmp,2)
      Ltmp(:,j) = Ltmp(:,j) - mean(Ltmp(:,j));
   end

   psi_tmp = zeros(2,M);
   corr_tmp = zeros(1,M);
   for j = 1:M
      [Q,V] = eig((Ltmp'*winv(:,j))*(winv(:,j)'*Ltmp),Ltmp'*Ltmp); % maximize normalized correlation
      [val,pos] = max(diag(V));
      if Q(:,pos)'*Ltmp'*winv(:,j) < 0, Q(:,pos) = -Q(:,pos); end % correct the sign
      psi_tmp(:,j) = Q(:,pos);

      L_std = sqrt(sum((Ltmp*Q(:,pos)).^2,1));
      if L_std < eps, L_std = eps; end
      corr_tmp(j) = (Q(:,pos)'*Ltmp'*winv(:,j))/(L_std*winv_std(j));
   end
   [val,pos] = max(corr_tmp);
   L(:,i) = Lp(:,:,i)*psi_tmp(:,pos);
end

% remove mean
for i=1:size(L,2)
   L(:,i) = L(:,i) - mean(L(:,i));
end

% compute standard deviation of each column
L_std = sqrt(sum(L.^2));
fL = find(L_std < eps); L_std(fL) = eps;

% compute absolute value of normalized corrrelation coefficient between columns
out = zeros(M,size(Lp,3));
for i = 1:size(Lp,3)
   for j = 1:M
      out(j,i) = (abs(L(:,i)'*winv(:,j))/(L_std(i)*winv_std(j)));
   end
end

if M > 1
   out = max(out);
end

% exponentiate
area_ratio = 0.1; % 0 <= area_ratio <= 1 (0.02)
thresh = 0.25; % 0 <= thresh <= 1

out = out - min(out);
out = out/max(out);
exponent = 0;
f = isnan(out); f = find(f == 0);
done = 0;
while ~done
   exponent = exponent + 1;
   f2 = find(out(f).^exponent > thresh);
   if length(f2)/length(f) < area_ratio
      done = 1;
   end
end
out = out.^exponent;
disp(' '); disp(['Exponentiating results. Exponent = ' int2str(exponent) '.'])

if max(max(out)) < 10^(-30)
   out = out*10^(-18)/max(max(out));
   disp('Note - output power has been increased for display purposes'), disp(' ')
end

weight = zeros(size(Lp));
weight(1,1,:) = sqrt(out); % this variable does not represent the weights of a beamformer; this was done to display the desired results using the existing beamformer code

output.W1=weight;
output.meg=[ones(N,1) zeros(N,NS-1)]; % changes data in order to display the desired results

return
