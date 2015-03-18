function nuts = nut_eeg_deref(nuts)
% NUT_EEG_DEREF  makes EEG data and lead potential reference free.
%   nuts = nut_eeg_deref(nuts)
% Reference: Van Veen et al. IEEE Trans Biomed Eng 1997; 44: 867

[U,dum] = svd(nuts.meg.data(:,:,1)');   % We only use first trial for SVD, since this does not seem to matter:
                                        % Using all trials or the covariance matrix produces exactly the same results
                                        % We could also use qr instead of svd, again does not seem to matter.
U(:,end)=[];
UT = U';

numnewchan=size(U,2);

for k=1:size(nuts.meg.data,3)
    nuts.meg.data(:,1:numnewchan,k) = (UT * nuts.meg.data(:,:,k)')';
end
nuts.meg.data(:,end,:)=[];

[nc,no,nv]=size(nuts.Lp);
if nc>size(nuts.meg.data,2)
    if nc==max(nuts.meg.goodchannels)
        nuts.Lp = nuts.Lp(nuts.meg.goodchannels,:,:);
    else
        error('Data does not correspond to leadfield.')
    end
end

for k=1:no
    Lo = permute(nuts.Lp(:,k,:),[1 3 2]);
    nuts.Lp(1:numnewchan,k,:) = reshape(UT * Lo, [numnewchan 1 nv]);
end
nuts.Lp(end,:,:)=[];

nuts.meg.reference='deref';
nuts.meg.referenceidx=[];

    
    