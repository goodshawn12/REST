function [weight,eta]=nut_Eigenspace_Scalar_Beamformer(Lp,data, flags) %---------------------------------------------------------
% Lp : lead field ( channel X 3 )
% InvRyy : Inverse data covariance matrix
% InvES : Inverse signal-component-data covariance matrix
% data.eta : (optional) optimum orientation for given voxel
% note: this method works for 3-component lead fields, but not tested for 2 component

dualstate = false;

% vectorized section -- we can compute Lp'*InvRyy and InvES*Lp for
% each voxel in one big matrix -----------
Lp2 = reshape(Lp,size(Lp,1),size(Lp,2)*size(Lp,3));

InvRyyLp = reshape(data.InvRyy*Lp2,size(Lp,1),size(Lp,2),size(Lp,3));
InvESLp = reshape(data.InvES*Lp2,size(Lp,1),size(Lp,2),size(Lp,3));

clear Lp2;
% end vectorized section ------------

weight=zeros(size(Lp,1),size(Lp,3));
if isfield(data,'eta')
    eta=data.eta;
else
    eta=zeros(size(Lp,3),size(Lp,2));
end

rnkLp = nut_rank(Lp(:,:,1));
switch(rnkLp)
    case 1
        disp('Lead field of rank 1 detected: scalar lead field');
    case 2
        disp('Lead field of rank 2 detected: single sphere model assumed');
    case 3
        disp('Lead field of rank 3 detected: fancy head model assumed');
    otherwise
        error(['Doooood, your lead field is cracked out! It''s rank is ' num2str(rnkLp) ' but needs to be either 2 or 3.']);
end

for i=1:size(Lp,3)
    % equivalent to: inv(Lp' * InvRyy* Lp)

    invJ = Lp(:,:,i)'*InvRyyLp(:,:,i);

    % We use pinv here to accomodate single-sphere head models (which
    % will be of rank 2 and fail with a true inverse).
    % For more sophisticated models, the rank will generally be 3, and
    % then pinv yields the same results as inv.
    % see nut_pinv for modifications to pinv.
    J=double(nut_pinv(invJ));

    % note InvRyy1 = InvRyy'
    if(0 & dualstate) % unnnnntested!!!
        % good for dual state? compute using Rall
        [v,d]=svd(invJall);
        %             [v,d]=svd(Jall*InvRallLp(:,:,i)'*InvRallLp(:,:,i));
    else  % original: compute using Ryy
        [v,d]=svd(invJ);
        %        [v,d]=svd(J*InvRyyLp(:,:,i)'*InvRyyLp(:,:,i));
    end
    eta(i,:) = v(:,rnkLp);
    L = Lp(:,:,i)*v(:,rnkLp);

    InvRyyL = data.InvRyy*L;
    InvESL = data.InvES*L;

    
    % note InvRyy1 = InvRyy' and InvES = InvES';
    weight(:,i)=InvESL*inv(L'*InvRyyL);
%     weight(:,:,i)=InvESLp(:,:,i)*J; % Nsensors*3 matrix

    if(flags.wn)
        % equivalent to: inv(Lp' * InvRyy* Lp) * (Lp' * InvRyy^2 * Lp) * inv(Lp' * InvRyy * Lp)
        % Lp' * InvRyy^2 * Lp = (Lp' * InvRyy) * (InvRyy*Lp) = (InvRyy*Lp)' * (InvRyy*Lp)
        omega=J*(InvRyyLp(:,:,i)' * InvRyyLp(:,:,i))*J;

        % real is needed for MATLAB acceleration --
        % potential complex output of sqrt operation disables acceleration!
        weight(:,i)=weight(:,i)/real(sqrt(omega(1,1)));
    end
    
    if(flags.progressbar)
        if(~mod(i,500))  % update progress every 500th voxel
            %         if(ndefaults.bf.viewwaitbar)
            %             waitbar(iv/nv,h);
            %         end
            disp(['Please wait. Eigenspace Scalar Beamformer has finished ' num2str(i) ' out of ' num2str(size(Lp,3)) ' voxels.']);
        end
    end
end
