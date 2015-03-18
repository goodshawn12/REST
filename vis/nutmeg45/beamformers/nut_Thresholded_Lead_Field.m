function [weight]=nut_Thresholded_Lead_Field(Lp,data, flags) %---------------------------------------------------------
% Lp : lead field ( channel X 3 )
% InvRzz : Inverse data covariance matrix
% InvES : Inverse signal-component-data covariance matrix
% Note that weight normalization is the only option here!

% global bolts nuts

% timept1 = dsearchn(nuts.meg.latency,bolts.timewindow(1));
% timept2 = dsearchn(nuts.meg.latency,bolts.timewindow(end));
% Rzz=bolts.meg(timept1:timept2,:)'*bolts.meg(timept1:timept2,:);

% vectorized section -- we can compute Lp'*InvRzz and InvES*Lp for
% each voxel in one big matrix -----------
Lp2 = reshape(Lp,size(Lp,1),size(Lp,2)*size(Lp,3));

InvRzzLp = reshape(data.InvRyy*Lp2,size(Lp,1),size(Lp,2),size(Lp,3));
% waitbar(0.25,bolts.barhandle)
% 
% if(bolts.stop)  % if cancel button on waitbar is pressed...
%     return
% end

% InvESLp = reshape(data.InvES*Lp2,size(Lp,1),size(Lp,2),size(Lp,3));
% waitbar(0.5,bolts.barhandle)

clear Lp2;
% end vectorized section ------------

weight=zeros(size(Lp,1),2,size(Lp,3));
for i=1:size(Lp,3)
%     if(bolts.stop)  % if cancel button on waitbar is pressed...
%         return
%     end

    [Ltemp1,Lindex1] = sort(abs(Lp(:,1,i)));
    [Ltemp2,Lindex2] = sort(abs(Lp(:,2,i)));
    %     Lthresh1 = Ltemp1(floor(.95*size(Lp,1)));
    %     Lthresh2 = Ltemp2(floor(.95*size(Lp,1)));

    Lsig = (floor(.8*size(Lp,1))):size(Lp,1);
    Lindex = union(Lindex1(Lsig),Lindex2(Lsig));
    L = Lp(Lindex,:,i);

    InvRyy = inv(data.Ryy(Lindex,Lindex));

    [u,q,v]=svd(data.Ryy(Lindex,Lindex));
    q=diag(q);
    qinv=zeros(size(q));
    qinv(flags.signalspace)=1./q(flags.signalspace);
    InvES=v*diag(qinv)*u';

    % equivalent to: inv(Lp' * InvRzz* Lp)
    %     J=inv(Lp(:,:,i)'*InvRzzLp(:,:,i));    % 2*2 matrix
    J = inv(L' * InvRyy * L);

    % note InvRzz1 = InvRzz' and InvES = InvES';
    weight(Lindex,:,i)= InvES' * L * J; % Nsensors*2 matrix

    if(flags.wn)
        % equivalent to: inv(Lp' * InvRzz* Lp) * (Lp' * InvRzz^2 * Lp) * inv(Lp' * InvRzz * Lp)
        % Lp' * InvRzz^2 * Lp = (Lp' * InvRzz) * (InvRzz*Lp) = (InvRzz*Lp)' * (InvRzz*Lp)
        %        omega=J*(InvRzzLp(:,:,i)' * InvRzzLp(:,:,i))*J;
        omega = J * (L' * InvRyy) * (InvRyy * L) * inv(L' * InvRyy *L);


        % real is needed for MATLAB acceleration --
        % potential complex output of sqrt operation disables acceleration!
        weight(Lindex,1,i)=weight(Lindex,1,i)/real(sqrt(omega(1,1)));
        weight(Lindex,2,i)=weight(Lindex,2,i)/real(sqrt(omega(2,2)));
    end

    if(~mod(i,500))  % update progress every 500th iteration
        disp(['Please wait. Thresholded Lead Field Beamformer has finished ' num2str(i) ' out of ' num2str(size(Lp,3)) ' voxels.']);
%         waitbar(0.5 + i/(2*size(Lp,3)),bolts.barhandle)
    end
end
% end
