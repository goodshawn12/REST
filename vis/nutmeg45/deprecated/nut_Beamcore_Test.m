function [weight]=beamcore_test(Lp,InvRzz, InvES) %---------------------------------------------------------
% Lp : lead field ( channel X 3 )
% InvRzz : Inverse data covariance matrix
% InvES : Inverse signal-component-data covariance matrix
% Note that weight normalization is the only option here!

% unvectorized (but matlab-optimized) version of beamcore
% remains here for demonstration purposes only

global bolts
weight=zeros(size(Lp));

for i=1:size(Lp,3)
    if(bolts.stop)  % if cancel button on waitbar is pressed...
        return
    end
    LpInvRzz = Lp(:,:,i)'*InvRzz;
    J=inv(LpInvRzz*Lp(:,:,i));    % 2*2 matrix
    weight(:,:,i)=(InvES'*Lp(:,:,i))*J(:,1:2); % Nsensors*2 matrix

    omega=J*(LpInvRzz * LpInvRzz')*J;

    weight(:,1,i)=weight(:,1,i)/real(sqrt(omega(1,1)));
    weight(:,2,i)=weight(:,2,i)/real(sqrt(omega(2,2)));

    if(~mod(i,500))  % update progress every 100th iteration
        if(ndefaults.bf.viewwaitbar)
            waitbar(i/size(Lp,3),bolts.barhandle)
        end
    end
end
return; %---------------------------------------------------------------------------------


