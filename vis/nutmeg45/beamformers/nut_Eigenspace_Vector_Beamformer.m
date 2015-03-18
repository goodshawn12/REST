function [weight]=nut_Eigenspace_Vector_Beamformer(Lp,data,flags) %---------------------------------------------------------
% Formally known as nut_Default_Beamformer
% Lp: lead field ( channel X orientations )
% data.InvRyy: data covariance matrix
% data.InvES : Inverse signal-component-data covariance matrix
% flags.progressbar

% global bolts
% global ndefaults


% stop=false; % keep track of cancel button hit, start off false
if(flags.progressbar)
    disp(['Please wait. Eigenspace Vector Beamformer has begun']);
end

% condwarning = true;
% 
% if(isunix) % suppress condition warning for knowledgeable users :)
%     [crap,whoami]=unix('whoami');
%     if(any(strcmp(whoami(1:(end-1)),{'sarang','zumer'})))
%         condwarning = false;
%     end
% end

% % cond_meaningless must match logcond_meaningless found in nut_beamforming_gui.m
% cond_meaningless = 1e19;
% if(condwarning && data.condRyy >= cond_meaningless)
%     warndlg('Your results will be meaningless. Get that condition number down! (Choose a bigger time window, remove lowpass filters, and/or discard channels.)')
%     %return
% end



% vectorized section -- we can compute Lp'*InvRzz and InvES*Lp for
% each voxel in one big matrix -----------
Lp2 = reshape(Lp,size(Lp,1),size(Lp,2)*size(Lp,3));

InvRyyLp = reshape(data.InvRyy*Lp2,size(Lp,1),size(Lp,2),size(Lp,3));

InvESLp = reshape(data.InvES*Lp2,size(Lp,1),size(Lp,2),size(Lp,3));

clear Lp2;
% end vectorized section ------------

weight=zeros(size(Lp));
for i=1:size(Lp,3)

    % equivalent to: inv(Lp' * data.Ryy* Lp)

    invJ = Lp(:,:,i)'*InvRyyLp(:,:,i);

    % We use pinv here to accomodate single-sphere head models (which
    % will be of rank 2 and fail with a true inverse).
    % For more sophisticated models, the rank will generally be 3, and
    % then pinv yields the same results as inv.
    J=double(nut_pinv(invJ));

    % note InvRzz1 = InvRzz' and InvES = InvES';
    weight(:,:,i)=InvESLp(:,:,i)*J; % Nsensors*3 matrix

    if(flags.wn)
        % equivalent to: inv(Lp' * InvRzz* Lp) * (Lp' * InvRzz^2 * Lp) * inv(Lp' * InvRzz * Lp)
        % Lp' * InvRzz^2 * Lp = (Lp' * InvRzz) * (InvRzz*Lp) = (InvRzz*Lp)' * (InvRzz*Lp)
        omega=J*(InvRyyLp(:,:,i)' * InvRyyLp(:,:,i))*J;

        % real is needed for MATLAB acceleration --
        % potential complex output of sqrt operation disables acceleration!
        for jj=1:size(weight,2)
            weight(:,jj,i)=weight(:,jj,i)/real(sqrt(omega(jj,jj)));
        end
        %          weight(:,2,i)=weight(:,2,i)/real(sqrt(omega(2,2)));
        %          weight(:,3,i)=weight(:,3,i)/real(sqrt(omega(3,3)));
    end

    if(flags.progressbar)
        if(~mod(i,500))  % update progress every 500th voxel
            %         if(ndefaults.bf.viewwaitbar)
            %             waitbar(iv/nv,h);
            %         end
            disp(['Please wait. Eigenspace Vector Beamformer has finished ' num2str(i) ' out of ' num2str(size(Lp,3)) ' voxels.']);
        end
    end

end
% end
