function [weight]=nut_LCMV_Vector_Beamformer(Lp,data, flags) %---------------------------------------------------------
% Lp : lead field ( channel X 3 x voxels)
% data.InvRyy : Inverse data covariance matrix
% flags.wn : 1 for weight normalization
% flags.progressbar : 1 to turn on progress report

%% provides weights for each vector direction with no orientation optimization

% global bolts

  
   
   % vectorized section -- we can compute Lp'*InvRzz and InvES*Lp for
   % each voxel in one big matrix -----------
   Lp2 = reshape(Lp,size(Lp,1),size(Lp,2)*size(Lp,3));
   
   InvRyyLp = reshape(data.InvRyy*Lp2,size(Lp,1),size(Lp,2),size(Lp,3));
   
   clear Lp2;
   % end vectorized section ------------
   
   weight=zeros(size(Lp));
   for i=1:size(Lp,3)
      
      % equivalent to: inv(Lp' * InvRzz* Lp) 
      J=inv(Lp(:,:,i)'*InvRyyLp(:,:,i));    % 2*2 matrix
      
      % note InvRzz1 = InvRzz'
      weight(:,:,i)=InvRyyLp(:,:,i)*J; % Nsensors*2 matrix
      
      if(flags.wn)
         % equivalent to: inv(Lp' * InvRzz* Lp) * (Lp' * InvRzz^2 * Lp) * inv(Lp' * InvRzz * Lp)
         % Lp' * InvRzz^2 * Lp = (Lp' * InvRzz) * (InvRzz*Lp) = (InvRzz*Lp)' * (InvRzz*Lp)
         omega=J*(InvRyyLp(:,:,i)' * InvRyyLp(:,:,i))*J;
         
         % real is needed for MATLAB acceleration --
         % potential complex output of sqrt operation disables acceleration!
         for jj=1:size(weight,2)
             weight(:,jj,i)=weight(:,jj,i)/real(sqrt(omega(jj,jj)));
         end
      end
      
          if(flags.progressbar)
      if(~mod(i,500))  % update progress every 500th iteration
            disp(['Please wait. LCMV Vector Beamformer has finished ' num2str(i) ' out of ' num2str(size(Lp,3)) ' voxels.']);
%          waitbar(0.5 + i/(2*size(Lp,3)),bolts.barhandle)
      end
          end
   end
% end
