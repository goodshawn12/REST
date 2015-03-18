function [weight]=nut_Beamspace(Lp,data, flags) %---------------------------------------------------------
% Lp : lead field ( channel X 3 )
% InvRzz : Inverse data covariance matrix
% InvES : Inverse signal-component-data covariance matrix
% Note that weight normalization is the only option here!
% based on Van Veen 

% global bolts

 
   
   % vectorized section -- we can compute Lp'*InvRzz and InvES*Lp for
   % each voxel in one big matrix -----------
   Lp2 = reshape(Lp,size(Lp,1),size(Lp,2)*size(Lp,3));
   G = Lp2*Lp2';
   
   [u,s,v]=svd(G);
   
   % take top 99% of eigenvalues
   [Tidx,b]=dsearchn(cumsum(diag(s))/sum(diag(s)),.9);
   Tidx
   T = u(:,1:Tidx);
   % T=u;
%    Rzz = bolts.params.Rzz1;
   
   InvTRyyT=inv(T'*data.Ryy*T);
   InvRyy = T*InvTRyyT*T';
   
   InvRyyLp = reshape(InvRyy*Lp2,size(Lp,1),size(Lp,2),size(Lp,3));

   InvESLp = reshape(data.InvES*Lp2,size(Lp,1),size(Lp,2),size(Lp,3));
   
   clear Lp2;
   % end vectorized section ------------
   
   weight=zeros(size(Lp,1),2,size(Lp,3));
   for i=1:size(Lp,3)
      
      % equivalent to: inv(Lp' * InvRzz* Lp) 
      J=inv(Lp(:,:,i)'*InvRyyLp(:,:,i));    % 2*2 matrix
      
      % note InvRzz1 = InvRzz' and InvES = InvES';
      weight(:,:,i)=InvESLp(:,:,i)*J; % Nsensors*2 matrix (original method)
      %weight(:,:,i)=InvRzzLp(:,:,i)*J; % Nsensors*2 matrix (alternative method)
      
      if(flags.wn)
         % equivalent to: inv(Lp' * InvRzz* Lp) * (Lp' * InvRzz^2 * Lp) * inv(Lp' * InvRzz * Lp)
         % Lp' * InvRzz^2 * Lp = (Lp' * InvRzz) * (InvRzz*Lp) = (InvRzz*Lp)' * (InvRzz*Lp)
         omega=J*(InvRyyLp(:,:,i)' * InvRyyLp(:,:,i))*J;
         
         % real is needed for MATLAB acceleration --
         % potential complex output of sqrt operation disables acceleration!
         weight(:,1,i)=weight(:,1,i)/real(sqrt(omega(1,1)));
         weight(:,2,i)=weight(:,2,i)/real(sqrt(omega(2,2)));
      end
  if(flags.progressbar)    
      if(~mod(i,500))  % update progress every 500th iteration
                      disp(['Please wait. Beamspace has finished ' num2str(i) ' out of ' num2str(size(Lp,3)) ' voxels.']);
%          waitbar(0.5 + i/(2*size(Lp,3)),bolts.barhandle)
      end
  end
   end
% end
