function [weight]=nut_Point_Suppression(Lp,data, flags) %---------------------------------------------------------
% HIGHLY EXPERIMENTAL!

% Lp : lead field ( channel X 2 )
% data.InvRyy : Inverse data covariance matrix
% data.InvES : Inverse signal-component-data covariance matrix
% flags.wn : 1 to enablee weight normalization, 0 to disable


% global nuts bolts

   disp('suppression baby!');
   
   prompt={'Enter coordinates of the point to be suppressed:'};
   name='Point Suppression';
   defaultanswer={'[0 30 40]'};
   answer=inputdlg(prompt,name,1,defaultanswer);
   r_c0 = str2num(answer{1});
   
   %%%%% suppression lead field
   % r_c0 = [30 30 0];   % suppressed point
   
   [LpC,r_c] = nut_compute_lead_field(2,r_c0);
   Lp = [Lp repmat(LpC,[1 1 size(Lp,3)])];
   
   % compute weights
   weight=zeros(size(Lp));
   
%    waitbar(0.5,bolts.barhandle)
   
   
   % vectorized section -- we can compute Lp'*InvRzz and InvES*Lp for
   % each voxel in one big matrix -----------
   Lp2 = reshape(Lp,size(Lp,1),size(Lp,2)*size(Lp,3));
   
   InvRzzLp = reshape(data.InvRyy*Lp2,size(Lp,1),size(Lp,2),size(Lp,3));
%    waitbar(0.25,bolts.barhandle)
   
%    if(bolts.stop)  % if cancel button on waitbar is pressed...
%       return
%    end
   
   InvESLp = reshape(data.InvES*Lp2,size(Lp,1),size(Lp,2),size(Lp,3));
%    waitbar(0.5,bolts.barhandle)
   
   clear Lp2;
   % end vectorized section ------------
   
   weight=zeros(size(Lp,1),2,size(Lp,3));
   for i=1:size(Lp,3)
%       if(bolts.stop)  % if cancel button on waitbar is pressed...
%          return
%       end
      
      % equivalent to: inv(Lp' * InvRzz* Lp) 
      J=inv(Lp(:,:,i)'*InvRzzLp(:,:,i));    % 2*2 matrix
      
      % note InvRzz1 = InvRzz' and InvES = InvES';
      weight(:,:,i)=InvESLp(:,:,i)*J(:,1:2); % Nsensors*2 matrix (original method)
      %weight(:,:,i)=InvRzzLp(:,:,i)*J; % Nsensors*2 matrix (alternative method)
      
      if(flags.wn)
         % equivalent to: inv(Lp' * InvRzz* Lp) * (Lp' * InvRzz^2 * Lp) * inv(Lp' * InvRzz * Lp)
         % Lp' * InvRzz^2 * Lp = (Lp' * InvRzz) * (InvRzz*Lp) = (InvRzz*Lp)' * (InvRzz*Lp)
         omega=J*(InvRzzLp(:,:,i)' * InvRzzLp(:,:,i))*J;
         
         % real is needed for MATLAB acceleration --
         % potential complex output of sqrt operation disables acceleration!
         weight(:,1,i)=weight(:,1,i)/real(sqrt(omega(1,1)));
         weight(:,2,i)=weight(:,2,i)/real(sqrt(omega(2,2)));
      end
      
      if(~mod(i,500))  % update progress every 500th iteration
           disp(['Please wait. Point Suppression has finished ' num2str(i) ' out of ' num2str(size(Lp,3)) ' voxels.']);
%          waitbar(0.5 + i/(2*size(Lp,3)),bolts.barhandle)
      end
   end
% end

return; %---------------------------------------------------------------------------------
