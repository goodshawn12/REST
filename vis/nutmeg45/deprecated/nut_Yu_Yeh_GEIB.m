function [weight]=nut_Yu_Yeh_GEIB(Lp,InvRzz, InvES, flags) %---------------------------------------------------------
% Lp : lead field ( channel X 3 )
% InvRzz : Inverse data covariance matrix
% InvES : Inverse signal-component-data covariance matrix
% Note that weight normalization is the only option here!

global bolts

% if nargin == 1
%    handles = Lp;
%    lh = [handles.nut_eigen_axes; get(handles.nut_eigen_axes,'Children')];
%    viewsingle_num = str2num(get(handles.nut_viewsingle_text,'String'));
%    if length(viewsingle_num) ~= 1 || viewsingle_num ~= 0 % don't make it visible unless it should be visible
%       lh = [lh; handles.nut_inset_axes; get(handles.nut_inset_axes,'Children')];
%    end
%    lh = [lh; handles.nut_signalspace_text; handles.nut_viewsingle_text; handles.text7'];
%    lh = [lh; handles.nut_logplot_box; handles.nut_cumulative_box];
%    set(lh,'Visible','on')
% else
   disp('Yu. Yeh.');
%    cn = 1  % force column normalization for now
%    wn = 0
%    warning('forcing column normalization!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');

   % vectorized section -- we can compute Lp'*InvRzz and InvES*Lp for
   % each voxel in one big matrix -----------
   Lp2 = reshape(Lp,size(Lp,1),size(Lp,2)*size(Lp,3));
   
   InvRzzLp = reshape(InvRzz*Lp2,size(Lp,1),size(Lp,2),size(Lp,3));
   
   waitbar(0.25,bolts.barhandle)
   
   if(bolts.stop)  % if cancel button on waitbar is pressed...
      return
   end
   
   Es = bolts.params.u(:,bolts.params.signalspace);
   EsEs = Es*Es';
   
   % InvESLp = reshape(InvES*Lp2,size(Lp,1),size(Lp,2),size(Lp,3));
   waitbar(0.5,bolts.barhandle)
   
   clear Lp2;
   % end vectorized section ------------
   
   weight=zeros(size(Lp,1),2,size(Lp,3));
   for i=1:size(Lp,3)
      if(bolts.stop)  % if cancel button on waitbar is pressed...
         return
      end
      
      if(flags.YYcn) % column normalization ( lead field normalization )
         Lp(:,:,i) = Lp(:,:,i)./repmat(nut_colnorm(Lp(:,:,i)),[size(Lp,1) 1]);
      end
      
      % equivalent to: inv(Lp' * InvRzz* Lp) 
      J=inv(Lp(:,:,i)'*InvRzzLp(:,:,i));    % 2*2 matrix
      
      w_cp = InvRzzLp(:,:,i)*J(:,1:2);
      
      % note InvRzz1 = InvRzz' and InvES = InvES';
      weight(:,:,i)=EsEs*w_cp; % Nsensors*2 matrix
      
      if(1)   % screwy results with this normalization!
         % equivalent to: inv(Lp' * InvRzz* Lp) * (Lp' * InvRzz^2 * Lp) * inv(Lp' * InvRzz * Lp)
         % Lp' * InvRzz^2 * Lp = (Lp' * InvRzz) * (Lp' * InvRzz)'
         
         % omega=Lp(:,:,i)'*weight(:,:,i);
         omega=Lp(:,:,i)'*InvRzz*w_cp;
         
         
         weight(:,1,i)=weight(:,1,i)/omega(1,1);
         weight(:,2,i)=weight(:,2,i)/omega(2,2);
      end
      
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
         waitbar(0.5 + i/(2*size(Lp,3)),bolts.barhandle)
      end
   end
% end
