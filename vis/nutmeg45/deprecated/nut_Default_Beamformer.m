function [weight]=nut_Default_Beamformer(Lp,InvRzz, InvES, flags) %---------------------------------------------------------
% Lp : lead field ( channel X 3 )
% InvRzz : Inverse data covariance matrix
% InvES : Inverse signal-component-data covariance matrix

global bolts ndefaults

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
   condwarning = true;
   
   if(isunix) % suppress condition warning for knowledgeable users :)
      [crap,whoami]=unix('whoami');
      if(any(strcmp(whoami(1:(end-1)),{'sarang','johannaz'})))
         condwarning = false;
      end
   end
   
   % cond_meaningless must match logcond_meaningless found in nut_beamforming_gui.m
   cond_meaningless = 1e19;
   if(condwarning && cond(bolts.params.Rzz1) >= cond_meaningless)
      warndlg('Your results will be meaningless. Get that condition number down! (Choose a bigger time window, remove lowpass filters, and/or discard channels.)')
      %return
   end
   
%    if(cn) % column normalization ( lead field normalization )
%       for i=1:size(Lp,3)
%          Lp(:,:,i) = Lp(:,:,i)./repmat(nut_colnorm(Lp(:,:,i)),[size(Lp,1) 1]);
%       end
%    end
   
   
   % vectorized section -- we can compute Lp'*InvRzz and InvES*Lp for
   % each voxel in one big matrix -----------
   Lp2 = reshape(Lp,size(Lp,1),size(Lp,2)*size(Lp,3));
   
   InvRzzLp = reshape(InvRzz*Lp2,size(Lp,1),size(Lp,2),size(Lp,3));
   if(ndefaults.bf.viewwaitbar)
       waitbar(0.25,bolts.barhandle)
   end

   if(bolts.stop)  % if cancel button on waitbar is pressed...
       return
   end

   InvESLp = reshape(InvES*Lp2,size(Lp,1),size(Lp,2),size(Lp,3));
   if(ndefaults.bf.viewwaitbar)
       waitbar(0.5,bolts.barhandle)
   end

   clear Lp2;
   % end vectorized section ------------
   
   weight=zeros(size(Lp));
   for i=1:size(Lp,3)
      if(bolts.stop)  % if cancel button on waitbar is pressed...
         return
      end
      
      % equivalent to: inv(Lp' * InvRzz* Lp) 
      
      invJ = Lp(:,:,i)'*InvRzzLp(:,:,i);
      
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
         omega=J*(InvRzzLp(:,:,i)' * InvRzzLp(:,:,i))*J;
         
         % real is needed for MATLAB acceleration --
         % potential complex output of sqrt operation disables acceleration!
         for jj=1:size(weight,2)
             weight(:,jj,i)=weight(:,jj,i)/real(sqrt(omega(jj,jj)));
         end
%          weight(:,2,i)=weight(:,2,i)/real(sqrt(omega(2,2)));
%          weight(:,3,i)=weight(:,3,i)/real(sqrt(omega(3,3)));
      end
      
      if(~mod(i,500))  % update progress every 500th iteration
          if(ndefaults.bf.viewwaitbar)
         waitbar(0.5 + i/(2*size(Lp,3)),bolts.barhandle)
          end
      end
   end
% end
