function [weight]=nut_Original_Beamformer(Lp,InvRzz, InvES, flags) %---------------------------------------------------------
% Lp : lead field ( channel X 3 )
% InvRzz : Inverse data covariance matrix
% InvES : Inverse signal-component-data covariance matrix
% wn : flag for weight normalization
% cn : flat for column normalization( lead field normalization )

% It's the O.B. -- the original beamcore. slower than dubya.

global bolts;


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
   weight = zeros(size(Lp));
   
   for i=1:size(Lp,3)
      if(bolts.stop) % kill waitbar and bust out of function if user hit cancel
         return
      end;
      
%       if(cn) % column normalization ( lead field normalization )
%          Lp(:,:,i) = Lp(:,:,i)./repmat(nut_colnorm(Lp(:,:,i)),[size(Lp,1) 1]);
%       end
%       
      J=inv(Lp(:,:,i)'*InvRzz*Lp(:,:,i));    % 2*2 matrix
      weight(:,:,i)=(InvES'*Lp(:,:,i))*J; % Nsensors*2 matrix
      
      if(flags.wn) % weight normalization
         omega=J*(Lp(:,:,i)'*(InvRzz^2)*Lp(:,:,i))*J;  % note: InvRzz*InvRzz is much faster than InvRzz^2
         weight(:,1,i)=weight(:,1,i)/sqrt(omega(1,1));
         weight(:,2,i)=weight(:,2,i)/sqrt(omega(2,2));
         weight(:,3,i)=weight(:,3,i)/sqrt(omega(3,3));
      end;
      
      if(~mod(i,50))  % update progress every 50th iteration
         waitbar(i/size(Lp,3),bolts.barhandle)
      end
      
   end
% end

return; %---------------------------------------------------------------------------------
