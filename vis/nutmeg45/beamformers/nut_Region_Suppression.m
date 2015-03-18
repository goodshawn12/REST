function [weight]=nut_Region_Suppression(Lp,data, flags) %---------------------------------------------------------
% HIGHLY EXPERIMENTAL!

% Lp : lead field ( channel X 2 )
% data.InvRyy : Inverse data covariance matrix
% data.InvES : Inverse signal-component-data covariance matrix
% flags.wn : 1 to enablee weight normalization, 0 to disable

% unvectorized (but matlab-optimized) version of beamcore

% cn=1
% wn=0
% warning('forcing column normalization!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');

% to do: remove global dependence!
global nuts bolts

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
   bolts.params.suppressiontype = 1;
   if(isfield(bolts,'signalspace_C'))
      signalspace_C = bolts.signalspace_C; %should make this flags.signalspace_C
   end
   
   
   disp('suppression baby!');
   
   if(0)
      % fit to desired MEG voxel grid, shift such that origin is (0,0,0)mm
      for i = 1:3
         VOIvoxels(:,i) = nuts.voxelsize(i)*round(VOIvoxels(:,i)/nuts.voxelsize(i)) + st.vols{1}.mat(i,4);
      end
      
      % remove duplicate voxels from fit, keep track of mapping with voxelindex
      [r_c crap voxelindex_sup] = unique(VOIvoxels,'rows');
      
      clear crap voxelindex_sup VOIvoxels lx1 lx2 ly1 ly2 voxels;
      
      % compute lead field of suppression region -- this should definitely be functionized
      [lx1c, ly1c] = nut_leadf(r_c,nuts.meg.sensorCoord(:,:,1),nuts.meg.sensorOrient,nuts.meg.lsc);
      [lx2c, ly2c] = nut_leadf(r_c,nuts.meg.sensorCoord(:,:,2),nuts.meg.sensorOrient,nuts.meg.lsc);
      
      if(isfield(nuts.meg,'refsensorpos'))
         [ref_lx1c, ref_ly1c] = nut_leadf(r_c,nuts.meg.refSensorCoord(:,:,1),ref_direction,nuts.meg.lsc);
         % some reference channels are plain magnetometers; their nuts.meg.refSensorCoord(:,:,2) = [0 0 0];
         % we should verify they result in ref_lx2 = ref_ly2 = 0 and make it so if it isn't
         [ref_lx2c, ref_ly2c] = nut_leadf(r_c,nuts.meg.refSensorCoord(:,:,2),ref_direction,nuts.meg.lsc);
      else  % no reference/gradient correction, so set values to 0
         [ref_lx1c,ref_lx2c,ref_ly1c,ref_ly2c, Gcoef] = deal(0);
      end
      
      Cx = lx1c-lx2c - Gcoef*(ref_lx1c-ref_lx2c); clear lx1c lx2c ref_lx1c ref_lx2c;
      Cy = ly1c-ly2c - Gcoef*(ref_ly1c-ref_ly2c); clear ly1c ly2c ref_ly1c ref_ly2c s v; pack;
   else
       if(1)
          load voxels_sup.mat;
%           [LpC,r_c] = nut_compute_lead_field(nuts.voxelsize(1),VOIvoxels);
          [LpC,r_c] = nut_compute_lead_field(2,VOIvoxels);
          clear VOIvoxels;
          Cx = squeeze(LpC(:,1,:));
          Cy = squeeze(LpC(:,2,:));
       else
           load Lp_sup.mat
       end
   end
   
   
   
   method = 'orth';
   % compute weights
    switch(method)
        case 'generalizedeig'
            Lx = squeeze(nuts.Lp(:,1,:));
            Ly = squeeze(nuts.Lp(:,2,:));
            switch 3
                case 1
                    [uc,U,V,s,Q] = gsvd([Cx Cy]',[Lx Ly]');
                case 2
                    [U,uc,V,Q,s] = gsvd([Lx Ly]',[Cx Cy]');
%                 [U,uc,V,Q,s] = gsvd(bolts.params.Rzz1,[Cx Cy]');
                %%%% wait -- maybe replace [Lx Ly] with Rzz???

                case 3
%                     [uc,s] = eig([Lx Ly]*[Lx Ly]', [Cx Cy]*[Cx Cy]');
                   A=inv([Lx Ly]*[Lx Ly]'+1e-6*eye(275)*max(eig([Lx Ly]*[Lx Ly]')))*( [Cx Cy]*[Cx Cy]');
%                    A=inv([Cx Cy]*[Cx Cy]'+1e-6*eye(275)*max(eig([Cx Cy]*[Cx Cy]')))*( [Lx Ly]*[Lx Ly]');
%                     A=inv([Lx Ly]*[Lx Ly]') * ([Cx Cy]*[Cx Cy]');
%                     A=inv([Lx Ly]*[Lx Ly]')*( [Cx Cy]*[Cx Cy]');
                    [uc,s]=eig(A);
                    
                    s = real(s);
                case 4
                    [uc,s] = eig([Cx Cy]*[Cx Cy]', [Lx Ly]*[Lx Ly]');
                    s = real(s);
            end
        otherwise
            [uc,s]=svd([Cx Cy],0);
%             [uc,s]=eig([Cx Cy]*[Cx Cy]');
%             s=real(s);
        %    [uc,s,v]=svd([Cx Cy],0);
    end

    clear Cx Cy

   weight=zeros(size(Lp));
   
   waitbar(0.5,bolts.barhandle)
   
   if(~exist('signalspace_C','var'))
      q = diag(s);
%       q = q(isfinite(q));
      qmax = min(50,length(q));
      eigfig=figure;subplot(1,3,1); plot(q(1:qmax)/sum(q)*100,'o-');
      subplot(1,3,2); semilogy(q(1:qmax),'o-');
      subplot(1,3,3); plot(cumsum(q(1:qmax))/sum(q)*100,'o-');
      xlabel('Eigenspectrum'); ylabel('% of Total');
      saveas(eigfig,'eigsup.fig');
      
      
      prompt   = 'What''s the signalspace of your suppressed region?';
      title    = 'Input signalspace';
      lines = 1;
      eig95 = dsearchn(cumsum(q(1:qmax))/sum(q)*100,95);
      def{1}   = ['1:' num2str(eig95)];
      answer   = inputdlg(prompt,title,lines,def);
      if (isempty(answer)) msgbox('What''s wrong with you?!');error('user is whack.'); end;
      signalspace_C = str2num(cell2mat(answer));
      close(eigfig);
   end
   
   bolts.params.signalspace_C = signalspace_C;   % retain in beam structure for future reference
   
   switch(method)
      case 'simple'
         Es = [bolts.params.u(:,bolts.params.signalspace) uc(nuts.meg.goodchannels,signalspace_C)];
         disp('ah, the simple life.');
      case 'qr'
         [Es,r] = qr([bolts.params.u(:,bolts.params.signalspace) uc(nuts.meg.goodchannels,signalspace_C)],0);
         disp('queer.');
      case 'orth'
         Es = orth([bolts.params.u(:,bolts.params.signalspace) uc(nuts.meg.goodchannels,signalspace_C)]);
         disp('orth orth!');
      case 'generalizedeig'
         Es = orth([bolts.params.u(:,bolts.params.signalspace) uc(nuts.meg.goodchannels,signalspace_C)]);
%           Es = [bolts.params.u(:,bolts.params.signalspace) uc(nuts.meg.goodchannels,signalspace_C)];
         disp('gen. eig.')
   end
   
   EsEs = Es*Es';
   
   % scale for augmentation of Lp
%    uc = uc(nuts.meg.goodchannels,signalspace_C)*s(signalspace_C,signalspace_C); %*v(:,signalspace_C)';
   
   switch(bolts.params.suppressiontype)
      case 1 % based on equation 24
         for i=1:size(Lp,3)   % loop through voxels
            if(bolts.stop)  % if cancel button on waitbar is pressed...
               return
            end
            LpC = [Lp(:,:,i) uc(nuts.meg.goodchannels,signalspace_C)];
            
            if(flags.RScn)  % column normalization ( lead field normalization )
               LpC = LpC./repmat(nut_colnorm(LpC),[size(LpC,1) 1]);
            end     % end column normalization
            
            InvRzzLpC = data.invRyy*LpC;
            
            % equivalent to: inv(Lp' * InvRzz* Lp) 
            J=inv(LpC'*InvRzzLpC);    % 2*2 matrix
            
            w_cp = InvRzzLpC*J(:,1:2);
            
            % note InvRzz1 = InvRzz' and InvES = InvES';
            weight(:,:,i)=EsEs*w_cp; % Nsensors*2 matrix
            
            
            % note InvRzz1 = InvRzz' and InvES = InvES';
            %             weight(:,:,i)=EsEs*InvRzzLp*J(:,1:2); % Nsensors*2 matrix
            
            if(0)   % screwy results with this normalization!
               % equivalent to: inv(Lp' * InvRzz* Lp) * (Lp' * InvRzz^2 * Lp) * inv(Lp' * InvRzz * Lp)
               % Lp' * InvRzz^2 * Lp = (Lp' * InvRzz) * (Lp' * InvRzz)'
               
               % omega=LpC'*weight(:,:,i);
               omega=Lp(:,:,i)'*data.invRyy*w_cp;
               
               weight(:,1,i)=weight(:,1,i)/omega(1,1);
               weight(:,2,i)=weight(:,2,i)/omega(2,2);
            end
            
            if(flags.wn)
               % equivalent to: inv(Lp' * InvRzz* Lp) * (Lp' * InvRzz^2 * Lp) * inv(Lp' * InvRzz * Lp)
               % Lp' * InvRzz^2 * Lp = (Lp' * InvRzz) * (InvRzz*Lp) = (InvRzz*Lp)' * (InvRzz*Lp)
               omega=J*(InvRzzLp' * InvRzzLp)*J;
               
               % real is needed for MATLAB acceleration --
               % potential complex output of sqrt operation disables acceleration!
               weight(:,1,i)=weight(:,1,i)/real(sqrt(omega(1,1)));
               weight(:,2,i)=weight(:,2,i)/real(sqrt(omega(2,2)));
            end
            
            if(~mod(i,500))  % update progress every 500th iteration
               waitbar(0.5 + i/(2*size(Lp,3)),bolts.barhandle)
            end
         end
         
      case 2  % based on equation 44
         for i=1:size(Lp,3)   % loop through voxels
            if(bolts.stop)  % if cancel button on waitbar is pressed...
               return
            end
            LpC = [Lp(:,:,i) uc(nuts.meg.goodchannels,signalspace_C)];
            
            if(flags.RScn) % column normalization ( lead field normalization )
               LpC = LpC./repmat(nut_colnorm(LpC),[size(LpC,1) 1]);
            end % end column normalization
            
            InvRzzLp = data.invRyy*LpC;
            
            % equivalent to: inv(Lp' * InvRzz* Lp) 
            J=inv(LpC'*InvRzzLp);    % 2*2 matrix
            
            % note InvRzz1 = InvRzz' and InvES = InvES';
            weight(:,:,i)=data.InvES*LpC*J(:,1:2); % Nsensors*2 matrix
            
            if(flags.wn)
               % equivalent to: inv(Lp' * InvRzz* Lp) * (Lp' * InvRzz^2 * Lp) * inv(Lp' * InvRzz * Lp)
               % Lp' * InvRzz^2 * Lp = (Lp' * InvRzz) * (InvRzz*Lp) = (InvRzz*Lp)' * (InvRzz*Lp)
               omega=J*(InvRzzLp' * InvRzzLp)*J;
               
               % real is needed for MATLAB acceleration --
               % potential complex output of sqrt operation disables acceleration!
               weight(:,1,i)=weight(:,1,i)/real(sqrt(omega(1,1)));
               weight(:,2,i)=weight(:,2,i)/real(sqrt(omega(2,2)));
            end
            
            if(~mod(i,500))  % update progress every 500th iteration
               waitbar(0.5 + i/(2*size(Lp,3)),bolts.barhandle)
            end
         end
   end
% end

return; %---------------------------------------------------------------------------------
