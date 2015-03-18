function [weight]=nut_agmnrug(Lp,data, flags) %---------------------------------------------------------
%function [S] = rens1(L,B,gamma)

%
% Nc : number of sensors
% Np : number of voxels
% Nt : number of time points
%
% L(Nc,Np,2) :lead field
% B(Nc,Nt)   :sensor data
% gamma=1e-6:relative regularization parameter
%
% Empirically, the optimum value of this gamma is given by
%        gamma=(mean(lambda_noise)/lambda_max)*10,
%  where
%        the noise-level eigenvalues of hRbb are denoted lambda_noise,
%        the maximum eigenvalue of hRbb is denoted lambda_max,
%    and hRbb is the sample covariance matrix obtained using hRbb=B*B'/Nt.
%  Here, if you are unsure about the exact threshold between the signal- and
%  the noise-level eigenvalues, just overestimate it.
%
% S(1:Np,1:Nt,2):source reconstruction results
%
% repeat : number of iteration

warning('this algo is not in a polished and tested state!')

repeat=8;

gamma = 1e-6;

L = permute(Lp,[1 3 2]);
B = data.meg';


Nc = size(L,1);
Np = size(L,2); 
Nt = size(B,2);
S = zeros(Np,Nt,3); % 

% lead field normalization
Ln = zeros(size(L));
for index=1:Np
    Ln(:,index,1) = L(:,index,1)/norm(squeeze(L(:,index,1)));
    Ln(:,index,2) = L(:,index,2)/norm(squeeze(L(:,index,2)));
    Ln(:,index,3) = L(:,index,3)/norm(squeeze(L(:,index,3)));
end

% RENS initial step
h = waitbar(0,'initial step');
G=squeeze(L(:,:,1))*squeeze(L(:,:,1))'+ squeeze(L(:,:,2))*squeeze(L(:,:,2))' + squeeze(L(:,:,3))*squeeze(L(:,:,3))';
InvG=inv(G+gamma*max(eig(G))*eye(size(G)));
for npix=1:Np
    lpf=squeeze( Ln(:,npix,1) );
    lpg=squeeze( Ln(:,npix,2) );
    lph=squeeze( Ln(:,npix,3) );
    Lp=[lpf,lpg,lph];
    WT=inv(Lp'*InvG*Lp)*Lp'*InvG;
    Srec=WT*B;
    S(npix,:,:)=Srec';
    waitbar(npix/Np,h)
end
close(h)
     
% RENS recursive steps
h = waitbar(0,'recursion process');
for nt =1:Nt
    for tempI=1:repeat
        diagI11 = sparse(Np,Np);
        diagI22 = sparse(Np,Np);
        diagI33 = sparse(Np,Np);
        for ip =1:Np
            diagI11(ip,ip) = S(ip,nt,1)^2;
            diagI22(ip,ip) = S(ip,nt,2)^2;
            diagI33(ip,ip) = S(ip,nt,3)^2;
        end
        G=squeeze(L(:,:,1))*diagI11*squeeze(L(:,:,1))'+ squeeze(L(:,:,2))*diagI22*squeeze(L(:,:,2))' + squeeze(L(:,:,3))*diagI33*squeeze(L(:,:,3))';
        InvG=inv(G+gamma*max(eig(G))*eye(size(G)));
        for npix=1:Np
            lpf=squeeze( Ln(:,npix,1) );
            lpg=squeeze( Ln(:,npix,2) );
            lph=squeeze( Ln(:,npix,3) );
            Lp=[lpf,lpg,lph];
            WT=inv(Lp'*InvG*Lp)*Lp'*InvG;
            Srec=WT*B(:,nt);
            S(npix,nt,:)=Srec';
        end
%         waitbar(nt/Nt,h,['recursion process:
%         iteration',num2str(tempI),'/',num2str(repeat)]);
       waitbar(nt/Nt,h,'recursion process');
    end
end
close(h)
weight = zeros(size(Lp));
save Stmp S
