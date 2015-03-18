function [L,eta] = nut_dipoleorientation(Lp,data,algo)
% [LpScal,eta] = nut_dipoleorientation(Lp,data,algo)
% finds optimal dipole orientation from data and makes leadfield scalar.
% LpScal    Scalarized leadfield
% eta       optimal dipole orientations
% Lp        vector leadfield
% data      must contain data.InvRyy, i.e., the inv covariance matrix

disp('Calculating optimal dipole orientation at each voxel...')

if nargin<3, algo='SAM'; end
Lp2 = reshape(Lp,size(Lp,1),size(Lp,2)*size(Lp,3));

switch algo
    case {'SAM' 'LCMV'}
        InvRyyLp = reshape(data.InvRyy*Lp2,size(Lp,1),size(Lp,2),size(Lp,3));
        clear Lp2;

        L = zeros(size(Lp,1),1,size(Lp,3));
        eta = zeros(size(Lp,3),3);
        for i=1:size(Lp,3)
            invJ = Lp(:,:,i)'*InvRyyLp(:,:,i);
            J = nut_pinv(invJ);
            [vvv,ddd] = eig(J*InvRyyLp(:,:,i)'*InvRyyLp(:,:,i));
            [jnk,mineig]=min(diag(ddd));
            eta(i,:)=vvv(:,mineig);
            L(:,1,i) = Lp(:,:,i)*eta(i,:)';
        end
        
    case 'MinNorm'
        G = Lp2*Lp2'; clear Lp2;
        InvG=inv(G+data.gamma*eye(size(G)));
        RGG = InvG * data.Ryy * InvG;
        L = zeros(size(Lp,1),1,size(Lp,3));
        eta = zeros(size(Lp,3),3);
        for i=1:size(Lp,3);
            [vvv,ddd] = eig(Lp(:,:,i)' * RGG * Lp(:,:,i));
            [jnk,maxeig]=max(diag(ddd));
            eta(i,:)=vvv(:,maxeig);
            L(:,1,i) = Lp(:,:,i) * eta(i,:)';
        end
        
    otherwise
        error('Dipole orientation for this inverse solution is not programmed.')
end
            
        