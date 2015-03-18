%%
% Performs a thin-plate spline warp based on Bookstein, 1989.
% @param Cp the coordinates of the original reference points.
% @param Cg the coordinates of the mapped reference points.
% @param XYZ the points to be warped.
% @return XYZp the warped points.
% @author Daniel D.E. Wong
%%
function XYZp = tps_warp(Cp,Cg,XYZ)
warning off

NPs = size(Cp,1);

Xp = Cp(:,1)';
Yp = Cp(:,2)';
Zp  = Cp(:,3)';
Xg = Cg(:,1)';
Yg = Cg(:,2)';
Zg = Cg(:,3)';

rXp = repmat(Xp(:),1,NPs);
rYp = repmat(Yp(:),1,NPs);
rZp = repmat(Zp(:),1,NPs);

wR = sqrt((rXp-rXp').^2 + (rYp-rYp').^2 + (rZp-rZp').^2);

wK = 2*(wR.^2).*log(wR+1e-20);
wP = [single(ones(NPs,1)) Xp(:) Yp(:) Zp(:)];
wL = [wK wP;wP' single(zeros(4,4))];
wY = [Xg(:) Yg(:) Zg(:); single(zeros(4,3))];
wW = inv(wL)*wY;

XYZp = single(zeros(size(XYZ)));
for coord = 1:3
    V = wW(:,coord)';
    w = V(1:size(Cp,1));
    a_vec = V(size(Cp,1)+1:end);
    XYZp(:,coord) = sum([ones(size(XYZ,1),1) XYZ].*repmat(a_vec,size(XYZ,1),1),2);
    for i = 1:size(Cp,1)
        r_i = sqrt(sum((repmat(Cp(i,:),size(XYZ,1),1)-XYZ).^2,2));
        wU_i = w(i).*r_i.^2.*log(r_i.^2);
        wU_i(find(isnan(wU_i))) = 0;
        XYZp(:,coord) = XYZp(:,coord)+wU_i;
    end
end