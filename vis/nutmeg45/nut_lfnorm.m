function Lp = nut_lfnorm(Lp);
% Lp_norm = nut_lfnorm(Lp);
% performs columnwise normalization of the lead field

for i=1:size(Lp,3)
    Lp(:,:,i) = Lp(:,:,i)./repmat(nut_colnorm(Lp(:,:,i)),[size(Lp,1) 1]);
end
