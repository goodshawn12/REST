function colnorm = nut_colnorm(A)
% norm = colnorm(A)
% treat each column of A as vector and take norm
colnorm = sqrt(sum(A.^2));
