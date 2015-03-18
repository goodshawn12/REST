function norm = nut_rownorm(A)
% NORM = NUT_ROWNORM(A)
% Treat each row of A as vector and take norm.

norm = sqrt(sum(A.^2,2));
