function C = nut_arrayadd(A,B)
% C = arrayadd  (A,B)
% add M x 1 array from each column of M x N matrix A
 
C = A+repmat(B,[size(A,1) 1]);
 
