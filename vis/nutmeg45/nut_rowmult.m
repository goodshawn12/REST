function C = nut_rowmult(A,B)
% C = NUT_ROWMULT(A,B)
%
% allows multiplication of Nx1 array A with Nx3 matrix B
% e.g., C = [A.*B(:,1) A.*B(:,2) A.*B(:,3)]

C = [A.*B(:,1) A.*B(:,2) A.*B(:,3)];
