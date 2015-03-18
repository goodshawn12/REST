function C = nut_rowmult_jz(A,B)
% C = NUT_ROWMULT(A,B)
%
% allows multiplication of Nx1 array A with NxM matrix B
% e.g., C = [A.*B(:,1) A.*B(:,2) A.*B(:,3)]

M=size(B,2);
for ii=1:M
    C(:,ii)=A.*B(:,ii);
end

