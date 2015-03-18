function C=nut_dotfast(A,B,dim)
% C = NUT_DOTFAST(A,B,DIM)
%
% faster dot product, treating each row (dim=1) or column (dim=2)
% as set of vectors to dot sequentially

C = sum(A.*B,dim);
