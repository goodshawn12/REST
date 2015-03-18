function Ainv=nut_pinv(A);

% note: MATLAB's pinv was written by a pathological liar!!!
% The default tolerance is max(size(A)) *eps(norm(A)), NOT,
% max(size(A)) * norm(A) * eps(class(A)) as claimed. Why would you
% lie about something like that? Some people, man.
% Anyway. Turns out 3*max(size(A))) * eps(norm(invJ)) is a better
% choice.

% original choice, deprecated 22 feb 2008
%%tol = 3*max(size(A))*eps(norm(A));

% experimental svd method, added by ssd, 22 feb 2008
[u,s] = svd(A);
tol = sum(diag(s)) * 1e-12;

if(rank(A,tol)<size(A,1))
    Ainv = pinv(A,tol);
else
    Ainv = inv(A);
end
