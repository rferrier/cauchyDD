function [Atilde] = closestDeficient (A, niter)
% THis function computes the closest rank 2 matrix for A
% input  : A      : squared matrix
%          niter  : nb of fixed point iterations
% output : Atilde : matrix

u1 = ones(size(A,1),1); u2 = ones(size(A,1),1);
v1 = ones(size(A,1),1); v2 = ones(size(A,1),1);
for i=1:niter
   u1 = A*v1/(v1'*v1);
   v1 = A'*u1/(u1'*u1);
end
Atest = A - u1*v1';
for i=1:niter
   u2 = Atest*v2/(v2'*v2);
   v2 = Atest'*u2/(u2'*u2);
end
Atilde = u1*v1'+u2*v2';