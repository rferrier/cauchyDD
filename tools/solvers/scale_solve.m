function x = scale_solve(A,b,n,m)
% This function solves Ax = b with scaled lines and columns of A (and b)
% Useful for Lagrange multipliers
% input : A : matrix to invert
%         b : (M)rhs
%         n : size of the unmodified part on A
%         m : multiplier

A(n+1:end,:) = m*A(n+1:end,:);
A(:,n+1:end) = m*A(:,n+1:end);

b(n+1:end,:)   = m*b(n+1:end,:);

x = A\b;

x(n+1:end,:) = m*x(n+1:end,:);

end
