function Rk = fatNull(M)
% This function computes the null space of a FAT matrix (code from Pierre Gosselet)

 epsi=1.e-14; % criterion for the kernel
 
 try % Handle strange error in LU (too many arguments)
    [Ll, Uu, Pp, Qq, Rr] = lu (M);  % P * (R \ M) * Q = L * U
 catch
    warning('Il LU, could not use Q and R matrices');
    [Ll, Uu, Pp] = lu (M);  % P * (R \ M) * Q = L * U
%     Rr = eye(size(M,1)); % Actually, it's useless
    Qq = eye(size(M,2));
 end
 
 
 
 if (norm(diag(Ll)-1,'inf')~=0), warning('diag L has 0 values'); end
 if (size(Uu,1)~=size(Uu,2)), warning('U is not squared'); end
 
 % Find the kernel of U 
 z1 = find(abs(diag(Uu))<=epsi);
 z1p = setdiff([1:size(Uu,1)]',z1);% z1p is invertible,
 
 % Schur complement on z1
 UU = Uu(z1p,z1p)\Uu(z1p,z1);
 U2 = Uu(z1,z1)-Uu(z1,z1p)*UU;
 N2 = null(full(U2));
 Rk  = zeros(size(Uu,1),size(N2,2));
 Rk(z1,:) = N2;
 Rk(z1p,:) = -UU*N2;
 Rk = Qq*Rk; % operate the permutations

end
