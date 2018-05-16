function [usol,indm2] = spd_cauchy(ndofs, nbloq1d, nbloq2d, K1d, K2d, fD, fN, Cu, niter, v, varargin)
% This function solves the Cauchy problem by Dual Steklov-Poincaré method

% input : ndofs    : nb of primal dofs in the problem
%         nbloq1d  : nb of dirichlet imposed dofs in the dual Dirichlet problem
%         nbloq2d  : nb of dirichlet imposed dofs in the dual Neumann problem
%         K1d      : dual Dirichlet matrix
%         K2d      : dual Neumann matrix
%         fD       : Dirichlet Rhs (on K1d) 
%         fN       : Neumann Rhs (on K2d)
%         Cu       : Cu'*u = u on the unknown
%         niter    : nb of iterations
%         v        : should the algorithm talk ?
%         varargin : if known, the reference force

% output : usol   : deplacement solution
%          fsol   : force solution

   nmbodes = 0;
   if numel(varargin)>0
      nbmodes = cell2mat(varargin(1));
   end

   kn = 0;
   if numel(varargin)>1
      kn = 1;
      fref = cell2mat(varargin(2));
   end

   % Restriction operator
   Cru = sparse(Cu*Cu');

   % Find the indices of the unknown nodes
   indexC = sum(Cu'); indexa = find(indexC);
   clear Cu;

   iszer = 0;

   Itere = zeros( ndofs, 2 );
   d     = zeros( ndofs, 2*(niter+1) );  % 2 directions per step
   Ad    = zeros( ndofs, 2*(niter+1) );
   Res   = zeros( ndofs, 2*(niter+1) );
   Zed   = zeros( ndofs, 2*(niter+1) );
   
   % Equivalent Saad matrices
   unsuralpha   = zeros( 2*(niter+1) );
   betasuralpha = zeros( 2*(niter+1) );
   etaeta       = zeros( 2*(niter+1) );
   
   % No need to recompute this one as it is 0
%   %% Perform A x0 :
%   % Solve 1
%   f1 = [Itere(:,1); zeros(nbloq1d,1)];
%   uin1 = K1d\f1;
%   u1i = uin1(1:ndofs,1);
%   u1 = Cru*u1i;
%   % Solve 2
%   f2 = [Itere(:,2); zeros(nbloq2d,1)];
%   uin2 = K2d\f2;
%   u2i = uin2(1:ndofs,1);
%   u2 = Cru*u2i;
%   %
%   Axz = -[u1,u2];
   
   %%%%
   %% Compute Rhs :
   % Solve 1
   uin1 = K1d\fD;
   u1i = uin1(1:ndofs);
   u1 = Cru*u1i;
   % Solve 2
   uin2 = K2d\fN;
   u2i = uin2(1:ndofs);
   u2 = Cru*u2i;
   %
   b = [u1-u2,u1+u2];
   if norm(u2)/norm(u1) < 1e-10
      iszer = 1;
      b(indexa,2) = norm(u1-u2)*randn( size(indexa,2), 1 ); % To avoid problems
   end
   %%
   Res(:,[1,2]) = b;% - Axz;

   Zed(:,[1,2]) = Res(:,[1,2]); % No preconditionner
   d(:,[1,2]) = Zed(:,[1,2]);
   
   residual(1) = norm( Res(indexa,1) );
   if kn == 1
      error(1)    = norm(Itere(indexa,1) - fref(indexa)) / ...
                                          norm(fref(indexa));
   end
%   regulari(1) = sqrt( (Itere(:,1)')* ... 
%                        regul( Itere(:,1) , nodes, boundary, 2) );
   %%
   V  = zeros(ndofs, 2*niter);
   AV = zeros(ndofs, 2*niter);
   MV = zeros(ndofs, 2*niter);
   H  = zeros(2*niter);
   
   num = [0,0;0,0]; % useless, but eta needs initialization #lazy
   den = [0,0;0,0];
   %%
   for iter = 1:niter

       %% Optimal step
       % Solve 1
       f1 = [d(:,[2*iter-1,2*iter]); zeros(nbloq1d,2)];
       uin1 = K1d\f1; u1i = uin1(1:ndofs,:);
       u1 = Cru*u1i; % Restrict on the unknown dofs
       % Solve 2
       f2 = [d(:,[2*iter-1,2*iter]); zeros(nbloq2d,2)];
       uin2 = K2d\f2; u2i = uin2(1:ndofs,:);
       u2 = Cru*u2i;
       %
       Ad(:,[2*iter-1,2*iter]) = u2-u1;

       denprec = den; numprec = num; % Store those ones
       den = d(indexa,[2*iter-1,2*iter])'*Ad(indexa,[2*iter-1,2*iter]);
       sqD = den^(1/2);
       d(:,[2*iter-1,2*iter]) = d(:,[2*iter-1,2*iter]) * inv(sqD);
       Ad(:,[2*iter-1,2*iter]) = Ad(:,[2*iter-1,2*iter]) * inv(sqD);
       num = Res(indexa,[2*iter-1,2*iter])'*Zed(indexa,[2*iter-1,2*iter]);
       num = sqD\num; % because of Zed and not d
       
       Itere = Itere + d(:,[2*iter-1,2*iter])*num;
       Res(:,[2*iter+1,2*iter+2]) = Res(:,[2*iter-1,2*iter]) - ...
                                       Ad(:,[2*iter-1,2*iter])*num;
       
       residual(iter+1) = norm( Res(indexa,2*iter+1) );
       if kn == 1
          error(iter+1)    = norm(Itere(indexa,1) - fref(indexa)) / ...
                                           norm(fref(indexa));
       end
%       regulari(iter+1) = sqrt( (Itere(:,1)')* ... 
%                             regul( Itere(:,1) , nodes, boundary, 2) );
       

       Zed(:,[2*iter+1,2*iter+2]) = Res(:,[2*iter+1,2*iter+2]); % No precond
       
       % Ritz variables (Saad++)
       unsuralpha( [2*iter-1,2*iter], [2*iter-1,2*iter] ) = ...
                                       (sqD*num)^(-1/2)*den*(sqD*num)^(-1/2);
       
       if iter > 1
          betasuralpha( [2*iter-3,2*iter-2], [2*iter-3,2*iter-2] ) = ...
                            (sqD*num)^(-1/2) * betaij'*betaij * (sqD*num)^(-1/2);
                            % use betaij from the previous iteration
                            
          etaeta( [2*iter-1,2*iter], [2*iter-1,2*iter] ) = ...
                  (denprec^(1/2)*numprec)^(-1/2) * denprec * ...
                  inv( denprec^(1/2)*numprec ) * (sqD*num)^(1/2);
       end
       
       % Reorthogonalize the residual (as we use it next), in sense of M
       for jter=1:iter-1
           betac = (Res(indexa,[2*jter-1,2*jter])'*Zed(indexa,[2*jter-1,2*jter])) \...
                   (Res(indexa,[2*jter-1,2*jter])'*Zed(indexa,[2*iter+1,2*iter+2])) ;
   
           Zed(:,[2*iter+1,2*iter+2]) = Zed(:,[2*iter+1,2*iter+2]) - ...
                                         Zed(:,[2*jter-1,2*jter]) * betac;
           Res(:,[2*iter+1,2*iter+2]) = Res(:,[2*iter+1,2*iter+2]) - ...
                                         Res(:,[2*jter-1,2*jter]) * betac;
       end
   
       %% Orthogonalization
       d(:,[2*iter+1,2*iter+2]) = Zed(:,[2*iter+1,2*iter+2]);
       for jter=iter:iter % Only the last one
           betaij = ( Ad(indexa,[2*jter-1,2*jter])'*d(indexa,[2*jter-1,2*jter]) ) \ ...
               ( Ad(indexa,[2*jter-1,2*jter])'*d(indexa,[2*iter+1,2*iter+2]) );
   
           d(:,[2*iter+1,2*iter+2]) = d(:,[2*iter+1,2*iter+2]) - ...
                                      d(:,[2*jter-1,2*jter]) * betaij;
       end
   
       %% The Ritz elements
       V(indexa,[2*iter-1,2*iter]) = (-1)^(iter-1) * Zed(indexa,[2*iter-1,2*iter]) * ...
         ((Res(indexa,[2*iter-1,2*iter])'*Zed(indexa,[2*iter-1,2*iter]))^(-1/2)) ;
                          
       delta  = unsuralpha( [2*iter-1,2*iter], [2*iter-1,2*iter] ) ;
       if iter > 1
          delta = delta + betasuralpha( [2*iter-3,2*iter-2], [2*iter-3,2*iter-2] );
       end
   
       eta = etaeta( [2*iter-1,2*iter], [2*iter-1,2*iter] ); % what a stupid variable name
       
       if iter > 1
          H( [2*iter-1,2*iter] , [2*iter-3,2*iter-2,2*iter-1,2*iter] ) = ... 
                                                                 [eta', delta];
          H( [2*iter-3,2*iter-2] , [2*iter-1,2*iter] ) = eta;
       else
          H( [2*iter-1,2*iter] , [2*iter-1,2*iter] ) = delta;
       end
   end
   
   % Compute eigenelems of the Hessenberg :
   [Q,Theta1] = eig(H);
   theta = diag(Theta1);
   % Sort it
   [theta,Ind] = sort(theta,'descend');
   Q = Q(:,Ind);
   Theta1 = Theta1(Ind,Ind);
   Y = V*Q;
   
   % Build residual and such
   regD = zeros(niter,1); resD = zeros(niter,1); bt = Y'*b;
   for i = 1:2*iter+1
      chiD   = inv(Theta1)*Y'*(b(:,1)); chiD(i:end) = 0;
      %ItereD = Y*chiD;
      %resD(i) = sqrt( sum( bt(i:end).^2) );  
      %regD(i) = sqrt( ItereD'*regul(ItereD, nodes, boundary, 2) );
   end
   
   % Compute the number of pertinent Ritz modes
   [indm2,pol] = findPicard2 (log10(abs(chiD)), ceil(niter/5), 1, 3);%findPicard2(log10(abs(chiD)), 3, 1);
   n = size(pol,1);
   t = 1:.05:2*niter; tt = zeros(n,20*(2*niter-1)+1);
   for j=1:n
      tt(j,:) = t.^(n-j);
   end
   px = pol'*tt;
   
   
   % Compute the solution
   chi = inv(Theta1)*Y'*( b(:,1) );
   if nbmodes > 0
      chi(nbmodes:end) = 0;
   else
      chi(indm2:end) = 0;
   end
   ItereR = Y*chi;

   if v == 1 % Outputs
      if kn == 1
         try
         figure;
         hold on
         semilogy(error,'Color','blue')
         semilogy(residual/residual(1),'Color','red')
         legend('error','residual')
         end
      else
         try
         figure;
         semilogy(residual/residual(1),'Color','red')
         legend('residual')
         end
      end

      %% L-curve :
      %loglog(residual(2:end),regulari(2:end));
      %figure
   
      try
      figure;
      hold on;
      plot(log10(theta),'Color','blue')
      plot(log10(abs(Y'*( b(:,1) ))),'Color','red')
      plot(log10(abs(chiD)),'Color','black')
      plot(t,px,'Color','cyan')
      legend( 'Ritz Values','RHS values','solution coefficients', ...
              'polynomial approximation' );
      end
   end

   %%%%%%%
   %% Final resolution
   f_D = fD; % Recover the Dirichlet Rhs where it is avaliable
   f_N = [ItereR; zeros(nbloq1d,1)]; % Use the reconstructed rhs elsewhere
   uin  = K1d\(f_N+f_D);

   usol = uin(1:ndofs,1);
end
