function [usol,indm2] = spd_mrhs(ndofs, nbloq1d, nbloq2d, K1d, K2d, fD, fN, Cu, niter, v, varargin)
% This function solves many Cauchy problem by MRHS Dual Steklov-Poincaré method

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

   % Compute the multiplier
   %mul = norm(K1d(1:ndofs,1:ndofs),'fro');

   % Find the indices of the unknown nodes
   indexC = sum(Cu'); indexa = find(indexC);
   clear Cu;

   % Recover the nb of test-cases
   ncase = size(fD,2);
   if size(fN,2) ~= ncase
      warning('sizes of Dirichlet and Neumann data inconsistent');
   end

   Itere = zeros( ndofs, ncase );
   d     = zeros( ndofs, ncase*(niter+1) );  % ncase directions per step
   Ad    = zeros( ndofs, ncase*(niter+1) );
   Res   = zeros( ndofs, ncase*(niter+1) );
   Zed   = zeros( ndofs, ncase*(niter+1) ); % Yep, we store all that (my name is Leak, Memory Leak)
   
   % Equivalent Saad matrices
   unsuralpha   = zeros( ncase*(niter+1) );
   betasuralpha = zeros( ncase*(niter+1) );
   etaeta       = zeros( ncase*(niter+1) );
   
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
   uin1 = K1d\fD;%scale_solve(K1d,fD,ndofs,mul);%
   u1i = uin1(1:ndofs,:);
   u1 = Cru*u1i;
   % Solve 2
   uin2 = K2d\fN;%scale_solve(K2d,fN,ndofs,mul);%
   u2i = uin2(1:ndofs,:);
   u2 = Cru*u2i;
   %
   b = u1-u2;

   %%
   Res(:,[1:ncase]) = b;% - Axz;

   Zed(:,[1:ncase]) = Res(:,[1:ncase]); % No preconditionner
   d(:,[1:ncase]) = Zed(:,[1:ncase]);
   
   residual(1) = norm( Res(indexa,[1:ncase]), 'fro' );
   if kn == 1
      error(1)    = norm(Itere(indexa,:) - fref(indexa,:)) / norm(fref(indexa,:));
   end
%   regulari(1) = sqrt( (Itere(:,1)')* ... 
%                        regul( Itere(:,1) , nodes, boundary, 2) );
   %%
   V  = zeros(ndofs, ncase*niter);
   AV = zeros(ndofs, ncase*niter);
   MV = zeros(ndofs, ncase*niter);
   H  = zeros(ncase*niter);
   
   num = zeros(ncase); % useless, but eta needs initialization #lazy
   den = zeros(ncase);
   %%
   for iter = 1:niter

       multindex   = [ncase*(iter-1)+1:ncase*iter];
       multindexp1 = multindex + ncase;
       multindexm1 = multindex - ncase;

       %% Optimal step
       % Solve 1
       f1 = [d(:,multindex); zeros(nbloq1d,ncase)];
       uin1 = K1d\f1;%scale_solve(K1d,f1,ndofs,mul);%
       u1i = uin1(1:ndofs,:);
       u1 = Cru*u1i; % Restrict on the unknown dofs
       % Solve 2
       f2 = [d(:,multindex);zeros(nbloq2d,ncase)];
       uin2 = K2d\f2;%scale_solve(K2d,f2,ndofs,mul);%
       u2i = uin2(1:ndofs,:);
       u2 = Cru*u2i;
       %
       Ad(:,multindex) = u2-u1;

       denprec = den; numprec = num; % Store those ones
       den = d(indexa,multindex)'*Ad(indexa,multindex);

       if rank(den) ~= size(den,1)
          warning('Your test-cases are of lower dimension as their number');
       end

       sqD = den^(1/2);
       d(:,multindex) = d(:,multindex) * inv(sqD);
       Ad(:,multindex) = Ad(:,multindex) * inv(sqD);
       num = Res(indexa,multindex)'*Zed(indexa,multindex);
       num = sqD\num; % because of Zed and not d
       
       Itere = Itere + d(:,multindex)*num;
       Res(:,multindexp1) = Res(:,multindex) - Ad(:,multindex)*num;
       
       residual(iter+1) = norm( Res(indexa,multindexp1), 'fro' );
       if kn == 1
          error(iter+1)    = norm(Itere(indexa,:) - fref(indexa,:),'fro') / ...
                                           norm(fref(indexa,:),'fro');
       end
%       regulari(iter+1) = sqrt( (Itere(:,1)')* ... 
%                             regul( Itere(:,1) , nodes, boundary, 2) );

       Zed(:,multindexp1) = Res(:,multindexp1); % No precond
       
       % Ritz variables (Saad++)
       unsuralpha( multindex, multindex ) = (sqD*num)^(-1/2)*den*(sqD*num)^(-1/2);
       
       if iter > 1
          betasuralpha( multindexm1, multindexm1 ) = ...
                            (sqD*num)^(-1/2) * betaij'*betaij * (sqD*num)^(-1/2);
                            % use betaij from the previous iteration
                            
          etaeta( multindex, multindex ) = ...
                  (denprec^(1/2)*numprec)^(-1/2) * denprec * ...
                  inv( denprec^(1/2)*numprec ) * (sqD*num)^(1/2);
       end
       
       % Reorthogonalize the residual (as we use it next), in sense of M
       for jter=1:iter-1
           multjndex = [ncase*(jter-1)+1:ncase*jter];
           betac = (Res(indexa,multjndex)'*Zed(indexa,multjndex)) \...
                   (Res(indexa,multjndex)'*Zed(indexa,multindexp1)) ;
   
           Zed(:,multindexp1) = Zed(:,multindexp1) - Zed(:,multjndex) * betac;
           Res(:,multindexp1) = Res(:,multindexp1) - Res(:,multjndex) * betac;
       end
   
       %% Orthogonalization
       d(:,multindexp1) = Zed(:,multindexp1);
       for jter=iter:iter % Only the last one
           multjndex = [ncase*(jter-1)+1:ncase*jter];
           betaij = ( Ad(indexa,multjndex)'*d(indexa,multjndex) ) \ ...
               ( Ad(indexa,multjndex)'*d(indexa,multindexp1) );
   
           d(:,multindexp1) = d(:,multindexp1) - d(:,multjndex) * betaij;
       end
   
       %% The Ritz elements
       V(indexa,multindex) = (-1)^(iter-1) * Zed(indexa,multindex) * ...
         ((Res(indexa,multindex)'*Zed(indexa,multindex))^(-1/2)) ;
                          
       delta  = unsuralpha( multindex, multindex ) ;
       if iter > 1
          delta = delta + betasuralpha( multindexm1, multindexm1 );
       end
   
       eta = etaeta( multindex, multindex ); % what a stupid variable name
       
       if iter > 1
          H( multindex , [multindexm1,multindex] ) = [eta', delta];
          H( multindexm1 , multindex ) = eta;
       else
          H( multindex , multindex ) = delta;
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
   for i = 1:ncase*iter+1
      for j=1:ncase
         chiD(:,j) = inv(Theta1)*Y'*b(:,j); chiD(i:end,j) = 0;
      end
   end
   
   % Compute the number of pertinent Ritz modes
   indm2 = zeros(ncase,1); px = zeros(20*(ncase*niter-1)+1,ncase);
   t = 1:.05:ncase*niter; 
   for i=1:ncase
      [indm2(i),pol] = findPicard2 (log10(abs(chiD(:,i))), ceil(ncase*niter/5), 1, 3);%findPicard2(log10(abs(chiD)), 3, 1);
      n = size(pol,1);
      tt = zeros(n,20*(ncase*niter-1)+1);
      for j=1:n
         tt(j,:) = t.^(n-j);
      end
      px(:,i) = pol'*tt;
   end
   
   % Compute the solution
   chi = inv(Theta1)*Y'*b;
   if nbmodes > 0
      chi(nbmodes:end,:) = 0;
   else
      for i=1:ncase
         chi(indm2(i):end,i) = 0;
      end
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
   
%      try
%      figure;
%      hold on;
%      plot(log10(theta),'Color','blue')
%      plot(log10(abs(Y'*( b(:,1) ))),'Color','red')
%      plot(log10(abs(chiD(:,1))),'Color','black')
%      plot(t,px(:,1),'Color','cyan')
%      legend( 'Ritz Values','RHS values','solution coefficients', ...
%              'polynomial approximation' );
%      end

      try
      figure;
      hold on;
      for i=1:ncase
         plot(log10(abs(chiD(:,i))),'Color','black')
         plot(t,px(:,i),'Color','cyan')
      end
      end
   end

   %%%%%%%
   %% Final resolution
   f_D = fD; % Recover the Dirichlet Rhs where it is avaliable
   f_N = [ItereR; zeros(nbloq1d,ncase)]; % Use the reconstructed rhs elsewhere
   uin  = K1d\(f_N+f_D);

   usol = uin(1:ndofs,:);
end
