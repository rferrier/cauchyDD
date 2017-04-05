% 05/01/2017
% Problème de la plaque trouée

close all;
clear all;

% Parameters
E        = 210000; % MPa : Young modulus
nu       = 0.3;    % Poisson ratio
fscalar  = 250;    % N.mm-2 : Loading on the plate
mat      = [0, E, nu];
niter    = 10;
precond  = 1;      %
mu       = 0.;     % Regularization parameter
ratio    = 1e-300; % Maximal ratio (for eigenfilter)
br       = .0;      % noise
brt      = 0;      % "translation" noise
epsilon  = 1e-1;   % Convergence criterion for ritz value
ntrunc   = 0;      % In case the algo finishes at niter
bestiter = 0;     % chooses the best iteration with the L-curve

% methods : 1-> SPD, 2-> KMF-O, 3-> KMF-fixed point
methods = [1];

% Boundary conditions
% first index  : index of the boundary
% second index : 1=x, 2=y
% third        : value
% [0,1,value] marks a dirichlet regularization therm on x
dirichlet = [ 3,1,0 ; 3,2,0 ];
%dirichlet = [ 0,3,0 ];
neumann   = [3,2,fscalar];

% First, import the mesh
[ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/plate_hole.msh' );
nnodes = size(nodes,1);

noises = load('./noises/noisehole0.mat'); % Particular noise vector
noise  = noises.noise;
%noise  = randn(2*nnodes,1);

% Then, build the stiffness matrix :
[K,C,nbloq,node2c,c2node] = Krig2 (nodes,elements,mat,order,boundary,dirichlet);
Kinter = K( 1:2*nnodes, 1:2*nnodes );

% RHS 
udi = zeros( nnodes, 1 );
ind = 2:2:2*nnodes;
udi(ind,1) = 3e-3;
f = dirichletRhs2( udi, 3, c2node, boundary, nnodes );
%f = loading(nbloq,nodes,boundary,neumann);

% Center of the pion / Radius
Rc = .997;%.997
xc = 0; yc = (Rc-1);

% Build Cc
[ node2b5, b2node5 ] = mapBound( 5, boundary, nnodes );
Cc = zeros( 2*nnodes+nbloq, size(b2node5) );
b = zeros( size(b2node5), 1 ); % RHS
for i=1:size(b2node5)
   no = b2node5(i);
   x = nodes(no,1); y = nodes(no,2);
   Dis = sqrt( (x-xc)^2 + (y-yc)^2 );
   Cc(2*no-1,i) = -(x-xc)/Dis; Cc(2*no,i) = -(y-yc)/Dis;  % U.N <= b
   b(i) = Dis - Rc;
end

% Solve the problem :
uin = statusIneq( K, Cc, b, f, 10 );

% Extract displacement :
uref = uin(1:2*nnodes,1);
fref = Kinter*uref;
urefb = ( 1 + br*noise ) .* uref;
%frefb = ( 1 + br*noise ) .* fref;
ui = reshape(uref,2,[])';  ux = ui(:,1);  uy = ui(:,2);

% Compute stress :
sigma = stress(uref,E,nu,nodes,elements,order,1,ntoelem);

% Output :
plotGMSH({ux,'U_x';uy,'U_y';uref,'U_vect';sigma,'stress'}, elements, nodes, 'field');

% Plot displacement on the interface :
index = 2*b2node5;
thetax = 0:2*pi/size(index,1):2*pi*(1-1/size(index,1));
%hold on
%set(gca, 'fontsize', 20);
%plot(thetax,uref(index,1),'linewidth',3);
%plot(thetax,uref(index-1,1),'Color','red','linewidth',3);
%legend('uy','ux')
%xlabel('angle(rad)')
%ylabel('U(mm)')
%figure
%hold on
%set(gca, 'fontsize', 20);
%plot(thetax,fref(index,1),'linewidth',3);
%plot(thetax,fref(index-1,1),'Color','red','linewidth',3);
%legend('fy','fx')
%xlabel('angle(rad)')
%ylabel('F(mm)')

indexxy = [index-1;index];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition of the operators
if find(methods==1)
   % First problem
   dirichlet1 = [5,1,0;5,2,0;
                 4,1,0;4,2,0;
                 3,1,0;3,2,0;
                 2,1,0;2,2,0;
                 1,1,0;1,2,0];
   
   [K1,C1,nbloq1,node2c1,c2node1] = Krig (nodes,elements,E,nu,order,boundary,dirichlet1);
   
   % Second problem
   dirichlet2 = [3,1,0;3,2,0;
                 5,1,0;5,2,0];
   neumann2   = []; % no imposed force
   neumann0   = [];
   [K2,C2,nbloq2,node2c2,c2node2] = Krig (nodes,elements,E,nu,order,boundary,dirichlet2);
   
   error    = zeros(niter+1,1);
   residual = zeros(niter+1,1);
   regulari = zeros(niter+1,1);
   
   %% Dual problems
   % First problem
   dirichlet1d = [4,1,0;4,2,0;
                  3,1,0;3,2,0;
                  2,1,0;2,2,0;
                  1,1,0;1,2,0];
   [K1d,C1d,nbloq1d,node2c1d,c2node1d] = Krig (nodes,elements,E,nu,order,boundary,dirichlet1d);
   
   % Second problem
   dirichlet2d = [3,1,0;3,2,0];
   [K2d,C2d,nbloq2d,node2c2d,c2node2d] = Krig (nodes,elements,E,nu,order,boundary,dirichlet2d);
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% Conjugate Gradient for the problem : (D20-D10) x = D1-D2
   Itere  = zeros( 2*nnodes, 1 );
   d      = zeros( 2*nnodes, niter+1 );
   Ad     = zeros( 2*nnodes, niter+1 );
   Res    = zeros( 2*nnodes, niter+1 );
   Zed    = zeros( 2*nnodes, niter+1 );
   alpha  = zeros( niter+1, 1 );
   beta   = zeros( niter+1, 1 );
   alpha2 = zeros( niter+1, 1 );
   
   %% Perform A x0 :
   % Solve 1
   f1 = [Itere; zeros(nbloq1d,1)];
   uin1 = K1d\f1;
   u1i = uin1(1:2*nnodes,1);
   u1 = keepField( u1i, 5, boundary );
   % Solve 2
   f2 = [Itere; zeros(nbloq2d,1)];
   uin2 = K2d\f2;
   u2i = uin2(1:2*nnodes,1);
   u2 = keepField( u2i, 5, boundary );
   % Regularization term
   Nu = regul(Itere, nodes, boundary, 5);
   %
   Axz = mu*Nu+u2-u1;
   %%%%
   %% Compute Rhs :
   % Solve 1
   %f11 = dirichletRhs2( urefb, 1, c2node1d, boundary, nnodes );
   %f12 = dirichletRhs2( urefb, 2, c2node1d, boundary, nnodes );
   %f14 = dirichletRhs2( urefb, 4, c2node1d, boundary, nnodes );
   f16 = dirichletRhs2( urefb, 6, c2node1d, boundary, nnodes ); %(6) = (1)+(2)+(4)
   f13 = dirichletRhs2( uref, 3, c2node1d, boundary, nnodes );
   uin1 = K1d \ assembleDirichlet( [ f16, f13 ] ); %assembleDirichlet( [ f11, f12, f14 ] );
   u1i = uin1(1:2*nnodes,1);
   u1 = keepField( u1i, 5, boundary );
   % Solve 2
   f13 = dirichletRhs2( uref, 3, c2node2d, boundary, nnodes );
   uin2 = K2d\f13;
   u2i = uin2(1:2*nnodes,1);
   u2 = keepField( u2i, 5, boundary );
   %
   b = u1-u2;
   %%
   Res(:,1) = b - Axz;
   
   if precond == 1
       % Solve 1
       f1 = dirichletRhs2( Res(:,1), 5, c2node1, boundary, nnodes );
       uin1 = K1\f1;
       lagr1 = uin1(2*nnodes+1:end,1);
       lamb1 = lagr2forces2( lagr1, c2node1, 5, boundary, nnodes );
   %    % Solve 2
   %    f2 = dirichletRhs2( Res(:,1)/2, 5, c2node2, boundary, nnodes );
   %    uin2 = K2\f2;
   %    lagr2 = uin2(2*nnodes+1:end,1);
   %    lamb2 = lagr2forces( lagr2, C2, 5, boundary );
   %    %
   %    Zed(:,1) = -lamb2/2+lamb1/2;
       Zed(:,1) = lamb1;
%       Zed(:,1) = Res(:,1);
   else
       Zed(:,1) = Res(:,1);
   end
   
   d(:,1) = Zed(:,1);
   
   residual(1) = norm(Res( indexxy,1));
   error(1)    = norm(Itere(indexxy) - fref(indexxy)) / norm(fref(indexxy));
   regulari(1) = sqrt( Itere'*regul(Itere, nodes, boundary, 5) );
   
   ritzval  = 0; % Last ritz value that converged
   oldtheta = 0;
   eta      = 0;
   getmeout = 0; % utility
   %V = zeros(2*nnodes, iter);
   %H = zeros(iter);
   %%
   for iter = 1:niter
       %% Optimal step
       
       % Solve 1
       f1 = [d(:,iter); zeros(nbloq1d,1)];
       uin1 = K1d\f1;
       u1i = uin1(1:2*nnodes,1);
       u1 = keepField( u1i, 5, boundary );
       % Solve 2
       f2 = [d(:,iter); zeros(nbloq2d,1)];
       uin2 = K2d\f2;
       u2i = uin2(1:2*nnodes,1);
       u2 = keepField( u2i, 5, boundary );
       % Regularization term
       Nu = regul(d(:,iter), nodes, boundary, 5);
       %
       Ad(:,iter) = mu*Nu+u2-u1;
       
       den = (d(indexxy,iter)'*Ad(indexxy,iter));
       d(:,iter) = d(:,iter)/sqrt(den); Ad(:,iter) = Ad(:,iter)/sqrt(den);
       num = Res(indexxy,iter)'*d(indexxy,iter);
       
       Itere         = Itere + d(:,iter)*num;%/den;
       Res(:,iter+1) = Res(:,iter) - Ad(:,iter)*num;%/den;
       
       residual(iter+1) = norm(Res(indexxy,iter+1));
       error(iter+1)    = norm(Itere(indexxy) - fref(indexxy)) / norm(fref(indexxy));
       regulari(iter+1) = sqrt( Itere'*regul(Itere, nodes, boundary, 5) );
       
       if precond == 1
           % Solve 1
           f1 = dirichletRhs2( Res(:,iter+1), 5, c2node1, boundary, nnodes );
           uin1 = K1\f1;
           lagr1 = uin1(2*nnodes+1:end,1);
           lamb1 = lagr2forces2( lagr1, c2node1, 5, boundary, nnodes );
   %        % Solve 2
   %        f2 = dirichletRhs2( Res(:,iter+1)/2, 5, c2node2, boundary, nnodes );
   %        uin2 = K2\f2;
   %        lagr2 = uin2(2*nnodes+1:end,1);
   %        lamb2 = lagr2forces( lagr2, C2, 5, boundary );
   %        %
   %        Zed(:,iter+1) = -lamb2/2+lamb1/2;
           Zed(:,iter+1) = lamb1;
           %Zed(:,iter+1) = Res(:,iter+1);
       else
           Zed(:,iter+1) = Res(:,iter+1);
       end
       
       % Needed values for the Ritz stuff
       alpha(iter) = num/sqrt(den);
       beta(iter)  = - Zed(indexxy,iter+1)'*Ad(indexxy,iter)/sqrt(den);
   %    alpha(iter) = Res(indexxy,iter)'*Res(indexxy,iter) / den;
   %    beta(iter)  = Zed(indexxy,iter+1)'*Res(indexxy,iter+1) /... 
   %                                (Zed(indexxy,iter)'*Res(indexxy,iter));
       
       % First Reorthogonalize the residual (as we use it next), in sense of M
       for jter=1:iter
           betac = Zed(indexxy,iter+1)'*Res(indexxy,jter) / (Zed(indexxy,jter)'*Res(indexxy,jter));
           Zed(:,iter+1) = Zed(:,iter+1) - Zed(:,jter) * betac;
           Res(:,iter+1) = Res(:,iter+1) - Res(:,jter) * betac;
       end
       
       %% Orthogonalization
       d(:,iter+1) = Zed(:,iter+1);
       
       for jter=iter:iter
           betaij = ( Zed(indexxy,iter+1)'*Ad(indexxy,jter) );%/...
               %( d(indexxy,jter)'*Ad(indexxy,jter) );
           d(:,iter+1) = d(:,iter+1) - d(:,jter) * betaij;
       end
       
       %% Ritz algo : find the Ritz elements
       % Build the matrices
       V(:,iter) = zeros(2*nnodes,1);
       V(indexxy,iter) = (-1)^(iter-1)*Zed(indexxy,iter)/(sqrt(Res(indexxy,iter)'*Zed(indexxy,iter)));
   
       etap   = eta;
       delta  = 1/alpha(iter);
       if iter > 1
          delta  = delta + beta(iter-1)/alpha(iter-1);
       end
       eta    = sqrt(beta(iter))/alpha(iter);
       
       if iter > 1
          H(iter,[iter-1,iter]) = [etap, delta];
          H(iter-1,iter)        = etap;
       else
          H(iter,iter) = delta;
       end
       
       % Compute eigenelems of the Hessenberg :
       [Q,Theta1] = eig(H);
       theta = diag(Theta1);
       % Sort it
       [theta,Ind] = sort(theta,'descend');
       Q = Q(:,Ind);
       Theta1 = Theta1(Ind,Ind);
       Y = V*Q;
       
       % See if the current one converged
       if abs(theta(ritzval+1) - oldtheta) < epsilon*oldtheta
          % increment oldtheta
          ritzval = ritzval+1;
          if size(theta,1) > ritzval
             oldtheta = theta(ritzval+1);
          else
             oldtheta = 0;
          end
          
          % Check small value / hold value
          if ritzval > 1
             if theta(ritzval) < ratio*theta(1)
                ntrunc = ritzval
                getmeout = 1;
                break;
             end
          end
       else
          oldtheta = theta(ritzval+1);
       end
   %    if getmeout == 1  % In case I use a while over there
   %       break;
   %    end
   end
   
   % Compute the solution
   chi = inv(Theta1)*Y'*b;
   if ntrunc > 0
      chi(ntrunc:end) = 0;
   end
   ItereR = Y*chi;
   
   regS = zeros(niter,1);
   resS = zeros(niter,1);
   %% Build the L-curve regul, ntrunc
   for i = 1:iter+1
      chiS   = inv(Theta1)*Y'*b; chiS(i:end) = 0;
      ItereS = Y*chiS;
      
      % Solve 1
      f1 = [ItereS; zeros(nbloq1d,1)];
      uin1 = K1d\f1;
      u1i = uin1(1:2*nnodes,1);
      u1 = keepField( u1i, 5, boundary );
      % Solve 2
      f2 = [ItereS; zeros(nbloq2d,1)];
      uin2 = K2d\f2;
      u2i = uin2(1:2*nnodes,1);
      u2 = keepField( u2i, 5, boundary );
      %
      AI = u2-u1;
      
      ResS = AI-b;
      resS(i) = norm(ResS);   
      regS(i) = sqrt( ItereS'*regul(ItereS, nodes, boundary, 5) );
      errS(i) = norm(ItereS(indexxy)-fref(indexxy)) / norm(fref(indexxy));
   end
   
   traD = 1 - sum(theta(1:ntrunc-1))/sum(theta);
   regD = zeros(niter,1); resD = zeros(niter,1); bt = Y'*b;
   for i = 1:iter+1
      chiD   = inv(Theta1)*Y'*b; chiD(i:end) = 0;
      ItereD = Y*chiD;
      resD(i) = sqrt( sum( bt(i:end).^2) );  
      regD(i) = sqrt( ItereD'*regul(ItereD, nodes, boundary, 5) );
   end
   
   [indm2,pol] = findPicard2(log10(abs(chiD)), 5, 1);
   n = size(pol,1);
   t = 1:.05:niter; tt = zeros(n,20*(niter-1)+1);
   for j=1:n
   tt(j,:) = t.^(n-j);
   end
   px = pol'*tt;
   
   figure;
   hold on;
   plot(log10(theta),'Color','blue')
   plot(log10(abs(Y'*b)),'Color','red')
   plot(log10(abs(chiD)),'Color','black')
   plot(t,px,'Color','cyan')
   legend('Ritz Values','RHS values','solution coefficients', ...
          'polynomial approxiamtion')
   %
%   figure;
%   hold on;
%   plot(Itere(2*b2node5-1),'Color','red')
%   plot(ItereR(2*b2node5-1),'Color','blue')
%   plot(fref(2*b2node5-1),'Color','green')
%   legend('brutal solution','filtred solution', 'reference')
   %
   figure;
   hold on;
   set(gca, 'fontsize', 20);
   %plot(Itere(2*b2node5),'Color','red')
   plot(thetax,ItereR(2*b2node5),'Color','blue','linewidth',3)
   plot(thetax,fref(2*b2node5),'Color','green','linewidth',3)
   %legend('brutal solution','filtred solution', 'reference')
   legend('SPD solution', 'reference')
   xlabel('angle'); ylabel('F');
   
   figure;
   hold on
   plot(error(2:end),'Color','blue')
   plot(residual(2:end),'Color','red')
   plot(errS(2:end),'Color','green')
   legend('error','residual','Ritz error')
   %L-curve :
%   figure;
%   set(gca, 'fontsize', 20);
%   hold on;
%   loglog(residual(2:iter+1),regulari(2:iter+1),'Color','red','-*','linewidth',3);
%   loglog(resS(2:iter+1),regS(2:iter+1),'-+','linewidth',3);
%   legend('L-curve','RL-curve')
%   xlabel('residual')
%   ylabel('H1 norm')
   %findCorner (residual(2:iter+1), regulari(2:iter+1), 3)
   %findCorner (resS(2:iter+1), regS(2:iter+1), 3)
   
   %%hold on;
   %%loglog(resS(2:iter+1),regS(2:iter+1),'Color','red','-*');
   %loglog(resD(2:iter),regD(2:iter),'-+');
   %%legend('RL-curve (natural basis)','RL-curve (diagonal basis)')
   %figure
   %findCorner (resD(2:iter), regD(2:iter), 3)
   %%%%%
   %% Final problem : compute u
   dirichlet = [4,1,0;4,2,0;
                3,1,0;3,2,0;
                2,1,0;2,2,0;
                1,1,0;1,2,0];
   
   [K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
   f5    = [keepField( Itere, 5, boundary ); zeros(nbloq,1)];
   %fdir1 = dirichletRhs(urefb, 1, C, boundary);
   %fdir2 = dirichletRhs(urefb, 2, C, boundary);
   fdir3 = dirichletRhs(uref, 3, C, boundary);
   %fdir4 = dirichletRhs(urefb, 4, C, boundary);
   fdir6 = dirichletRhs(urefb, 6, C, boundary);
   usoli = K \ (assembleDirichlet( [ fdir6, fdir3 ] ) + f5);
   
   usol = usoli(1:2*nnodes,1);
   fsol = Kinter*usol;
   
   % With the reduced solution
   dirichlet = [4,1,0;4,2,0;
                3,1,0;3,2,0;
                2,1,0;2,2,0;
                1,1,0;1,2,0];
   
   [K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
   f5    = [keepField( ItereR, 5, boundary ); zeros(nbloq,1)];
   %fdir1 = dirichletRhs(urefb, 1, C, boundary);
   %fdir2 = dirichletRhs(urefb, 2, C, boundary);
   fdir3 = dirichletRhs(uref, 3, C, boundary);
   %fdir4 = dirichletRhs(urefb, 4, C, boundary);
   fdir6 = dirichletRhs(urefb, 6, C, boundary);
   usoli = K \ (assembleDirichlet( [ fdir6, fdir3 ] ) + f5);
   
   usolR = usoli(1:2*nnodes,1);
   fsolR = Kinter*usolR;
   
%   figure;
%   hold on;
%   plot(usol(2*b2node5-1),'Color','red')
%   plot(usolR(2*b2node5-1),'Color','blue')
%   plot(uref(2*b2node5-1),'Color','green')
%   legend('brutal solution','filtred solution','reference')
   
%   figure;
%   hold on;
%   set(gca, 'fontsize', 20);
%   %plot(thetax,usol(2*b2node5),'Color','red')
%   plot(thetax,usolR(2*b2node5),'Color','blue','linewidth',3)
%   plot(thetax,uref(2*b2node5),'Color','green','linewidth',3)
%   legend('SPD solution','reference')
%   %legend('brutal solution','filtred solution','reference')
%   xlabel('angle'); ylabel('U');
   
   total_error = norm(usol-uref)/norm(uref);
   % Compute stress :
   sigma = stress(usol,E,nu,nodes,elements,order,1,ntoelem);
   plotGMSH({usol,'U_vect';sigma,'stress'}, elements, nodes, 'solution');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KMF Prequisites
if size(find(methods==2),1) == 1 || size(find(methods==3),1) == 1
   %% Definition of the operators
   % DN problem
   dirichlet1 = [5,1,0;5,2,0;
                 3,1,0;3,2,0];
   neumann1   = [];
   neumann0   = [];
   [K1,C1,nbloq1,node2c1,c2node1] = ...
                   Krig (nodes,elements,E,nu,order,boundary,dirichlet1);
   
   % ND problem
   dirichlet2 = [6,1,0;6,2,0;
                 3,1,0;3,2,0];
   neumann2   = [];% [3,1,lagr1; 3,2,lagr1] is managed by lagr2forces
   [K2,C2,nbloq2,node2c2,c2node2] = ...
                   Krig (nodes,elements,E,nu,order,boundary,dirichlet2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KMF + Orthodir
if find(methods==2)
   
   error   = zeros(niter+1,1);
   residual = zeros(niter+1,1);
   regulari = zeros(niter+1,1);

   %% ORTHODIR for the problem : (1 - D0 o S0 + mu*N) x = DoS(0), with N, a 
   % regulalization term
   Itere = zeros( 2*nnodes, niter+1 ); 
   p     = zeros( 2*nnodes, niter+1 );
   q     = zeros( 2*nnodes, niter+1 );
   Delta = zeros( niter+1, 1 );
   Res   = zeros( 2*nnodes, niter+1 );

   Itere(2:2:2*nnodes,1) = 3e-3;
   
   %% Perform A x0 :
   % Solve DN
   fdir = dirichletRhs(Itere(:,1), 5, C1, boundary);
   f1 = fdir;
   uin1 = K1\f1;
   lagr1 = uin1(2*nnodes+1:end,1);
   % Solve ND
   fr = lagr2forces( lagr1, C1, 5, boundary );
   fr = [ fr; zeros(size(C2,2),1) ];
   f2 = fr;% + fdir1 + fdir2;
   uin2 = K2\f2;
   % Regularization term
   Nu = regul(Itere(:,1), nodes, boundary, 5);
   
   atimesItere = mu*Nu + Itere(:,1) - uin2(1:2*nnodes,1);
   %%%%

   %% Write Rhs :
   % Solve DN
   fdir = dirichletRhs(uref, 3, C1, boundary); %uref
   f1 = fdir;
   uin1 = K1\f1;
   u1 = uin1(1:2*nnodes,1);
   lagr1 = uin1(2*nnodes+1:end,1);
   % Solve ND
   fr = lagr2forces( lagr1, C1, 5, boundary );
   fr = [ fr; zeros(size(C2,2),1) ]; % Give to fr the right size
   fdir1 = dirichletRhs(urefb, 6, C2, boundary);
   fdir2 = dirichletRhs(uref, 3, C2, boundary);
   f2 = fr + assembleDirichlet( [fdir1,fdir2] );
   uin2 = K2\f2;
   b = uin2(1:2*nnodes,1);
   %%%% 
   %%
   Res(:,1) = b - atimesItere;
   p(:,1) = Res(:,1);

   residual(1) = norm(Res(indexxy,1));
   error(1)    = norm(Itere(indexxy,1)-uref(indexxy)) / norm(uref(indexxy));
   regulari(1) = sqrt( Itere(:,1)'*regul(Itere(:,1), nodes, boundary, 5) );

   %% Perform Q1 = A P1 :
   % Solve DN
   fdir = dirichletRhs(p(:,1), 5, C1, boundary);
   f1 = fdir;
   uin1 = K1\f1;
   lagr1 = uin1(2*nnodes+1:end,1);
   % Solve ND
   fr = lagr2forces( lagr1, C1, 5, boundary );
   fr = [ fr; zeros(size(C2,2),1) ];
   f2 = fr;
   uin2 = K2\f2;
   % Regularization term
   Nu = regul(p(:,1), nodes, boundary, 5);
   
   q(:,1) = mu*Nu + p(:,1) - uin2(1:2*nnodes,1);
   %%%%
   %%
   for iter = 1:niter

      Delta(iter,1)   = q(:,iter)'*q(:,iter);
      gammai          = q(:,iter)'*Res(:,iter);
      alphai          = gammai/Delta(iter,1);

      Itere(:,iter+1) = Itere(:,iter) + p(:,iter)*alphai;
      Res(:,iter+1)   = Res(:,iter) - q(:,iter)*alphai;

      residual(iter+1) = norm(Res(indexxy,iter+1));
      error(iter+1)    = norm(Itere(indexxy,iter+1)-uref(indexxy)) / norm(uref(indexxy));
      regulari(iter+1) = sqrt( Itere(:,iter+1)'*regul(Itere(:,iter+1), nodes, boundary, 5) );

      %% Perform Ari = A*Res
      % Solve DN
      fdir = dirichletRhs(Res(:,iter+1), 5, C1, boundary);
      f1 = fdir;
      uin1 = K1\f1;
      lagr1 = uin1(2*nnodes+1:end,1);
      % Solve ND
      fr = lagr2forces( lagr1, C1, 5, boundary );
      fr = [ fr; zeros(size(C2,2),1) ];
      f2 = fr;% + fdir1 + fdir2;
      uin2 = K2\f2;
      % Regularization term
      Nu = regul(Res(:,iter+1), nodes, boundary, 5);
      
      Ari = mu*Nu + Res(:,iter+1) - uin2(1:2*nnodes,1);
      %%%%

      %% Orthogonalization
      p(:,iter+1) = Res(:,iter+1);
      q(:,iter+1) = Ari;

      for jter=1:iter
         phiij  = q(:,jter)'*Ari;
         betaij = phiij/Delta(jter,1);
         p(:,iter+1) = p(:,iter+1) - p(:,jter) * betaij;
         q(:,iter+1) = q(:,iter+1) - q(:,jter) * betaij;
         %myps( q(:,iter+1), q(:,jter), Kinter, boundary, M, nodes )
      end

   end

   figure
   hold on
   set(gca, 'fontsize', 15);
%   set(gca,'ylim',[-11 1])
   plot(log10(error(2:end)),'Color','blue')
   plot(log10(residual(2:end)/residual(1)),'Color','red')
   legend('error (log)','residual (log)')

    
   % L-curve
   figure
   loglog(residual(2:end),regulari(2:end));

   %% Automatic stop (regularization tools, broken on Octave) :
   %[stopHere,~,~] = l_corner(residual,regulari,1:1:niter+1);
   if bestiter == 1
      stopHere = findCorner(residual(2:end),regulari(2:end)) + 1;
   else
      stopHere = niter+1;
   end
   
   % Plot chosen Itere
   figure
   hold on;
   set(gca, 'fontsize', 20);
   plot(uref(index));
   plot(Itere(index,stopHere),'Color','red');
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KMF + fixed point
if find(methods==3)
    % init :
    u1    = uref-uref;
    u2    = u1;
    u2(2:2:2*nnodes,1) = 3e-3;

    error1   = zeros(niter,1);
    error2   = zeros(niter,1);
    residual = zeros(niter,1);
    regulari = zeros(niter,1);

    for iter = 1:niter
        % Solve DN
        f1 = dirichletRhs2(u2, 5, c2node1, boundary, nnodes);
        fdir = dirichletRhs2(uref, 3, c2node1, boundary, nnodes);
        uin1 = K1\(f1+fdir);
        u1 = uin1(1:2*nnodes,1);
        lagr1 = uin1(2*nnodes+1:end,1);
        error1(iter) = norm(u1(indexxy)-uref(indexxy)) / norm(uref(indexxy));

        % Solve ND
        fr = lagr2forces2( lagr1, c2node1, 5, boundary, nnodes );
        fr = [ fr; zeros(size(C2,2),1) ]; % Give to fr the right size
        fdir1 = dirichletRhs2(urefb, 6, c2node2, boundary, nnodes);
        fdir2 = dirichletRhs2(uref, 3, c2node2, boundary, nnodes);
        f2 = fr + assembleDirichlet( [fdir1,fdir2] );
        uin2 = K2\f2;
        u2 = uin2(1:2*nnodes,1);
        
        error2(iter) = norm(u2(indexxy)-uref(indexxy)) / norm(uref(indexxy));
        if iter>1  
           residual(iter) = norm(u1(indexxy)-u2(indexxy)) /...
                            sqrt ( norm(u1(indexxy))*norm(u2(indexxy)) );
        end
        regulari(iter) = sqrt(u2'*( regul(u2, nodes, boundary, 5) ));

    end
    
    figure
    hold on
    set(gca, 'fontsize', 20);
%    set(gca,'ylim',[-11 1])
    plot(log10(error1(2:end)),'Color','blue','linewidth',3)
    plot(log10(error2(2:end)),'Color','black','linewidth',3)
    plot(log10(residual(2:end)),'Color','red','linewidth',3)
    legend('error1 (log)','error2 (log)','residual (log)')

%    % L-curve
%    figure
%    loglog(residual(2:end),regulari(2:end));
   
%    % Plot chosen Itere
%    figure
%    hold on;
%    set(gca, 'fontsize', 20);
%    plot(uref(index),'linewidth',3);
%    plot(u2(index),'Color','red','linewidth',3);
    
    % Last problem
    dirichlet = [5,1,0;5,2,0;
                 4,1,0;4,2,0;
                 3,1,0;3,2,0;
                 2,1,0;2,2,0;
                 1,1,0;1,2,0];

    [K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
    f5    = dirichletRhs(u2, 5, C, boundary );
    fdir3 = dirichletRhs(uref, 3, C, boundary);
    fdir6 = dirichletRhs(urefb, 6, C, boundary);
    usoli = K \ (assembleDirichlet( [ fdir6, fdir3 ] ) + f5);

    usolR = usoli(1:2*nnodes,1);
    fsolR = Kinter*usolR;
    
    figure;
    hold on;
    plot(usolR(2*b2node5),'Color','blue')
    plot(uref(2*b2node5),'Color','green')
    legend('solution','reference')
    
    figure;
    hold on;
    plot(fsolR(2*b2node5),'Color','blue')
    plot(fref(2*b2node5),'Color','green')
    legend('solution','reference')
    
end