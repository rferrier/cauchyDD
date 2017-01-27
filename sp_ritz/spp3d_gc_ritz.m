% 26/01/2017
% Algo Steklov-Poincaré primal avec Gradient Conjugué, cas 3D

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 70000;  % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 1;      % N.mm-1 : Loading on the plate
niter   = 50;
precond = 0;      % 1 : Use a dual precond, 2 : use H1/2 precond, 3 : use gradient precond
mu      = 0.;     % Regularization parameter
ratio   = 5e-200; % Maximal ratio (for eigenfilter)
br      = 0.;    % noise
brt     = 0;      % "translation" noise
epsilon = 1e-200; % Convergence criterion for ritz value
ntrunc  = 20;      % In case the algo finishes at niter
inhomog = 0;      % inhomogeneous medium

% methods : 1-> SPP, 2-> SPD
methods = [2];

%if inhomog == 2  % load previously stored matrix
%   mat = [2, E, nu, .1, 1];
%   Kinter = load('./noises/stocrig1.mat'); Kinter = Kinter.Kinter;
%elseif inhomog == 1  % /!\ IN THE INHOMOGENEOUS CASE, ALL THE SIGMAS ARE WRONG
%   mat = [2, E, nu, .1, 1];
%else
%   mat = [0, E, nu];
%end
mat = [0, E, nu];

% Boundary conditions
% first index  : index of the boundary
% second index : 1=x, 2=y
% third        : value
% [0,1,value] marks a dirichlet regularization therm on x
dirichlet = [3,1,0; 3,2,0 ; 3,3,0
             4,1,0; 4,2,0 ; 4,3,0
             5,1,0; 5,2,0 ; 5,3,0
             6,1,0; 6,2,0 ; 6,3,0 ];
neumann   = [7,3,-fscalar];

% Import the mesh
[ nodes,elements,ntoelem,boundary,order ] = readmesh3D( 'meshes/plate3d_charge2.msh' );
nnodes = size(nodes,1);

%noises = load('./noises/noise3d.mat'); % Particular noise vector
%noise  = noises.noise;
noise = randn(3*nnodes,1);

% suppress the aretes:
boundaryp1 = suppressBound( boundary, 'extreme', 2, nodes,1e-5 );
boundaryp2 = boundaryp1;
%
%[~, b2node2] = mapBound(2, boundaryp1, nnodes);
%: = [2*b2node2-1; 2*b2node2];

% Then, build the stiffness matrix :
[K,C,nbloq] = Krig3D (nodes,elements,mat,order,boundary,dirichlet);
%if inhomog == 2
%   K(1:3*nnodes, 1:3*nnodes) = Kinter;
%else
%   Kinter = K(1:3*nnodes, 1:3*nnodes);
%end
Kinter = K(1:3*nnodes, 1:3*nnodes);
%M      = mass_mat(nodes, elements);

% The right hand side :
f = loading3D(nbloq,nodes,boundary,neumann);

% Solve the problem :
uin = K\f;

% Extract displacement and Lagrange multiplicators :
uref = uin(1:3*nnodes,1);
fref = Kinter*uref;
lagr = uin(3*nnodes+1:end,1);

urefb = ( 1 + brt + br*noise ) .* uref;

sigma = stress3D(uref,mat,nodes,elements,order,1,ntoelem);
plotGMSH3D({uref,'Vect_U';sigma,'stress'}, elements, nodes, 'output/reference');

% Extract the boundary 2
bounindex = find( boundary(:,1)==2 );
%indexa = [3*bounindex-2;3*bounindex-1;3*bounindex];
patch2    = boundary(bounindex,2:end);

% find the index of the upper nodes
Ckk = Cbound ( nodes, [2,1,0;2,2,0;2,3,0], boundaryp1 );
indexC= sum(Ckk'); indexa = find(indexC);

% H1 norm
Kreg = H1 (nodes, patch2, order); Kreg = Kreg(indexa,indexa);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition of the operators

% First problem
dirichlet1 = [6,1,0;6,2,0;6,3,0;
              5,1,0;5,2,0;5,3,0;
              4,1,0;4,2,0;4,3,0;
              3,1,0;3,2,0;3,3,0;
              2,1,0;2,2,0;2,3,0;
              1,1,0;1,2,0;1,3,0];

[K1,C1,nbloq1,node2c1,c2node1] = Krig3D(nodes,elements,mat,order,boundaryp1,dirichlet1);
if inhomog >= 1  % Because of the random stuff
   K1(1:3*nnodes, 1:3*nnodes) = Kinter;
end
C12 = Cbound ( nodes, [2,1,0;2,2,0;2,3,0], boundaryp1 ); % C12'u = u on 2 (without the aretes)
C11 = Cbound ( nodes, [1,1,0;1,2,0;1,3,0], boundary );

% Second problem
dirichlet2 = [6,1,0;6,2,0;6,3,0;
              5,1,0;5,2,0;5,3,0;
              4,1,0;4,2,0;4,3,0;
              3,1,0;3,2,0;3,3,0;
              2,1,0;2,2,0;2,3,0];
neumann2   = [7,3,-fscalar];
neumann0   = [];
[K2,C2,nbloq2,node2c2,c2node2] = Krig3D(nodes,elements,mat,order,boundaryp2,dirichlet2);
if inhomog >= 1
   K2(1:3*nnodes, 1:3*nnodes) = Kinter;
end
C22 = Cbound ( nodes, [2,1,0;2,2,0;2,3,0], boundaryp1 ); % OK, I admit it's the same as C12
C21 = Cbound ( nodes, [1,1,0;1,2,0;1,3,0], boundary );

%% Dual problems
% First problem
dirichlet1d = [6,1,0;6,2,0;6,3,0;
               5,1,0;5,2,0;5,3,0;
               4,1,0;4,2,0;4,3,0;
               3,1,0;3,2,0;3,3,0;
               1,1,0;1,2,0;1,3,0];
[K1d,C1d,nbloq1d] = Krig3D(nodes,elements,mat,order,boundary,dirichlet1d);
if inhomog >= 1
   K1d(1:3*nnodes, 1:3*nnodes) = Kinter;
end

% Second problem
dirichlet2d = [6,1,0;6,2,0;6,3,0;
               5,1,0;5,2,0;5,3,0;
               4,1,0;4,2,0;4,3,0;
               3,1,0;3,2,0;3,3,0];
[K2d,C2d,nbloq2d] = Krig3D(nodes,elements,mat,order,boundary,dirichlet2d);
if inhomog >= 1
   K2d(1:3*nnodes, 1:3*nnodes) = Kinter;
end

%% Anti-cancellation trick
K1r = K1; K2r = K2; K1dr = K1d; K2dr = K2d;
%K1(indexa,indexa) = 0;
%K2(indexa,indexa) = 0;
%K1d(indexa,indexa) = 0;
%K2d(indexa,indexa) = 0;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Conjugate Gradient for the problem : (S10-S20) x = S2-S1
if find(methods==1)

   % First, plot the reference
   ret = patch('Faces',patch2,'Vertices',nodes, ... 
               'FaceVertexCData',uref(3:3:end),'FaceColor','interp');
   colorbar;

   error    = zeros(niter+1,1);
   residual = zeros(niter+1,1);
   regulari = zeros(niter+1,1);
   
   Itere = zeros( 3*nnodes, 1 );
   d     = zeros( 3*nnodes, niter+1 );
   Ad    = zeros( 3*nnodes, niter+1 );
   Res   = zeros( 3*nnodes, niter+1 );
   Zed   = zeros( 3*nnodes, niter+1 );
   alpha = zeros( niter+1, 1 );
   beta  = zeros( niter+1, 1 );
   alpha2 = zeros( niter+1, 1 );
   %H12   = H1demi(size(Res,1), nodes, boundary, 2 );
   %Mgr   = Mgrad(size(Res,1), nodes, boundary, 2 );
   
   %% Perform A x0 :
   % Solve 1
   f1 = [ zeros(3*nnodes,1) ; C1'*C12*C12'*Itere ]; % u=0 outside of 2
   uin1 = K1\f1;
   lagr1 = uin1(3*nnodes+1:end,1);
   lamb1 = -C12*C12'*C1*lagr1;    % restriction of lagr to 2
   % Solve 2
   f2 = [ zeros(3*nnodes,1) ; C2'*C22*C22'*Itere ];
   uin2 = K2\f2;
   lagr2 = uin2(3*nnodes+1:end,1);
   lamb2 = -C22*C22'*C2*lagr2;
   % Regularization term
   %Nu = regul(Itere, nodes, boundary, 2);
   %
   Axz = lamb1-lamb2;%mu*Nu+lamb1-lamb2;
   %%%%
   %% Compute Rhs :
   % Solve 1
   f1 = [ zeros(3*nnodes,1) ; C1'*C11*C11'*urefb ];
   uin1 = K1\f1;
   lagr1 = uin1(3*nnodes+1:end,1);
   lamb1 = -C12*C12'*C1*lagr1;
   % Solve 2 (zero 'cause no loading except on 7)
   f2 = loading3D(nbloq2,nodes,boundary,neumann0);
   uin2 = K2\f2;
   lagr2 = uin2(3*nnodes+1:end,1);
   lamb2 = -C22*C22'*C2*lagr2;
   %
   b = lamb2-lamb1;
   %%
   Res(:,1) = b - Axz;
   
   if precond == 1
       % Solve 1
       f1 = [Res(:,1); zeros(nbloq1d,1)];
       uin1 = K1dr\f1;
       u1i = uin1(1:3*nnodes,1);
       u1 = C12*C12'*u1i;
       Zed(:,1) = u1;
   elseif precond == 2 % Not working
       Zed(index,1) = 1/E*H12(index,index)\Res(index,1);
   elseif precond == 3
       Zed(index,1) = 1/E*Mgr(index,index)\Res(index,1);
   else
       Zed(:,1) = Res(:,1);
   end
   
   d(:,1) = Zed(:,1);
   
   residual(1) = norm(Res( :,1));
   error(1)    = norm(Itere(indexa) - uref(indexa)) / norm(uref(indexa));
   regulari(1) = sqrt(Itere(indexa)'*Kreg*Itere(indexa));%norm(Itere);
   
   ritzval  = 0; % Last ritz value that converged
   oldtheta = 0;
   eta      = 0;
   getmeout = 0; % utility
   %V = zeros(3*nnodes, iter);
   %H = zeros(iter);
   %%
   for iter = 1:niter
       %% Optimal step
       
       % Solve 1
       f1 = [ zeros(3*nnodes,1) ; C1'*C12*C12'*d(:,iter) ]; % u=0 outside of 2
       uin1 = K1\f1;
       lagr1 = uin1(3*nnodes+1:end,1);
       lamb1 = -C12*C12'*C1*lagr1;    % restriction of lagr1 to 2
       % Solve 2
       f2 = [ zeros(3*nnodes,1) ; C2'*C22*C22'*d(:,iter) ];
       uin2 = K2\f2;
       lagr2 = uin2(3*nnodes+1:end,1);
       lamb2 = -C22*C22'*C2*lagr2;
       %
        Ad(:,iter) = lamb1-lamb2;
       
       den = (d(:,iter)'*Ad(:,iter));
       d(:,iter) = d(:,iter)/sqrt(den); Ad(:,iter) = Ad(:,iter)/sqrt(den);
       num = Res(:,iter)'*d(:,iter);
       
       Itere         = Itere + d(:,iter)*num;%/den;
       Res(:,iter+1) = Res(:,iter) - Ad(:,iter)*num;%/den;
       
       residual(iter+1) = norm(Res(:,iter+1));
       error(iter+1)    = norm(Itere(indexa) - uref(indexa)) / norm(uref(indexa));
       regulari(iter+1) = sqrt(Itere(indexa)'*Kreg*Itere(indexa));%norm(Itere);
       
       if precond == 1
           % Solve 1
           f1 = [Res(:,iter+1); zeros(nbloq1d,1)];
           uin1 = K1dr\f1;
           u1i = uin1(1:3*nnodes,1);
           u1 = C12*C12'*u1i;
           Zed(:,iter+1) = u1;
       elseif precond == 2
           Zed(index,iter+1) = 1/E*H12(index,index)\Res(index,iter+1);
       elseif precond == 3
           Zed(index,iter+1) = 1/E*Mgr(index,index)\Res(index,iter+1);
       else
           Zed(:,iter+1) = Res(:,iter+1);
       end
       
       % Needed values for the Ritz stuff
       alpha(iter) = num/sqrt(den);
       alpha2(iter) = num;
       beta(iter)  = - Zed(:,iter+1)'*Ad(:,iter)/sqrt(den);
       
       % First Reorthogonalize the residual (as we use it next), in sense of M
       for jter=1:iter-1
           betac = Zed(:,iter+1)'*Res(:,jter) / (Zed(:,jter)'*Res(:,jter));
           Zed(:,iter+1) = Zed(:,iter+1) - Zed(:,jter) * betac;
       end
       
       %% Orthogonalization
       d(:,iter+1) = Zed(:,iter+1);
       
       for jter=iter:iter % No need to reorthogonalize (see above)
           betaij = ( Zed(:,iter+1)'*Ad(:,jter) );%/...
               %( d(:,jter)'*Ad(:,jter) );
           d(:,iter+1) = d(:,iter+1) - d(:,jter) * betaij;
   
       end
       
       %% Ritz algo : find the Ritz elements
       % Build the matrices
       V(:,iter) = zeros(3*nnodes,1);
       V(:,iter) = (-1)^(iter-1)*Zed(:,iter)/(sqrt(Res(:,iter)'*Zed(:,iter)));
       %norm(Zed(:,iter))
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
   
   %%for i=1:niter
   %   plot(V(2*b2node2-1,20:25))
   %   figure;
   %%end
   
   regS = zeros(niter,1);
   resS = zeros(niter,1);
   %% Build the L-curve regul, ntrunc
   for i = 1:iter+1
      chiS   = inv(Theta1)*Y'*b; chiS(i:end) = 0;
      ItereS = Y*chiS;
      
      % Solve 1
      f1 = [ zeros(3*nnodes,1) ; C1'*C12*C12'*ItereS ]; % u=0 outside of 2
      uin1 = K1\f1;
      lagr1 = uin1(3*nnodes+1:end,1);
      lamb1 = -C12*C12'*C1*lagr1;    % restriction of lagr to 2
      % Solve 2
      f2 = [ zeros(3*nnodes,1) ; C2'*C22*C22'*ItereS ];
      uin2 = K2\f2;
      lagr2 = uin2(3*nnodes+1:end,1);
      lamb2 = -C22*C22'*C2*lagr2;
      %
      AI = lamb1-lamb2;
      
      ResS = AI-b;
      resS(i) = norm(ResS);   
      regS(i) = sqrt(ItereS(indexa)'*Kreg*ItereS(indexa));%norm(ItereS);
      errS(i) = norm(ItereS(indexa)-uref(indexa))/norm(uref(indexa));
   end
   
   % Residual in the diagonal base :
   %resD = zeros(iter,1);  regD = zeros(iter,1);  bt = Y'*b;
   %for i=1:iter
   %   resD(i) = sqrt( sum( bt(i:end).^2) );
   %   regD(i) = sqrt( sum((theta(1:i) .* chit(1:i)).^2) );
   %end
   traD = 1 - sum(theta(1:ntrunc-1))/sum(theta);
   regD = zeros(niter,1); resD = zeros(niter,1); bt = Y'*b;
   for i = 1:iter+1
      chiD   = inv(Theta1)*Y'*b; chiD(i:end) = 0;
      ItereD = Y*chiD;
      %AI = Y*Theta1*Y'*ItereD;
      %ResD = AI-b;
      %resD(i) = sqrt(norm(ResD));  %Problem with the preconditionner
      resD(i) = sqrt( sum( bt(i:end).^2) );  
      regD(i) = sqrt(ItereD(indexa)'*Kreg*ItereD(indexa));%norm(ItereD);
   end
   
   %figure;
   %hold on;
   %plot(log10(theta),'Color','blue')
   %plot(log10(abs(Y'*b)),'Color','red')
   %plot(log10(abs(chi)),'Color','black')
   %legend('Ritz Values','RHS values','solution coefficients')
   %
%   errorI = (Itere(3:3:end)-uref(3:3:end))/norm(uref(3:3:end));
%   nerrorI = norm(Itere(indexa)-uref(indexa))/norm(uref(indexa));
%   figure;
%   ret = patch('Faces',patch2,'Vertices',nodes, ... 
%               'FaceVertexCData',errorI,'FaceColor','interp');
%   colorbar;
%   %
%   errorR = (ItereR(3:3:end)-uref(3:3:end))/norm(uref(3:3:end));
%   nerrorR = norm(ItereR(indexa)-uref(indexa))/norm(uref(indexa));
%   figure;
%   ret = patch('Faces',patch2,'Vertices',nodes, ... 
%               'FaceVertexCData',errorR,'FaceColor','interp');
%   colorbar;
   %
   figure;
   ret = patch('Faces',patch2,'Vertices',nodes, ... 
               'FaceVertexCData',Itere(3:3:end),'FaceColor','interp');
   colorbar;
   %
   figure;
   ret = patch('Faces',patch2,'Vertices',nodes, ... 
               'FaceVertexCData',ItereR(3:3:end),'FaceColor','interp');
   colorbar;
   %
   %hold on;
   %plot(Itere(2*b2node2-1),'Color','red')
   %plot(ItereR(2*b2node2-1),'Color','blue')
   %plot(uref(2*b2node2-1),'Color','green')
   %legend('brutal solution','filtred solution', 'reference')
   %figure;
   %
   %hold on;
   %plot(Itere(2*b2node2),'Color','red')
   %plot(ItereR(2*b2node2),'Color','blue')
   %plot(uref(2*b2node2),'Color','green')
   %legend('brutal solution','filtred solution', 'reference')
   %figure;
   
   figure;
   hold on
   %plot(log10(error(2:end)),'Color','blue')
   plot(error(2:end),'Color','blue')
   plot(residual(2:end),'Color','red')
   plot(errS(2:end),'Color','green')
   legend('error','residual','Ritz error')
   
   %L-curve :
   figure;
   hold on;
   loglog(residual(2:iter+1),regulari(2:iter+1),'Color','red','-*');
   loglog(resS(2:iter+1),regS(2:iter+1),'-+');
   legend('L-curve','RL-curve')
   %legend('RL-curve')
   %findCorner (residual(2:iter+1), regulari(2:iter+1), 3)
   %findCorner (resS(2:iter+1), regS(2:iter+1), 3)
   
   %%hold on;
   %%loglog(resS(2:iter+1),regS(2:iter+1),'Color','red','-*');
   %loglog(resD(2:iter),regD(2:iter),'-+');
   %%legend('RL-curve (natural basis)','RL-curve (diagonal basis)')
   %figure
   %findCorner (resD(2:iter), regD(2:iter), 3)
   
   %%%%%
   %% Final problem : compute u ( not necessary for now)
   %dirichlet = [4,1,0;4,2,0;
   %             3,1,0;3,2,0;
   %             1,1,0;1,2,0;
   %             2,1,0;2,2,0];
   %neumann   = [];
   %[K,C,nbloq] = Krig2 (nodes,elements,mat,order,boundary,dirichlet);
   %if inhomog >= 1
   %   K(1:3*nnodes, 1:3*nnodes) = Kinter;
   %end
   %fdir2 = dirichletRhs(Itere, 2, C, boundary);
   %fdir4 = dirichletRhs(urefb, 4, C, boundary);
   %usoli = K \ (fdir4 + fdir2);
   %
   %usol = usoli(1:3*nnodes,1);
   %fsol = Kinter*usol;
   %
   %% With the reduced solution
   %dirichlet = [4,1,0;4,2,0;
   %             3,1,0;3,2,0;
   %             1,1,0;1,2,0;
   %             2,1,0;2,2,0];
   %neumann   = [];
   %[K,C,nbloq] = Krig2 (nodes,elements,mat,order,boundary,dirichlet);
   %if inhomog >= 1
   %   K(1:3*nnodes, 1:3*nnodes) = Kinter;
   %end
   %fdir2 = dirichletRhs(ItereR, 2, C, boundary);
   %fdir4 = dirichletRhs(urefb, 4, C, boundary);
   %usoli = K \ (fdir4 + fdir2);
   %
   %usolR = usoli(1:3*nnodes,1);
   %fsolR = Kinter*usolR;
   %
   %hold on;
   %plot(fsol(2*b2node2-1),'Color','red')
   %plot(fsolR(2*b2node2-1),'Color','blue')
   %plot(f(2*b2node2-1),'Color','green')
   %legend('brutal solution','filtred solution','reference')
   %
   %total_error = norm(usol-uref)/norm(uref);
   %% Compute stress :
   %sigma = stress(usol,E,nu,nodes,elements,order,1,ntoelem);
   %plotGMSH({usol,'U_vect';sigma,'stress'}, elements, nodes, 'solution');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dual SP method
if methods == 2

   % First, plot the reference
   ret = patch('Faces',patch2,'Vertices',nodes, ... 
               'FaceVertexCData',f(3:3:end),'FaceColor','interp');
   colorbar;

   error    = zeros(niter+1,1);
   residual = zeros(niter+1,1);
   regulari = zeros(niter+1,1);
   
   Itere = zeros( 3*nnodes, 1 );
   d     = zeros( 3*nnodes, niter+1 );
   Ad    = zeros( 3*nnodes, niter+1 );
   Res   = zeros( 3*nnodes, niter+1 );
   Zed   = zeros( 3*nnodes, niter+1 );
   alpha = zeros( niter+1, 1 );
   beta  = zeros( niter+1, 1 );
   alpha2 = zeros( niter+1, 1 );
   %H12   = H1demi(size(Res,1), nodes, boundary, 2 );
   %Mgr   = Mgrad(size(Res,1), nodes, boundary, 2 );
   
   %% Perform A x0 :
   % Solve 1
   f1 = [Itere; zeros(nbloq1d,1)];
   uin1 = K1d\f1;
   u1i = uin1(1:3*nnodes,1);
   u1 = C12*C12'*u1i;
   % Solve 2
   f2 = [Itere; zeros(nbloq2d,1)];
   uin2 = K2d\f2;
   u2i = uin2(1:3*nnodes,1);
   u2 = C22*C22'*u2i;
   %
   Axz = u2-u1;
   
   %%%%
   %% Compute Rhs :
   % Solve 1
   f1 = [ zeros(3*nnodes,1) ; C1d'*C11*C11'*urefb ];
   uin1 = K1d\f1;
   u1i = uin1(1:3*nnodes);
   u1 = C12*C12'*u1i;
   % Solve 2 (zero 'cause no loading except on 7)
   f2 = loading3D(nbloq2d,nodes,boundary,neumann0);
   uin2 = K2d\f2;
   u2i = uin2(1:3*nnodes);
   u2 = C22*C22'*u2i;
   %
   b = u1-u2;
   %%
   Res(:,1) = b - Axz;
   
   if precond == 1
       f1 = [ zeros(3*nnodes,1) ; C1'*C12*C12'*Res(:,1) ];
       uin1 = K1\f1;
       lagr1 = uin1(3*nnodes+1:end,1);
       lamb1 = -C12*C12'*C1*lagr1;
       Zed(:,1) = lamb1;
   elseif precond == 2 % Not working
       Zed(index,1) = E*H12(index,index)\Res(index,1);
   elseif precond == 3
       Zed(index,1) = E*Mgr(index,index)\Res(index,1);
   else
       Zed(:,1) = Res(:,1);
   end
   
   d(:,1) = Zed(:,1);
   
   residual(1) = norm(Res( :,1));
   error(1)    = norm(Itere(indexa) - fref(indexa)) / norm(fref(indexa));
   %regulari(1) = sqrt( Itere'*regul(Itere, nodes, boundary, 2) );
   regulari(1) = norm(Itere);
   
   ritzval  = 0; % Last ritz value that converged
   oldtheta = 0;
   eta      = 0;
   getmeout = 0; % utility
   %V = zeros(3*nnodes, iter);
   %H = zeros(iter);
   %%
   for iter = 1:niter
       %% Optimal step
       % Solve 1
       f1 = [d(:,iter); zeros(nbloq1d,1)];
       uin1 = K1d\f1;
       u1i = uin1(1:3*nnodes,1);
       u1 = C12*C12'*u1i;
       % Solve 2
       f2 = [d(:,iter); zeros(nbloq2d,1)];
       uin2 = K2d\f2;
       u2i = uin2(1:3*nnodes,1);
       u2 = C12*C12'*u2i;
       %
       Ad(:,iter) = u2-u1;
       
       den = (d(:,iter)'*Ad(:,iter));
       d(:,iter) = d(:,iter)/sqrt(den); Ad(:,iter) = Ad(:,iter)/sqrt(den);
       num = Res(:,iter)'*d(:,iter);
       
       Itere         = Itere + d(:,iter)*num;%/den;
       Res(:,iter+1) = Res(:,iter) - Ad(:,iter)*num;%/den;
       
       residual(iter+1) = norm(Res(:,iter+1));
       error(iter+1)    = norm(Itere(indexa) - fref(indexa)) / norm(fref(indexa));
   %    regulari(iter+1) = sqrt( Itere'*regul(Itere, nodes, boundary, 2) );
       regulari(iter+1) = norm(Itere);
       
       if precond == 1
          f1 = [ zeros(3*nnodes,1) ; C1'*C12*C12'*Res(:,iter+1) ];
          uin1 = K1\f1;
          lagr1 = uin1(3*nnodes+1:end,1);
          lamb1 = -C12*C12'*C1*lagr1;
          Zed(:,iter+1) = lamb1;
       elseif precond == 2
           Zed(index,iter+1) = 1/E*H12(index,index)\Res(index,iter+1);
       elseif precond == 3
           Zed(index,iter+1) = 1/E*Mgr(index,index)\Res(index,iter+1);
       else
           Zed(:,iter+1) = Res(:,iter+1);
       end
       
       % Needed values for the Ritz stuff
       alpha(iter) = num/sqrt(den);
       alpha2(iter) = num;
       beta(iter)  = - Zed(:,iter+1)'*Ad(:,iter)/sqrt(den);
       
       % First Reorthogonalize the residual (as we use it next), in sense of M
       for jter=1:iter-1
           betac = Zed(:,iter+1)'*Res(:,jter) / (Zed(:,jter)'*Res(:,jter));
           Zed(:,iter+1) = Zed(:,iter+1) - Zed(:,jter) * betac;
       end
       
       %% Orthogonalization
       d(:,iter+1) = Zed(:,iter+1);
       
       for jter=iter:iter % No need to reorthogonalize (see above)
           betaij = ( Zed(:,iter+1)'*Ad(:,jter) );%/...
               %( d(:,jter)'*Ad(:,jter) );
           d(:,iter+1) = d(:,iter+1) - d(:,jter) * betaij;
   
       end
       
       %% Ritz algo : find the Ritz elements
       % Build the matrices
       V(:,iter) = zeros(3*nnodes,1);
       V(:,iter) = (-1)^(iter-1)*Zed(:,iter)/(sqrt(Res(:,iter)'*Zed(:,iter)));
       %norm(Zed(:,iter))
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
   
   %%for i=1:niter
   %   plot(V(2*b2node2-1,20:25))
   %   figure;
   %%end
   
   regS = zeros(niter,1);
   resS = zeros(niter,1);
   %% Build the L-curve regul, ntrunc
   for i = 1:iter+1
      chiS   = inv(Theta1)*Y'*b; chiS(i:end) = 0;
      ItereS = Y*chiS;
      
      % Solve 1
      f1 = [ItereS; zeros(nbloq1d,1)];
      uin1 = K1d\f1;
      u1i = uin1(1:3*nnodes,1);
      u1 = C12*C12'*u1i;
      % Solve 2
      f2 = [ItereS; zeros(nbloq2d,1)];
      uin2 = K2d\f2;
      u2i = uin2(1:3*nnodes,1);
      u2 = C12*C12'*u2i;
      %
      AI = u2-u1;
      
      ResS = AI-b;
      resS(i) = norm(ResS);   
   %   regS(i) = sqrt( ItereS'*regul(ItereS, nodes, boundary, 2) );
      regS(i) = norm(ItereS);
      errS(i) = norm(ItereS(indexa)-fref(indexa))/norm(fref(indexa));
   end
   
   % Residual in the diagonal base :
   %resD = zeros(iter,1);  regD = zeros(iter,1);  bt = Y'*b;
   %for i=1:iter
   %   resD(i) = sqrt( sum( bt(i:end).^2) );
   %   regD(i) = sqrt( sum((theta(1:i) .* chit(1:i)).^2) );
   %end
   traD = 1 - sum(theta(1:ntrunc-1))/sum(theta);
   regD = zeros(niter,1); resD = zeros(niter,1); bt = Y'*b;
   for i = 1:iter+1
      chiD   = inv(Theta1)*Y'*b; chiD(i:end) = 0;
      ItereD = Y*chiD;
      %AI = Y*Theta1*Y'*ItereD;
      %ResD = AI-b;
      %resD(i) = sqrt(norm(ResD));  %Problem with the preconditionner
      resD(i) = sqrt( sum( bt(i:end).^2) );  
   %   regD(i) = sqrt( ItereD'*regul(ItereD, nodes, boundary, 2) );
      regD(i) = norm(ItereD);
   end
   
   %figure;
   %hold on;
   %plot(log10(theta),'Color','blue')
   %plot(log10(abs(Y'*b)),'Color','red')
   %plot(log10(abs(chi)),'Color','black')
   %legend('Ritz Values','RHS values','solution coefficients')
   %
%   errorI = (Itere(3:3:end)-uref(3:3:end))/norm(uref(3:3:end));
%   nerrorI = norm(Itere(indexa)-uref(indexa))/norm(uref(indexa));
%   figure;
%   ret = patch('Faces',patch2,'Vertices',nodes, ... 
%               'FaceVertexCData',errorI,'FaceColor','interp');
%   colorbar;
%   %
%   errorR = (ItereR(3:3:end)-uref(3:3:end))/norm(uref(3:3:end));
%   nerrorR = norm(ItereR(indexa)-uref(indexa))/norm(uref(indexa));
%   figure;
%   ret = patch('Faces',patch2,'Vertices',nodes, ... 
%               'FaceVertexCData',errorR,'FaceColor','interp');
%   colorbar;
   %
   figure;
   ret = patch('Faces',patch2,'Vertices',nodes, ... 
               'FaceVertexCData',Itere(3:3:end),'FaceColor','interp');
   colorbar;
   %
   figure;
   ret = patch('Faces',patch2,'Vertices',nodes, ... 
               'FaceVertexCData',ItereR(3:3:end),'FaceColor','interp');
   colorbar;
   %
   %hold on;
   %plot(Itere(2*b2node2-1),'Color','red')
   %plot(ItereR(2*b2node2-1),'Color','blue')
   %plot(uref(2*b2node2-1),'Color','green')
   %legend('brutal solution','filtred solution', 'reference')
   %figure;
   %
   %hold on;
   %plot(Itere(2*b2node2),'Color','red')
   %plot(ItereR(2*b2node2),'Color','blue')
   %plot(uref(2*b2node2),'Color','green')
   %legend('brutal solution','filtred solution', 'reference')
   %figure;
   
   figure;
   hold on
   %plot(log10(error(2:end)),'Color','blue')
   plot(error(2:end),'Color','blue')
   plot(residual(2:end),'Color','red')
   plot(errS(2:end),'Color','green')
   legend('error','residual','Ritz error')
%   legend('error','residual')
   
   %L-curve :
   figure;
   hold on;
   loglog(residual(2:iter+1),regulari(2:iter+1),'Color','red','-*');
   loglog(resS(2:iter+1),regS(2:iter+1),'-+');
   legend('L-curve','RL-curve')
   %legend('RL-curve')
   %findCorner (residual(2:iter+1), regulari(2:iter+1), 3)
   %findCorner (resS(2:iter+1), regS(2:iter+1), 3)
   
   %%hold on;
   %%loglog(resS(2:iter+1),regS(2:iter+1),'Color','red','-*');
   %loglog(resD(2:iter),regD(2:iter),'-+');
   %%legend('RL-curve (natural basis)','RL-curve (diagonal basis)')
   %figure
   %findCorner (resD(2:iter), regD(2:iter), 3)
   
   %%%%%
   %% Final problem : compute u ( not necessary for now)
   %dirichlet = [4,1,0;4,2,0;
   %             3,1,0;3,2,0;
   %             1,1,0;1,2,0;
   %             2,1,0;2,2,0];
   %neumann   = [];
   %[K,C,nbloq] = Krig2 (nodes,elements,mat,order,boundary,dirichlet);
   %if inhomog >= 1
   %   K(1:3*nnodes, 1:3*nnodes) = Kinter;
   %end
   %fdir2 = dirichletRhs(Itere, 2, C, boundary);
   %fdir4 = dirichletRhs(urefb, 4, C, boundary);
   %usoli = K \ (fdir4 + fdir2);
   %
   %usol = usoli(1:3*nnodes,1);
   %fsol = Kinter*usol;
   %
   %% With the reduced solution
   %dirichlet = [4,1,0;4,2,0;
   %             3,1,0;3,2,0;
   %             1,1,0;1,2,0;
   %             2,1,0;2,2,0];
   %neumann   = [];
   %[K,C,nbloq] = Krig2 (nodes,elements,mat,order,boundary,dirichlet);
   %if inhomog >= 1
   %   K(1:3*nnodes, 1:3*nnodes) = Kinter;
   %end
   %fdir2 = dirichletRhs(ItereR, 2, C, boundary);
   %fdir4 = dirichletRhs(urefb, 4, C, boundary);
   %usoli = K \ (fdir4 + fdir2);
   %
   %usolR = usoli(1:3*nnodes,1);
   %fsolR = Kinter*usolR;
   %
   %hold on;
   %plot(fsol(2*b2node2-1),'Color','red')
   %plot(fsolR(2*b2node2-1),'Color','blue')
   %plot(f(2*b2node2-1),'Color','green')
   %legend('brutal solution','filtred solution','reference')
   %
   %total_error = norm(usol-uref)/norm(uref);
   %% Compute stress :
   %sigma = stress(usol,E,nu,nodes,elements,order,1,ntoelem);
   %plotGMSH({usol,'U_vect';sigma,'stress'}, elements, nodes, 'solution');

end