% 26/01/2017
% Algo Steklov-Poincaré primal avec Gradient Conjugué, cas 3D

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 70000;  % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 1;      % N.mm-1 : Loading on the plate
niter   = 25;
precond = 0;      % 1 : Use a dual precond, 2 : use H1/2 precond, 3 : use gradient precond
mu      = 0.;     % Regularization parameter
ratio   = 5e-200; % Maximal ratio (for eigenfilter)
br      = 0.;    % noise
brt     = 0;      % "translation" noise
epsilon = 1e-200; % Convergence criterion for ritz value
ntrunc  = 0;      % In case the algo finishes at niter
inhomog = 0;      % inhomogeneous medium
assembl = 0;      % Shall I assemble the operators ?

% methods : 1-> SPP, 2-> SPD, 3-> SPD Block
methods = [3];

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
neumann   = [7,3,-fscalar, ];

% Import the mesh
[ nodes,elements,ntoelem,boundary,order ] = readmesh3D( 'meshes/plate3d_charge2.msh' );
nnodes = size(nodes,1);

noises = load('./noises/noise3d0.mat'); % Particular noise vector
noise  = noises.noise;
%noise = randn(3*nnodes,1);

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

%% Debug : get the index of
if assembl == 1
   boundaryp3 = boundaryp1 = suppressBound( boundaryp1, 'extreme', 1, nodes,1e-5 );
   C11 = Cbound ( nodes, [1,1,0;1,2,0;1,3,0], boundaryp3 );
   indexC11= sum(C11'); index1 = find(indexC11)';
   index2 = indexa';
   C3456 = Cbound ( nodes, ...
                    [3,1,0;3,2,0;3,3,0 ; 4,1,0;4,2,0;4,3,0 ;...
                     5,1,0;5,2,0;5,3,0 ; 6,1,0;6,2,0;6,3,0], boundaryp3 );
   indexC3456= sum(C3456'); index3456 = find(indexC3456)';
   
   Krr = Kinter(index1, index1);
   Kgg = Kinter( index2, index2 );
   Kjj = Kinter;
   Kjj([index1;index2;index3456],:) = [];
   Kjj(:,[index1;index2;index3456]) = [];
   
   Krj = Kinter( index1, : );
   Krj(:,[index1;index2;index3456]) = [];
   
   Kgj = Kinter( index2, : );
   Kgj( :, [index1;index2;index3456] ) = [];
   
   Ikjj = inv(Kjj);
   SR = full(Krr - Krj*Ikjj*Krj'); Crm = Krj*Ikjj*Kgj';
   K1k = Crm'*inv(SR)*Crm;
   
   ur = uin(index1); um = uin(index2);
   tr = f(index1); tm = f(index2);
   b = Kgj*Ikjj*Krj'*(ur - SR\tr);
   
   Crmp = Crm*inv(Kgg-Kgj*Ikjj*Kgj');
   umtilde = Crm'*((Crm*K1k*Crm')\(Crm*K1k*um));
   ur2 = SR\Crm*um;
   umtildep = Crmp'*((Crmp*K1k*Crmp')\(Crmp*K1k*um));
   umorth = um-umtilde; umorthp = um-umtildep;
%   norm(umtilde)
   errlimu = norm(umorth)/norm(um); errlimup = norm(umorthp)/norm(um);
   %%
   Ikgg = inv(Kgg);
   L = Kjj-Kgj'*Ikgg*Kgj; Il = inv(L);
   SR2 = Krr-Krj*Il*Krj'; Crm2 = Krj*Il*Kgj'*Ikgg;
   K2k = Crm2'*inv(SR2)*Crm2;
   
   tmtilde = Crm2'*((Crm2*K2k*Crm2')\(Crm2*K2k*tm));
   tmorth = tm-tmtilde;
   umdtilde = (Kgg-Kgj*Ikjj*Kgj')\tmtilde+Crm2'*ur; umdorth = um-umdtilde;
%   norm(tmtilde)
   errlimt = norm(tmorth)/norm(tm); errlimud = norm(umdorth)/norm(um);
end

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
neumannd   = [1,3,fscalar];
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
   
   [indm2,pol] = findPicard2(log10(abs(chiD)), 10, 1);
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
          'polynomial interpolation')
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
%   figure;
%   ret = patch('Faces',patch2,'Vertices',nodes, ... 
%               'FaceVertexCData',Itere(3:3:end),'FaceColor','interp');
%   colorbar;
%   %
   figure;
   ret = patch('Faces',patch2,'Vertices',nodes, ... 
               'FaceVertexCData',ItereR(3:3:end),'FaceColor','interp');
   colorbar;
   
   % Plot the stuff on a line x = 4
%   Yma = max(nodes(:,2)); step = Yma/100; Ymi = min(nodes(:,2)); Zma = max(nodes(:,3));
%   nodeploty = Ymi:step:Yma; nodeplotx = 4*ones(1,101); nodeplotz = Zma*ones(1,101);
%   nodeplot = [nodeplotx;nodeploty;nodeplotz]';
%   uplo = passMesh3D(nodes, elements, nodeplot, [], [Itere,ItereR,uref]);
%   figure;
%   hold on;
%   plot(nodeploty, uplo(3:3:end,1),'Color','red');
%   plot(nodeploty, uplo(3:3:end,2),'Color','blue');
%   plot(nodeploty, uplo(3:3:end,3),'Color','green');
%   legend('brutal solution','filtred solution', 'reference')
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
%   plot(error(2:end),'Color','blue')
%   plot(residual(2:end),'Color','red')
   plot(errS(2:end),'Color','green')
%   legend('error','residual','Ritz error')
   legend('Ritz error')
   
   %L-curve :
   figure;
   hold on;
   set(gca, 'fontsize', 20);
   loglog(residual(2:iter+1),regulari(2:iter+1),'Color','red','-*','linewidth',3);
   loglog(resS(2:iter+1),regS(2:iter+1),'-+','linewidth',3);
   legend('L-curve','RL-curve')
   xlabel('Residual')
   ylabel('H1 norm')
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

%   % First, plot the reference
%   ret = patch('Faces',patch2,'Vertices',nodes, ... 
%               'FaceVertexCData',f(3:3:end),'FaceColor','interp');
               
   % First, plot the reference (sigma normal)
   ret = patch('Faces',patch2,'Vertices',nodes, ... 
               'FaceVertexCData',sigma(3:6:end),'FaceColor','interp');
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
   
   % Find the best Ritz value
   [indm2,pol] = findPicard2(log10(abs(chiD)), 10, 1);
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
          'polynomial interpolation')
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
%   figure;
%   ret = patch('Faces',patch2,'Vertices',nodes, ... 
%               'FaceVertexCData',Itere(3:3:end),'FaceColor','interp');
%   colorbar;
%   %
%   figure;
%   ret = patch('Faces',patch2,'Vertices',nodes, ... 
%               'FaceVertexCData',ItereR(3:3:end),'FaceColor','interp');
%   colorbar;
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
%   legend('Ritz error')
%   legend('error','residual')
   
   %L-curve :
   figure;
   hold on;
   set(gca, 'fontsize', 20);
   loglog(residual(2:iter+1),regulari(2:iter+1),'Color','red','-*','linewidth',3);
   loglog(resS(2:iter+1),regS(2:iter+1),'-+','linewidth',3);
   legend('L-curve','RL-curve')
   xlabel('Residual')
   ylabel('H1 norm')
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
   %% Final resolution
   f_D = [ zeros(3*nnodes,1) ; C1d'*C1d*C1d'*uref ];
   f_N = [ItereR; zeros(nbloq1d,1)];
   uin  = K1d\(f_N+f_D);

   usol = uin(1:3*nnodes,1);
   fsol = Kinter*usol;
   
   sigmasol = stress3D(usol,mat,nodes,elements,order,1,ntoelem);
   plotGMSH3D({usol,'Vect_U';sigmasol,'stress'}, elements, nodes, 'output/solution');
   
   % plot the solution (sigma normal)
   figure;
   ret = patch('Faces',patch2,'Vertices',nodes, ... 
               'FaceVertexCData',sigmasol(3:6:end),'FaceColor','interp');
   colorbar;
   
   erroru = norm(usol(indexa) - uref(indexa)) / norm(uref(indexa));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Conjugate Gradient for the problem : (D10-D20) x = D2-D1
if methods==3
%   % First, plot the reference (sigma normal)
%   ret = patch('Faces',patch2,'Vertices',nodes, ... 
%               'FaceVertexCData',sigma(3:6:end),'FaceColor','interp');

   % First, plot the reference
   ret = patch('Faces',patch2,'Vertices',nodes, ... 
               'FaceVertexCData',uref(3:3:end),'FaceColor','interp');
   colorbar;
   
   Itere = zeros( 3*nnodes, 2 );
   d     = zeros( 3*nnodes, 2*(niter+1) );  % 2 directions per step
   Ad    = zeros( 3*nnodes, 2*(niter+1) );
   Res   = zeros( 3*nnodes, 2*(niter+1) );
   Zed   = zeros( 3*nnodes, 2*(niter+1) );
   AZed  = zeros( 3*nnodes, 2*(niter+1) );
   H12   = H1demi(size(Res,1), nodes, boundary, 2 );
   Mgr   = Mgrad(size(Res,1), nodes, boundary, 2 );
   alpha = zeros( 2*(niter+1) );
   beta  = zeros( 2*(niter+1) );
   
   % Equivalent Saad matrices
   unsuralpha   = zeros( 2*(niter+1) );
   betasuralpha = zeros( 2*(niter+1) );
   etaeta       = zeros( 2*(niter+1) );
   
   %% Perform A x0 :
   % Solve 1
   f1 = [Itere(:,1); zeros(nbloq1d,1)];
   uin1 = K1d\f1;
   u1i = uin1(1:3*nnodes,1);
   u1 = C12*C12'*u1i;
   % Solve 2
   f2 = [Itere(:,2); zeros(nbloq2d,1)];
   uin2 = K2d\f2;
   u2i = uin2(1:3*nnodes,1);
   u2 = C22*C22'*u2i;
   %
   Axz = -[u1,u2];
   
   %%%%
   %% Compute Rhs :
   % Solve 1
   f1 = [ zeros(3*nnodes,1) ; C1d'*C11*C11'*urefb ];
   uin1 = K1d\f1;
   u1i = uin1(1:3*nnodes);
   u1 = C12*C12'*u1i;
   % Solve 2 (zero 'cause no loading except on 7)
%   f2 = loading3D(nbloq2d,nodes,boundary,neumann0);
%   f2 = loading3D(nbloq2d,nodes,boundary,neumannd); % Dummy rhs
   f2 = zeros( 3*nnodes+nbloq2d, 1 );
   f2(indexa,1) = randn( size(indexa,2), 1 );
   uin2 = K2d\f2;
   u2i = uin2(1:3*nnodes);
   u2 = C22*C22'*u2i;
%   u2 = zeros( 3*nnodes, 1 );
%   u2(indexa,1) = norm(u1)*randn( size(indexa,2), 1 );
   %
   b = [u1,u2]; %!!!!!!!!!! Problem because ZERO
   %%
   Res(:,[1,2]) = b - Axz;
   
   if precond == 1
       f11 = [ zeros(3*nnodes,1) ; C1'*C12*C12'*Res(:,1) ];
       f12 = [ zeros(3*nnodes,1) ; C1'*C12*C12'*Res(:,2) ];
       uin1 = K1\[f11,f12];
       lagr1 = uin1(3*nnodes+1:end,:);
       lamb1 = -C12*C12'*C1*lagr1;
       Zed(:,[1,2]) = lamb1;
   elseif precond == 2 % Not working
       % Solve 2 (Sn)
       f21 = dirichletRhs2( Res(:,1), 2, c2node2, boundaryp2, nnodes );
       f22 = dirichletRhs2( Res(:,2), 2, c2node2, boundaryp2, nnodes );
       uin2 = K2\[f21,f22];
       lagr2 = uin2(2*nnodes+1:end,:);
       lamb21 = lagr2forces2( lagr2(:,1), c2node2, 2, boundaryp2, nnodes );
       lamb22 = lagr2forces2( lagr2(:,2), c2node2, 2, boundaryp2, nnodes );
       %
       Zed(:,[1,2]) = [lamb21,lamb22];
   else
       Zed(:,[1,2]) = Res(:,[1,2]);
   end
   
   d(:,[1,2]) = Zed(:,[1,2]);
   
   residual(1) = norm( Res(indexa,1)-Res(indexa,2) );
   error(1)    = norm(Itere(indexa,1) - fref(indexa)) / ...
                                       norm(fref(indexa));
   regulari(1) = sqrt( (Itere(:,1)')* ... 
                        regul( Itere(:,1) , nodes, boundary, 2) );
   %%
   V  = zeros(3*nnodes, 2*niter);
   AV = zeros(3*nnodes, 2*niter);
   MV = zeros(3*nnodes, 2*niter);
   H  = zeros(2*niter);
   
   num = [0,0;0,0]; % useless, but eta needs initialization #lazy
   den = [0,0;0,0];
   %%
   for iter = 1:niter

       %% Optimal step
       % Solve 1
       f1 = [d(:,[2*iter-1,2*iter]); zeros(nbloq1d,2)];
       uin1 = K1d\f1;
       u1i = uin1(1:3*nnodes,:);
       u1 = C12*C12'*u1i;
       % Solve 2
       f2 = [d(:,[2*iter-1,2*iter]); zeros(nbloq2d,2)];
       uin2 = K2d\f2;
       u2i = uin2(1:3*nnodes,:);
       u2 = C12*C12'*u2i;
       %
       Ad(:,[2*iter-1,2*iter]) = u2-u1;
   
       denprec = den; numprec = num; % Store those ones
       den = d(indexa,[2*iter-1,2*iter])'*Ad(indexa,[2*iter-1,2*iter]);
       sqD = den^(1/2);
       d(:,[2*iter-1,2*iter]) = d(:,[2*iter-1,2*iter]) * inv(sqD);
       Ad(:,[2*iter-1,2*iter]) = Ad(:,[2*iter-1,2*iter]) * inv(sqD);
       num = Res(indexa,[2*iter-1,2*iter])'*Zed(indexa,[2*iter-1,2*iter]);
       num = sqD\num; % because of Zed and not d
       
   %    Itere = Itere + d(:,[2*iter-1,2*iter])*alpha;
   %    Res(:,[2*iter+1,2*iter+2]) = Res(:,[2*iter-1,2*iter]) - ...
   %                                    Ad(:,[2*iter-1,2*iter])*alpha;
       Itere = Itere + d(:,[2*iter-1,2*iter])*num;
       Res(:,[2*iter+1,2*iter+2]) = Res(:,[2*iter-1,2*iter]) - ...
                                       Ad(:,[2*iter-1,2*iter])*num;
       
       residual(iter+1) = norm( Res(indexa,2*iter+1) );
       error(iter+1)    = norm(Itere(indexa,1) - fref(indexa)) / ...
                                       norm(fref(indexa));
       regulari(iter+1) = sqrt( (Itere(:,1)')* ... 
                             regul( Itere(:,1) , nodes, boundary, 2) );
       
       if precond == 1
          f11 = [ zeros(3*nnodes,1) ; C1'*C12*C12'*Res(:,2*iter+1) ];
          f12 = [ zeros(3*nnodes,1) ; C1'*C12*C12'*Res(:,2*iter+2) ];
          uin1 = K1\[f11,f12];
          lagr1 = uin1(3*nnodes+1:end,:);
          lamb1 = -C12*C12'*C1*lagr1;
          Zed(:,[2*iter+1,2*iter+2]) = lamb1;
       elseif precond == 2 % Not working
          % Solve 2 (Sn)
          f21 = dirichletRhs2( Res(:,2*iter+1), 2, c2node2, boundaryp2, nnodes );
          f22 = dirichletRhs2( Res(:,2*iter+2), 2, c2node2, boundaryp2, nnodes );
          uin2 = K2\[f21,f22];
          lagr2 = uin2(2*nnodes+1:end,:);
          lamb21 = lagr2forces2( lagr2(:,1), c2node2, 2, boundaryp2, nnodes );
          lamb22 = lagr2forces2( lagr2(:,2), c2node2, 2, boundaryp2, nnodes );
          %
          Zed(:,[2*iter+1,2*iter+2]) = [lamb21,lamb22];
       else
           Zed(:,[2*iter+1,2*iter+2]) = Res(:,[2*iter+1,2*iter+2]);
       end
       
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
       
       % First Reorthogonalize the residual (as we use it next), in sense of M
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
       for jter=iter:iter
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
   
   % Compute the solution
   chi = inv(Theta1)*Y'*( b(:,1) );
   if ntrunc > 0
      chi(ntrunc:end) = 0;
   end
   ItereR = Y*chi;
   
   % Build residual and such
   regD = zeros(niter,1); resD = zeros(niter,1); bt = Y'*b;
   for i = 1:2*iter+1
      chiD   = inv(Theta1)*Y'*(b(:,1)); chiD(i:end) = 0;
      ItereD = Y*chiD;
      resD(i) = sqrt( sum( bt(i:end).^2) );  
      regD(i) = sqrt( ItereD'*regul(ItereD, nodes, boundary, 2) );
   end
   
   % Compute the number of pertinent Ritz modes
   [indm2,pol] = findPicard2(log10(abs(chiD)), 3, 1);
   n = size(pol,1);
   t = 1:.05:2*niter; tt = zeros(n,20*(2*niter-1)+1);
   for j=1:n
      tt(j,:) = t.^(n-j);
   end
   px = pol'*tt;
   
   figure;
   hold on
   plot(error,'Color','blue')
   plot(residual,'Color','red')
   legend('error','residual')
   %% L-curve :
   %loglog(residual(2:end),regulari(2:end));
   %figure
   
   figure;
   hold on;
   plot(log10(theta),'Color','blue')
   plot(log10(abs(Y'*( b(:,1) ))),'Color','red')
   plot(log10(abs(chiD)),'Color','black')
   plot(t,px,'Color','cyan')
   legend( 'Ritz Values','RHS values','solution coefficients', ...
           'polynomial approximation' )

   %%%%%%%
   %% Final resolution
   f_D = [ zeros(3*nnodes,1) ; C1d'*C1d*C1d'*uref ];
   f_N = [ItereR; zeros(nbloq1d,1)];
   uin  = K1d\(f_N+f_D);

   usol = uin(1:3*nnodes,1);
   fsol = Kinter*usol;
   
   sigmasol = stress3D(usol,mat,nodes,elements,order,1,ntoelem);
   plotGMSH3D({usol,'Vect_U';sigmasol,'stress'}, elements, nodes, 'output/solution');
   relerror = abs(usol-uref)/max(abs(uref(indexa)));
   plotGMSH3D({relerror,'Vect_U'}, elements, nodes, 'output/error');
   
   % plot the solution (sigma normal)
   figure;
%   ret = patch('Faces',patch2,'Vertices',nodes, ... 
%               'FaceVertexCData',sigmasol(3:6:end),'FaceColor','interp');
%   colorbar;
   maxuref = max(abs(uref(3:3:end)));
   ret = patch('Faces',patch2,'Vertices',nodes, ... 
               'FaceVertexCData',abs(uref(3:3:end)-usol(3:3:end))/maxuref, ...
               'FaceColor','interp');
   colorbar;
   
   erroru = norm(usol(indexa) - uref(indexa)) / norm(uref(indexa));
           
end