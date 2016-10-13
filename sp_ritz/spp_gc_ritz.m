% 08/09/2016
% Algo Steklov-Poincaré primal avec Gradient Conjugué

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
ratio   = 5e-200;    % Maximal ratio (for eigenfilter)
br      = .0;      % noise
brt     = 0;      % "translation" noise
epsilon = 1e-1;   % Convergence criterion for ritz value

% Boundary conditions
% first index  : index of the boundary
% second index : 1=x, 2=y
% third        : value
% [0,1,value] marks a dirichlet regularization therm on x
dirichlet = [1,1,0; 1,2,0 ;
             3,1,0; 3,2,0 ];
neumann   = [2,1,fscalar,0,fscalar;
             4,1,fscalar,0,-fscalar];

% Import the mesh
[ nodes,elements,ntoelem,boundary,order ] = readmesh( 'meshes/plate.msh' );
nnodes = size(nodes,1);

noises = load('./noises/noise1.mat'); % Particular noise vector
noise  = noises.bruit1;
%noise = randn(2*nnodes,1);

% find the nodes in the corners and suppress the element :
xmax = max(nodes(:,1));
xmin = min(nodes(:,1));
ymax = max(nodes(:,2));
ymin = min(nodes(:,2));
no1  = findNode(xmin, ymin, nodes, 1e-5);
no2  = findNode(xmax, ymin, nodes, 1e-5);
no3  = findNode(xmax, ymax, nodes, 1e-5);
no4  = findNode(xmin, ymax, nodes, 1e-5);

boundaryp1 = suppressBound( boundary, no2, 2 );
boundaryp1 = suppressBound( boundaryp1, no3, 2 );
boundaryp2 = boundaryp1;

[~, b2node2] = mapBound(2, boundaryp1, nnodes);
indexa = [2*b2node2-1; 2*b2node2];

% Then, build the stiffness matrix :
[K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
Kinter = K(1:2*nnodes, 1:2*nnodes);
M      = mass_mat(nodes, elements);
[node2b4, b2node4] = mapBound(4, boundaryp1, nnodes);
[node2b3, b2node3] = mapBound(3, boundaryp1, nnodes);
[node2b1, b2node1] = mapBound(1, boundaryp1, nnodes);
[node2b2, b2node2] = mapBound(2, boundaryp1, nnodes);
index = [2*b2node2-1;2*b2node2];
% The right hand side :
f = loading(nbloq,nodes,boundary,neumann);

% Solve the problem :
uin = K\f;

% Extract displacement and Lagrange multiplicators :
uref = uin(1:2*nnodes,1);
lagr = uin(2*nnodes+1:end,1);

urefb = ( 1 + brt + br*noise ) .* uref;

sigma = stress(uref,E,nu,nodes,elements,order,1,ntoelem);
plotGMSH({uref,'Vect_U';sigma,'stress'}, elements, nodes, 'reference');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition of the operators

% First problem
dirichlet1 = [4,1,0;4,2,0;
              3,1,0;3,2,0;
              2,1,0;2,2,0;
              1,1,0;1,2,0];

[K1,C1,nbloq1,node2c1,c2node1] = Krig (nodes,elements,E,nu,order,boundaryp1,dirichlet1);

% Second problem
dirichlet2 = [1,1,0;1,2,0;
              2,1,0;2,2,0;
              3,1,0;3,2,0];
neumann2   = [4,1,fscalar,0,-fscalar];
neumann0   = [];
[K2,C2,nbloq2,node2c2,c2node2] = Krig (nodes,elements,E,nu,order,boundaryp2,dirichlet2);

error    = zeros(niter+1,1);
residual = zeros(niter+1,1);
regulari = zeros(niter+1,1);

%% Dual problems
% First problem
dirichlet1d = [4,1,0;4,2,0;
               3,1,0;3,2,0;
               1,1,0;1,2,0];
[K1d,C1d,nbloq1d] = Krig (nodes,elements,E,nu,order,boundary,dirichlet1d);

% Second problem
dirichlet2d = [1,1,0;1,2,0;
               3,1,0;3,2,0];
[K2d,C2d,nbloq2d] = Krig (nodes,elements,E,nu,order,boundary,dirichlet2d);

%% Anti-cancellation trick
K1r = K1; K2r = K2; K1dr = K1d; K2dr = K2d;
%K1(indexa,indexa) = 0;
%K2(indexa,indexa) = 0;
%K1d(indexa,indexa) = 0;
%K2d(indexa,indexa) = 0;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Conjugate Gradient for the problem : (S10-S20) x = S2-S1
Itere = zeros( 2*nnodes, 1 );
d     = zeros( 2*nnodes, niter+1 );
Ad    = zeros( 2*nnodes, niter+1 );
Res   = zeros( 2*nnodes, niter+1 );
Zed   = zeros( 2*nnodes, niter+1 );
alpha = zeros( niter+1, 1 );
beta  = zeros( niter+1, 1 );
alpha2 = zeros( niter+1, 1 );
ntrunc = 0;  % In case the algo finishes at niter
H12   = H1demi(size(Res,1), nodes, boundary, 2 );
Mgr   = Mgrad(size(Res,1), nodes, boundary, 2 );

%% Perform A x0 :
% Solve 1
f1 = dirichletRhs2( Itere, 2, c2node1, boundaryp1, nnodes );
uin1 = K1\f1;
lagr1 = uin1(2*nnodes+1:end,1);
lamb1 = lagr2forces( lagr1, C1, 2, boundaryp1 );
% Solve 2
f2 = dirichletRhs2( Itere, 2, c2node2, boundaryp2, nnodes );
uin2 = K2\f2;
lagr2 = uin2(2*nnodes+1:end,1);
lamb2 = lagr2forces( lagr2, C2, 2, boundaryp2 );
% Regularization term
Nu = regul(Itere, nodes, boundary, 2);
%
Axz = mu*Nu+lamb1-lamb2;
%%%%
%% Compute Rhs :
% Solve 1
f1 = dirichletRhs2( urefb, 4, c2node1, boundary, nnodes );
uin1 = K1\f1;
lagr1 = uin1(2*nnodes+1:end,1);
lamb1 = lagr2forces2( lagr1, c2node1, 2, boundaryp1, nnodes );
% Solve 2
f2 = loading(nbloq2,nodes,boundary,neumann2);
uin2 = K2\f2;
lagr2 = uin2(2*nnodes+1:end,1);
lamb2 = lagr2forces2( lagr2, c2node2, 2, boundaryp2, nnodes );
%
b = lamb2-lamb1;
%plot(b(2*b2node2-1));
%figure;
%%
Res(:,1) = b - Axz;

if precond == 1
    % Solve 1
    f1 = [Res(:,1); zeros(nbloq1d,1)];
    uin1 = K1dr\f1;
    u1i = uin1(1:2*nnodes,1);
    u1 = keepField( u1i, 2, boundaryp1 );
    Zed(:,1) = u1;
elseif precond == 2
    Zed(index,1) = 1/E*H12(index,index)\Res(index,1);
elseif precond == 3
    Zed(index,1) = 1/E*Mgr(index,index)\Res(index,1);
else
    Zed(:,1) = Res(:,1);
end

d(:,1) = Zed(:,1);

residual(1) = sqrt( norm(Res( indexa,1)));
error(1)    = sqrt( norm(Itere(indexa) - uref(indexa)) / norm(uref(indexa)));
regulari(1) = sqrt( Itere'*regul(Itere, nodes, boundary, 2) );

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
    rhs1 = d(:,iter);
    f1 = dirichletRhs(rhs1, 2, C1, boundaryp1);
    uin1 = K1\f1;
    lagr1 = uin1(2*nnodes+1:end,1);
    lamb1 = lagr2forces( lagr1, C1, 2, boundaryp1 );
    % Solve 2
    rhs2 = d(:,iter);
    f2 = dirichletRhs(rhs2, 2, C2, boundaryp2);
    uin2 = K2\f2;
    lagr2 = uin2(2*nnodes+1:end,1);
    lamb2 = lagr2forces( lagr2, C2, 2, boundaryp2 );
    % Regularization term
    Nu = regul(d(:,iter), nodes, boundaryp2, 2);
    %
    Ad(:,iter) = mu*Nu+lamb1-lamb2;
    
    den = (d(indexa,iter)'*Ad(indexa,iter));
    d(:,iter) = d(:,iter)/sqrt(den); Ad(:,iter) = Ad(:,iter)/sqrt(den);
    num = Res(indexa,iter)'*d(indexa,iter);
    
    Itere         = Itere + d(:,iter)*num;%/den;
    Res(:,iter+1) = Res(:,iter) - Ad(:,iter)*num;%/den;
    
    residual(iter+1) = sqrt( norm(Res(indexa,iter+1)));
    error(iter+1)    = sqrt( norm(Itere(indexa) - uref(indexa)) / norm(uref(indexa)));
    regulari(iter+1) = sqrt( Itere'*regul(Itere, nodes, boundary, 2) );
    
    if precond == 1
        % Solve 1
        f1 = [Res(:,iter+1); zeros(nbloq1d,1)];
        uin1 = K1dr\f1;
        u1i = uin1(1:2*nnodes,1);
        u1 = keepField( u1i, 2, boundaryp1 );
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
    beta(iter)  = - Zed(indexa,iter+1)'*Ad(indexa,iter)/sqrt(den);
    
    % First Reorthogonalize the residual (as we use it next), in sense of M
    for jter=1:iter-1
        betac = Zed(indexa,iter+1)'*Res(indexa,jter) / (Zed(indexa,jter)'*Res(indexa,jter));
        Zed(:,iter+1) = Zed(:,iter+1) - Zed(:,jter) * betac;
    end
    
    %% Orthogonalization
    d(:,iter+1) = Zed(:,iter+1);
    
    for jter=iter:iter % No need to reorthogonalize (see above)
        betaij = ( Zed(indexa,iter+1)'*Ad(indexa,jter) );%/...
            %( d(indexa,jter)'*Ad(indexa,jter) );
        d(:,iter+1) = d(:,iter+1) - d(:,jter) * betaij;

    end
    
    %% Ritz algo : find the Ritz elements
    % Build the matrices
    V(:,iter) = zeros(2*nnodes,1);
    V(indexa,iter) = (-1)^(iter-1)*Zed(indexa,iter)/(sqrt(Res(indexa,iter)'*Zed(indexa,iter)));
    %norm(Zed(indexa,iter))
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
   rhs1 = ItereS;
   f1 = dirichletRhs(rhs1, 2, C1, boundaryp1);
   uin1 = K1\f1;
   lagr1 = uin1(2*nnodes+1:end,1);
   lamb1 = lagr2forces( lagr1, C1, 2, boundaryp1 );
   % Solve 2
   rhs2 = ItereS;
   f2 = dirichletRhs(rhs2, 2, C2, boundaryp2);
   uin2 = K2\f2;
   lagr2 = uin2(2*nnodes+1:end,1);
   lamb2 = lagr2forces( lagr2, C2, 2, boundaryp2 );
   % Regularization term
   AI = lamb1-lamb2;
   
   ResS = AI-b;
   resS(i) = sqrt(norm(ResS));   
   regS(i) = sqrt( ItereS'*regul(ItereS, nodes, boundary, 2) );
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
   regD(i) = sqrt( ItereD'*regul(ItereD, nodes, boundary, 2) );
end

hold on;
plot(log10(theta),'Color','blue')
plot(log10(abs(Y'*b)),'Color','red')
plot(log10(abs(chi)),'Color','black')
legend('Ritz Values','RHS values','solution coefficients')
figure;
%
hold on;
plot(Itere(2*b2node2-1),'Color','red')
plot(ItereR(2*b2node2-1),'Color','blue')
plot(uref(2*b2node2-1),'Color','green')
legend('brutal solution','filtred solution', 'reference')
figure;

hold on;
plot(Itere(2*b2node2),'Color','red')
plot(ItereR(2*b2node2),'Color','blue')
plot(uref(2*b2node2),'Color','green')
legend('brutal solution','filtred solution', 'reference')
figure;

%hold on
%plot(log10(error(2:end)),'Color','blue')
%plot(log10(residual(2:end)),'Color','red')
%legend('error (log)','residual (log)')
%figure;
%L-curve :
hold on;
loglog(residual(2:iter+1),regulari(2:iter+1),'Color','red','-*');
%figure
loglog(resS(2:iter+1),regS(2:iter+1),'-+');
legend('L-curve','RL-curve')
figure
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
             1,1,0;1,2,0;
             2,1,0;2,2,0];
neumann   = [];
[K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
fdir2 = dirichletRhs(Itere, 2, C, boundary);
fdir4 = dirichletRhs(urefb, 4, C, boundary);
usoli = K \ (fdir4 + fdir2);

usol = usoli(1:2*nnodes,1);
fsol = Kinter*usol;

% With the reduced solution
dirichlet = [4,1,0;4,2,0;
             3,1,0;3,2,0;
             1,1,0;1,2,0;
             2,1,0;2,2,0];
neumann   = [];
[K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
fdir2 = dirichletRhs(ItereR, 2, C, boundary);
fdir4 = dirichletRhs(urefb, 4, C, boundary);
usoli = K \ (fdir4 + fdir2);

usolR = usoli(1:2*nnodes,1);
fsolR = Kinter*usolR;

hold on;
plot(fsol(2*b2node2-1),'Color','red')
plot(fsolR(2*b2node2-1),'Color','blue')
plot(f(2*b2node2-1),'Color','green')
legend('brutal solution','filtred solution','reference')

total_error = norm(usol-uref)/norm(uref);
% Compute stress :
sigma = stress(usol,E,nu,nodes,elements,order,1,ntoelem);
plotGMSH({usol,'U_vect';sigma,'stress'}, elements, nodes, 'solution');