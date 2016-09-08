% 06/06/2016
% Algo Steklov-Poincar� primal avec Gradient Conjugu�

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 70000;  % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 1;      % N.mm-1 : Loading on the plate
niter   = 28;
precond = 1;      % 1 : Use a dual precond
mu      = 0.;    % Regularization parameter
ratio   = 1e7;    % Maximal ratio (for eigenfilter)
br      = 0;     % noise
brt     = 0;   % "translation" noise

% Boundary conditions
% first index  : index of the boundary
% second index : 1=x, 2=y
% third        : value
% [0,1,value] marks a dirichlet regularization therm on x
dirichlet = [1,1,0; 1,2,0 ;
             3,1,0; 3,2,0 ];
neumann   = [2,1,fscalar;
             4,1,fscalar];

% Import the mesh
[ nodes,elements,ntoelem,boundary,order ] = readmesh( 'meshes/plate.msh' );
nnodes = size(nodes,1);

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

urefb = ( 1 + brt + br*randn(2*nnodes,1) ) .* uref;

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
neumann2   = [4,1,fscalar];
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
% K1(indexa,indexa) = 0;
% K2(indexa,indexa) = 0;
% K1d(indexa,indexa) = 0;
% K2d(indexa,indexa) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Conjugate Gradient for the problem : (S10-S20) x = S2-S1
Itere = zeros( 2*nnodes, 1 );
d     = zeros( 2*nnodes, niter+1 );
Ad    = zeros( 2*nnodes, niter+1 );
Res   = zeros( 2*nnodes, niter+1 );
Zed   = zeros( 2*nnodes, niter+1 );
alpha = zeros( niter+1, 1 );
beta  = zeros( niter+1, 1 );
%H12   = H1demi(size(Res,1), nodes, boundary, 2 );
%Mgr   = Mgrad(size(Res,1), nodes, boundary, 2 );

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
b = -lamb1+lamb2;

%%
Res(:,1) = b - Axz;

if precond == 1
    % Solve 1
    f1 = [Res(:,1); zeros(nbloq1d,1)];
    uin1 = K1d\f1;
    u1i = uin1(1:2*nnodes,1);
    u1 = keepField( u1i, 2, boundaryp1 );
    Zed(:,1) = u1;
else
    Zed(:,1) = Res(:,1);
end

d(:,1) = Zed(:,1);

residual(1) = sqrt( norm(Res( indexa,1)));
error(1)    = sqrt( norm(Itere(indexa) - uref(indexa)) / norm(uref(indexa)));
regulari(1) = sqrt( Itere'*regul(Itere, nodes, boundary, 2) );

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
    
    Itere         = Itere + d(:,iter)*num;
    Res(:,iter+1) = Res(:,iter) - Ad(:,iter)*num;
    
    residual(iter+1) = sqrt( norm(Res(indexa,iter+1)));
    error(iter+1)    = sqrt( norm(Itere(indexa) - uref(indexa)) / norm(uref(indexa)));
    regulari(iter+1) = sqrt( Itere'*regul(Itere, nodes, boundary, 2) );
    
    if precond == 1
        % Solve 1
        f1 = [Res(:,iter+1); zeros(nbloq1d,1)];
        uin1 = K1d\f1;
        u1i = uin1(1:2*nnodes,1);
        u1 = keepField( u1i, 2, boundaryp1 );
        Zed(:,iter+1) = u1;
    else
        Zed(:,iter+1) = Res(:,iter+1);
    end
    
    % Needed values for the Ritz stuff
    alpha(iter)   = num;
    beta(iter+1)  = -(Zed(indexa,iter+1)'*Ad(indexa,iter));
    
    %% Orthogonalization
    d(:,iter+1) = Zed(:,iter+1);
    
    for jter=1:iter
        betaij = ( Zed(indexa,iter+1)'*Ad(indexa,jter) );%/...
            %( d(indexa,jter)'*Ad(indexa,jter) );
        d(:,iter+1) = d(:,iter+1) - d(:,jter) * betaij;
    end

%     betaij = ( Zed(indexa,iter+1)'*Ad(indexa,iter) )/...
%         ( d(indexa,iter)'*Ad(indexa,iter) );
%     d(:,iter+1) = d(:,iter+1) - d(:,iter) * betaij;
    
end

% The last one (for Ritz analysis)
% Solve 1
rhs1 = d(:,niter+1);
f1 = dirichletRhs(rhs1, 2, C1, boundaryp1);
uin1 = K1\f1;
lagr1 = uin1(2*nnodes+1:end,1);
lamb1 = lagr2forces( lagr1, C1, 2, boundaryp1 );
% Solve 2
rhs2 = d(:,niter+1);
f2 = dirichletRhs(rhs2, 2, C2, boundaryp2);
uin2 = K2\f2;
lagr2 = uin2(2*nnodes+1:end,1);
lamb2 = lagr2forces( lagr2, C2, 2, boundaryp2 );
% Regularization term
Nu = regul(d(:,niter+1), nodes, boundaryp2, 2);
%
Ad(:,niter+1) = mu*Nu+lamb1-lamb2;

alpha(niter+1) = Res(indexa,niter+1)'*d(indexa,niter+1)/...
                 sqrt(d(indexa,niter+1)'*Ad(indexa,niter+1));

%% Ritz algo : find the Ritz elements
V = zeros(2*nnodes, niter+1);
H = zeros(niter+1, niter+1);

% Build the matrices
V(:,1) = Zed(:,1)/(sqrt(Res(:,1)'*Zed(:,1)));
delta  = 1/alpha(1);
eta    = sqrt(beta(1))/alpha(1);
H(1,[1,2]) = [delta, eta];

for i=1:niter-1
   V(:,i+1) = (-1)^i*Zed(:,i+1)/(sqrt(Res(:,i+1)'*Zed(:,i+1)));
   etap = eta;
   delta = 1/alpha(i+1) + beta(i)/alpha(i);
   eta = sqrt(beta(i+1))/alpha(i+1);
   H(i+1,[i,i+1,i+2]) = [etap,delta,eta];
end

V(:,niter+1) = (-1)^niter*Zed(:,niter+1)/(sqrt(Res(:,niter+1)'*Zed(:,niter+1)));
delta = 1/alpha(niter+1) + beta(niter)/alpha(niter);
H(niter+1,[niter,niter+1]) = [eta,delta];

% Compute eigenelems of the Hessenberg :
[Q,Theta1] = eig(H);
theta = diag(Theta1);
Y = V*Q;

%% Determine the truncature
%ntrunc = 0;
%for i=1:niter
%   if theta(niter+1-i) < ratio*theta(niter+1)  % /!\ Order must be growing
%      ntrunc = niter+1-i
%      break;
%   end
%end
%
%% Actually do the truncature
%chi = Y'*Itere;
%if ntrunc > 0
%   chi(1:ntrunc) = 0;
%end

% Determine the truncature (growing)
ntrunc = 0;
for i=2:niter+1
   if theta(i) > ratio*theta(1)  % /!\ Order must be growing
      ntrunc = i
      break;
   end
end

% Actually do the truncature
chi = Y'*Itere;
if ntrunc > 0
   chi(ntrunc:end) = 0;
   if precond == 1  % Remember Y'MY = 1
      % Apply M = S1 = D1^-1
      rhs1 = Y*chi;
      f1 = dirichletRhs(rhs1, 2, C1, boundaryp1);
      uin1 = K1\f1;
      lagr1 = uin1(2*nnodes+1:end,1);
      lamb1 = lagr2forces( lagr1, C1, 2, boundaryp1 );
      ItereR = lamb1;
   else
      ItereR = Y*chi;
   end
else
   ItereR = Itere;
end

hold on;
plot(Itere(2*b2node2-1),'Color','red')
plot(ItereR(2*b2node2-1),'Color','blue')
plot(uref(2*b2node2-1),'Color','green')
legend('brutal solution','filtred solution', 'reference')
figure;

hold on
plot(log10(error(2:end)),'Color','blue')
plot(log10(residual(2:end)),'Color','red')
legend('error (log)','residual (log)')
figure;
% L-curve :
loglog(residual(2:end),regulari(2:end));
figure
%%%%
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
legend('brutal solution','filtred solution')

total_error = norm(usol-uref)/norm(uref);
% Compute stress :
sigma = stress(usol,E,nu,nodes,elements,order,1,ntoelem);
plotGMSH({usol,'U_vect';sigma,'stress'}, elements, nodes, 'solution');