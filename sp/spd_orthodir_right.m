% 18/10/2016
% Algo Steklov-Poincar� dual avec Orthodir et précond à droite

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 70000;  % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 1;      % N.mm-1 : Loading on the plate
precond = 1;
niter   = 20;
mu      = 100;    % Regularization parameter
br      = 0.;     % noise

% Boundary conditions
% first index  : index of the boundary
% second index : 1=x, 2=y
% third        : value
% [0,1,value] marks a dirichlet regularization therm on x
dirichlet = [4,1,0;
             4,2,0];
neumann   = [1,2,fscalar;
             2,1,fscalar;
             3,2,fscalar];

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

boundaryp1 = suppressBound( boundary, [no3;no4], 3 );
%boundaryp2 = suppressBound( boundary, [no4], 3 );
boundaryp2 = boundaryp1;
%boundary = suppressBound( boundary, [no3;no4], 3 );
% Suppress no2 from 1 for a Schur complement computation
boundarym = suppressBound( boundary, [no2], 2 );
boundarym = suppressBound( boundarym, [no1], 1 );

% Then, build the stiffness matrix :
[K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
Kinter = K(1:2*nnodes, 1:2*nnodes);
M      = mass_mat(nodes, elements);
[node2b4, b2node4] = mapBound(4, boundaryp1, nnodes);
[node2b3, b2node3] = mapBound(3, boundaryp1, nnodes);
[node2b1, b2node1] = mapBound(1, boundarym, nnodes);
[node2b2, b2node2] = mapBound(2, boundarym, nnodes);
b2node12 = [b2node1;b2node2];
% The right hand side :
f = loading(nbloq,nodes,boundary,neumann);

% Solve the problem :
uin = K\f;

% Extract displacement and Lagrange multiplicators :
uref = uin(1:2*nnodes,1);
lagr = uin(2*nnodes+1:end,1);

urefb = ( 1 + br*randn(2*nnodes,1) ) .* uref;

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
dirichlet2 = [4,1,0;4,2,0;
              3,1,0;3,2,0];
neumann2   = [1,2,fscalar;
              2,1,fscalar];
neumann0   = [];
[K2,C2,nbloq2,node2c2,c2node2] = Krig (nodes,elements,E,nu,order,boundaryp2,dirichlet2);

error    = zeros(niter+1,1);
residual = zeros(niter+1,1);
regulari = zeros(niter+1,1);

%% Dual problems
% First problem
dirichlet1d = [4,1,0;4,2,0;
               2,1,0;2,2,0;
               1,1,0;1,2,0];
[K1d,C1d,nbloq1d] = Krig (nodes,elements,E,nu,order,boundary,dirichlet1d);

% Second problem
dirichlet2d = [4,1,0;4,2,0];
[K2d,C2d,nbloq2d] = Krig (nodes,elements,E,nu,order,boundary,dirichlet2d);

%% Anti-cancellation trick
index = [2*b2node3-1; 2*b2node3];
K1(index,index) = 0;
K2(index,index) = 0;
K1d(index,index) = 0;
K2d(index,index) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ORTHODIR for the problem : (D20-D10)(S10-S20) x = D2-D1
Itere = zeros( 2*nnodes, 1 );
p     = zeros( 2*nnodes, niter+1 );
q     = zeros( 2*nnodes, niter+1 );
Delta = zeros( niter+1, 1 );
Res   = zeros( 2*nnodes, niter+1 );
Zed   = zeros( 2*nnodes, niter+1 );

%% Perform A x0 :
% Solve 1
f1 = dirichletRhs2( Itere, 3, c2node1, boundaryp1, nnodes );
uin1 = K1\f1;
lagr1 = uin1(2*nnodes+1:end,1);
lamb1 = lagr2forces( lagr1, C1, 3, boundaryp1 );
% Solve 2
f2 = dirichletRhs2( Itere, 3, c2node2, boundaryp2, nnodes );
uin2 = K2\f2;
lagr2 = uin2(2*nnodes+1:end,1);
lamb2 = lagr2forces( lagr2, C2, 3, boundaryp2 );
%
Axz = lamb1-lamb2;
% Solve 1 (second computation)
f1 = [Axz; zeros(nbloq1d,1)];
uin1 = K1d\f1;
u1i = uin1(1:2*nnodes,1);
u1 = keepField( u1i, 3, boundaryp1 );
% Solve 2
f2 = [Axz; zeros(nbloq2d,1)];
uin2 = K2d\f2;
u2i = uin2(1:2*nnodes,1);
u2 = keepField( u2i, 3, boundaryp2 );
% Regularization term
Nu = regul(Itere, nodes, boundaryp2, 3);
Axz = mu*Nu + u1-u2;
%%%%
%% Compute Rhs :
% Solve 1
f11 = dirichletRhs(urefb, 1, C1d, boundary);
f12 = dirichletRhs(urefb, 2, C1d, boundary);
f1 = max(f11, f12);
uin1 = K1d\f1;
u1i = uin1(1:2*nnodes,1);
u1 = keepField( u1i, 3, boundary );
% Solve 2
f2 = loading(nbloq2d,nodes,boundary,neumann2);
uin2 = K2d\f2;
u2i = uin2(1:2*nnodes,1);
u2 = keepField( u2i, 3, boundary );
%
b = -u1+u2;
% hold on;
% urb = reshape(b(map1,1),2,[])';
% plot(urb(:,1), 'Color', 'red');
% urb = reshape(btot,2,[])';
% plot(urb(:,1), 'Color', 'blue');
% figure;
%%%%
Res(:,1) = b - Axz;
% urb = reshape(Res(map1,1),2,[])';
% plot(urb(:,1));
% figure;
Zed(:,1) = Res(:,1);

p(:,1) = Zed(:,1);

%residual(1) = sqrt( myps( Zed(:,1), Zed(:,1), Kinter, boundaryp1, M, nodes ) );
residual(1) = sqrt( myps( Res(:,1), Res(:,1), Kinter, boundary, M, nodes ) );
error(1)    = sqrt( myps( Itere - uref, Itere - uref, Kinter, boundary, M, nodes )...
                 / myps( uref, uref, Kinter, boundary, M, nodes ) );
regulari(1) = sqrt( Itere'*regul(Itere, nodes, boundary, 3) );

%% Perform Q1 = A P1 :
% Solve 1
f1 = dirichletRhs(p(:,1), 3, C1, boundaryp1);
uin1 = K1\f1;
lagr1 = uin1(2*nnodes+1:end,1);
lamb1 = lagr2forces( lagr1, C1, 3, boundaryp1 );
% Solve 2
f2 = dirichletRhs(p(:,1), 3, C2, boundaryp2);
uin2 = K2\f2;
lagr2 = uin2(2*nnodes+1:end,1);
lamb2 = lagr2forces( lagr2, C2, 3, boundaryp2 );
%
q(:,1) = lamb1-lamb2;
% Solve 1 (second computation)
f1 = [q(:,1); zeros(nbloq1d,1)];
uin1 = K1d\f1;
u1i = uin1(1:2*nnodes,1);
u1 = keepField( u1i, 3, boundaryp1 );
% Solve 2
f2 = [q(:,1); zeros(nbloq2d,1)];
uin2 = K2d\f2;
u2i = uin2(1:2*nnodes,1);
u2 = keepField( u2i, 3, boundaryp2 );
% Regularization term
Nu = regul(p(:,1), nodes, boundary, 3);
Ari = mu*Nu + u1-u2;
%%%%
%%
for iter = 1:niter
    
    Delta(iter,1) = myps( q(:,iter), q(:,iter), Kinter, boundary, M, nodes ); %q(:,iter)'*q(:,iter);
    gammai        = myps( q(:,iter), Res(:,iter), Kinter, boundary, M, nodes ); %q(:,iter)'*Res;
    alphai        = gammai/Delta(iter,1);
    
    Itere         = Itere + p(:,iter)*alphai;
    Res(:,iter+1) = Res(:,iter) - q(:,iter)*alphai;
    
    Zed(:,iter+1) = Res(:,iter+1);
    
    %residual(iter+1) = sqrt( myps( Zed(:,iter+1), Zed(:,iter+1), Kinter, boundaryp1, M, nodes ) );
    residual(iter+1) = sqrt( myps( Res(:,iter+1), Res(:,iter+1), Kinter, boundary, M, nodes ) );
    error(iter+1)    = sqrt( myps( Itere - uref, Itere - uref, Kinter, boundary, M, nodes )...
                     / myps( uref, uref, Kinter, boundary, M, nodes ) );
    regulari(iter+1) = sqrt( Itere'*regul(Itere, nodes, boundary, 3) );

    %% Perform Ari = A*Res
    % Solve 1
    rhs1 = Zed(:,iter+1);
    f1 = dirichletRhs(rhs1, 3, C1, boundaryp1);
    uin1 = K1\f1;
    lagr1 = uin1(2*nnodes+1:end,1);
    lamb1 = lagr2forces( lagr1, C1, 3, boundaryp1 );
    % Solve 2
    rhs2 = Zed(:,iter+1);
    f2 = dirichletRhs(rhs2, 3, C2, boundaryp2);
    uin2 = K2\f2;
    lagr2 = uin2(2*nnodes+1:end,1);
    lamb2 = lagr2forces( lagr2, C2, 3, boundaryp2 );
    %
    Ari = lamb1-lamb2;
    % Solve 1 (second computation)
    f1 = [Ari; zeros(nbloq1d,1)];
    uin1 = K1d\f1;
    u1i = uin1(1:2*nnodes,1);
    u1 = keepField( u1i, 3, boundaryp1 );
    % Solve 2
    f2 = [Ari; zeros(nbloq2d,1)];
    uin2 = K2d\f2;
    u2i = uin2(1:2*nnodes,1);
    u2 = keepField( u2i, 3, boundaryp2 );
    % Regularization term
    Nu = regul(Zed(:,iter+1), nodes, boundaryp2, 3);
    Ari = mu*Nu + u1-u2;
    
    %% Orthogonalization
    p(:,iter+1) = Zed(:,iter+1);
    q(:,iter+1) = Ari;
    
    for jter=1:iter
        phiij  = myps( q(:,jter), Ari, Kinter, boundary, M, nodes ); %q(:,jter)'*Ari;
        betaij = phiij/Delta(jter,1);
        p(:,iter+1) = p(:,iter+1) - p(:,jter) * betaij;
        q(:,iter+1) = q(:,iter+1) - q(:,jter) * betaij;
    end
    
end

hold on
% plot(error,'Color','blue')
% plot(residual,'Color','red')
plot(log10(error),'Color','blue')
plot(log10(residual),'Color','red')
legend('error (log)','residual (log)')
figure;
% L-curve :
loglog(residual,regulari);
figure
%%%%
%% Final problem : compute u
% DN problem
dirichlet = [4,1,0;4,2,0;
             3,1,0;3,2,0;
             1,1,0;1,2,0;
             2,1,0;2,2,0];
neumann   = [];
[K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
fdir1 = dirichletRhs(urefb, 1, C, boundary);
fdir2 = dirichletRhs(urefb, 2, C, boundary);
fdir3 = dirichletRhs(Itere, 3, C, boundary);
usoli = K \ assembleDirichlet( [fdir1+fdir3,fdir2] );%( max(fdir1, max(fdir2, fdir3)) );

usol = usoli(1:2*nnodes,1);
fsol = Kinter*usol;

plot(fsol(2*b2node3))

total_error = norm(usol-uref)/norm(uref);
% Compute stress :
sigma = stress(usol,E,nu,nodes,elements,order,1,ntoelem);
plotGMSH({usol,'U_vect';sigma,'stress'}, elements, nodes, 'solution');