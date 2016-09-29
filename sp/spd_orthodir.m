% 06/04/2016
% Algo Steklov-Poincaré dual avec Orthodir

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 70000;  % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 1;      % N.mm-1 : Loading on the plate
niter   = 20;
precond = 0;      % Use a primal precond
mu      = 0./E;   % Regularization parameter
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
no3  = findNode(xmax, ymax, nodes, 1e-5);
no4  = findNode(xmin, ymax, nodes, 1e-5);
boundaryp1 = suppressBound( boundary, [no3;no4], 3 );
%boundaryp2 = suppressBound( boundary, [no4], 3 );
boundaryp2 = boundaryp1;
%boundary = suppressBound( boundary, [no3;no4], 3 );

% Then, build the stiffness matrix :
[K,C,nbloq,node2c,c2node] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
Kinter = K(1:2*nnodes, 1:2*nnodes);
M      = mass_mat(nodes, elements);
[node2b3, b2node3] = mapBound(3, boundaryp1, nnodes);
% The right hand side :
f = loading(nbloq,nodes,boundary,neumann);

% Solve the problem :
uin = K\f;

% Extract displacement and Lagrange multiplicators :
uref = uin(1:2*nnodes,1);
lagr = uin(2*nnodes+1:end,1);
fref = keepField( f(1:2*nnodes,1), 3, boundary );

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

[K1,C1,nbloq1] = Krig (nodes,elements,E,nu,order,boundaryp1,dirichlet1);

% Second problem
dirichlet2 = [4,1,0;4,2,0;
              3,1,0;3,2,0];
neumann2   = [1,2,fscalar;
              2,1,fscalar];
neumann0   = [];
[K2,C2,nbloq2] = Krig (nodes,elements,E,nu,order,boundaryp2,dirichlet2);

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

%% ORTHODIR for the problem : (S1-S2) x = b
Itere = zeros( 2*nnodes, 1 );
p     = zeros( 2*nnodes, niter+1 );
q     = zeros( 2*nnodes, niter+1 );
Delta = zeros( niter+1, 1 );
Res   = zeros( 2*nnodes, niter+1 );
Zed   = zeros( 2*nnodes, niter+1 );

%% Perform A x0 :
% Solve 1
rhs1 = [ Itere ; zeros(nbloq1d,1) ];
uin1 = K1d\rhs1;
u1i = uin1(1:2*nnodes,1);
u1 = keepField( u1i, 3, boundary );
% Solve 2
rhs2 = [ Itere ; zeros(nbloq2d,1) ];
uin2 = K2d\rhs2;
u2i = uin2(1:2*nnodes,1);
u2 = keepField( u2i, 3, boundary );
% Regularization term
Nu = regul(Itere, nodes, boundary, 3);
%
Axz = mu*Nu+u1-u2;
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
%%%%
%%
Res(:,1) = b - Axz;
% urb = reshape(Res(map1,1),2,[])';
% plot(urb(:,1));
% figure;

if precond == 1
    % Solve 1
    rhs1 = Res(:,1)/2;
    f1 = dirichletRhs(rhs1, 3, C1, boundaryp1);
    uin1 = K1\f1;
    lagr1 = uin1(2*nnodes+1:end,1);
    lamb1 = lagr2forces( lagr1, C1, 3, boundary );
    % Solve 2
    rhs2 = Res(:,1)/2;
    f2 = dirichletRhs(rhs2, 3, C2, boundaryp2);
    uin2 = K2\f2;
    lagr2 = uin2(2*nnodes+1:end,1);
    lamb2 = lagr2forces( lagr2, C2, 3, boundary );
    %
    Zed(:,1) = lamb1/2-lamb2/2;
else
    Zed(:,1) = Res(:,1);
end

p(:,1) = Zed(:,1);

%residual(1) = sqrt( myps( Zed(:,1), Zed(:,1), Kinter, boundaryp1, M, nodes ) );
residual(1) = sqrt( myps( Res(:,1), Res(:,1), Kinter, boundaryp1, M, nodes ) );
error(1)    = sqrt( myps( Itere - fref, Itere - fref, Kinter, boundaryp1, M, nodes )...
                 / myps( fref, fref, Kinter, boundaryp1, M, nodes ) );
regulari(1) = sqrt( Itere'*regul(Itere, nodes, boundaryp1, 3) );

%% Perform Q1 = A P1 :
% Solve 1
rhs1 = [ p(:,1) ; zeros(nbloq1d,1) ];
uin1 = K1d\rhs1;
u1i = uin1(1:2*nnodes,1);
u1 = keepField( u1i, 3, boundary );
% Solve 2
rhs2 = [ p(:,1) ; zeros(nbloq2d,1) ];
uin2 = K2d\rhs2;
u2i = uin2(1:2*nnodes,1);
u2 = keepField( u2i, 3, boundary );
% Regularization term
Nu = regul(p(:,1), nodes, boundary, 3);
%
q(:,1) = mu*Nu+u1-u2;
%%%%
%%
for iter = 1:niter
    
    Delta(iter,1) = myps( q(:,iter), q(:,iter), Kinter, boundaryp1, M, nodes ); %q(:,iter)'*q(:,iter);
    gammai        = myps( q(:,iter), Res(:,iter), Kinter, boundaryp1, M, nodes ); %q(:,iter)'*Res;
    alphai        = gammai/Delta(iter,1);
    
    Itere         = Itere + p(:,iter)*alphai;
    Res(:,iter+1) = Res(:,iter) - q(:,iter)*alphai;
    
    if precond == 1
        % Solve 1
        rhs1 = Res(:,iter+1)/2;
        f1 = dirichletRhs(rhs1, 3, C1, boundaryp1);
        uin1 = K1\f1;
        lagr1 = uin1(2*nnodes+1:end,1);
        lamb1 = lagr2forces( lagr1, C1, 3, boundary );
        % Solve 2
        rhs2 = Res(:,iter+1)/2;
        f2 = dirichletRhs(rhs2, 3, C2, boundaryp2);
        uin2 = K2\f2;
        lagr2 = uin2(2*nnodes+1:end,1);
        lamb2 = lagr2forces( lagr2, C2, 3, boundary );
        %
        Zed(:,iter+1) = lamb1/2-lamb2/2;
    else
        Zed(:,iter+1) = Res(:,iter+1);
    end
    
    %residual(iter+1) = sqrt( myps( Zed(:,iter+1), Zed(:,iter+1), Kinter, boundaryp1, M, nodes ) );
    residual(iter+1) = sqrt( myps( Res(:,iter+1), Res(:,iter+1), Kinter, boundaryp1, M, nodes ) );
    error(iter+1)    = sqrt( myps( Itere - fref, Itere - fref, Kinter, boundaryp1, M, nodes )...
                     / myps( fref, fref, Kinter, boundaryp1, M, nodes ) );
    regulari(iter+1) = sqrt( Itere'*regul(Itere, nodes, boundaryp1, 3) );

    %% Perform Ari = A*Res
    % Solve 1
    rhs1 = [ Zed(:,iter+1) ; zeros(nbloq1d,1) ];
    uin1 = K1d\rhs1;
    u1i = uin1(1:2*nnodes,1);
    u1 = keepField( u1i, 3, boundary );
    % Solve 2
    rhs2 = [ Zed(:,iter+1) ; zeros(nbloq2d,1) ];
    uin2 = K2d\rhs2;
    u2i = uin2(1:2*nnodes,1);
    u2 = keepField( u2i, 3, boundary );
    % Regularization term
    Nu = regul(Zed(:,iter+1), nodes, boundary, 3);
    %
    Ari = mu*Nu+u1-u2;

    %% Orthogonalization
    p(:,iter+1) = Zed(:,iter+1);
    q(:,iter+1) = Ari;
    
    for jter=1:iter
        phiij  = myps( q(:,jter), Ari, Kinter, boundaryp1, M, nodes ); %q(:,jter)'*Ari;
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
%%%%
%% Final problem : compute u
% DN problem
dirichlet = [4,1,0;4,2,0;
             1,1,0;1,2,0;
             2,1,0;2,2,0];
neumann   = [];
[K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
fdir1 = dirichletRhs(urefb, 1, C, boundary);
fdir2 = dirichletRhs(urefb, 2, C, boundary);
f3 = [ Itere ; zeros(nbloq,1) ];
usoli = K \ ( max(fdir1, fdir2) + f3 );
usol = usoli(1:2*nnodes,1);

total_error = norm(usol-uref)/norm(uref);
% Compute stress :
sigma1 = stress(usol,E,nu,nodes,elements,order,1,ntoelem);

plotGMSH({usol,'U_vect';sigma1,'stress'}, elements, nodes, 'solution');
