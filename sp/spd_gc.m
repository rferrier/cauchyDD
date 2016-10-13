% 06/06/2016
% Algo Steklov-Poincar� primal avec Gradient Conjugu�

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 70000;  % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 1;      % N.mm-1 : Loading on the plate
niter   = 20;
precond = 0;      % 1 : Use a dual precond, 2 : use the H1/2 precond
mu      = 0.;    % Regularization parameter
br      = 0.;     % noise

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

fref = Kinter*uref;

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
[K1d,C1d,nbloq1d,node2c1d,c2node1d] = Krig (nodes,elements,E,nu,order,boundary,dirichlet1d);

% Second problem
dirichlet2d = [1,1,0;1,2,0;
               3,1,0;3,2,0];
[K2d,C2d,nbloq2d,node2c2d,c2node2d] = Krig (nodes,elements,E,nu,order,boundary,dirichlet2d);

%% Anti-cancellation trick
% K1(indexa,indexa) = 0;
% K2(indexa,indexa) = 0;
% K1d(indexa,indexa) = 0;
% K2d(indexa,indexa) = 0;

%% Conjugate Gradient for the problem : (D20-D10) x = D1-D2
Itere = zeros( 2*nnodes, 1 );
d     = zeros( 2*nnodes, niter+1 );
Ad    = zeros( 2*nnodes, niter+1 );
Res   = zeros( 2*nnodes, niter+1 );
Zed   = zeros( 2*nnodes, niter+1 );
H12   = H1demi(size(Res,1), nodes, boundary, 2 );

%% Perform A x0 :
% Solve 1
f1 = [Itere; zeros(nbloq1d,1)];
uin1 = K1d\f1;
u1i = uin1(1:2*nnodes,1);
u1 = keepField( u1i, 2, boundaryp1 );
% Solve 2
f2 = [Itere; zeros(nbloq2d,1)];
uin2 = K2d\f2;
u2i = uin2(1:2*nnodes,1);
u2 = keepField( u2i, 2, boundaryp2 );
% Regularization term
Nu = regul(Itere, nodes, boundary, 2);
%
Axz = mu*Nu+u2-u1;
%%%%
%% Compute Rhs :
% Solve 1
f1 = dirichletRhs2( urefb, 4, c2node1d, boundary, nnodes );
uin1 = K1d\f1;
u1i = uin1(1:2*nnodes,1);
u1 = keepField( u1i, 2, boundaryp1 );
% Solve 2
f2 = loading(nbloq2d,nodes,boundary,neumann2);
uin2 = K2d\f2;
u2i = uin2(1:2*nnodes,1);
u2 = keepField( u2i, 2, boundaryp2 );
%
b = u1-u2;

%%
Res(:,1) = b - Axz;

if precond == 1
    % Solve 1
    f1 = dirichletRhs2( Res(:,1)/2, 2, c2node1, boundaryp1, nnodes );
    uin1 = K1\f1;
    lagr1 = uin1(2*nnodes+1:end,1);
    lamb1 = lagr2forces( lagr1, C1, 2, boundaryp1 );
    % Solve 2
    f2 = dirichletRhs2( Res(:,1)/2, 2, c2node2, boundaryp2, nnodes );
    uin2 = K2\f2;
    lagr2 = uin2(2*nnodes+1:end,1);
    lamb2 = lagr2forces( lagr2, C2, 2, boundaryp2 );
    %
    Zed(:,1) = -lamb1/2+lamb2/2;
elseif precond == 2
    Zed(index,1) = E*H12(index,index)*Res(index,1);
else
    Zed(:,1) = Res(:,1);
end

d(:,1) = Zed(:,1);

residual(1) = norm(Res( indexa,1));
error(1)    = norm(Itere(indexa) - fref(indexa) / norm(fref(indexa));
regulari(1) = sqrt( Itere'*regul(Itere, nodes, boundary, 2) );

%%
for iter = 1:niter
    %% Optimal step
    
    % Solve 1
    f1 = [d(:,iter); zeros(nbloq1d,1)];
    uin1 = K1d\f1;
    u1i = uin1(1:2*nnodes,1);
    u1 = keepField( u1i, 2, boundaryp1 );
    % Solve 2
    f2 = [d(:,iter); zeros(nbloq2d,1)];
    uin2 = K2d\f2;
    u2i = uin2(1:2*nnodes,1);
    u2 = keepField( u2i, 2, boundaryp2 );
    % Regularization term
    Nu = regul(d(:,iter), nodes, boundary, 2);
    %
    Ad(:,iter) = mu*Nu+u2-u1;
    
    den = (d(indexa,iter)'*Ad(indexa,iter));
    d(:,iter) = d(:,iter)/sqrt(den); Ad(:,iter) = Ad(:,iter)/sqrt(den);
    num = Res(indexa,iter)'*d(indexa,iter);

    Itere         = Itere + d(:,iter)*num;
    Res(:,iter+1) = Res(:,iter) - Ad(:,iter)*num;
    
    residual(iter+1) = norm(Res(indexa,iter+1));
    error(iter+1)    = norm(Itere(indexa) - fref(indexa)) / norm(fref(indexa));
    regulari(iter+1) = sqrt( Itere'*regul(Itere, nodes, boundary, 2) );
    
    if precond == 1
        % Solve 1
        f1 = dirichletRhs2( Res(:,iter+1)/2, 2, c2node1, boundaryp1, nnodes );
        uin1 = K1\f1;
        lagr1 = uin1(2*nnodes+1:end,1);
        lamb1 = lagr2forces( lagr1, C1, 2, boundaryp1 );
        % Solve 2
        f2 = dirichletRhs2( Res(:,iter+1)/2, 2, c2node2, boundaryp2, nnodes );
        uin2 = K2\f2;
        lagr2 = uin2(2*nnodes+1:end,1);
        lamb2 = lagr2forces( lagr2, C2, 2, boundaryp2 );
        %
        Zed(:,iter+1) = -lamb1/2+lamb2/2;
    elseif precond == 2
        Zed(index,iter+1) = E*H12(index,index)*Res(index,iter+1);
    else
        Zed(:,iter+1) = Res(:,iter+1);
    end
    
    %% Orthogonalization
    d(:,iter+1) = Zed(:,iter+1);
    
    for jter=1:iter
        betaij = ( Zed(indexa,iter+1)'*Ad(indexa,jter) );
        d(:,iter+1) = d(:,iter+1) - d(:,jter) * betaij;
    end    
end

hold on
plot(log10(error),'Color','blue')
plot(log10(residual),'Color','red')
legend('error (log)','residual (log)')
figure;
% L-curve :
loglog(residual(2:end),regulari(2:end));
figure
%%%%
%% Final problem : compute u
% DN problem
dirichlet = [4,1,0;4,2,0;
             3,1,0;3,2,0;
             1,1,0;1,2,0];
neumann   = [];
[K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
f2    = [keepField( Itere, 2, boundary ); zeros(nbloq,1)];
fdir4 = dirichletRhs(urefb, 4, C, boundary);
usoli = K \ (fdir4 + f2);

usol = usoli(1:2*nnodes,1);
fsol = Kinter*usol;

plot(fsol(2*b2node2-1))

total_error = norm(usol-uref)/norm(uref);
% Compute stress :
sigma = stress(usol,E,nu,nodes,elements,order,1,ntoelem);
plotGMSH({usol,'U_vect';sigma,'stress'}, elements, nodes, 'solution');