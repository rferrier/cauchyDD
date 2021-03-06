% 31/03/2017
% Algo Steklov-Poincaré primal avec Gradient Conjugué, anti-cancellation

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 70000;  % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 1;      % N.mm-1 : Loading on the plate
niter   = 10;
precond = 1;      % 1 : Use a dual precond, 2 : use regul precond
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
[ nodes,elements,ntoelem,boundary,order ] = readmesh( 'meshes/platee.msh' );
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
boundaryp1 = suppressBound( boundaryp1, no1, 4 );
boundaryp1 = suppressBound( boundaryp1, no4, 4 );
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
fref = Kinter*uref;
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
[K2d,C2d,nbloq2d,node2c1d,c2node1d] = Krig (nodes,elements,E,nu,order,boundary,dirichlet2d);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition of the special anti-cancellation operators
dirichletm2r = [4,1,0;4,2,0;
                3,1,0;3,2,0;
                2,1,0;2,2,0;
                1,1,0;1,2,0];
[Km2r,Cm2r,nbloqm2r,node2cm2r,c2nodem2r] = Krig (nodes,elements,E,nu,order,boundaryp1,dirichletm2r);

dirichletr2m = [4,1,0;4,2,0;
                3,1,0;3,2,0;
                2,1,0;2,2,0;
                1,1,0;1,2,0];
[Kr2m,Cr2m,nbloqr2m,node2cr2m,c2noder2m] = Krig (nodes,elements,E,nu,order,boundaryp1,dirichletr2m);

dirichletZ = [3,1,0;3,2,0;
              2,1,0;2,2,0;
              1,1,0;1,2,0];
[KZ,CZ,nbloqZ,node2cZ,c2nodeZ] = Krig (nodes,elements,E,nu,order,boundaryp1,dirichletZ);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Conjugate Gradient for the problem : (S10-S20) x = S2-S1
Itere = zeros( 2*nnodes, 1 );
usolC = zeros( 2*nnodes, 1 ); % Reconstructed total solution
d     = zeros( 2*nnodes, niter+1 );
Ad    = zeros( 2*nnodes, niter+1 );
Res   = zeros( 2*nnodes, niter+1 );
Zed   = zeros( 2*nnodes, niter+1 );
H12   = H1demi(size(Res,1), nodes, boundary, 2 );
Mgr   = Mgrad(size(Res,1), nodes, boundary, 2 );

%% Perform A x0 :
%% Solve 1
%f1 = dirichletRhs2( Itere, 2, c2node1, boundaryp1, nnodes );
%uin1 = K1\f1;
%lagr1 = uin1(2*nnodes+1:end,1);
%lamb1 = lagr2forces( lagr1, C1, 2, boundaryp1 );
%% Solve 2
%f2 = dirichletRhs2( Itere, 2, c2node2, boundaryp2, nnodes );
%uin2 = K2\f2;
%lagr2 = uin2(2*nnodes+1:end,1);
%lamb2 = lagr2forces( lagr2, C2, 2, boundaryp2 );
%%
%Axz = lamb1-lamb2;
% First resolution
f1 = dirichletRhs2( Itere, 2, c2nodem2r, boundaryp1, nnodes );
uin1 = Km2r\f1;
lagr1 = uin1(2*nnodes+1:end,1);
lamb1 = lagr2forces2( lagr1, c2nodem2r, 4, boundaryp1, nnodes );
% Z-1
uin1 = KZ \ [lamb1 ; zeros(nbloqZ,1)];
u1i = uin1(1:2*nnodes,1);
u1 = keepField( u1i, 4, boundaryp1 );
% Last resolution
f1 = dirichletRhs2( u1, 4, c2noder2m, boundaryp1, nnodes );
uin1 = Kr2m\f1;
lagr1 = uin1(2*nnodes+1:end,1);
lamb1 = lagr2forces2( lagr1, c2noder2m, 2, boundaryp1, nnodes );
%
Axz = lamb1;

%%%%
% Compute Rhs :
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

%% Z-1*tr % not so useful
%tr = loading(nbloqZ,nodes,boundary,neumann2);
%uin1 = KZ\tr;
%u1i = uin1(1:2*nnodes,1);
%u1 = keepField( u1i, 4, boundaryp1 );
%% r2m resolution
%f1 = dirichletRhs2( urefb-u1, 4, c2node1, boundaryp1, nnodes );
%uin1 = K1\f1;
%lagr1 = uin1(2*nnodes+1:end,1);
%lamb1 = lagr2forces2( lagr1, c2node1, 2, boundaryp1, nnodes );
%%
%b = -lamb1;
%%%
Res(:,1) = b - Axz;

if precond == 1
    % Solve 1
    f1 = [Res(:,1)/2; zeros(nbloq1d,1)];
    uin1 = K1d\f1;
    u1i = uin1(1:2*nnodes,1);
    u1 = keepField( u1i, 2, boundaryp1 );
%    % Solve 2
%    f2 = [Res(:,1)/2; zeros(nbloq2d,1)];
%    uin2 = K2d\f2;
%    u2i = uin2(1:2*nnodes,1);
%    u2 = keepField( u2i, 2, boundaryp2 );
    %
    Zed(:,1) = u1;%/2-u2/2;
elseif precond == 2
    Zed(index,1) = 1/E*H12(index,index)\Res(index,1);
elseif precond == 3
    Zed(index,1) = 1/E*Mgr(index,index)\Res(index,1);
else
    Zed(:,1) = Res(:,1);
end

d(:,1) = Zed(:,1);

residual(1) = norm(Res( indexa,1));
error(1)    = norm(Itere(indexa) - uref(indexa)) / norm(uref(indexa));
regulari(1) = sqrt( Itere'*regul(Itere, nodes, boundary, 2) );

%%
for iter = 1:niter
    %% Optimal step
    
%    % Solve 1
%    rhs1 = d(:,iter);
%    f1 = dirichletRhs2(rhs1, 2, c2node1, boundaryp1, nnodes);
%    uin1 = K1\f1;
%    lagr1 = uin1(2*nnodes+1:end,1);
%    lamb1 = lagr2forces2( lagr1, c2node1, 2, boundaryp1, nnodes );
%    % Solve 2
%    rhs2 = d(:,iter);
%    f2 = dirichletRhs2(rhs2, 2, c2node2, boundaryp2, nnodes);
%    uin2 = K2\f2;
%    lagr2 = uin2(2*nnodes+1:end,1);
%    lamb2 = lagr2forces2( lagr2, c2node2, 2, boundaryp2, nnodes );
%    %
%    Ad(:,iter) = lamb1-lamb2;
    
    % First resolution
    f1 = dirichletRhs2( d(:,iter), 2, c2nodem2r, boundaryp1, nnodes );
    uin1 = Km2r\f1;
    lagr1 = uin1(2*nnodes+1:end,1);
    lamb1 = lagr2forces2( lagr1, c2nodem2r, 4, boundaryp1, nnodes );
    % Z-1
    uin1 = KZ \ [lamb1 ; zeros(nbloqZ,1)];
    u1i = uin1(1:2*nnodes,1);
    u1 = keepField( u1i, 4, boundaryp1 );
    % Last resolution
    f1 = dirichletRhs2( u1, 4, c2noder2m, boundaryp1, nnodes );
    uin1 = Kr2m\f1;
    lagr1 = uin1(2*nnodes+1:end,1);
    lamb1 = lagr2forces2( lagr1, c2noder2m, 2, boundaryp1, nnodes );
    %
    Ad(:,iter) = lamb1;
    
    den = (d(indexa,iter)'*Ad(indexa,iter));
    d(:,iter) = d(:,iter)/sqrt(den); Ad(:,iter) = Ad(:,iter)/sqrt(den);
    num = Res(indexa,iter)'*d(indexa,iter);
    Itere         = Itere + d(:,iter)*num;
    Res(:,iter+1) = Res(:,iter) - Ad(:,iter)*num;
    
    residual(iter+1) = norm(Res(indexa,iter+1));
    error(iter+1)    = norm(Itere(indexa) - uref(indexa)) / norm(uref(indexa));
    regulari(iter+1) = sqrt( Itere'*regul(Itere, nodes, boundary, 2) );
    
    if precond == 1
        % Solve 1
        f1 = [Res(:,iter+1)/2; zeros(nbloq1d,1)];
        uin1 = K1d\f1;
        u1i = uin1(1:2*nnodes,1);
        u1 = keepField( u1i, 2, boundaryp1 );
%        % Solve 2
%        f2 = [Res(:,iter+1)/2; zeros(nbloq2d,1)];
%        uin2 = K2d\f2;
%        u2i = uin2(1:2*nnodes,1);
%        u2 = keepField( u2i, 2, boundaryp2 );
%        %
        Zed(:,iter+1) = u1;%/2-u2/2;
    elseif precond == 2
        Zed(index,iter+1) = 1/E*H12(index,index)\Res(index,iter+1);
    elseif precond == 3
        Zed(index,iter+1) = 1/E*Mgr(index,index)\Res(index,iter+1);
    else
        Zed(:,iter+1) = Res(:,iter+1);
    end
    
    %% Orthogonalization
    d(:,iter+1) = Zed(:,iter+1);
    
    for jter=1:iter
        betaij = ( Zed(indexa,iter+1)'*Ad(indexa,jter) );%/...
            %( d(indexa,jter)'*Ad(indexa,jter) );
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
             1,1,0;1,2,0;
             2,1,0;2,2,0];
neumann   = [];
[K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
fdir2 = dirichletRhs(Itere, 2, C, boundary);
fdir4 = dirichletRhs(urefb, 4, C, boundary);
usoli = K \ (fdir4 + fdir2);

usol = usoli(1:2*nnodes,1);
fsol = Kinter*usol;

hold on;
plot(usol(2*b2node2-1),'Color','red');
plot(uref(2*b2node2-1),'Color','blue');

figure;
hold on;
plot(fsol(2*b2node2-1),'Color','red');
plot(fref(2*b2node2-1),'Color','blue');

total_error = norm(usol-uref)/norm(uref);

% Compute stress :
sigma = stress(usol,E,nu,nodes,elements,order,1,ntoelem);
plotGMSH({usol,'U_vect';sigma,'stress'}, elements, nodes, 'solution');