%21/03/2016
%Algo KMF avec Orthodir

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 70000;      % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 1;      % N.mm-1 : Loading on the plate
niter   = 20;
br      = 0.0;      % bruit

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

% Then, build the stiffness matrix :
[K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
Kinter = K(1:2*nnodes, 1:2*nnodes);
% The right hand side :
f = loading(nbloq,nodes,boundary,neumann);

% Solve the problem :
uin = K\f;

% Extract displacement and Lagrange multiplicators :
uref = uin(1:2*nnodes,1);
lagr = uin(2*nnodes+1:end,1);
urefb = ( 1 + br*randn(2*nnodes,1) ) .* uref;
lagrb = ( 1 + br*randn(nbloq,1) ) .* lagr;

% Post-pro :
sigma = stress(uref,E,nu,nodes,elements,order,1,ntoelem);
plotGMSH({uref,'U_vect';sigma,'stress'}, elements, nodes, 'reference');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition of the operators

% DN problem
dirichlet1 = [4,1,0;
              4,2,0;
              3,1,0;
              3,2,0];
neumann1   = [1,2,fscalar;
              2,1,fscalar];
neumann0   = [];
[K1,C1,nbloq1] = Krig (nodes,elements,E,nu,order,boundary,dirichlet1);

% ND problem
dirichlet2 = [4,1,0;
              4,2,0;
              1,1,0;
              1,2,0;
              2,1,0;
              2,2,0];
neumann2   = [];% [3,1,lagr1; 3,2,lagr1]
                % is managed by lagr2forces
[K2,C2,nbloq2] = Krig (nodes,elements,E,nu,order,boundary,dirichlet2);

error   = zeros(niter+1,1);
residual = zeros(niter+1,1);
regulari = zeros(niter+1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ORTHODIR for the problem : (1 - D0 o S0) x = DoS(0)
Itere = zeros( 2*nnodes, 1 );
p     = zeros( 2*nnodes, niter+1 );
q     = zeros( 2*nnodes, niter+1 );
Delta = zeros( niter+1, 1 );
Res   = zeros( 2*nnodes, niter+1 );

%% Perform A x0 :
% Solve DN
f1in = loading(nbloq1,nodes,boundary,neumann0);
fdir = dirichletRhs(Itere, 3, C1, boundary);
f1 = f1in + fdir;

uin1 = K1\f1;
u1 = uin1(1:2*nnodes,1);
lagr1 = uin1(2*nnodes+1:end,1);

% Solve ND
fr = lagr2forces( lagr1, C1, 3, boundary );
fr = [ fr; zeros(size(C2,2),1) ]; % Give to fr the right size (because Lagrange)
%fdir1 = dirichletRhs(uref, 1, C2, boundary);
%fdir2 = dirichletRhs(uref, 2, C2, boundary);
f2 = fr;% + fdir1 + fdir2;

uin2 = K2\f2;
atimesItere = Itere - uin2(1:2*nnodes,1);
%%%%

%% Write Rhs :
% Solve DN
f1in = loading(nbloq1,nodes,boundary,neumann1);
fdir = dirichletRhs(zeros( 2*nnodes,1 ), 3, C1, boundary);
f1 = f1in + fdir;

uin1 = K1\f1;
u1 = uin1(1:2*nnodes,1);
lagr1 = uin1(2*nnodes+1:end,1);

% Solve ND
fr = lagr2forces( lagr1, C1, 3, boundary );
fr = [ fr; zeros(size(C2,2),1) ]; % Give to fr the right size
fdir1 = dirichletRhs(urefb, 1, C2, boundary);
fdir2 = dirichletRhs(urefb, 2, C2, boundary);
f2 = fr + assembleDirichlet( [fdir1,fdir2] );

uin2 = K2\f2;
b = uin2(1:2*nnodes,1);
%%%%

%%
Res(:,1) = b - atimesItere;
p(:,1) = Res(:,1);

residual(1) = sqrt( myps( Res(:,1), Res(:,1), Kinter, boundary, C1, nodes ) );
error(1)    = sqrt( myps( Itere - uref, Itere - uref, Kinter, boundary, C1, nodes )...
                 / myps( uref, uref, Kinter, boundary, C1, nodes ) );
regulari(1) = Itere'*regul(Itere, nodes, boundary, 3);

%% Perform Q1 = A P1 :
% Solve DN
f1in = loading(nbloq1,nodes,boundary,neumann0);
fdir = dirichletRhs(p(:,1), 3, C1, boundary);
f1 = f1in + fdir;

uin1 = K1\f1;
u1 = uin1(1:2*nnodes,1);
lagr1 = uin1(2*nnodes+1:end,1);

% Solve ND
fr1 = lagr2forces( lagr1, C1, 3, boundary );
fr = [ fr1; zeros(size(C2,2),1) ]; % Give to fr the right size
%fdir1 = dirichletRhs(uref, 1, C2, boundary); Those guys have nothing to do here : we are in a linear problem
%fdir2 = dirichletRhs(uref, 2, C2, boundary);
f2 = fr;% + fdir1 + fdir2;

uin2 = K2\f2;
q(:,1) = p(:,1) - uin2(1:2*nnodes,1);
%%%%
%%
for iter = 1:niter
    
    Delta(iter,1) = myps( q(:,iter), q(:,iter), Kinter, boundary, C1, nodes ); %q(:,iter)'*q(:,iter);
    gammai        = myps( q(:,iter), Res(:,iter), Kinter, boundary, C1, nodes ); %q(:,iter)'*Res;
    alphai        = gammai/Delta(iter,1);
    
    Itere         = Itere + p(:,iter)*alphai;
    Res(:,iter+1) = Res(:,iter) - q(:,iter)*alphai;
    
    residual(iter+1) = sqrt( myps( Res(:,iter+1), Res(:,iter+1), Kinter, boundary, C1, nodes ) );
    error(iter+1)    = sqrt( myps( Itere - uref, Itere - uref, Kinter, boundary, C1, nodes )...
                     / myps( uref, uref, Kinter, boundary, C1, nodes ) );
    regulari(iter+1) = Itere'*regul(Itere, nodes, boundary, 3);

    %% Perform Ari = A*Res
    % Solve DN
    f1in = loading(nbloq1,nodes,boundary,neumann0);
    fdir = dirichletRhs(Res(:,iter+1), 3, C1, boundary);
    f1 = f1in + fdir;

    uin1 = K1\f1;
    u1 = uin1(1:2*nnodes,1);
    lagr1 = uin1(2*nnodes+1:end,1);
    
    % Solve ND
    fr1 = lagr2forces( lagr1, C1, 3, boundary );
    fr = [ fr1; zeros(size(C2,2),1) ]; % Give to fr the right size
    %fdir1 = dirichletRhs(uref, 1, C2, boundary);
    %fdir2 = dirichletRhs(uref, 2, C2, boundary);
    f2 = fr;% + fdir1 + fdir2;

    uin2 = K2\f2;
    Ari = Res(:,iter+1) - uin2(1:2*nnodes,1);
    
    %% Orthogonalization
    p(:,iter+1) = Res(:,iter+1);
    q(:,iter+1) = Ari;
    
    for jter=1:iter
        phiij  = myps( q(:,jter), Ari, Kinter, boundary, C1, nodes ); %q(:,jter)'*Ari;
        betaij = phiij/Delta(jter,1);
        p(:,iter+1) = p(:,iter+1) - p(:,jter) * betaij;
        q(:,iter+1) = q(:,iter+1) - q(:,jter) * betaij;
        %myps( q(:,iter+1), q(:,jter), Kinter, boundary, C1, nodes )
    end
    
end

hold on
plot(log10(error),'Color','blue')
plot(log10(residual),'Color','red')
figure;
% L-curve
loglog(regulari, residual);
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
usoli = K \ assembleDirichlet( [fdir1,fdir2,fdir3] );
usol = usoli(1:2*nnodes,1);

total_error = norm(usol-uref)/norm(uref);

% Post-pro :
sigma = stress(usol,E,nu,nodes,elements,order,1,ntoelem);
plotGMSH({usol,'U_vect';sigma,'stress'}, elements, nodes, 'solution');

plotGMSH({(usol-uref)/norm(uref),'Error'}, elements, nodes, 'error');