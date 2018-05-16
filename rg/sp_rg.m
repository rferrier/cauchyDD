%21/11/16
%Détection de fissure 1D plane par Cauchy puis écart à la réciprocité

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E          = 210000; % MPa : Young modulus
nu         = 0.3;    % Poisson ratio
fscalar    = 250;    % N.mm-1 : Loading on the plate
mat        = [0, E, nu];
mu         = 1e-2;%1e-5;%5e-3;     % Regularization coef
%mu         = 3;     % Regularization coef
dolcurve   = 0;      % Do a L-curve or not
niter1     = 30;
niter2     = 20;

usefourier = 1;
usepolys   = 0;

nbase = 2; % Number of Fourier basis functions
ordp = 6;  % Number of Polynomial basis functions

useorder = 1; % Order of the FE computation

% Boundary conditions
% first index  : index of the boundary
% second index : 1=x, 2=y
% third        : value
% [0,1,value] marks a dirichlet regularization therm on x
dirichlet  = [1,1,0 ; 1,2,0];
neumann1   = [3,2,fscalar];
neumann2   = [2,1,fscalar ; 7,1,fscalar ; 4,1,-fscalar ; 6,1,-fscalar];

% First, import the mesh
if useorder == 1
   [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_sp/plate_totr2.msh' );
elseif useorder == 2
   error('T6 unimplemented');
end

nnodes = size(nodes,1);

% mapBounds
[ node2b5, b2node5 ]   = mapBound( 5, boundary, nnodes );
[ node2b6, b2node6 ]   = mapBound( 6, boundary, nnodes );
[ node2b7, b2node7 ]   = mapBound( 7, boundary, nnodes );
[ node2b2, b2node2 ]   = mapBound( 2, boundary, nnodes );
[ node2b3, b2node3 ]   = mapBound( 3, boundary, nnodes );
[ node2b4, b2node4 ]   = mapBound( 4, boundary, nnodes );
[ node2b8, b2node8 ]   = mapBound( 8, boundary, nnodes );
[ node2b9, b2node9 ]   = mapBound( 9, boundary, nnodes );
index5 = [ 2*b2node5-1 ; 2*b2node5 ];
index6 = [ 2*b2node6-1 ; 2*b2node6 ];
index7 = [ 2*b2node7-1 ; 2*b2node7 ];

% Then, build the stiffness matrix :
[K,C,nbloq,node2c,c2node] = Krig2 (nodes,elements,mat,order,boundary,dirichlet);
Kinter = K( 1:2*nnodes, 1:2*nnodes );

f1  = loading(nbloq,nodes,boundary,neumann1);
uin = K\f1;
u1 = uin(1:2*nnodes,1);
f1s = Kinter*u1; u1s = u1;

f2  = loading(nbloq,nodes,boundary,neumann2);
uin = K\f2;
u2 = uin(1:2*nnodes,1);
f2s = Kinter*u2; u2s = u2;

%% Cut the mesh in order to get only the top to compute the reaction forces
%line = [ 0 ; 20/6 ; 1 ; 20/6];
%[ nodes_d, elements_d, boundary_d, map_d,...
%  nodes_u, elements_u, boundary_u, map_u, newbound ] =...
%                            cutMesh (nodes, elements, boundary, line);
%[K_u,C_u,nbloq_u,node2c_u,c2node_u] = Krig2 (nodes_u,elements_u,mat,order,boundary_u,dirichlet);
%nnu = size(nodes_u,1);
%Kinter_u = K_u( 1:2*nnu, 1:2*nnu );
%
%[ b1to2, b2to1 ] = superNodes( nodes_u, nodes, 1e-6 );
%
%u1u = zeros(2*size(nodes_u,1),1); u2u = zeros(2*size(nodes_u,1),1);
%u1u([(1:2:2*nnu-1)';(2:2:2*nnu)']) = u1s([2*b1to2-1;2*b1to2]); 
%u2u([(1:2:2*nnu-1)';(2:2:2*nnu)']) = u2s([2*b1to2-1;2*b1to2]);
%
%%u1u([2*map_u-1;2*map_u]) = u1s([(1:2:2*nnodes-1)';(2:2:2*nnodes)']); % Why doesn't that shit work ?
%%u2u([2*map_u-1;2*map_u]) = u2s([(1:2:2*nnodes-1)';(2:2:2*nnodes)']);
%
%f1u = Kinter_u*u1u; f2u = Kinter_u*u2u;
%f1su = zeros(2*nnodes); f2su = zeros(2*nnodes);
%f1su([2*b1to2-1;2*b1to2]) = f1u([(1:2:2*nnu-1)';(2:2:2*nnu)']);
%f2su([2*b1to2-1;2*b1to2]) = f2u([(1:2:2*nnu-1)';(2:2:2*nnu)']);

%f1su([1:2:2*nnodes-1;2:2:2*nnodes]) = f1u([2*map_u-1;2*map_u]);
%f2su([1:2:2*nnodes-1;2:2:2*nnodes]) = f2u([2*map_u-1;2*map_u]);

%figure; hold on;
%ret  = patch('Faces',elements_u(:,1:3),'Vertices',nodes_u,'FaceAlpha',0);
%axis('equal');

%% Display both fields
ui = reshape(u1,2,[])';  ux = ui(:,1);  uy = ui(:,2);
sigma = stress(u1,E,nu,nodes,elements,order,1,ntoelem);
plotGMSH({ux,'U_x';uy,'U_y';u1,'U_vect';sigma,'stress'}, elements, nodes, 'reference1');
%
ui = reshape(u2,2,[])';  ux = ui(:,1);  uy = ui(:,2);
sigma = stress(u2,E,nu,nodes,elements,order,1,ntoelem);
plotGMSH({ux,'U_x';uy,'U_y';u2,'U_vect';sigma,'stress'}, elements, nodes, 'reference2');

uref1 = u1; uref2 = u2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First step : solve the Cauchy problem

% Import the data
if useorder == 1
   [ nodesd,elementsd,ntoelemd,boundaryd,order] = readmesh( 'meshes/rg_sp/plate_downr2.msh' );
elseif useorder == 2
   error('T6 unimplemented');
end
nnodesd = size(nodesd,1);

[ b1to2, b2to1 ] = superNodes( nodes, nodesd, 1e-3 );

% Reference reactions on 5
%f1sus = zeros(nnodesd,1); f2sus = zeros(nnodesd,1);
%f1sus([(1:2:2*nnu-1)';(2:2:2*nnu)']) = f1su([2*b1to2-1;2*b1to2]); 
%f2sus([(1:2:2*nnu-1)';(2:2:2*nnu)']) = f2su([2*b1to2-1;2*b1to2]);

[ node2b5d, b2node5d ]   = mapBound( 5, boundaryd, nnodesd );
[ node2b6d, b2node6d ]   = mapBound( 6, boundaryd, nnodesd );
[ node2b7d, b2node7d ]   = mapBound( 7, boundaryd, nnodesd );
index5d = [2*b2node5d-1;2*b2node5d];
index6d = [2*b2node6d-1;2*b2node6d];
index7d = [2*b2node7d-1;2*b2node7d];

uref1 = zeros(2*nnodesd,1);
%uref1(index5d) = u1(index5); uref1(index6d) = u1(index6); uref1(index7d) = u1(index7);
uref2 = zeros(nnodesd,1);
%uref2(index5d) = u2(index5);uref2(index6d) = u2(index6); uref2(index7d) = u2(index7);
uref1([1:2:2*nnodesd-1,2:2:2*nnodesd]) = u1( [2*b2to1-1;2*b2to1] );
uref2([1:2:2*nnodesd-1,2:2:2*nnodesd]) = u2( [2*b2to1-1;2*b2to1] );

% DEBUG
%plotGMSH({uref1,'U_vect'}, elementsd, nodesd, 'uref');

%% Define the stiffness
% First problem
dirichlet1d = [1,1,0;1,2,0;
               7,1,0;7,2,0;
               6,1,0;6,2,0];
[K1d,C1d,nbloq1d,node2c1d,c2node1d] = Krig (nodesd,elementsd,E,nu,order,boundaryd,dirichlet1d);
Kinter1d = K1d(1:2*nnodesd, 1:2*nnodesd);

% Second problem
dirichlet2d = [1,1,0;1,2,0];
[K2d,C2d,nbloq2d,node2c2d,c2node2d] = Krig (nodesd,elementsd,E,nu,order,boundaryd,dirichlet2d);

% The loading (it's 0 for the first problem)
neumann1d   = [];
neumann2d   = [7,1,fscalar ; 6,1,-fscalar];
fref1 = loading(nbloq2d,nodesd,boundaryd,neumann1d); % fref1 refers to the first load case
fref2 = loading(nbloq2d,nodesd,boundaryd,neumann2d); % nbloq2d refers to the second problem

% The refercence solution
fsr1 = Kinter1d*uref1 - fref1(1:2*nnodesd) ;
fsr2 = Kinter1d*uref2 - fref2(1:2*nnodesd) ;
fsr1(~index5d) = 0; fsr2(~index5d) = 0;
fsr1(index7d) = 0;  fsr2(index7d) = 0;  % Because of the corners
fsr1(index6d) = 0;  fsr2(index6d) = 0;
%fsr1(~[2*index5d-1,2*index5d]) = 0; fsr2(~[2*index5d-1,2*index5d]) = 0;
%fsr1([2*index7d-1,2*index7d]) = 0;  fsr2([2*index7d-1,2*index7d]) = 0;  % Because of the corners
%fsr1([2*index6d-1,2*index6d]) = 0;  fsr2([2*index6d-1,2*index6d]) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ORTHODIR for the first problem : (S1-S2) x = b
Itere = zeros( 2*nnodesd, 1 );
p     = zeros( 2*nnodesd, niter1+1 );
q     = zeros( 2*nnodesd, niter1+1 );
Delta = zeros( niter1+1, 1 );
Res   = zeros( 2*nnodesd, niter1+1 );
Zed   = zeros( 2*nnodesd, niter1+1 );

%% Perform A x0 :
% Solve 1
rhs1 = [ Itere ; zeros(nbloq1d,1) ];
uin1 = K1d\rhs1;
u1i = uin1(1:2*nnodesd,1);
u1 = keepField( u1i, 5, boundaryd );
% Solve 2
rhs2 = [ Itere ; zeros(nbloq2d,1) ];
uin2 = K2d\rhs2;
u2i = uin2(1:2*nnodesd,1);
u2 = keepField( u2i, 5, boundaryd );
%
Axz = u1-u2;
%%%%
%% Compute Rhs :
% Solve 1
f11 = dirichletRhs2(uref1, 6, c2node1d, boundaryd, nnodesd);
f12 = dirichletRhs2(uref1, 7, c2node1d, boundaryd, nnodesd);
f1 = f11+f12;%assembleDirichlet( [f11, f12 ] );
uin1 = K1d\f1;
u1i = uin1(1:2*nnodesd,1);
u1 = keepField( u1i, 5, boundaryd );
% Solve 2 (zero)
uin2 = K2d\fref1;
u2i = uin2(1:2*nnodesd,1);
u2 = keepField( u2i, 5, boundaryd );
%
b = -u1+u2;
%%%%
%%
Res(:,1) = b - Axz;

Zed(:,1) = Res(:,1); % No precond
p(:,1) = Zed(:,1);

residual(1) = sqrt( Res(:,1)'*Res(:,1) );
error(1)    = sqrt( (Itere(index5d) - fsr1(index5d))'*...
                       (Itere(index5d) - fsr1(index5d) ) / ...
                       (fsr1(index5d)'*fsr1(index5d)) );
regulari(1) = sqrt( Itere'*regul(Itere, nodes, boundaryd, 5) );

%% Perform Q1 = A P1 :
% Solve 1
rhs1 = [ p(:,1) ; zeros(nbloq1d,1) ];
uin1 = K1d\rhs1;
u1i = uin1(1:2*nnodesd,1);
u1 = keepField( u1i, 5, boundaryd );
% Solve 2
rhs2 = [ p(:,1) ; zeros(nbloq2d,1) ];
uin2 = K2d\rhs2;
u2i = uin2(1:2*nnodesd,1);
u2 = keepField( u2i, 5, boundaryd );
%
q(:,1) = u1-u2;
%%%%
%%
for iter = 1:niter1
    
    Delta(iter,1) = q(:,iter)'*q(:,iter);
    gammai        = q(:,iter)'*Res(:,iter);
    alphai        = gammai/Delta(iter,1);
    
    Itere         = Itere + p(:,iter)*alphai;
    Res(:,iter+1) = Res(:,iter) - q(:,iter)*alphai;
    
    Zed(:,iter+1) = Res(:,iter+1);

    residual(iter+1) = sqrt( (Res(:,iter+1))'*Res(:,iter+1) ) ;
    error(iter+1)    = sqrt( (Itere(index5d) - fsr1(index5d))'*...
                       (Itere(index5d) - fsr1(index5d) ) / ...
                       (fsr1(index5d)'*fsr1(index5d)) );
    regulari(iter+1) = sqrt( Itere'*regul(Itere, nodes, boundaryd, 5) );

    %% Perform Ari = A*Res
    % Solve 1
    rhs1 = [ Zed(:,iter+1) ; zeros(nbloq1d,1) ];
    uin1 = K1d\rhs1;
    u1i = uin1(1:2*nnodesd,1);
    u1 = keepField( u1i, 5, boundaryd );
    % Solve 2
    rhs2 = [ Zed(:,iter+1) ; zeros(nbloq2d,1) ];
    uin2 = K2d\rhs2;
    u2i = uin2(1:2*nnodesd,1);
    u2 = keepField( u2i, 5, boundaryd );
    %
    Ari = u1-u2;

    %% Orthogonalization
    p(:,iter+1) = Zed(:,iter+1);
    q(:,iter+1) = Ari;
    
    for jter=1:iter
        phiij  = q(:,jter)'*Ari;
        betaij = phiij/Delta(jter,1);
        p(:,iter+1) = p(:,iter+1) - p(:,jter) * betaij;
        q(:,iter+1) = q(:,iter+1) - q(:,jter) * betaij;
    end 
end

%% Plot convergence
%figure; hold on;
%plot(error,'Color','blue')
%plot(residual,'Color','red')
%legend('error','residual')

% Last resolution
f3 = [ Itere ; zeros(nbloq2d,1) ]; fsol1 = Itere;
usoli = K2d \ ( fref1 + f3 );
usol1 = usoli(1:2*nnodesd,1);

figure; hold on;
plot(usoli(2*b2node5d));
plot(uref1(2*b2node5d),'Color','red');
legend('solution','reference')

figure; hold on;
plot(Itere(2*b2node5d));
plot(fsr1(2*b2node5d),'Color','red');
legend('solution','reference')

% Compute stress :
sigma1 = stress(usol1,E,nu,nodesd,elementsd,order,1,ntoelemd);
plotGMSH({usol1,'U_vect';sigma1,'stress'}, elementsd, nodesd, 'solution_SP_1');

Itere1 = Itere;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ORTHODIR for the second problem : (S1-S2) x = b
Itere = zeros( 2*nnodesd, 1 );
p     = zeros( 2*nnodesd, niter2+1 );
q     = zeros( 2*nnodesd, niter2+1 );
Delta = zeros( niter2+1, 1 );
Res   = zeros( 2*nnodesd, niter2+1 );
Zed   = zeros( 2*nnodesd, niter2+1 );

%% Perform A x0 :
% Solve 1
rhs1 = [ Itere ; zeros(nbloq1d,1) ];
uin1 = K1d\rhs1;
u1i = uin1(1:2*nnodesd,1);
u1 = keepField( u1i, 5, boundaryd );
% Solve 2
rhs2 = [ Itere ; zeros(nbloq2d,1) ];
uin2 = K2d\rhs2;
u2i = uin2(1:2*nnodesd,1);
u2 = keepField( u2i, 5, boundaryd );
%
Axz = u1-u2;
%%%%
%% Compute Rhs :
% Solve 1
f11 = dirichletRhs2(uref2, 6, c2node1d, boundaryd, nnodesd);
f12 = dirichletRhs2(uref2, 7, c2node1d, boundaryd, nnodesd);
f1 = f11+f12;%assembleDirichlet( [f11, f12 ] );
uin1 = K1d\f1;
u1i = uin1(1:2*nnodesd,1);
u1 = keepField( u1i, 5, boundaryd );
% Solve 2 (NAUTE zero)
uin2 = K2d\fref2;
u2i = uin2(1:2*nnodesd,1);
u2 = keepField( u2i, 5, boundaryd );
%
b = -u1+u2;
%%%%
%%
Res(:,1) = b - Axz;

Zed(:,1) = Res(:,1); % No precond
p(:,1) = Zed(:,1);

residual(1) = sqrt( Res(:,1)'*Res(:,1) );
error(1)    = sqrt( (Itere(index5d) - fsr2(index5d))'*...
                       (Itere(index5d) - fsr2(index5d) ) / ...
                       (fsr2(index5d)'*fsr2(index5d)) );
regulari(1) = sqrt( Itere'*regul(Itere, nodes, boundaryd, 5) );

%% Perform Q1 = A P1 :
% Solve 1
rhs1 = [ p(:,1) ; zeros(nbloq1d,1) ];
uin1 = K1d\rhs1;
u1i = uin1(1:2*nnodesd,1);
u1 = keepField( u1i, 5, boundaryd );
% Solve 2
rhs2 = [ p(:,1) ; zeros(nbloq2d,1) ];
uin2 = K2d\rhs2;
u2i = uin2(1:2*nnodesd,1);
u2 = keepField( u2i, 5, boundaryd );
%
q(:,1) = u1-u2;
%%%%
%%
for iter = 1:niter2
    
    Delta(iter,1) = q(:,iter)'*q(:,iter);
    gammai        = q(:,iter)'*Res(:,iter);
    alphai        = gammai/Delta(iter,1);
    
    Itere         = Itere + p(:,iter)*alphai;
    Res(:,iter+1) = Res(:,iter) - q(:,iter)*alphai;
    
    Zed(:,iter+1) = Res(:,iter+1);

    residual(iter+1) = sqrt( (Res(:,iter+1))'*Res(:,iter+1) ) ;
    error(iter+1)    = sqrt( (Itere(index5d) - fsr2(index5d))'*...
                       (Itere(index5d) - fsr2(index5d) ) / ...
                       (fsr2(index5d)'*fsr2(index5d)) );
    regulari(iter+1) = sqrt( Itere'*regul(Itere, nodes, boundaryd, 5) );

    %% Perform Ari = A*Res
    % Solve 1
    rhs1 = [ Zed(:,iter+1) ; zeros(nbloq1d,1) ];
    uin1 = K1d\rhs1;
    u1i = uin1(1:2*nnodesd,1);
    u1 = keepField( u1i, 5, boundaryd );
    % Solve 2
    rhs2 = [ Zed(:,iter+1) ; zeros(nbloq2d,1) ];
    uin2 = K2d\rhs2;
    u2i = uin2(1:2*nnodesd,1);
    u2 = keepField( u2i, 5, boundaryd );
    %
    Ari = u1-u2;

    %% Orthogonalization
    p(:,iter+1) = Zed(:,iter+1);
    q(:,iter+1) = Ari;
    
    for jter=1:iter
        phiij  = q(:,jter)'*Ari;
        betaij = phiij/Delta(jter,1);
        p(:,iter+1) = p(:,iter+1) - p(:,jter) * betaij;
        q(:,iter+1) = q(:,iter+1) - q(:,jter) * betaij;
    end 
end

%% Plot convergence
%figure; hold on;
%plot(error,'Color','blue')
%plot(residual,'Color','red')
%legend('error','residual')

% Last resolution
f3 = [ Itere ; zeros(nbloq2d,1) ]; fsol2 = Itere;
usoli = K2d \ ( fref2 + f3 );
usol2 = usoli(1:2*nnodesd,1);
% Compute stress :
sigma2 = stress(usol2,E,nu,nodesd,elementsd,order,1,ntoelemd);
plotGMSH({usol2,'U_vect';sigma2,'stress'}, elementsd, nodesd, 'solution_SP_2');

figure; hold on;
plot(usoli(2*b2node5d));
plot(uref2(2*b2node5d),'Color','red');
legend('solution','reference')

figure; hold on;
plot(Itere(2*b2node5d));
plot(fsr2(2*b2node5d),'Color','red');
legend('solution','reference')

% Build u1, f1, u2 and f2 on the total domain.
u1 = u1s; u2 = u2s; f1 = f1s; f2 = f2s;

u1(index5) = usol1(index5d); u2(index5) = usol2(index5d);
f1(index5) = -fsol1(index5d); f2(index5) = -fsol2(index5d);

% Node 5 and 6 : add exterior loading
index56 = [9,10,11,12];
f1(index56) = f1(index56) + f1s(index56); % Because of the corner
f2(index56) = f2(index56) + f2s(index56); %

% Compute some errors due to the Cauchy solving
uerror1 = ( u1(index5)-u1s(index5) )'*( u1(index5)-u1s(index5) ) /...
           ( u1s(index5)'*u1s(index5) );
uerror2 = ( u2(index5)-u2s(index5) )'*( u2(index5)-u2s(index5) ) /...
           ( u2s(index5)'*u2s(index5) );
ferror1 = ( f1(index5)-f1s(index5) )'*( f1(index5)-f1s(index5) ) /...
           ( f1s(index5)'*f1s(index5) );
ferror2 = ( f2(index5)-f2s(index5) )'*( f2(index5)-f2s(index5) ) /...
           ( f2s(index5)'*f2s(index5) );

%plotGMSH({u1,'U_vect';f1,'F_vect'}, elements, nodes, 'uref');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now, solve the RG stuff
% Import the uncracked domain
if useorder == 1
   [ nodes2,elements2,ntoelem2,boundary2,order2] = readmesh( 'meshes/rg_sp/plate_upnr.msh' );
elseif useorder == 2
   error('How did you get so far ?');
end

nnodes2 = size(nodes2,1);
[K2,C2,nbloq2,node2c2,c2node2] = Krig2 (nodes2,elements2,mat,order,boundary2,dirichlet);
Kinter2 = K2( 1:2*nnodes2, 1:2*nnodes2 );

% Mapbounds
[ node2b52, b2node52 ] = mapBound( 5, boundary2, nnodes2 );
[ node2b22, b2node22 ] = mapBound( 2, boundary2, nnodes2 );
[ node2b32, b2node32 ] = mapBound( 3, boundary2, nnodes2 );
[ node2b42, b2node42 ] = mapBound( 4, boundary2, nnodes2 );

indexbound  = [2*b2node5-1 ; 2*b2node5 ; 2*b2node2-1 ; 2*b2node2 ;...
               2*b2node3-1 ; 2*b2node3 ; 2*b2node4-1 ; 2*b2node4];

 % Sort the nodes (in case order > 1)
[~,order8]  = sort( nodes( b2node8, 1 ) );
[~,order9]  = sort( nodes( b2node9, 1 ) );

icrack8x    = [ 2*b2node8(order8)-1 ];
icrack9x    = [ 2*b2node9(order9)-1 ];
icrack8y    = [ 2*b2node8(order8) ];
icrack9y    = [ 2*b2node9(order9) ];

indexbound2 = [2*b2node52-1 ; 2*b2node52 ; 2*b2node22-1 ; 2*b2node22 ;...
               2*b2node32-1 ; 2*b2node32 ; 2*b2node42-1 ; 2*b2node42];
% Pass f and u on the uncracked mesh
ur1 = zeros( 2*nnodes2, 1 ); fr1 = zeros( 2*nnodes2, 1 );
ur1(indexbound2) = u1(indexbound);
fr1(indexbound2) = f1(indexbound);
% Same for the second one
ur2 = zeros( 2*nnodes2, 1 ); fr2 = zeros( 2*nnodes2, 1 );
ur2(indexbound2) = u2(indexbound);
fr2(indexbound2) = f2(indexbound);

%plotGMSH({ur2,'U_vect';fr2,'F_vect'}, elements2, nodes2, 'uref');

% Debug : compute the displacement gap
Mc = bMass_mat (nodes, boundary, 5);
onfi = ones(2*nnodes,1);
udepx = u1; udepx(icrack8x) = u1(icrack9x);
intx = onfi'*Mc*(u1-udepx);
udepy = u1; udepy(icrack8y) = u1(icrack9y);
inty = onfi'*Mc*(u1-udepy);
intt = sqrt(intx^2+inty^2);
%
% Debug : compute n
x7 = nodes(7,1); y7 = nodes(7,2); x8 = nodes(8,1); y8 = nodes(8,2);
n1 = 1; n2 = -n1*(x7-x8)/(y7-y8);
non = sqrt(n1^2+n2^2);
n1 = n1/non; n2 = n2/non;
CteR = x8*n1+y8*n2;

%%% Preliminary stuff : find the volumic elements corresponding to the boundaries
%nboun2 = size(boundary2,1); nelem2 = size(elements2,1);
%boun2vol2 = zeros( nboun2, 1 ); extnorm2 = zeros( nboun2, 2 );
%for i=1:nboun2
%   % Volumic element
%   no1 = boundary2(i,2); no2 = boundary2(i,3); % only with 2 nodes even if order > 1
%   cand1 = rem( find(elements2==no1),nelem2 ); % find gives line + column*size
%   cand2 = rem( find(elements2==no2),nelem2 );
%   boun2vol2(i) = intersect(cand1, cand2); % If everything went well, there is only one
%   
%   % Exterior normal
%   elt = boun2vol2(i); no3 = setdiff( elements2( elt, 1:3 ), [no1,no2]);
%   x1 = nodes2(no1,1); y1 = nodes2(no1,2);
%   x2 = nodes2(no2,1); y2 = nodes2(no2,2);
%   x3 = nodes2(no3,1); y3 = nodes2(no3,2);
%   extnorm2(i,:) = [ y1-y2 , -(x1-x2) ];
%   extnorm2(i,:) = extnorm2(i,:)/norm(extnorm2(i,:));
%   if (x3-x2)*(y1-y2) - (y3-y2)*(x1-x2) > 0 % Check that the normal is exterior
%      extnorm2(i,:) = -extnorm2(i,:);
%   end
%end
%
%nboun1 = size(boundary,1); nelem1 = size(elements,1);
%boun2vol1 = zeros( nboun1, 1 ); extnorm1 = zeros( nboun1, 2 );
%sigr1 = zeros( nboun1, 3*order ); sigr2 = zeros( nboun1, 3*order );
%urr1  = zeros( nboun1, 2+2*order ); urr2 = zeros( nboun1, 2+2*order );
%for i=1:nboun1
%   % Volumic element
%   no1 = boundary(i,2); no2 = boundary(i,3); % only with 2 nodes even if order > 1
%   cand1 = rem( find(elements==no1),nelem1 ); % find gives line + column*size
%   cand2 = rem( find(elements==no2),nelem1 );
%   boun2vol1(i) = intersect(cand1, cand2); % If everything went well, there is only one
%   
%   % Exterior normal
%   elt = boun2vol1(i); no3 = setdiff( elements( elt, 1:3 ), [no1,no2]);
%   x1 = nodes(no1,1); y1 = nodes(no1,2);
%   x2 = nodes(no2,1); y2 = nodes(no2,2);
%   x3 = nodes(no3,1); y3 = nodes(no3,2);
%   extnorm1(i,:) = [ y1-y2 , -(x1-x2) ];
%   extnorm1(i,:) = extnorm1(i,:)/norm(extnorm1(i,:));
%   if (x3-x2)*(y1-y2) - (y3-y2)*(x1-x2) > 0 % Check that the normal is exterior
%      extnorm1(i,:) = -extnorm1(i,:);
%   end
%   
%   % ur
%   urr1(i,1:4) = u1( [2*no1-1,2*no1,2*no2-1,2*no2] );
%   urr2(i,1:4) = u2( [2*no1-1,2*no1,2*no2-1,2*no2] );
%   if order == 2
%      no4 = boundary(i,4);
%      urr1(i,5:6) = u1( [2*no4-1,2*no4] );
%      urr2(i,5:6) = u2( [2*no4-1,2*no4] );
%   end
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First determine the crack's line.

% Compute and apply v fields (for plane constraint)
v1 = zeros(2*nnodes2, 1);
v2 = zeros(2*nnodes2, 1);
v3 = zeros(2*nnodes2, 1);
for i=1:nnodes2
   x = nodes2(i,1);
   y = nodes2(i,2);
   v1(2*i-1) = 1/E*x;
   v2(2*i-1) = -nu/E*x;
   v3(2*i-1) = (1+nu)/(2*E)*y;
   v1(2*i)   = -nu/E*y;
   v2(2*i)   = 1/E*y;
   v3(2*i)   = (1+nu)/(2*E)*x;
end
%
f1 = Kinter2*v1;
f2 = Kinter2*v2;
f3 = Kinter2*v3;

% Clean redondant stuff in indexbound2
i = 1;
while i <= size(indexbound2,1)
   if find( indexbound2(i) == indexbound2(1:i-1) )
      indexbound2(i) = [];
      i = i-1;
   end
   i = i+1;
end

R11 = (fr1(indexbound2)'*v1(indexbound2) - f1(indexbound2)'*ur1(indexbound2));
R21 = (fr1(indexbound2)'*v2(indexbound2) - f2(indexbound2)'*ur1(indexbound2));
R31 = (fr1(indexbound2)'*v3(indexbound2) - f3(indexbound2)'*ur1(indexbound2));
%[R11s, R31s ; R31s, R21s]

R12 = (fr2(indexbound2)'*v1(indexbound2) - f1(indexbound2)'*ur2(indexbound2));
R22 = (fr2(indexbound2)'*v2(indexbound2) - f2(indexbound2)'*ur2(indexbound2));
R32 = (fr2(indexbound2)'*v3(indexbound2) - f3(indexbound2)'*ur2(indexbound2));

% Normalize R1
%R11 = R1d; R21 = R2d; R31 = R3d;%DEBUG
normR1 = sqrt( 2*(R11^2+R21^2+2*R31^2) - (R11+R21)^2 );
R1b1 = R11/normR1; R2b1 = R21/normR1; R3b1 = R31/normR1;
%[R1b1, R3b1 ; R3b1, R2b1]

lam11 = (1+R1b1+R2b1)/2; lam21 = -(1-R1b1-R2b1)/2; R1 = [R11,R31;R31,R21];
[phi1,Lambda1] = eig( [R1b1,R3b1;R3b1,R2b1] );
dir11 = [ sqrt( abs(Lambda1(1,1)) ) ; sqrt( abs(Lambda1(2,2)) ) ];
dir21 = [ sign(Lambda1(1,1)) * sqrt( abs(Lambda1(1,1)) ) ;...
          sign(Lambda1(2,2)) * sqrt( abs(Lambda1(2,2)) ) ];

dir31 = [ sign(Lambda1(1,1)) * sqrt( abs(Lambda1(1,1)) ) ;...
          sqrt( abs(Lambda1(2,2)) ) ];
dir41 = [ sqrt( abs(Lambda1(1,1)) ) ;...
          sign(Lambda1(2,2)) * sqrt( abs(Lambda1(2,2)) ) ];

dir11 = phi1*dir11; dir11 = dir11/norm(dir11);
dir21 = phi1*dir21; dir21 = dir21/norm(dir21); % Normal candidates (1)
dir31 = phi1*dir31; dir31 = dir31/norm(dir31);
dir41 = phi1*dir41; dir41 = dir41/norm(dir41);

% Normalize R2
%R12 = R1d2; R22 = R2d2; R32 = R3d2;%DEBUG
normR2 = sqrt( 2*(R12^2+R22^2+2*R32^2) - (R12+R22)^2 );
R1b2 = R12/normR2; R2b2 = R22/normR2; R3b2 = R32/normR2;
%[R12, R32 ; R32, R22]

R2 = [R12,R32;R32,R22];
[phi2,Lambda2] = eig( [R1b2,R3b2;R3b2,R2b2] );
dir12 = [ sqrt( abs(Lambda2(1,1)) ) ; sqrt( abs(Lambda2(2,2)) ) ];
dir22 = [ sign(Lambda2(1,1)) * sqrt( abs(Lambda2(1,1)) ) ;...
          sign(Lambda2(2,2)) * sqrt( abs(Lambda2(2,2)) ) ];

dir12 = phi2*dir12; dir12 = dir12/norm(dir12);
dir22 = phi2*dir22; dir22 = dir22/norm(dir22); % Normal candidates (2)
dir32 = [ sign(Lambda2(1,1)) * sqrt( abs(Lambda2(1,1)) ) ;...
          sqrt( abs(Lambda2(2,2)) ) ];
dir42 = [ sqrt( abs(Lambda2(1,1)) ) ;...
          sign(Lambda2(2,2)) * sqrt( abs(Lambda2(2,2)) ) ];

% Find the real normal (closest candidates)
dist11 = norm(dir11-dir12); dist12 = norm(dir11-dir22);
dist21 = norm(dir21-dir12); dist22 = norm(dir21-dir22);
dist11m = norm(dir11+dir12); dist12m = norm(dir11+dir22);
dist21m = norm(dir21+dir12); dist22m = norm(dir21+dir22);
mindist = min([dist11,dist12,dist21,dist22,dist11m,dist12m,dist21m,dist22m]);

if dist11 == mindist
   normal = .5*(dir11+dir12);
elseif dist22 == mindist
   normal = .5*(dir21+dir22);
elseif dist12 == mindist
   normal = .5*(dir11+dir22);
elseif dist21 == mindist
   normal = .5*(dir21+dir12);
elseif dist11m == mindist
   normal = .5*(dir11-dir12);
elseif dist22m == mindist
   normal = .5*(dir21-dir22);
elseif dist12m == mindist
   normal = .5*(dir11-dir22);
elseif dist21m == mindist
   normal = .5*(dir21-dir12);
end

% Build the base-change matrix : [x;y] = Q*[X;Y], [X;Y] = Q'*[x;y]
Q = [ normal(2), normal(1) ; - normal(1), normal(2) ];

% Plot the normal
figure
hold on;
ret = patch('Faces',elements2(:,1:3),'Vertices',nodes2,'FaceAlpha',0);
x6 = nodes(6,1); y6 = nodes(6,2); x5 = nodes(5,1); y5 = nodes(5,2); 
xc = .5*(max(nodes(:,1)) + min(nodes(:,1)));
yc = .5*(max(nodes(:,2)) + min(nodes(:,2)));

x2 = xc + 3*normal(1); y2 = yc + 3*normal(2);
xn = xc + 3*n1; yn = yc + 3*n2;
plot( [xc,x2], [yc,y2] ,'Color', 'red', 'LineWidth',3);
plot( [xc,xn], [yc,yn] ,'Color', 'cyan', 'LineWidth',2);
plot( [x7,x8], [y7,y8] ,'Color', 'magenta', 'LineWidth',3);
axis('equal');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine the constant

% First, find the minimal point (in order to have Cte < 0)
norep = Q'*nodes'; K = min(norep(2,:));

vt = zeros(2*nnodes2, 1);
Xs = zeros(nnodes2, 1);
Ys = zeros(nnodes2, 1);
for i=1:nnodes2
   x = nodes2(i,1);
   y = nodes2(i,2);
   % Change base (for coordiantes)
   ixigrec = Q'*[x;y]; X = ixigrec(1); Y = ixigrec(2)-K;
   Xs(i) = X; Ys(i) = Y;

   vloc = [ -X^2/(2*E) + (2+nu)*Y^2/(2*E) ; nu*X*Y/E ];
   % Change base (for vector), in the other direction
   vxy = Q*vloc; vt(2*i-1) = vxy(1); vt(2*i) = vxy(2);
end

% Norm of [[ut]] (case 1)
normT  = sqrt( abs( 2*(R11^2+R21^2+2*R31^2) - 2*(R11+R21)^2 ) );
normT2 = sqrt( abs( 2*(R12^2+R22^2+2*R32^2) - 2*(R12+R22)^2 ) );

ft = Kinter2*vt;
Rt = (fr1(indexbound2)'*vt(indexbound2) - ft(indexbound2)'*ur1(indexbound2));
Rt2 = (fr2(indexbound2)'*vt(indexbound2) - ft(indexbound2)'*ur2(indexbound2));

Cte  = min( Rt/normT, -Rt/normT) - K;       % Select the negative one
Cte2 = min( Rt2/normT2, -Rt2/normT2 ) - K;  %  /!\ The sign depends on the test case

% Plot the crack, and its estimated lines (there are Y +/- Cte)
figure
hold on;
ret = patch('Faces',elements2(:,1:3),'Vertices',nodes2,'FaceAlpha',0);
x7 = nodes(7,1); y6 = nodes(7,2); x8 = nodes(8,1); y8 = nodes(8,2); 

Vp1 = [20;-Cte]; Vp2 = [-20;-Cte];
Vm1 = [20;-Cte2]; Vm2 = [-20;-Cte2];
vp1 = Q*Vp1; vm1 = Q*Vm1; vp2 = Q*Vp2; vm2 = Q*Vm2;

plot( [vp1(1), vp2(1)], [vp1(2), vp2(2)] ,'Color', 'red', 'LineWidth',1);
plot( [vm1(1), vm2(1)], [vm1(2), vm2(2)] ,'Color', 'green', 'LineWidth',1);
plot( [x8,x7], [y8,y7] ,'Color', 'magenta', 'LineWidth',2);
axis('equal');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% And now, the crack itself

% Compute L (provisionnal formula assuming that something<Cte<0
tangent = [normal(2) ; -normal(1)]; b = max(nodes2(:,2)) - min(nodes2(:,2));
a = tangent(1)/tangent(2)*b; L = sqrt(a^2+b^2);

% Convenient values for plot
udepx  = u2(icrack8x)-u2(icrack9x);
udepy  = u2(icrack8y)-u2(icrack9y);
Xx     = nodes(b2node8,1); Yy = nodes(b2node8,2);
ubase  = Q'*[udepx,udepy]'; ubase = ubase';
XY     = Q'*[Xx,Yy]'; XY = XY';
[newX, orderX] = sort(XY(:,1));          % In case order > 1, need to sort
xy5    = Q'*[ nodes(5,1) ; nodes(5,2) ];
offset = Q(2,2)/Q(2,1)*CteR + xy5(1);   % /!\ Only in case rectangle

left  = offset : (-offset+newX(1))/(size(newX,1)) : newX(1) ;
right = newX(end) : (offset + L - newX(end))/(size(newX,1)) : offset + L ;
newXo = newX;
newX  = [left(1:end-1)';
         newX;
         right(2:end)']; % Add a few points

if usefourier == 1
   Rp     = zeros(nbase+1,1);
   Rm     = zeros(nbase+1,1);
   lambda = zeros(nbase+1,1);
   fourn  = zeros(nbase+1,1);
   fournm = zeros(nbase+1,1);
   akan   = zeros(nbase,1);
   bkan   = zeros(nbase,1);
   
   for kp=2:nbase+1
      k = kp-1;
      vp = zeros(2*nnodes2,1);
      vm = zeros(2*nnodes2,1);
      lambda(kp) = 2*k*pi/L;
      for sx = [1,-1]  % sx=-1 is not really used, but Debug stuff
         lambdae = sx*lambda(kp);
         for i=1:nnodes2
            x = nodes2(i,1);
            y = nodes2(i,2);
            % Change base (for coordinates)
            ixigrec = Q'*[x;y]; X = ixigrec(1); Y = ixigrec(2)+CteR;
            Xs(i) = X; Ys(i) = Y;
            
            v1 = -I*lambda(kp)*exp(-I*lambdae*X)* ...
                               ( exp(lambda(kp)*Y)+exp(-lambda(kp)*Y) );
            v2 = lambda(kp)*exp(-I*lambdae*X)* ...
                               ( exp(lambda(kp)*Y)-exp(-lambda(kp)*Y) );
            vloc = [ v1 ; v2 ];
            % Change base (for vector), in the other direction
            vxy = Q*vloc; vp(2*i-1) = vxy(1); vp(2*i) = vxy(2);
         end
         fp = Kinter2*vp;
         
         % Fourier coefficient
         if sx == 1
            %Rp(kp) = (fr1(indexbound2)'*vp(indexbound2) - fp(indexbound2)'*ur1(indexbound2));
            Rp(kp) = (fr2(indexbound2)'*vp(indexbound2) - fp(indexbound2)'*ur2(indexbound2));
            if dolcurve == 0 % No L-curve stuff
               fourn(kp) = - 1/(1+mu*k^2) *(1+nu)/(2*E*L*lambda(kp)^2)*Rp(kp);
               akan(k) = 2*real(fourn(kp));
               bkan(k) = 2*imag(fourn(kp));
            end
         else
            Rpm = (fr2(indexbound2)'*vp(indexbound2) - fp(indexbound2)'*ur2(indexbound2));
            fournm(kp) = -(1+nu)/(2*E*L*lambda(kp)^2)*Rpm;
         end
      end
   end
   
   % The constant term
   fourn(1) = -(R12+R22)/L;
   
   % Invert the operator
   if dolcurve == 1
      listmu = -1:.1:3;
      resid  = zeros( size(listmu) );
      regno  = zeros( size(listmu) );
      rhsk   = zeros(nbase+1,1);
      i = 1;
      for lnmu1 = listmu
         mu1 = 10^lnmu1;
         for kp=2:nbase+1
            %lhsk(kp)  = - (2*E*L*lambda(kp)^2)/(1+nu);
            rhsk(kp)  = - (1+nu)/(2*E*L*lambda(kp)^2)*Rp(kp);
            fourn(kp) = 1/(1+mu1*(kp-1)^2) * rhsk(kp);
            akan(k)   = 2*real(fourn(kp));
            bkan(k)   = 2*imag(fourn(kp));
            resid(i)  = resid(i) + abs(fourn(kp) - rhsk(kp))^2;
            regno(i)  = regno(i) + kp^2 * abs(fourn(kp))^2;
         end
         i = i+1;
      end
      figure;
      loglog(resid,regno);
      xlabel('residual (log)')
      ylabel('norm (log)')
   end
   
   % Plot the reference normal displacement (first test case)
   solu = fourn(1) + sum ( [0*newX' ;...  % Hack in case there is 1 element only
              akan.*cos(lambda(2:end)*newX') + bkan.*sin(lambda(2:end)*newX') ] );
   solu = solu';
   
   figure
   hold on;
   plot(newX, [0*newXo;ubase(:,2);0*newXo])
   plot(newX, solu, 'Color', 'red')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mean square polynoms
if usepolys == 1
   % First, build the polynomial test functions.
   coef = zeros(ordp+2, ordp+1);
   coef( 1:2, 1 ) = [-(1-nu^2)/(nu*E) ; 0];
   
   for k = 1:ordp
      Rhsco = zeros(k+1,1);
      Lhsco = zeros(k+1);
      
      Axx = zeros(k-1,2);
      Axy = zeros(k-1,2);
      Ayy = zeros(k-1,2);
      Bxx = zeros(k-1,2);
      Bxy = zeros(k-1,2);
      Byy = zeros(k-1,2);
      
      azero = -(1-nu^2)/(nu*E*(k+1));
   
      for i=0:floor( (k-1)/2 )
         Lhsco(2*i+1,2*i+1) = (k+1-2*i)*(k-2*i)*E/(1-nu^2);  % a_i
         Lhsco(2*i+1,2*i+2) = (2*i+1)*(k-2*i) * ( nu*E/(1-nu^2) + E/(2*(1+nu)) ); % b_i
         Lhsco(2*i+1,2*i+3) = (2*i+1)*(2*i+2)*E/(2*(1+nu));  % a_{i+1}
      end
      
      for i=1:floor( k/2 )
         Lhsco(2*i,2*i)   = (k-2*i+2)*(k-2*i+1)*E/(2*(1+nu)) ; % b_{i-1}
         Lhsco(2*i,2*i+1) = 2*i*(k+1-2*i)*( E/(2*(1+nu)) + nu*E/(1-nu^2) ) ; % a_i
         Lhsco(2*i,2*i+2) = 2*i*(2*i+1)*E/(1-nu^2) ; %b_i
      end
   
      C = [ eye(2) , zeros(2,k)];
      Lhsco = [ Lhsco ; C ]; % find an unique solution
      Rhsco = [ Rhsco ; azero ; 0 ];
      
      Lhsco(size(Lhsco,1)-2,:) = [];  % 'cause square matrix is life
      Rhsco(size(Rhsco,1)-2,:) = [];
      
      coef( 1:k+2, k+1 ) = Lhsco\Rhsco;
   end
   
   % Place zeros in coef
   coefa = coef;
   coefb = coef;
   for i=1:size(coefa,1)
      if mod(i,2) == 0
         coefa(i,:) = 0;
      end
      if mod(i,2) == 1
         coefb(i,:) = 0;
      end
   end
   
   %% Compute the RG
   Rhs = zeros(ordp+1,1);
   vpa = zeros(2*nnodes2, 1);
   for k=0:ordp
      for i=1:nnodes2
         x = nodes2(i,1);
         y = nodes2(i,2);
         ixigrec = Q'*[x;y]; X = (ixigrec(1)-offset)/L; Y = (ixigrec(2)+CteR)/L; %% /!\ CteR
         
         % Build X^k*Y^j
         GROX = zeros(ordp+2,1);
         for j = 0:k+1
            GROX(j+1) = X^(k+1-j)*Y^j;
         end
         
         vloc = [ coefa(:,k+1)'*GROX ; coefb(:,k+1)'*GROX ];
         vxy = Q*vloc; vpa(2*i-1) = vxy(1); vpa(2*i) = vxy(2);
      end
      fpa = Kinter2*vpa;
      %Rhs(k+1) = (fr1'*vpa - fpa'*ur1);
      Rhs(k+1) = (fr2(indexbound2)'*vpa(indexbound2) - fpa(indexbound2)'*ur2(indexbound2));
   end
   
   L1 = offset; L2 = offset+L;
   L1 = 0; L2 = 1;
   for i=0:ordp
      for j=0:ordp
         ord = i+j+1;
         Lhs(i+1,j+1) = (L2^ord - L1^ord)/ord;
         if i>1 && j>1
            Lhs(i+1,j+1) = Lhs(i+1,j+1) + mu*i*j/(i+j-1)*...
                                          (L2^(i+j-1) - L1^(i+j-1));
         end
      end
   end

      % Invert the operator
   if dolcurve == 1
      listmu = -6:.1:-2;
      resid  = zeros( size(listmu) );
      regno  = zeros( size(listmu) );
      rhsk   = zeros(nbase+1,1);
      index = 1;
      for lnmu1 = listmu
         Regop  = Lhs-Lhs;
         mu1 = 10^lnmu1;
         
         % Regularize
         for i=0:ordp
            for j=0:ordp
               if i>1 && j>1
                  Regop(i+1,j+1) = Regop(i+1,j+1) + mu1*i*j/(i+j-1)*...
                                                   (L2^(i+j-1) - L1^(i+j-1));
               end
            end
         end
         
         % Invert
         McCoef = (Lhs+Regop)\Rhs;
         resid(index)  = norm(Lhs*McCoef - Rhs)^2; % Actually, it's resid^2
         regno(index)  = McCoef'*(Regop/mu1)*McCoef;
         index = index+1;
      end
      figure;
      plot(resid,regno);
      xlabel('residual')
      ylabel('norm')
      %res10 = 10.^resid'; reg10 = 10.^regno';
      mu = 10^listmu( findCorner (resid', regno',2,0) );
   end
   
   % Regularize
   for i=0:ordp
      for j=0:ordp
         if i>1 && j>1
            Lhs(i+1,j+1) = Lhs(i+1,j+1) + mu*i*j/(i+j-1)*...
                                          (L2^(i+j-1) - L1^(i+j-1));
         end
      end
   end
   McCoef = Lhs\Rhs;
   nbase = size(McCoef,1);
   
   % Plot the result
   solref = [0*newXo ; ubase(:,2) ; 0*newXo];
   McCoefref = polyfit( newX, solref, nbase-1 );
   McCoefref = McCoefref(end:-1:1);
   
   solu = zeros(size(newX,1),1); solpref = zeros(size(newX,1),1);
   for i=1:nbase
      solu = solu + McCoef(i)*( (newX-offset)./L).^(i-1);
   end
   for i=1:nbase
      solpref = solpref + McCoefref(i)*newX.^(i-1);
   end
   
   figure
   hold on;
   plot(newX, solref)
   plot(newX, solu, 'Color', 'red')
   %plot(newX, solpref, 'Color', 'green')
   
   % Check the MS :
   Sref = sum((solpref-solref).^2) / sum(solref.^2);
   Smc  = sum((solu-solref).^2) / sum(solref.^2);
end
