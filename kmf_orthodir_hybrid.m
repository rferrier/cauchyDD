%04/05/2016
%Algo KMF hybride avec Orthodir

close all;
clear all;

addpath(genpath('./tools'))
addpath(genpath('./regu'))

% Parameters
E        = 70000;      % MPa : Young modulus
nu       = 0.3;    % Poisson ratio
fscalar  = 1;      % N.mm-1 : Loading on the plate
niter    = 10;
br       = 0.;      % noise
mu       = 0;      % Regularization parameter

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

[node2t, t2node] = mapBound(3, boundary, nnodes);

% Then, build the stiffness matrix :
[K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
Kinter = K(1:2*nnodes, 1:2*nnodes);
M      = mass_mat(nodes, elements);
% The right hand side :
f = loading(nbloq,nodes,boundary,neumann);

% Solve the problem :
uin = K\f;

% Extract the index of the boundary
[ node2b3, b2node3 ] = mapBound( 3, boundary, nnodes );
index = 2*b2node3-1;

% Extract on the suraboundant surfaces
[ node2b1, b2node1 ] = mapBound( 1, boundary, nnodes );
[ node2b2, b2node2 ] = mapBound( 2, boundary, nnodes );
index12 = [2*b2node1; 2*b2node1-1; 2*b2node2; 2*b2node2-1];
index1 = [2*b2node1; 2*b2node1-1];
index2 = [2*b2node2; 2*b2node2-1];

% Extract displacement and Lagrange multiplicators :
uref = uin(1:2*nnodes,1);
lagr = uin(2*nnodes+1:end,1);
%urefb = uref;
%urefb(index1) = ( 1 + br*randn( size(index1,1) ,1) ) .* uref(index1);
urefb = ( 1 + br*randn(2*nnodes,1) ) .* uref;
lagrb = ( 1 + br*randn(nbloq,1) ) .* lagr;

% Post-pro :
sigma = stress(uref,E,nu,nodes,elements,order,1,ntoelem);
ui = reshape(uref,2,[])';  ux = ui(:,1);  uy = ui(:,2);
plotGMSH({ux,'U_x';uy,'U_y';uref,'U_vect';sigma,'stress'}, elements, nodes, 'reference');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition of the operators

% DN problem
dirichlet11 = [4,1,0;4,2,0;
               1,1,0;
               2,1,0;
               3,2,0];
neumann11   = [1,2,fscalar];

dirichlet10 = [4,1,0;4,2,0;
               3,1,0;3,2,0];
neumann10   = [1,2,fscalar;
               2,1,fscalar];

[K11,C11,nbloq11,node2c11,c2node11] = Krig (nodes,elements,E,nu,order,boundary,dirichlet11);
[K10,C10,nbloq10,node2c10,c2node10] = Krig (nodes,elements,E,nu,order,boundary,dirichlet10);

% ND problem
dirichlet21 = [4,1,0;4,2,0;
               1,2,0;
               2,2,0;
               3,1,0];
neumann21   = [2,1,fscalar];

dirichlet20 = [4,1,0;4,2,0;
               1,1,0;1,2,0;
               2,1,0;2,2,0];
neumann20   = [];

[K21,C21,nbloq21,node2c21,c2node21] = Krig (nodes,elements,E,nu,order,boundary,dirichlet21);
[K20,C20,nbloq20,node2c20,c2node20] = Krig (nodes,elements,E,nu,order,boundary,dirichlet20);

error   = zeros(niter+1,1);
residual = zeros(niter+1,1);
regulari = zeros(niter+1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ORTHODIR for the problem : 
% (1 - D0 o S0 o DH0 o SH0 + mu*N) x = DoSoDHoSH(0), with N, a 
% regulalization term
Itere = zeros( 2*nnodes, niter+1 );
p     = zeros( 2*nnodes, niter+1 );
q     = zeros( 2*nnodes, niter+1 );
Delta = zeros( niter+1, 1 );
Res   = zeros( 2*nnodes, niter+1 );

%% Perform A x0 :
% Solve DN
fdir = dirichletRhs2(Itere(:,1), 3, c2node10, boundary, nnodes);
f1 = fdir;
uin1 = K10\f1;
u1 = uin1(1:2*nnodes,1);

% Solve ND
fr = [Kinter*u1;zeros(nbloq20,1)];
f2 = fr;
uin2 = K20\f2;
u2 = uin2(1:2*nnodes,1);

% Solve DN Hybrid
% fr = [Kinter*Itere(:,1);zeros(nbloq11,1)];
% fdir = dirichletRhs2(Itere(:,1), 3, c2node11, boundary, nnodes, 2);
fr = [Kinter*u2;zeros(nbloq11,1)];
fdir = dirichletRhs2(u2, 3, c2node11, boundary, nnodes, 2);
f1 = fdir + fr;
uin1 = K11\f1;
u1 = uin1(1:2*nnodes,1);

% Solve ND Hybrid
fr = [Kinter*u1;zeros(nbloq21,1)];
fdir = dirichletRhs2(u1, 3, c2node21, boundary, nnodes, 1);
f1 = fdir + fr;
uin1 = K21\f1;
u2 = uin1(1:2*nnodes,1);

% Regularization term
Nu = regul(Itere(:,1), nodes, boundary, 3);
atimesItere = mu*Nu + Itere(:,1) - u2;
%%%%

%% Write Rhs :
% Solve DN
f1in = loading(nbloq10,nodes,boundary,neumann10);
f1 = f1in;

uin1 = K10\f1;
u1 = uin1(1:2*nnodes,1);

% Solve ND
fr = [Kinter*u1-f1in(1:2*nnodes) ; zeros(nbloq20,1)];
fdir1 = dirichletRhs(urefb, 1, C20, boundary);
fdir2 = dirichletRhs(urefb, 2, C20, boundary);
f2 = fr + assembleDirichlet( [fdir1,fdir2] );

uin2 = K20\f2;
u2 = uin2(1:2*nnodes,1);

% Solve DN Hybrid
fr = [Kinter*u2 ; zeros(nbloq11,1)];
f1in = loading(nbloq11,nodes,boundary,neumann11);
fdir = dirichletRhs2(u2, 3, c2node11, boundary, nnodes, 2);
fdir1 = dirichletRhs2(urefb, 1, c2node11, boundary, nnodes, 1);
fdir2 = dirichletRhs2(urefb, 2, c2node11, boundary, nnodes, 1);
f1 = fr + f1in + assembleDirichlet( [fdir,fdir1,fdir2] );

uin1 = K11\f1;
u1 = uin1(1:2*nnodes,1);

% Solve ND Hybrid
fr = [Kinter*u1 - f1in(1:2*nnodes) ; zeros(nbloq21,1)];
f2in = loading(nbloq21,nodes,boundary,neumann21);
fdir = dirichletRhs2(u1, 3, c2node21, boundary, nnodes, 1);
fdir1 = dirichletRhs2(urefb, 1, c2node21, boundary, nnodes, 2);
fdir2 = dirichletRhs2(urefb, 2, c2node21, boundary, nnodes, 2);
f2 = fr + f2in + assembleDirichlet( [fdir,fdir1,fdir2] );

uin1 = K21\f2;
u2 = uin1(1:2*nnodes,1);

b = u2;
%%%%

%%
Res(:,1) = b - atimesItere;
p(:,1) = Res(:,1);

residual(1) = sqrt( myps( Res(:,1), Res(:,1), Kinter, boundary, M, nodes ) );
error(1)    = sqrt( myps( Itere(:,1) - uref, Itere(:,1) - uref, Kinter, boundary, M, nodes )...
                 / myps( uref, uref, Kinter, boundary, M, nodes ) );
regulari(1) = sqrt( Itere(:,1)'*regul(Itere(:,1), nodes, boundary, 3) );

%% Perform Q1 = A P1 :
% Solve DN
fdir = dirichletRhs(p(:,1), 3, C10, boundary);
f1 = fdir;
uin1 = K10\f1;
u1 = uin1(1:2*nnodes,1);

% Solve ND
fr = [Kinter*u1;zeros(nbloq20,1)];
f2 = fr;
uin2 = K20\f2;
u2 = uin2(1:2*nnodes,1);

% Solve DN Hybrid
% fr = [Kinter*p(:,1);zeros(nbloq11,1)];
fr = [Kinter*u2;zeros(nbloq11,1)];
% fdir = dirichletRhs2(p(:,1), 3, c2node11, boundary, nnodes, 2);
fdir = dirichletRhs2(u2, 3, c2node11, boundary, nnodes, 2);
f1 = fdir + fr;
uin1 = K11\f1;
u1 = uin1(1:2*nnodes,1);

% Solve ND Hybrid
fr = [Kinter*u1;zeros(nbloq21,1)];
fdir = dirichletRhs2(u1, 3, c2node21, boundary, nnodes, 1);
f1 = fdir + fr;
uin1 = K21\f1;
u2 = uin1(1:2*nnodes,1);

% Regularization term
Nu = regul(p(:,1), nodes, boundary, 3);

q(:,1) = mu*Nu + p(:,1) - u2;
%%%%
%%
for iter = 1:niter
    
    Delta(iter,1)   = myps( q(:,iter), q(:,iter), Kinter, boundary, M, nodes ); %q(:,iter)'*q(:,iter);
    gammai          = myps( q(:,iter), Res(:,iter), Kinter, boundary, M, nodes ); %q(:,iter)'*Res;
    alphai          = gammai/Delta(iter,1);
    
    Itere(:,iter+1) = Itere(:,iter) + p(:,iter)*alphai;
    Res(:,iter+1)   = Res(:,iter) - q(:,iter)*alphai;
    
    residual(iter+1) = sqrt( myps( Res(:,iter+1), Res(:,iter+1), Kinter, boundary, M, nodes ) );
    error(iter+1)    = sqrt( myps( Itere(:,iter+1) - uref, Itere(:,iter+1) - uref, Kinter, boundary, M, nodes )...
                     / myps( uref, uref, Kinter, boundary, M, nodes ) );
    regulari(iter+1) = sqrt( Itere(:,iter+1)'*regul(Itere(:,iter+1), nodes, boundary, 3) );

    %% Perform Ari = A*Res
    % Solve DN
    fdir = dirichletRhs(Res(:,iter+1), 3, C10, boundary);
    f1 = fdir;
    uin1 = K10\f1;
    u1 = uin1(1:2*nnodes,1);

    % Solve ND
    fr = [Kinter*u1;zeros(nbloq20,1)];
    f2 = fr;
    uin2 = K20\f2;
    u2 = uin2(1:2*nnodes,1);

    % Solve DN Hybrid
%     fr = [Kinter*Res(:,iter+1);zeros(nbloq11,1)];
%     fdir = dirichletRhs2(Res(:,iter+1), 3, c2node11, boundary, nnodes, 2);
    fr = [Kinter*u2;zeros(nbloq11,1)];
    fdir = dirichletRhs2(u2, 3, c2node11, boundary, nnodes, 2);
    f1 = fdir + fr;
    uin1 = K11\f1;
    u1 = uin1(1:2*nnodes,1);

    % Solve ND Hybrid
    fr = [Kinter*u1;zeros(nbloq21,1)];
    fdir = dirichletRhs2(u1, 3, c2node21, boundary, nnodes, 1);
    f1 = fdir + fr;
    uin1 = K21\f1;
    u2 = uin1(1:2*nnodes,1);

    % Regularization term
    Nu = regul(Res(:,iter+1), nodes, boundary, 3);
    Ari = mu*Nu + Res(:,iter+1) - u2;
    
    %% Orthogonalization
    p(:,iter+1) = Res(:,iter+1);
    q(:,iter+1) = Ari;
    
    for jter=1:iter
        phiij  = myps( q(:,jter), Ari, Kinter, boundary, M, nodes ); %q(:,jter)'*Ari;
        betaij = phiij/Delta(jter,1);
        p(:,iter+1) = p(:,iter+1) - p(:,jter) * betaij;
        q(:,iter+1) = q(:,iter+1) - q(:,jter) * betaij;
        %myps( q(:,iter+1), q(:,jter), Kinter, boundary, M, nodes )
    end
    
end

hold on
plot(log10(error),'Color','blue')
plot(log10(residual),'Color','red')
legend('error (log)','residual (log)')
figure
% L-curve
loglog(residual,regulari);

%% Automatic stop (regularization tools) :
% [stopHere,~,~] = l_corner(residual,regulari,1:1:niter+1);
% % Plot chosen Itere
% hold on;
% plot(uref(index));
% plot(Itere(index,stopHere),'Color','red');
% figure
%%%%
stopHere = niter+1;
%% Final problem : compute u
% DN problem
dirichlet = [4,1,0;4,2,0;
             3,1,0;3,2,0];
neumann   = [];
[K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
f1in = loading(nbloq,nodes,boundary,neumann);
fdir3 = dirichletRhs(Itere(:,stopHere), 3, C, boundary);
usoli = K \ ( assembleDirichlet(fdir3) + f1in );
usol = usoli(1:2*nnodes,1);

%% Boundary error => mesh adaptation

% hold on
% plot( abs( urefb(index12)-usol(index12) ) / norm( urefb(index12) ) )
% plot( abs( urefb(index12)-uref(index12) ) / norm( uref(index12) ), 'Color', 'red' )

total_error = norm(usol-uref)/norm(uref);

% Post-pro :
sigma = stress(usol,E,nu,nodes,elements,order,1,ntoelem);
ui = reshape(usol,2,[])';  ux = ui(:,1);  uy = ui(:,2);
plotGMSH({ux,'U_x';uy,'U_y';usol,'U_vect';sigma,'stress'}, elements, nodes, 'solution');

plotGMSH({(usol-uref)/norm(uref),'Error'}, elements, nodes, 'error');