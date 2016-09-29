%06/04/2016
%Algo KMF Robin avec Orthodir

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 70000;  % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 1;      % N.mm-1 : Loading on the plate
niter   = 20;
br      = 0.0;      % bruit
robin1  = 1e5*E;%-7*E;%-1e1*E;      % Robin parameters
robin2  = -1/7*E;%-1/7*E;%-1e-1*E;
mu      = 0.0;   % Regularization parameter

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

% problem 1
dirichlet1 = [4,1,0;4,2,0];
neumann1   = [1,2,fscalar;
              2,1,fscalar];
[K1,C1,nbloq1,node2c1,c2node1] = Krig (nodes,elements,E,nu,order,boundary,dirichlet1);

f1in = loading(nbloq1,nodes,boundary,neumann1);

% problem 2
dirichlet2 = [4,1,0;4,2,0;
              1,1,0;1,2,0;
              2,1,0;2,2,0];

[K2,C2,nbloq2,node2c2,c2node2] = Krig (nodes,elements,E,nu,order,boundary,dirichlet2);
fdir1 = dirichletRhs2(urefb, 1, c2node2, boundary, nnodes); % Dirichlet BCs
fdir2 = dirichletRhs2(urefb, 2, c2node2, boundary, nnodes);
fdir  = assembleDirichlet( [fdir1,fdir2] );

error   = zeros(niter+1,1);
residual = zeros(niter+1,1);
regulari = zeros(niter+1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ORTHODIR for the problem : (1 - D0 o S0) [u;f] = DoS(0)
Itere = zeros( 4*nnodes, 1 ); % 2 iterates (u and f)
p     = zeros( 4*nnodes, niter+1 );
q     = zeros( 4*nnodes, niter+1 );
Delta = zeros( niter+1, 1 );
Res   = zeros( 4*nnodes, niter+1 );

%% Perform A x0 :
u2 = Itere(1:2*nnodes);    % Split Itere
f2 = Itere(2*nnodes+1:end);
% Solve 1
kuplusf = u2 + f2/robin1;
[Kp1, fro1] = robinRHS( nbloq1, nodes, boundary, kuplusf, robin1, 3 );
K1e = K1 + Kp1;
uin1 = K1e\fro1;
u1 = uin1(1:2*nnodes,1);
f1 = Kinter*u1;  % Get the forces at the boundary
% Solve 2   
kuplusf = u1 + f1/robin2;
[Kp2, fro2] = robinRHS( nbloq2, nodes, boundary, kuplusf, robin2, 3 );
K2e = K2 + Kp2;
uin2 = K2e\fro2;
u2 = uin2(1:2*nnodes,1);
f2 = Kinter*u2; % Get the forces at the boundary
Nu = regul(Itere(1:2*nnodes), nodes, boundary, 3); % Regularization
Nf = regul(Itere(2*nnodes+1:end), nodes, boundary, 3); % Regularization
atimesItere = mu * [Nu;Nf/E] + Itere - [u2;f2];
%%%%

%% Write Rhs :
% Solve 1
uin1 = K1e\f1in;
u1 = uin1(1:2*nnodes,1);
f1 = Kinter*u1 - f1in(1:2*nnodes,1); % Get the forces at the boundary
% Solve 2   
kuplusf = u1 + f1/robin2;
[Kp2, fro2] = robinRHS( nbloq2, nodes, boundary, kuplusf, robin2, 3 );
K2e = K2 + Kp2;
uin2 = K2e\(fdir+fro2);
u2 = uin2(1:2*nnodes,1);
f2 = Kinter*u2; % Get the forces at the boundary
b = [u2;f2];
%%%%

%%
Res(:,1) = b - atimesItere;
p(:,1) = Res(:,1);

residual(1) = sqrt( myps( Res(:,1), Res(:,1), Kinter, boundary, C1, nodes ) );
error(1)    = sqrt( myps( Itere(1:2*nnodes) - uref, Itere(1:2*nnodes) - uref, Kinter, boundary, C1, nodes )...
                 / myps( uref, uref, Kinter, boundary, C1, nodes ) );
regulari(1) = Itere(1:2*nnodes)'*regul(Itere(1:2*nnodes), nodes, boundary, 3);

%% Perform Q1 = A P1 :
u2 = p(1:2*nnodes,1);    % Split p
f2 = p(2*nnodes+1:end,1);
% Solve 1
kuplusf = u2 + f2/robin1;
[Kp1, fro1] = robinRHS( nbloq1, nodes, boundary, kuplusf, robin1, 3 );
K1e = K1 + Kp1;
uin1 = K1e\fro1;
u1 = uin1(1:2*nnodes,1);
f1 = Kinter*u1;% - f1in(1:2*nnodes,1); % Get the forces at the boundary
% Solve 2   
kuplusf = u1 + f1/robin2;
[Kp2, fro2] = robinRHS( nbloq2, nodes, boundary, kuplusf, robin2, 3 );
K2e = K2 + Kp2;
uin2 = K2e\fro2;
u2 = uin2(1:2*nnodes,1);
f2 = Kinter*u2;% Get the forces at the boundary
Nu = regul(p(1:2*nnodes,1), nodes, boundary, 3); % Regularization
Nf = regul(p(2*nnodes+1:end,1), nodes, boundary, 3); % Regularization
q(:,1) = mu * [Nu;Nf/E] + p(:,1) - [u2;f2];
%%%%
%%
for iter = 1:niter
    
    Delta(iter,1) = myps( q(:,iter), q(:,iter), Kinter, boundary, C1, nodes ); %q(:,iter)'*q(:,iter);
    gammai        = myps( q(:,iter), Res(:,iter), Kinter, boundary, C1, nodes ); %q(:,iter)'*Res;
%     Delta(iter,1) = myps( q(1:2*nnodes,iter), q(1:2*nnodes,iter), Kinter, boundary, C1, nodes ); %q(:,iter)'*q(:,iter);
%     gammai        = myps( q(1:2*nnodes,iter), Res(1:2*nnodes,iter), Kinter, boundary, C1, nodes ); %q(:,iter)'*Res;
    alphai        = gammai/Delta(iter,1);
    
    Itere         = Itere + p(:,iter)*alphai;
    Res(:,iter+1) = Res(:,iter) - q(:,iter)*alphai;
    
    residual(iter+1) = sqrt( myps( Res(:,iter+1), Res(:,iter+1), Kinter, boundary, C1, nodes ) );
    error(iter+1)    = sqrt( myps( Itere(1:2*nnodes) - uref, Itere(1:2*nnodes) - uref, Kinter, boundary, C1, nodes )...
                     / myps( uref, uref, Kinter, boundary, C1, nodes ) );
    regulari(iter+1) = Itere(1:2*nnodes)'*regul(Itere(1:2*nnodes), nodes, boundary, 3);

    %% Perform Ari = A*Res
    % Solve DN
    u2 = Res(1:2*nnodes,iter+1);    % Split Res
    f2 = Res(2*nnodes+1:end,iter+1);
    % Solve 1
    kuplusf = u2 + f2/robin1;
    [Kp1, fro1] = robinRHS( nbloq1, nodes, boundary, kuplusf, robin1, 3 );
    K1e = K1 + Kp1;
    uin1 = K1e\fro1;
    u1 = uin1(1:2*nnodes,1);
    f1 = Kinter*u1;% - f1in(1:2*nnodes,1); % Get the forces at the boundary
    % Solve 2   
    kuplusf = u1 + f1/robin2;
    [Kp2, fro2] = robinRHS( nbloq2, nodes, boundary, kuplusf, robin2, 3 );
    K2e = K2 + Kp2;
    uin2 = K2e\fro2;
    u2 = uin2(1:2*nnodes,1);
    f2 = Kinter*u2; % Get the forces at the boundary
    Nu = regul(Res(1:2*nnodes,iter+1), nodes, boundary, 3); % Regularization
    Nf = regul(Res(2*nnodes+1:end,iter+1), nodes, boundary, 3); % Regularization
    Ari = mu * [Nu;Nf/E] + Res(:,iter+1) - [u2;f2];
    
    %% Orthogonalization
    p(:,iter+1) = Res(:,iter+1);
    q(:,iter+1) = Ari;
    
    for jter=1:iter
        phiij  = myps( q(:,jter), Ari, Kinter, boundary, C1, nodes ); %q(:,jter)'*Ari;
%         phiij  = myps( q(1:2*nnodes,jter), Ari(1:2*nnodes), Kinter, boundary, C1, nodes ); %q(:,jter)'*Ari;
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
loglog(residual,regulari);
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
fdir3 = dirichletRhs(Itere(1:2*nnodes), 3, C, boundary);
usoli = K \ assembleDirichlet( [fdir1,fdir2,fdir3] );
usol = usoli(1:2*nnodes,1);

total_error = norm(usol-uref)/norm(uref);

% Post-pro :
sigma = stress(usol,E,nu,nodes,elements,order,1,ntoelem);
plotGMSH({usol,'U_vect';sigma,'stress'}, elements, nodes, 'solution');

plotGMSH({(usol-uref)/norm(uref),'Error'}, elements, nodes, 'error');