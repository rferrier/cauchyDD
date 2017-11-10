% 09/11/2017
% Problème de tuyau avec Helmholtz, KMF Orthodir

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 70000;   % MPa : Young modulus
nu      = 0.3;     % Poisson ratio
rho     = 2700e-9;  % kg.mm-3 : Volumic mass
fscalar = 1;       % N.mm-1 : Loading on the plate
br      = 0.;      % noise
omega   = 0;%5e8;      % s-2 : square of the pulsation (dynamic case)
niter   = 20;
mu      = 0;      % Regularization parameter

% Methods : 1=KMF, 2=KMF Orthodir, 3=KMF Robin, 4=SPP, 5=SPD,
% 6=SPD flottant, 7=SPD flottant constraint, 8=evanescent regu
% 9=SPP GC Ritz, 10=SPD GC Ritz
% 100=KMF-R+ERC

% Boundary conditions
% first index  : index of the boundary
% second index : 1=x, 2=y
% third        : value
% [0,1,value] marks a dirichlet regularization therm on x
dirichlet   = [2,1,0; 2,2,0];

% Import the meshes
[ nodes,elements,ntoelem,boundary,order ] =...
    readmesh( 'meshes/tube_fat.msh' );
nnodes = size(nodes,1);

%noises = load('./noises/noisetube1.mat'); % Particular noise vector
%noise  = noises.bruit1;
noise  = randn(2*nnodes,1);

% patch('Faces',elements,'Vertices',nodes,'FaceAlpha',0);
% figure

% Then, build the stiffness matrix :
[K,C,nbloq,node2c,c2node] =...
    Krig (nodes,elements,E,nu,order,boundary,dirichlet,1);
Kinter = K(1:2*nnodes, 1:2*nnodes);
M0 = mass_mat(nodes, elements);
M0 = rho*M0;
M = [ M0 , zeros(2*nnodes,nbloq) ; zeros(2*nnodes,nbloq)' , zeros(nbloq) ];
[ node2b1, b2node1 ] = mapBound( 1, boundary, nnodes );
[ node2b2, b2node2 ] = mapBound( 2, boundary, nnodes );
[ node2b3, b2node3 ] = mapBound( 3, boundary, nnodes );
[ node2b9, b2node9 ] = mapBound( 9, boundary, nnodes );

% The right hand side :
f = pressureLoad( nbloq, nodes, boundary, fscalar*[10*sin(pi/6),0;-1,0], 8 );

% Solve the problem :
uin = (K-omega*M)\f;

% Extract displacement and Lagrange multiplicators :
uref = uin(1:2*nnodes,1);
fref = (Kinter-omega*M0)*uref;
lagr = uin(2*nnodes+1:end,1);
urefb = ( 1 + br*noise ) .* uref;
frefb = ( 1 + br*noise ) .* fref;

% Compute stress :
sigma = stress(uref,E,nu,nodes,elements,order,1,ntoelem,1);
% Output :
plotGMSH({uref,'U_vect';sigma,'stress'}, elements, nodes, 'output/reference');

% Plot displacement on the interface :
%index = 2*[b2node1;b2node2;b2node3];
index = 2*b2node3;
thetax = 0:2*pi/size(index,1):2*pi*(1-1/size(index,1));
hold on
set(gca, 'fontsize', 15);
plot(thetax,uref(index,1));
plot(thetax,uref(index-1,1),'Color','red');
legend('uy','ux')
xlabel('angle(rad)')
figure
hold on
set(gca, 'fontsize', 15);
plot(thetax,fref(index,1));
plot(thetax,fref(index-1,1),'Color','red');
legend('fy','fx')
xlabel('angle(rad)')

indexxy = [index-1;index];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KMF prequisites
% DN problem
dirichlet1 = [2,1,0;2,2,0;
              3,1,0;3,2,0];

[K1,C1,nbloq1,node2c1,c2node1] =...
    Krig (nodes,elements,E,nu,order,boundary,dirichlet1,1);
M1 = [ M0 , zeros(2*nnodes,nbloq1) ; zeros(2*nnodes,nbloq1)' , zeros(nbloq1) ];
% ND problem
dirichlet2 = [2,1,0;2,2,0;
              1,1,0;1,2,0];

[K2,C2,nbloq2,node2c2,c2node2] =...
    Krig (nodes,elements,E,nu,order,boundary,dirichlet2,1);
M2 = [ M0 , zeros(2*nnodes,nbloq2) ; zeros(2*nnodes,nbloq2)' , zeros(nbloq2) ];
% Dirichlet loading :
f2  = dirichletRhs2( urefb, 1, c2node2, boundary, nnodes );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KMF Orthodir 10 iterations

% Init
Itere    = zeros( 2*nnodes, 1 );
p        = zeros( 2*nnodes, niter+1 );
q        = zeros( 2*nnodes, niter+1 );
Delta    = zeros( niter+1, 1 );
Res      = zeros( 2*nnodes, niter+1 );
error    = zeros(niter+1,1);
residual = zeros(niter+1,1);
regulari = zeros(niter+1,1);
 
%% Compute residual
% Ax0 : Solve DN
f1 = dirichletRhs2(Itere, 3, c2node1, boundary, nnodes );
uin1 = (K1-omega*M1)\f1;
u1 = uin1(1:2*nnodes,1);
% Keep only the forces at the bottom boundary
fri = (Kinter-omega*M0)*u1;
% Solve ND
fr = [ fri ; zeros(size(C2,2),1) ]; % Give to fr the right size
uin2 = (K2-omega*M2)\fr;
%
Nu = regul(Itere, nodes, boundary, 3);
atimesItere = mu*Nu + Itere - uin2(1:2*nnodes,1);
%
% RHS : Solve DN (this one is useless because 0)
f1 = zeros(2*nnodes+nbloq1,1);
uin1 = (K1-omega*M1)\f1;
u1 = uin1(1:2*nnodes,1);
% Keep only the forces at the bottom boundary
fri = (Kinter-omega*M0)*u1;
% Solve ND
fr = [ fri ; zeros(size(C2,2),1) ]; % Give to fr the right size
uin2 = (K2-omega*M2)\(fr+f2);
%
b = uin2(1:2*nnodes,1);
%
Res(:,1) = b - atimesItere;
p(:,1) = Res(:,1);

residual(1) = 1; %norm(Res(indexxy,1));
error(1)    = norm( Itere(indexxy) - uref(indexxy )) / norm(uref(indexxy));
regulari(1) = sqrt(Itere'*( regul(Itere, nodes, boundary, 3)));

 
%% Compute Q1 = AP1
% Solve DN
f1 = dirichletRhs2(p(:,1), 3, c2node1, boundary, nnodes );
uin1 = (K1-omega*M1)\f1;
u1 = uin1(1:2*nnodes,1);
% Keep only the forces at the bottom boundary
fri = (Kinter-omega*M0)*u1;
% Solve ND
fr = [ fri ; zeros(size(C2,2),1) ]; % Give to fr the right size
uin2 = (K2-omega*M2)\fr;
%
Nu = regul(p(:,1), nodes, boundary, 3);
q(:,1) = mu*Nu + p(:,1) - uin2(1:2*nnodes,1);
 
for iter = 1:niter

    Delta(iter,1) = norm(q(indexxy,iter))^2;
    gammai        = q(indexxy,iter)'*Res(indexxy,iter);
    alphai        = gammai/Delta(iter,1);

    Itere         = Itere + p(:,iter)*alphai;
    Res(:,iter+1) = Res(:,iter) - q(:,iter)*alphai;

    residual(iter+1) = norm(Res(indexxy,iter+1)) / norm(Res(indexxy,1));
    error(iter+1)    = norm( Itere(indexxy) - uref(indexxy)) / norm(uref(indexxy));
    regulari(iter+1) = sqrt(Itere'*( regul(Itere, nodes, boundary, 3)));

    %% Perform Ari = A*Res
    % Solve DN
    f1 = dirichletRhs2(Res(:,iter+1), 3, c2node1, boundary, nnodes );
    uin1 = (K1-omega*M1)\f1;
    u1 = uin1(1:2*nnodes,1);
    % Keep only the forces at the bottom boundary
    fri = (Kinter-omega*M0)*u1; % only difference is on dirichlet boundary
%    fri = lagr2forces2( uin1(2*nnodes+1:end), c2node1, 3, boundary, nnodes );
    % Solve ND
    fr = [ fri ; zeros(size(C2,2),1) ]; % Give to fr the right size
    uin2 = (K2-omega*M2)\fr;
    %
    Nu = regul(Res(:,iter+1), nodes, boundary, 3);
    Ari = mu*Nu + Res(:,iter+1) - uin2(1:2*nnodes,1);

    %% Orthogonalization
    p(:,iter+1) = Res(:,iter+1);
    q(:,iter+1) = Ari;

     for jter=1:iter
        phiij  = q(indexxy,jter)'*Ari(indexxy);
        betaij = phiij/Delta(jter,1);
        p(:,iter+1) = p(:,iter+1) - p(:,jter) * betaij;
        q(:,iter+1) = q(:,iter+1) - q(:,jter) * betaij;
    end

end
figure
hold on
set(gca, 'fontsize', 15);
plot(log10(error),'Color','blue')
plot(log10(residual),'Color','red')
legend('error (log)','residual (log)')
%% L-curve
%figure
%loglog(residual,regulari);
% Output : terminal computation
dirichlet = [1,1,0;1,2,0;
             2,1,0;2,2,0;
             3,1,0;3,2,0];
[K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet,1);
M = [ M0 , zeros(2*nnodes,nbloq) ; zeros(2*nnodes,nbloq)' , zeros(nbloq) ];
fdir1 = dirichletRhs(urefb, 1, C, boundary);
fdir3 = dirichletRhs(Itere, 3, C, boundary);
usoli = (K-omega*M) \ ( fdir1 + fdir3 );
usol = usoli(1:2*nnodes,1);
 
efe = (Kinter-omega*M0)*usol;
figure
hold on
set(gca, 'fontsize', 15);
plot(thetax,usol(index,1));
plot(thetax,usol(index-1,1), 'Color', 'red');
legend('uy','ux')
xlabel('angle(rad)')
figure
hold on
set(gca, 'fontsize', 15);
plot(thetax,efe(index,1));
plot(thetax,efe(index-1,1), 'Color', 'red');
legend('fy','fx')
xlabel('angle(rad)')
 
total_error = norm(uref-usol)/norm(uref);
total_errorf = norm(fref-efe)/norm(fref);
 
% Compute stress :
sigma = stress(usol,E,nu,nodes,elements,order,1,ntoelem,1);
% Output :
plotGMSH({usol,'U_vect';sigma,'stress'}, elements, nodes, 'output/field2');