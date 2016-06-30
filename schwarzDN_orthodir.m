%24/03/2016
%Algorithme de Schwarz sans recouvrement (ND) accéléré avec Orthodir

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 70000;      % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 1;      % N.mm-1 : Loading on the plate
niter   = 15;

dirichlea = [4,1,0;4,2,0;
             1,1,0;1,2,0;
             2,1,0;2,2,0];
dirichleb = [4,1,0;4,2,0;
             3,1,0;3,2,0;
             2,1,0;2,2,0];
neumann   = [];

% First, import the mesh
[ nodea,elementa,ntoelea,boundara,ordea ] = readmesh( 'meshes/plate.msh' );
nodeb = translateMesh( nodea, 0, 10 );
nnodea = size(nodea,1);

% Then, build the stiffness matrix :
[Ka,Ca,nbloqa,node2ca,c2nodea] =...
    Krig (nodea,elementa,E,nu,ordea,boundara,dirichlea);
[Kb,Cb,nbloqb,node2cb,c2nodeb] =...
    Krig (nodeb,elementa,E,nu,ordea,boundara,dirichleb);
%Kab = penalEq( Ka, Kb, E*1e9, nodea, nodeb, 1e-5 );
[ Kab, naj] = lagrEq( Ka, Kb, nodea, nodeb, 1e-5,c2nodea,c2nodeb);
%Kab = [Ka,zeros(size(Ka,1), size(Kb,1)); zeros(size(Kb,1), size(Ka,1)),Kb];

% The right hand side :
fa = volumicLoad( nbloqa, nodea, elementa, 1, fscalar );
fb = volumicLoad( nbloqb, nodeb, elementa, 1, fscalar );

% Solve the problem :
uinab = Kab\[fa;fb;zeros(naj,1)];

% Extract displacement :
u = uinab(1:2*nnodea,1);

% Compute stress :
sigma = stress(u,E,nu,nodea,elementa,ordea,1,ntoelea);

% Output :
% plotNodes(u,elements,nodes); : TODO debug on Matlab r>2013
plotGMSH({u,'U_Vect';sigma,'stress'}, elementa, nodea(:,[1,2]), 'ref_field');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Schwarz ND algorithm

% We work on an half-mesh
[ nodes1,elements,ntoelem,boundary,order ] = readmesh( 'meshes/plate.msh' );
nnodes = size(nodes1,1);

[node2b1, b2node1] = mapBound(1, boundary, nnodes);
[node2b3, b2node3] = mapBound(3, boundary, nnodes);

% DN problem
dirichlet1 = [4,1,0;4,2,0;
              1,1,0;1,2,0;
              2,1,0;2,2,0;
              3,1,0;3,2,0];
neumann1   = [];
[K1,C1,nbloq1] = Krig (nodes1,elements,E,nu,order,boundary,dirichlet1);
f1in = volumicLoad( nbloq1, nodes1, elements, 1, fscalar );
Kinter = K1(1:2*nnodes, 1:2*nnodes);
% ND problem
nodes2 = translateMesh(nodes1, 0, 10); % Move the mesh2

dirichlet2 = [4,1,0;4,2,0;
              3,1,0;3,2,0;
              2,1,0;2,2,0];
neumann2   = [];% [1,1,lagr1; 1,2,lagr1]
                % is managed by lagr2forces
[K2,C2,nbloq2] = Krig (nodes2,elements,E,nu,order,boundary,dirichlet2);
f2in = volumicLoad( nbloq2, nodes2, elements, 1, fscalar );

%%
error = zeros(niter+1,1);
residual = zeros(niter+1,1);

% Ref solution
urni = projectField(u, nodea, nodes2, 1e-5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ORTHODIR for the problem : (1 - D0 o S0) x = DoS(0)
Itere = zeros( 2*nnodes, 1 );       % Itere lives on 2
p     = zeros( 2*nnodes, niter+1 );
q     = zeros( 2*nnodes, niter+1 );
Delta = zeros( niter+1, 1 );
Res   = zeros( 2*nnodes, niter+1 );

%% Perform A x0 :
% Solve DN
ui2 = projectField( Itere, nodes2, nodes1, 1e-5 );
f1 = dirichletRhs(ui2, 3, C1, boundary);
uin1 = K1\f1;
lagr1 = uin1(2*nnodes+1:end,1);
% Solve ND
fri = -lagr2forces( lagr1, C1, 3, boundary );
% fribis = keepField(Kinter*u1, 3, boundary);
% norm(fribis-fri)
frii = projectField( fri, nodes1, nodes2, 1e-5 );
f2 = [ frii; zeros(size(C2,2),1) ]; % Give to fr the right size
uin2 = K2\f2;
%
atimesItere = Itere - uin2(1:2*nnodes,1);
%%%%

%% Compute Rhs :
% Solve DN
f1 = f1in;
uin1 = K1\f1;
lagr1 = uin1(2*nnodes+1:end,1);
% Solve ND
fri = -lagr2forces( lagr1, C1, 3, boundary );
frii = projectField( fri, nodes1, nodes2, 1e-5 );
fr = [ frii; zeros(size(C2,2),1) ]; % Give to fr the right size
f2 = f2in + fr;
uin2 = K2\f2;
%
b = uin2(1:2*nnodes,1);
%%%%
%%
Res(:,1) = b - atimesItere;
p(:,1) = Res(:,1);

residual(1) = sqrt( norm( keepField( Res(:,1), 1, boundary ) ) );
error(1)    = sqrt( norm( keepField( Itere(:,1), 1, boundary )  - ...
   keepField( urni, 1, boundary ) ) / ...
   norm( keepField( urni, 1, boundary ) ) );

%% Perform Q1 = A P1 :
% Solve DN
ui2 = projectField( p(:,1), nodes2, nodes1, 1e-5 );
f1 = dirichletRhs(ui2, 3, C1, boundary);
uin1 = K1\f1;
lagr1 = uin1(2*nnodes+1:end,1);
% Solve ND
fri = -lagr2forces( lagr1, C1, 3, boundary );
frii = projectField( fri, nodes1, nodes2, 1e-5 );
f2 = [ frii; zeros(size(C2,2),1) ]; % Give to fr the right size
uin2 = K2\f2;
%
q(:,1) = p(:,1) - uin2(1:2*nnodes,1);
%%%%
%%
for iter = 1:niter
    
    Delta(iter,1) = norm(keepField( q(:,iter), 1, boundary ))^2; %q(:,iter)'*q(:,iter);
    gammai        = transpose(keepField( q(:,iter), 1, boundary )) *...
                    keepField( Res(:,iter), 1, boundary ); %q(:,iter)'*Res;
    alphai        = gammai/Delta(iter,1);
    
    Itere         = Itere + p(:,iter)*alphai;
    Res(:,iter+1) = Res(:,iter) - q(:,iter)*alphai;
    
    residual(iter+1) = sqrt( norm(keepField( Res(:,iter+1), 1, boundary )) );
    error(iter+1)    = sqrt( norm( keepField( Itere, 1, boundary )  - ...
       keepField( urni, 1, boundary ) ) / ...
       norm( keepField( urni, 1, boundary ) ) );

    %% Perform Ari = A*Res
    % Solve DN
    ui2 = projectField( Res(:,iter+1), nodes2, nodes1, 1e-5 );
    f1 = dirichletRhs(ui2, 3, C1, boundary);
    uin1 = K1\f1;
    lagr1 = uin1(2*nnodes+1:end,1);
    % Solve ND
    fri = -lagr2forces( lagr1, C1, 3, boundary );
    frii = projectField( fri, nodes1, nodes2, 1e-5 );
    f2 = [ frii; zeros(size(C2,2),1) ]; % Give to fr the right size
    uin2 = K2\f2;
    %
    Ari = Res(:,iter+1) - uin2(1:2*nnodes,1);
    
    %% Orthogonalization
    p(:,iter+1) = Res(:,iter+1);
    q(:,iter+1) = Ari;
    
    for jter=1:iter
        phiij  = transpose(keepField( q(:,jter), 1, boundary )) *...
                    keepField( Ari, 1, boundary ); %q(:,jter)'*Ari;
        betaij = phiij/Delta(jter,1);
        p(:,iter+1) = p(:,iter+1) - p(:,jter) * betaij;
        q(:,iter+1) = q(:,iter+1) - q(:,jter) * betaij;
        %transpose(keepField( q(:,iter+1), 1, boundary )) * keepField( q(:,jter), 1, boundary )
    end
    
end
%%%%
%% Post process : compute the solution of the Dirichlet problem.
% Solve lower
ui2 = projectField( Itere, nodes2, nodes1, 1e-5 );
f1 = dirichletRhs(ui2, 3, C1, boundary) + f1in;
uin1 = K1\f1;
u1 = uin1(1:2*nnodes,1);

% Solve upper
% this is the matrix with Dirichlet BCs :
[K2,C2,nbloq2] = Krig (nodes2,elements,E,nu,order,boundary,dirichlet1);
f2in = volumicLoad( nbloq2, nodes2, elements, 1, fscalar );

f2 = dirichletRhs(Itere, 1, C2, boundary) + f2in;
uin2 = K2\f2;    
u2 = uin2(1:2*nnodes,1);

sigma1 = stress(u1,E,nu,nodes1,elements,order,1,ntoelem);
sigma2 = stress(u2,E,nu,nodes2,elements,order,1,ntoelem);

plotGMSH({u1,'U_Vect';sigma1,'stress'}, elements, nodes1(:,[1,2]), 'field1');
plotGMSH({u2,'U_Vect';sigma2,'stress'}, elements, nodes1(:,[1,2]), 'field2');
%%%%
%%
hold on;
set(gca, 'fontsize', 15);
set(gca,'ylim',[-7 0])
plot(log10(error(2:end)),'Color','blue')
plot(log10(residual(2:end)/residual(1)),'Color','red')
legend('error (log)','residual (log)')