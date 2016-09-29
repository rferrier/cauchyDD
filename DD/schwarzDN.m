%15/03/2016
%Code FEM 2D contraintes planes

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 70000;      % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 1;      % N.mm-1 : Loading on the plate
niter   = 15;
relax   = 0.5;    % Relaxation parameter

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
u  = uinab(1:2*nnodea,1);
ua = uinab(1:2*nnodea,1);
ub = uinab(2*nnodea+nbloqa+1:4*nnodea+nbloqa,1);

% Compute stress :
sigmaa = stress(ua,E,nu,nodea,elementa,ordea,1,ntoelea);
sigmab = stress(ub,E,nu,nodeb,elementa,ordea,1,ntoelea);

% Output :
% plotNodes(u,elements,nodes); : TODO debug on Matlab r>2013
ui = reshape(ua,2,[])';  ux = ui(:,1);  uy = ui(:,2);
plotGMSH({ux,'U_x';uy,'U_y';ua,'U_vect';sigmaa,'stress'}, elementa, nodea, 'ref_fielda');
ui = reshape(ub,2,[])';  ux = ui(:,1);  uy = ui(:,2);
plotGMSH({ux,'U_x';uy,'U_y';ub,'U_vect';sigmab,'stress'}, elementa, nodea, 'ref_fieldb');
% hold on;
% patch('Faces',elementa(:,1:3),'Vertices',nodea,'FaceAlpha',0);
% patch('Faces',elementa(:,1:3),'Vertices',nodeb,'FaceAlpha',0);
% %axis([0,5,0,20])
% axis equal
% figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Schwarz ND algorithm

% We work on an half-mesh
[ nodes1,elements,ntoelem,boundary,order ] = readmesh( 'meshes/plate.msh' );
nnodes = size(nodes1,1);

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

error1   = zeros(niter,1);
error2   = zeros(niter,1);
residual = zeros(niter,1);
u1       = zeros(2*nnodes,1);
u2       = zeros(2*nnodes,1);
lagr1    = zeros(nbloq1,1);

% Reference field
urn = dirichletRhs(u, 3, C1, boundary);

for iter = 1:niter
    % Solve DN
    ui2 = projectField( u2, nodes2, nodes1, 1e-5 );
    fdir = dirichletRhs(ui2, 3, C1, boundary);
    f1 = f1in + fdir;

    uin1 = K1\f1;
    u1 = relax*u1 + (1-relax)*uin1(1:2*nnodes,1);
    lagr1 = relax*lagr1 + (1-relax)*uin1(2*nnodes+1:end,1);
    
    % Computation of the error :
    u1n = dirichletRhs(uin1, 3, C1, boundary);
    %urn = dirichletRhs(u, 3, C1, boundary);
    error1(iter) = sqrt( (u1n-urn)'*(u1n-urn)/(urn'*urn) );
    
    % Solve ND
    fri = -lagr2forces( lagr1, C1, 3, boundary );
    %fri2 = -Kinter*u1 + f1in(1:2*nnodes,1);
    frii = projectField( fri, nodes1, nodes2, 1e-5 );
    %frii2 = projectField( fri2, nodes1, nodes2, 1e-5 );
    %plotOnBound( frii, 1, boundary, nnodes, 1 );
    fr = [ frii; zeros(size(C2,2),1) ]; % Give to fr the right size

    f2 = f2in + fr;
    %plotOnBound( f2, 1, boundary, nnodes, 1 );
    uin2 = K2\f2;
    u2 = uin2(1:2*nnodes,1);
    %plotOnBound( u2, 1, boundary, nnodes, 1 );
    
    % Computation of the error :
    u2ni = projectField(u2, nodes2, nodes1, 1e-5);
    u2n = dirichletRhs(u2ni, 3, C1, boundary);
    %urni = projectField(u, nodes, nodes1, .5);
    %urn = dirichletRhs(urni, 3, C1, boundary);
    error2(iter) = sqrt( (u2n-urn)'*(u2n-urn)/(urn'*urn) );
    residual(iter) = sqrt( (u2n-u1n)'*(u2n-u1n)/sqrt((u2n'*u2n)*(u1n'*u1n)) );
end
%%
sigma1 = stress(u1,E,nu,nodes1,elements,order,1,ntoelem);
sigma2 = stress(u2,E,nu,nodes2,elements,order,1,ntoelem);

plotGMSH({u1,'U_vect';sigma1,'stress'}, elements, nodes1(:,[1,2]), 'field1');
plotGMSH({u2,'U_vect';sigma2,'stress'}, elements, nodes1(:,[1,2]), 'field2');

hold on;
set(gca, 'fontsize', 15);
set(gca,'ylim',[-16 0])
plot(log10(error1),'Color','black')
plot(log10(error2),'Color','blue')
plot(log10(residual),'Color','red')
legend('error1 (log)','error2 (log)','residual (log)')