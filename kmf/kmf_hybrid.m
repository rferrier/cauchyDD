%03/05/2016
%Algo KMF hybride

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E        = 70000;  % MPa : Young modulus
nu       = 0.3;    % Poisson ratio
fscalar  = 1;      % N.mm-1 : Loading on the plate
niter    = 10;
br       = 0.;     % noise
hybr_pb1 = 1;      % use hybrid conditions on the suraboundant boundary
hybr_pb2 = 1;
hybr_so1 = 1;      % use hybrid conditions on the unknown boundary
hybr_so2 = 1;

% Boundary conditions
% first index  : index of the boundary
% second index : 1=x, 2=y
% third        : value
% [0,1,value] marks a dirichlet regularization therm on x
dirichlet = [4,1,0;4,2,0];

neumann   = [1,2,fscalar;
             2,1,fscalar;
             3,2,fscalar];

% Import the mesh
[ nodes,elements,ntoelem,boundary,order ] = readmesh( 'meshes/plate.msh' );
nnodes = size(nodes,1);

% Then, build the stiffness matrix :
[K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
Kinter = K(1:2*nnodes, 1:2*nnodes);
M      = mass_mat(nodes, elements);

% The right hand side :
f = loading(nbloq,nodes,boundary,neumann);

% Solve the problem :
uin = K\f;

% Extract displacement and Lagrange multiplicators :
uref = uin(1:2*nnodes,1);
lagr = uin(2*nnodes+1:end,1);
urefb = ( 1 + br*randn(2*nnodes,1) ) .* uref;
lagrb = ( 1 + br*randn(nbloq,1) ) .* lagr;
fref = f( 1:2*nnodes,1 ); % Reaction forces

% Compute stress :
sigma = stress(uref,E,nu,nodes,elements,order,1,ntoelem);
% Output :
plotGMSH({uref,'U_vect';sigma,'stress'}, elements, nodes(:,[1,2]), 'reference');
% patch('Faces',elements(:,1:3),'Vertices',nodes,'FaceAlpha',0);
% axis equal
% figure

% Suppress redoundant corners
% xmax = max(nodes(:,1));
% xmin = min(nodes(:,1));
% ymax = max(nodes(:,2));
% no3  = findNode(xmax, ymax, nodes, 1e-5);
% no4  = findNode(xmin, ymax, nodes, 1e-5);
% boundaryd = suppressBound( boundary, [no3;no4], 3 );
boundaryd = boundary;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init :
u1    = uref-uref;
u2    = u1;
fri2   = u1;

% DN problem
if hybr_pb1 == 1 && hybr_so1 == 1
    dirichlet1 = [4,1,0;4,2,0;
                  1,1,0;
                  2,1,0;
                  3,2,0];
    neumann1   = [1,2,fscalar];
elseif hybr_pb1 == 1 && hybr_so1 == 0
    dirichlet1 = [4,1,0;4,2,0;
                  1,1,0;
                  2,1,0;
                  3,1,0;3,2,0];
    neumann1   = [1,2,fscalar];
elseif hybr_pb1 == 0 && hybr_so1 == 1
    dirichlet1 = [4,1,0;4,2,0;
                  3,2,0];
    neumann1   = [1,2,fscalar;
                  2,1,fscalar];
else % Everybody is 0
    dirichlet1 = [4,1,0;4,2,0;
                  3,1,0;3,2,0];
    neumann1   = [1,2,fscalar;
                  2,1,fscalar];
end

[K1,C1,nbloq1,node2c1,c2node1] = Krig (nodes,elements,E,nu,order,boundaryd,dirichlet1);

% ND problem
if hybr_pb2 == 1 && hybr_so2 == 1
    dirichlet2 = [4,1,0;4,2,0;
                  1,2,0;
                  2,2,0;
                  3,1,0];
    neumann2   = [2,1,fscalar];
elseif hybr_pb2 == 1 && hybr_so2 == 0
    dirichlet2 = [4,1,0;4,2,0;
                  1,1,0;
                  2,1,0];
    neumann2   = [2,1,fscalar];
elseif hybr_pb2 == 0 && hybr_so2 == 1
    dirichlet2 = [4,1,0;4,2,0;
                  1,1,0;1,2,0;
                  2,1,0;2,2,0;
                  3,1,0];
    neumann2   = [];
else % Everybody is 0
    dirichlet2 = [4,1,0;4,2,0;
                  1,1,0;1,2,0;
                  2,1,0;2,2,0];
    neumann2   = [];
end

[K2,C2,nbloq2,node2c2,c2node2] = Krig (nodes,elements,E,nu,order,boundaryd,dirichlet2);

error1   = zeros(niter,1);
error2   = zeros(niter,1);
residual = zeros(niter,1);
regulari = zeros(niter,1);

for iter = 1:niter
    % Solve DN
    fr1 = [ fri2; zeros(size(C1,2),1) ]; % Give to fr the right size
    fdir11 = dirichletRhs2(urefb, 1, c2node1, boundaryd, nnodes);
    fdir21 = dirichletRhs2(urefb, 2, c2node1, boundaryd, nnodes);
    f1in = loading(nbloq1,nodes,boundaryd,neumann1);
    fdir31 = dirichletRhs2(u2, 3, c2node1, boundaryd, nnodes );
    f1 = f1in + fr1 + assembleDirichlet( [fdir11,fdir31,fdir21] );

    uin1 = K1\f1;
    u1 = uin1(1:2*nnodes,1);
    %lagr1 = uin1(2*nnodes+1:end,1);
    fri1 = keepField( Kinter*u1-f1in(1:2*nnodes), 3, boundaryd );
    %fri1 = lagr2forces2( lagr1, c2node1, 3, boundary, nnodes );
    
    error1(iter) = sqrt( myps(u1-uref,u1-uref,Kinter,boundaryd,M,nodes)/...
        myps(uref,uref,Kinter,boundaryd,M,nodes) );
    
    % Solve ND
    fr2 = [ fri1; zeros(size(C2,2),1) ]; % Give to fr the right size
    fdir12 = dirichletRhs2(urefb, 1, c2node2, boundaryd, nnodes);
    fdir22 = dirichletRhs2(urefb, 2, c2node2, boundaryd, nnodes);
    f2in = loading(nbloq2,nodes,boundaryd,neumann2);
    fdir32  = dirichletRhs2(u1, 3, c2node2, boundaryd, nnodes);
    % Assembly the dirichlet RHS (ensure a redundant BC isn't imposed 2
    % times)
    f2 = f2in + fr2 + assembleDirichlet( [fdir12,fdir32,fdir22] );

    uin2 = K2\f2;
    u2 = uin2(1:2*nnodes,1);
    %lagr2 = uin1(2*nnodes+1:end,1);
    fri2 = keepField( Kinter*u2-f2in(1:2*nnodes), 3, boundaryd );
    %fri2 = lagr2forces2( lagr2, c2node2, 3, boundary, nnodes );
    
    error2(iter) = sqrt( myps(u2-uref,u2-uref,Kinter,boundary,M,nodes)/...
        myps(uref,uref,Kinter,boundary,M,nodes) );
    
    residual(iter) = sqrt( myps(u1-u2,u1-u2,Kinter,boundary,M,nodes)/...
                     sqrt(myps(u1,u1,Kinter,boundary,M,nodes)*...
                     myps(u2,u2,Kinter,boundary,M,nodes)) );
                 
    regulari(iter) = sqrt(u2'*regul(u2, nodes, boundary, 3));
end

hold on;
plot(log10(error1),'Color','black')
plot(log10(error2),'Color','blue')
plot(log10(residual),'Color','red')
legend('error1 (log)','error2 (log)','residual (log)')
figure
% L-curve
loglog(residual,regulari);
% figure
% plot(regulari.*residual)

% Compute stress :
sigma = stress(u2,E,nu,nodes,elements,order,1,ntoelem);
% Output :
plotGMSH({u2,'U_vect';sigma,'stress'}, elements, nodes(:,[1,2]), 'field2');