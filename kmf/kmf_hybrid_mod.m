%03/05/2016
%Algo KMF hybride à plus de 2 pas (4 en l'occurence (ou ocurrence, j'ai un doute))

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E          = 70000;  % MPa : Young modulus
nu         = 0.3;    % Poisson ratio
fscalar    = 1;      % N.mm-1 : Loading on the plate
niter      = 100;  % 200, pas plus (après, ça remonte)
br         = 0.;     % noise
fullhybrid = 1;      % Wether to use full hybrid problems

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
u1      = uref-uref;
u21     = u1;
fri21   = u1;

% DN problem
if fullhybrid == 1
    dirichlet11 = [4,1,0;4,2,0;
                   1,1,0;
                   2,1,0;
                   3,2,0];
    neumann11   = [1,2,fscalar];
else
    dirichlet11 = [4,1,0;4,2,0;
                   3,2,0];
    neumann11   = [1,2,fscalar;
                   2,1,fscalar];
end

dirichlet10 = [4,1,0;4,2,0;
               3,1,0;3,2,0];
neumann10   = [1,2,fscalar;
               2,1,fscalar];

[K11,C11,nbloq11,node2c11,c2node11] = Krig (nodes,elements,E,nu,order,boundaryd,dirichlet11);
[K10,C10,nbloq10,node2c10,c2node10] = Krig (nodes,elements,E,nu,order,boundaryd,dirichlet10);

% ND problem
if fullhybrid == 1
    dirichlet21 = [4,1,0;4,2,0;
                   1,2,0;
                   2,2,0;
                   3,1,0];
    neumann21   = [2,1,fscalar];
else
    dirichlet21 = [4,1,0;4,2,0;
                   1,1,0;1,2,0;
                   2,1,0;2,2,0;
                   3,1,0];
    neumann21   = [];
end

dirichlet20 = [4,1,0;4,2,0;
               1,1,0;1,2,0;
               2,1,0;2,2,0];
neumann20   = [];

[K21,C21,nbloq21,node2c21,c2node21] = Krig (nodes,elements,E,nu,order,boundaryd,dirichlet21);
[K20,C20,nbloq20,node2c20,c2node20] = Krig (nodes,elements,E,nu,order,boundaryd,dirichlet20);

error1   = zeros(niter,1);
error2   = zeros(niter,1);
residual = zeros(niter,1);
regulari = zeros(niter,1);

for iter = 1:niter
    % Solve DN
%     fri10 = fri21; u10 = u21;
    fr1 = [ fri21; zeros(size(C10,2),1) ]; % Give to fr the right size
    fdir11 = dirichletRhs2(urefb, 1, c2node10, boundaryd, nnodes);
    fdir21 = dirichletRhs2(urefb, 2, c2node10, boundaryd, nnodes);
    fret = fdir21;
    f1in = loading(nbloq10,nodes,boundaryd,neumann10);
    fdir31 = dirichletRhs2(u21, 3, c2node10, boundaryd, nnodes );
    freta = fdir31;
    f1 = f1in + fr1 + assembleDirichlet( [fdir11,fdir31,fdir21] );

    uin1 = K10\f1;
    u10 = uin1(1:2*nnodes,1);
    fri10 = keepField( Kinter*u10-f1in(1:2*nnodes), 3, boundaryd );
    
    error1(iter) = sqrt( myps(u10-uref,u10-uref,Kinter,boundaryd,M,nodes)/...
        myps(uref,uref,Kinter,boundaryd,M,nodes) );
    
    % Solve ND
%     fri20 = fri11; u20 = u11;
    fr2 = [ fri10; zeros(size(C20,2),1) ]; % Give to fr the right size
    fdir12 = dirichletRhs2(urefb, 1, c2node20, boundaryd, nnodes);
    fdir22 = dirichletRhs2(urefb, 2, c2node20, boundaryd, nnodes);
    f2in = loading(nbloq20,nodes,boundaryd,neumann20);
    fdir32  = dirichletRhs2(u10, 3, c2node20, boundaryd, nnodes );
    f2 = f2in + fr2 + assembleDirichlet( [fdir12,fdir32,fdir22] );

    uin2 = K20\f2;
    u20 = uin2(1:2*nnodes,1);
    fri20 = keepField( Kinter*u20-f2in(1:2*nnodes), 3, boundaryd );
    
    error2(iter) = sqrt( myps(u20-uref,u20-uref,Kinter,boundary,M,nodes)/...
        myps(uref,uref,Kinter,boundary,M,nodes) );
    
    residual(iter) = sqrt( myps(u10-u20,u10-u20,Kinter,boundary,M,nodes)/...
                     sqrt(myps(u10,u10,Kinter,boundary,M,nodes)*...
                     myps(u20,u20,Kinter,boundary,M,nodes)) );
                 
    regulari(iter) = sqrt(u20'*regul(u20, nodes, boundary, 3));
    
    % Solve DN Hybrid
%     fri11 = fri10; u11 = u10;
    fr1 = [ fri20; zeros(size(C11,2),1) ]; % Give to fr the right size
    fdir11 = dirichletRhs2(urefb, 1, c2node11, boundaryd, nnodes);
    fdir21 = dirichletRhs2(urefb, 2, c2node11, boundaryd, nnodes);
    f1in = loading(nbloq11,nodes,boundaryd,neumann11);
    fdir31 = dirichletRhs2(u20, 3, c2node11, boundaryd, nnodes );
    f1 = f1in + fr1 + assembleDirichlet( [fdir11,fdir31,fdir21] );

    uin1 = K11\f1;
    u11 = uin1(1:2*nnodes,1);
    fri11 = keepField( Kinter*u11-f1in(1:2*nnodes), 3, boundaryd );
    
%     error1(iter) = sqrt( myps(u11-uref,u11-uref,Kinter,boundaryd,M,nodes)/...
%         myps(uref,uref,Kinter,boundaryd,M,nodes) );

    % Solve ND Hybrid
%     fri21 = fri20; u21 = u20;
    fr2 = [ fri11; zeros(size(C21,2),1) ]; % Give to fr the right size
    fdir12 = dirichletRhs2(urefb, 1, c2node21, boundaryd, nnodes);
    fdir22 = dirichletRhs2(urefb, 2, c2node21, boundaryd, nnodes);
    f2in = loading(nbloq21,nodes,boundaryd,neumann21);
    fdir32  = dirichletRhs2(u11, 3, c2node21, boundaryd, nnodes );
    f2 = f2in + fr2 + assembleDirichlet( [fdir12,fdir32,fdir22] );

    uin2 = K21\f2;
    u21 = uin2(1:2*nnodes,1);
    fri21 = keepField( Kinter*u21-f2in(1:2*nnodes), 3, boundaryd );
%     
%    error2(iter) = sqrt( myps(u21-uref,u21-uref,Kinter,boundary,M,nodes)/...
%         myps(uref,uref,Kinter,boundary,M,nodes) );
%     
%     residual(iter) = sqrt( myps(u10-u21,u10-u21,Kinter,boundary,M,nodes)/...
%                      sqrt(myps(u10,u10,Kinter,boundary,M,nodes)*...
%                      myps(u21,u21,Kinter,boundary,M,nodes)) );
%                  
%     regulari(iter) = sqrt(u21'*regul(u21, nodes, boundary, 3));
    
end

hold on;
plot(log10(error1),'Color','black')
plot(log10(error2),'Color','blue')
plot(log10(residual),'Color','red')
legend('error1','error2','residual')
figure
% L-curve
loglog(residual,regulari);
% figure
% plot(regulari.*residual)

% Compute stress :
sigma = stress(u20,E,nu,nodes,elements,order,1,ntoelem);
% Output :
plotGMSH({u20,'U_vect';sigma,'stress'}, elements, nodes(:,[1,2]), 'field2');