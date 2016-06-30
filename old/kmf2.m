%21/03/2016
%Algo KMF

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 1;      % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 1;      % N.mm-1 : Loading on the plate
niter   = 40;

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
[ nodes,elements,ntoelem,boundary,order ] = readmesh( 'meshes/plate2.msh' );
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

% Extract reaction forces on the known boundaries :
% fr = lagr2forces( lagr, C, [1,2], boundary );
% fr1 = fr(:,1);
% fr2 = fr(:,2);
% Extract displacement
%ur = extrbound(uref, [1,2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init :
u1 = uref-uref;
u2 = u1;

% DN problem
dirichlet1 = [4,1,0;
              4,2,0;
              3,1,0;
              3,2,0];
neumann1   = [1,2,fscalar;
              2,1,fscalar];
[K1,C1,nbloq1,node2c1,c2node1] = Krig (nodes,elements,E,nu,order,boundary,dirichlet1);

% ND problem
dirichlet2 = [4,1,0;
              4,2,0;
              1,1,0;
              1,2,0;
              2,1,0;
              2,2,0];
neumann2   = [];% [3,1,lagr1; 3,2,lagr1]
                % is managed by lagr2forces
[K2,C2,nbloq2,node2c2,c2node2] = Krig (nodes,elements,E,nu,order,boundary,dirichlet2);

error1   = zeros(niter,1);
error2   = zeros(niter,1);
residual = zeros(niter,1);
for iter = 1:niter
    % Solve DN
    f1in = loading(nbloq1,nodes,boundary,neumann1);
    %fdir = dirichletRhs(u2, 3, C1, boundary);
    fdir = dirichletRhs2(u2, 3, c2node1, boundary, nnodes );
    f1 = f1in + fdir;

    uin1 = K1\f1;
    u1 = uin1(1:2*nnodes,1);
    lagr1 = uin1(2*nnodes+1:end,1);
    
    error1(iter) = sqrt( myps(u1-uref,u1-uref,Kinter,boundary,C1)/...
        myps(uref,uref,Kinter,boundary,C1) );
    
    % Solve ND
    %fri = lagr2forces( lagr1, C1, 3, boundary );
    fri = lagr2forces2( lagr1, c2node1, 3, boundary, nnodes );
    fr = [ fri; zeros(size(C2,2),1) ]; % Give to fr the right size
    %fdir1 = dirichletRhs(uref, 1, C2, boundary);
    %fdir2 = dirichletRhs(uref, 2, C2, boundary);
    fdir1 = dirichletRhs2(uref, 1, c2node2, boundary, nnodes);
    fdir2 = dirichletRhs2(uref, 2, c2node2, boundary, nnodes);
    % Assembly the dirichlet RHS (ensure a redundant CL isn't imposed 2
    % times
    f2 = fr + max(fdir1,fdir2);

    uin2 = K2\f2;
    u2 = uin2(1:2*nnodes,1);
    
    error2(iter) = sqrt( myps(u2-uref,u2-uref,Kinter,boundary,C1)/...
        myps(uref,uref,Kinter,boundary,C1) );
    
    residual(iter) = sqrt( myps(u1-u2,u1-u2,Kinter,boundary,C1)/...
                     sqrt(myps(u1,u1,Kinter,boundary,C1)*...
                     myps(u2,u2,Kinter,boundary,C1)) );
end

hold on;
% plot(error1,'Color','black')
% plot(error2,'Color','blue')
% plot(residual,'Color','red')
plot(log10(error1),'Color','black')
plot(log10(error2),'Color','blue')
plot(log10(residual),'Color','red')

% Compute stress :
% sigma = stress(u,E,nu,nodes,elements,order,1,ntoelem);

% Output :
% plotNodes(u,elements,nodes); : TODO debug on Matlab r>2013
% plotGMSH({u,sigma}, elements, nodes(:,[1,2]), 'reference');