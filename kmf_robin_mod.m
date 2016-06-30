%03/05/2016
% Algo KMF avec CL de Robin à plus de 2 pas

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E         = 70000;      % MPa : Young modulus
nu        = 0.3;    % Poisson ratio
fscalar   = 1;      % N.mm-1 : Loading on the plate
niter     = 100;
br        = 0.;     % noise
robin1    = [1e5;1e-5]*E;
robin2    = [1]*E;     % Robin parameters

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Obtain an approxiamtion of the Schur complements

xmax = max(nodes(:,1));
xmin = min(nodes(:,1));
ymax = max(nodes(:,2));
no3  = findNode(xmax, ymax, nodes, 1e-5);
no4  = findNode(xmin, ymax, nodes, 1e-5);
boundaryd = suppressBound( boundary, [no3;no4], 3 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init : get null vectors
uprec = uref-uref;
u2    = uprec;
fprec = uprec;
f2    = fprec;

% problem 1
dirichlet1 = [4,1,0;4,2,0];
neumann1   = [1,2,fscalar;
              2,1,fscalar];
[K1,C1,nbloq1,node2c1,c2node1] = Krig (nodes,elements,E,nu,order,boundary,dirichlet1);

f1in = loading(nbloq1,nodes,boundary,neumann1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fixed point algorithm

error1   = zeros(niter,1);
error2   = zeros(niter,1);
residual = zeros(niter,1);
regulari = zeros(niter,1);
for iter = 1:niter
    % Solve 1
    for i=1:size(robin1)
        unsrob = E^2/robin1(i); % Robin parameter of the known boundary
        kuplusf = uprec + fprec/robin1(i);
        kuplusf1 = urefb + f1in(1:2*nnodes,1)/unsrob;
        [Kp1, fro1] = robinRHS( nbloq1, nodes, boundary, kuplusf1, unsrob, 1 );
        [Kp2, fro2] = robinRHS( nbloq1, nodes, boundary, kuplusf1, unsrob, 2 );
        [Kp3, fro3] = robinRHS( nbloq1, nodes, boundary, kuplusf, robin1(i), 3 );
        f1 = fro1 + fro2 + fro3;
        K1e = K1 + Kp1 + Kp2 + Kp3;

        uin1 = K1e\f1;
        uprec = uin1(1:2*nnodes,1);
        fprec = Kinter*uprec - f1in(1:2*nnodes,1); % Get the forces at the boundary
    end
    u1 = uprec;
    f1 = fprec;
    error1(iter) = sqrt( myps(u1-uref,u1-uref,Kinter,boundary,M,nodes)/...
        myps(uref,uref,Kinter,boundary,M,nodes) );
    
    % Solve 2
    for i=1:size(robin2)
        unsrob = E^2/robin2(i); % Robin parameter of the known boundary
        kuplusf = uprec + fprec/robin2(i);
        kuplusf1 = urefb + f1in(1:2*nnodes,1)/unsrob;
        [Kp1, fro1] = robinRHS( nbloq1, nodes, boundary, kuplusf1, unsrob, 1 );
        [Kp2, fro2] = robinRHS( nbloq1, nodes, boundary, kuplusf1, unsrob, 2 );
        [Kp3, fro3] = robinRHS( nbloq1, nodes, boundary, kuplusf, robin2(i), 3 );
        f2 = fro1 + fro2 + fro3;
        K2e = K1 + Kp1 + Kp2 + Kp3;

        uin2 = K2e\f2;
        uprec = uin2(1:2*nnodes,1);
        fprec = Kinter*uprec - f1in(1:2*nnodes,1); % Get the forces at the boundary
    end
    u2 = uprec;
    f2 = fprec;
    error2(iter) = sqrt( myps(u2-uref,u2-uref,Kinter,boundary,M,nodes)/...
        myps(uref,uref,Kinter,boundary,M,nodes) );
    
    residual(iter) = sqrt( myps(u1-u2,u1-u2,Kinter,boundary,M,nodes)/...
                     sqrt(myps(u1,u1,Kinter,boundary,M,nodes)*...
                     myps(u2,u2,Kinter,boundary,M,nodes)) );
                 
    regulari(iter) = sqrt(uprec'*regul(uprec, nodes, boundary, 3));
end

hold on;
plot(log10(error1),'Color','black')
plot(log10(error2),'Color','blue')
plot(log10(residual),'Color','red')
figure;
% L-curve
loglog(residual,regulari);