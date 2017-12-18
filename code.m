%15/03/2016
%Code FEM 2D/3D

close all;
clear all;

addpath(genpath('./tools'))

%% Parameters
%E       = 200000; % MPa : Young modulus
%nu      = 0.3;    % Poisson ratio
%fscalar = 250;    % N.mm-2 : Loading on the plate
%mat = [0, E, nu];
%
%% Boundary conditions
%% first index  : index of the boundary
%% second index : 1=x, 2=y, 3=z
%% third        : value
%% [0,1,value] marks a dirichlet regularization therm on x
%%dirichlet = [1,1,0 ; 1,2,0 ; 1,3,0];
%dirichlet = [0,1,0 ; 0,2,0 ; 0,3,0 ; 0,4,0 ; 0,5,0 ; 0,6,0];
%%dirichlet = [2,1,0 ; 2,2,0 ; 1,3,0 ; 2,3,0];
%neumann   = [ 2,3,fscalar ; 1,3,-fscalar ];
%
%% First, import the mesh
%[ nodes,elements,ntoelem,boundary,order] = readmesh3D( 'meshes/rg3dm/platem_c.msh' );
%%[ nodes,elements,ntoelem,boundary,order] = readmesh3D( 'meshes/tetra10/plate3d.msh' );
%nnodes = size(nodes,1);
%
%% Then, build the stiffness matrix :
%[K,C,nbloq,node2c,c2node] = Krig3D (nodes,elements,mat,order,boundary,dirichlet);
%% The right hand side :
%f  = loading3D(nbloq,nodes,boundary,neumann);
%%f = volumicLoad( 0, nodes, elements, 2, -fscalar );
%%udir = ones( 3*nnodes, 1 );
%%udi = keepField3D( udir, 2, boundary, 3 );
%%f = [ zeros(3*nnodes,1) ; C'*udi ];
%
%% Solve the problem :
%uin = K\f;
%% Extract displacement :
%u = uin(1:3*nnodes,1);
%ui = reshape(u,3,[])';  ux = ui(:,1);  uy = ui(:,2);  uz = ui(:,3);
%
%% Compute stress :
%sigma = stress3D(u,mat,nodes,elements,order,1,ntoelem);
%
%% Output :
%% plotNodes(u,elements,nodes); : TODO debug on Matlab r>2013
%plotGMSH3D({ux,'U_x';uy,'U_y';uz,'U_z';u,'U_vect';sigma,'stress'}, elements, nodes, 'output/linear_field');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
%E       = 210000;  % MPa : Young modulus
%nu      = 0.3;     % Poisson ratio
%fscalar = 250;     % N.mm-2 : Loading on the plate
%rho     = 7500e-9; % kg.mm-3 : volumic mass
%%mat = [2, E, nu, 0.1, 1];
%mat     = [0, E, nu];
%dt      = 2e-5;      % s : time discrteization parameter
%
%% Boundary conditions
%% first index  : index of the boundary
%% second index : 1=x, 2=y
%% third        : value
%% [0,1,value] marks a dirichlet regularization therm on x
%dirichlet = [ 0,1,0 ; 1,2,0 ];
%neumann   = [ 3,2,fscalar ];
%
%% First, import the mesh
%[ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/plate.msh' );
%%[ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/t6/plate.msh' );
%nnodes = size(nodes,1);
%
%% Then, build the stiffness matrix :
%[K,C,nbloq,node2c,c2node] = Krig2 (nodes,elements,mat,order,boundary,dirichlet);
%Kinter = K(1:2*nnodes,1:2*nnodes);
%M = mass_mat(nodes, elements);
%M = rho*M;
%M = [ M , zeros(2*nnodes,nbloq) ; zeros(2*nnodes,nbloq)' , zeros(nbloq) ];
%C = zeros(size(M));
%%[ node2b, b2node ] = mapBound( 1, boundary, nnodes );
%% The right hand side :
%f  = loading(nbloq,nodes,boundary,neumann);
%T = 1:1:25;  fa = f*T/T(end);
%
%u0 = zeros(2*nnodes+nbloq,1); v0 = zeros(2*nnodes+nbloq,1); a0 = zeros(2*nnodes+nbloq,1);
%
%%f = volumicLoad( nbloq, nodes, elements, 2, fscalar );
%%udir = ones( 2*nnodes );
%%udi = keepField( udir, 4, boundary, 2 );
%%f = [ zeros(2*nnodes) ; C'*udi ];
%
%uin = Newmark (M, C, K, fa, u0, v0, a0, dt, .25, .5);
%
%% Extract displacement :
%u = uin(1:2*nnodes,:);
%%ux = zeros(nnodes,size(u,2)); uy = zeros(nnodes,size(u,2));
%ux = u(1:2:end-1,:); uy = u(2:2:end,:); 
%
%% Compute stress :
%sigma = stress(u,E,nu,nodes,elements,order,1,ntoelem);
%
%% Output :
%% plotNodes(u,elements,nodes); : TODO debug on Matlab r>2013
%plotGMSH({ux,'U_x';uy,'U_y';u,'U_vect';sigma,'stress'}, elements, nodes, 'output/linear_field');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
E       = 210000;  % MPa : Young modulus
nu      = 0.3;     % Poisson ratio
fscalar = 250;     % N.mm-2 : Loading on the plate
rho     = 7500e-9; % kg.mm-3 : volumic mass
%mat = [2, E, nu, 0.1, 1];
mat     = [0, E, nu];
omega   = 1e10;      % s-2 : square of the pulsation

% Boundary conditions
% first index  : index of the boundary
% second index : 1=x, 2=y
% third        : value
% [0,1,value] marks a dirichlet regularization therm on x
dirichlet = [ 0,1,0 ; 0,2,0 ; 0,3,0 ];
neumann   = [ 1,2,-fscalar ; 3,2,fscalar ];

% First, import the mesh
%[ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/cvg_mesh/plate02.msh' );
%[ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/t6/plate.msh' );
%[ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/plateer.msh' );
[ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/plate_hole_smlsml.msh' );
nnodes = size(nodes,1);

% Then, build the stiffness matrix :
[K,C,nbloq,node2c,c2node] = Krig2 (nodes,elements,mat,order,boundary,dirichlet);
Kinter = K(1:2*nnodes,1:2*nnodes);
M = mass_mat(nodes, elements);
M = rho*M;
M = [ M , zeros(2*nnodes,nbloq) ; zeros(2*nnodes,nbloq)' , zeros(nbloq) ];
%[ node2b, b2node ] = mapBound( 1, boundary, nnodes );
% The right hand side :
f  = loading(nbloq,nodes,boundary,neumann);
%f = volumicLoad( nbloq, nodes, elements, 2, fscalar );
%udir = ones( 2*nnodes );
%udi = keepField( udir, 4, boundary, 2 );
%f = [ zeros(2*nnodes) ; C'*udi ];

%% Cause corner pauses problem (as usual)
%boundary = suppressBound( boundary, [1], 4 );

% Solve the problem :
uin = (K-omega*M)\f;
uin = K\f;

% Extract displacement :
u = uin(1:2*nnodes,1);
ui = reshape(u,2,[])';  ux = ui(:,1);  uy = ui(:,2);

% Compute stress :
sigma = stress(u,E,nu,nodes,elements,order,1,ntoelem);

% Output :
% plotNodes(u,elements,nodes); : TODO debug on Matlab r>2013
plotGMSH({ux,'U_x';uy,'U_y';u,'U_vect';sigma,'stress'}, elements, nodes, 'output/linear_field');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Parameters
%E       = 200000; % MPa : Young modulus
%nu      = 0.3;    % Poisson ratio
%fscalar = 250;    % N.mm-1 : Loading on the plate
%
%% Boundary conditions
%% first index  : index of the boundary
%% second index : 1=x, 2=y
%% third        : value
%% [0,1,value] marks a dirichlet regularization therm on x
%%dirichlet = [0,1,0 ; 0,2,0 ; 0,3,0];
%%dirichlet = [1,1,0,0 ; 1,2,0,0 ; 4,1,-.0005,1];
%dirichlet = [1,1,0,0 ; 1,2,0,0 ; 3,2,0,-1];
%neumann   = [4,1,fscalar];
%
%% First, import the mesh
%[ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/plate.msh' );
%nnodes = size(nodes,1);
%
%%[node2to, to2node] = mapBound(2, boundary, nnodes);
%
%% Then, build the stiffness matrix :
%%[K,C,nbloq,node2c,c2node] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
%
%% The right hand side :
%%f  = loading(nbloq,nodes,boundary,neumann);
%f  = loading(0,nodes,boundary,neumann);
%%f = volumicLoad( 0, nodes, elements, 2, -fscalar );
%
%mat = [0, E, nu];
%eta = .3;
%
%% Cause corner pauses problem (as usual)
%boundary = suppressBound( boundary, [1], 4 );
%
%% Solve the problem :
%uin = statusContact( nodes, elements, mat, 1, boundary, dirichlet, eta, 10, ntoelem, f, 1 );
%
%%uin = K\f;
%% uin = uin1+uin2;
%
%% Extract displacement :
%u = uin(1:2*nnodes,1);
%ui = reshape(u,2,[])';  ux = ui(:,1);  uy = ui(:,2);
%
%% Compute stress :
%sigma = stress(u,E,nu,nodes,elements,order,1,ntoelem);
%
%% Output :
%% plotNodes(u,elements,nodes); : TODO debug on Matlab r>2013
%plotGMSH({ux,'U_x';uy,'U_y';u,'U_vect';sigma,'stress'}, elements, nodes, 'linear_field');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Boundary conditions
% % first index  : index of the boundary
% % second index : 1=x, 2=y
% % third        : value
% % [0,1,value] marks a dirichlet regularization therm on x
% dirichlet = [0,1,0;1,2,0];
%              %3,1,0;3,2,0];
% neumann   = [3,2,1];
% 
% % First, import the mesh
% [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/plate.msh' );
% nnodes = size(nodes,1);
% 
% [node2to, to2node] = mapBound(2, boundary, nnodes);
% 
% % Then, build the stiffness matrix :
% [K,C,nbloq,node2c,c2node] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
% 
% % The right hand side :
% udi = zeros( nnodes, 1 );
% ind = 2:2:2*nnodes;
% udi(ind,1) = 1e-4;
% %f = dirichletRhs2( udi, 3, c2node, boundary, nnodes );
% f = loading(nbloq,nodes,boundary,neumann);
% %farray = [1,0,0,0;1,0,0,0;1,0,0,0;1,0,0,0]; %f(x,y) = 1+y+y^2+y^3
% %f = volumicLoad( nbloq, nodes, elements, 2, 1 );
% 
% % Solve the problem :
% uin = K\f;
% 
% % Extract displacement :
% u = uin(1:2*nnodes,1);
% 
% % Compute stress :
% sigma = stress(u,E,nu,nodes,elements,order,1,ntoelem);
% 
% % Output :
% % plotNodes(u,elements,nodes); : TODO debug on Matlab r>2013
% plotGMSH({u,'U_Vect';sigma,'stress'}, elements, nodes(:,[1,2]), 'linear_field');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Non linear problem : hyperelastic
% alpha = 1e10;   % Nonlinearity parameter
% 
% % Solve the non-linear problem
% uin = solveHyp(nodes,elements,E,nu,order,boundary,dirichlet,...
%     alpha,1e-5,50,ntoelem,0,f);
% 
% % Extract displacement :
% unl = uin(1:2*nnodes,1);
% 
% % Compute stress :
% sigmanl = stressHyp(unl,E,nu,nodes,elements,order,1,ntoelem,alpha,1);
% 
% % Output :
% % plotNodes(u,elements,nodes); : TODO debug on Matlab r>2013
% plotGMSH({unl,'U_Vect';sigmanl,'stress'}, elements, nodes, 'NL_field');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visco-elastic problem : Kelvin-Voigt+instant
% tau   = 0.1;        % Caracteristic time
% alpha = 0;          % Parallel stifness
% T     = 0.2;        % Total duration
% nt    = 20;         % number of time intervals
% nts2  = floor(nt/2);  % The half of nt
% 
% global dotu
% global dotf
% 
% % assembly f
% time = 0:1:nt;
% fass = f * [0:1:nts2-1 , (nts2-1)*ones(1,nt-nts2+1)] / (nts2-1);
% %fass = f * sin(2*pi*time/nt);
% %fass(:,1) = 0;
% 
% % Solve the evolution problem
% uinv = solveKVI(nodes,elements,E,nu,order,boundary,dirichlet,...
%     tau,alpha,ntoelem,fass,T,nt);
% 
% % Extract displacement :
% uv = uinv(1:2*nnodes,:);
% 
% ind = 8;%68
% ind = size(uin,1)-1;
% hold on
% %plot( uin(ind)*ones(1,nt+1)*(1+alpha)/alpha, 'Color', 'red' );
% plot( uin(ind)*ones(1,nt+1), 'Color', 'red' );
% plot( uinv(ind,:) );
% %plot( fass(end,:), 'Color', 'black' );
% plot( dotu(ind,:)/10, 'Color', 'black' );
% %plot( fass(68,:)/E*30, 'Color', 'black' );
% 
% % C2la MR2 KI PU
% 
% % Compute stress /!\ it's false right now because of the time-dependant stiffness :
% sigmav = zeros(3*nnodes,nt+1);
% for i=1:nt+1
%     sigmav(:,i) = stress(uv(:,i),E,nu,nodes,elements,order,1,ntoelem);
% end
% 
% % Output :
% % plotNodes(u,elements,nodes); : TODO debug on Matlab r>2013
% plotGMSH({uv,'U_Vect';sigmav,'stress'}, elements, nodes, 'viscosity_field');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
%E       = 210000; % MPa : Young modulus
%nu      = 0.3;    % Poisson ratio
%fscalar = 250;    % N.mm-2 : Loading on the plate
%mat = [0, E, nu];
%
%% Boundary conditions
%% first index  : index of the boundary
%% second index : 1=x, 2=y
%% third        : value
%% [0,1,value] marks a dirichlet regularization therm on x
%dirichlet = [ 3,1,0 ; 3,2,0 ];
%%dirichlet = [ 0,3,0 ];
%neumann   = [3,2,fscalar];
%
%% First, import the mesh
%[ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/plate_hole.msh' );
%nnodes = size(nodes,1);
%
%% Then, build the stiffness matrix :
%[K,C,nbloq,node2c,c2node] = Krig2 (nodes,elements,mat,order,boundary,dirichlet);
%
%% RHS 
%udi = zeros( nnodes, 1 );
%ind = 2:2:2*nnodes;
%udi(ind,1) = 3e-3;
%f = dirichletRhs2( udi, 3, c2node, boundary, nnodes );
%%f = loading(nbloq,nodes,boundary,neumann);
%
%% Center of the pion / Radius
%Rc = .99;
%xc = 0; yc = Rc-1;
%
%% Build Cc
%[ node2b, b2node ] = mapBound( 5, boundary, nnodes );
%Cc = zeros( 2*nnodes+nbloq, size(b2node) );
%b = zeros( size(b2node), 1 ); % RHS
%for i=1:size(b2node)
%   no = b2node(i);
%   x = nodes(no,1); y = nodes(no,2);
%   Dis = sqrt( (x-xc)^2 + (y-yc)^2 );
%   Cc(2*no-1,i) = -(x-xc)/Dis; Cc(2*no,i) = -(y-yc)/Dis;  % U.N <= b
%   b(i) = Dis - Rc;
%end
%
%% Solve the problem :
%%uin = K\f;
%uin = statusIneq( K, Cc, b, f, 10 );
%%Ktot = [ K, Cc ; Cc', zeros(size(b2node,1)) ]; ftot = [f;b];
%%uin = Ktot \ ftot;
%
%% Extract displacement :
%u = uin(1:2*nnodes,1);
%ui = reshape(u,2,[])';  ux = ui(:,1);  uy = ui(:,2);
%
%% Compute stress :
%sigma = stress(u,E,nu,nodes,elements,order,1,ntoelem);
%
%% Output :
%plotGMSH({ux,'U_x';uy,'U_y';u,'U_vect';sigma,'stress'}, elements, nodes, 'field');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
%E       = 200000; % MPa : Young modulus
%nu      = 0.3;    % Poisson ratio
%Slim    = 250;    % MPa : elasticity limit
%alpha   = 30000;  % MPa : Isotropic hardening coefficient 1000000
%H       = 0000;  % MPa : Cinematic hardening (not implemented
%fscalar = 330;    % N.mm-2 : Loading on the plate
%mat = [10, E, nu, Slim, alpha, H]; % Elasto-plastic with linear isotropic hardening

% Boundary conditions
% first index  : index of the boundary
% second index : 1=x, 2=y, 3=z
% third        : value
% [0,1,value] marks a dirichlet regularization therm on x
%dirichlet = [ 1,1,0 ; 1,2,0 ];
%dirichlet = [ 0,1,0 ; 1,2,0 ];
%dirichlet = [ 1,1,0 ; 1,2,0 ];% 2,1,0 ; 4,1,0 ];
%dirichlet = [ 1,1,0 ; 1,2,0 ; 3,1,0 ; 3,2,0 ];
%%dirichlet = [ 0,1,0 ; 1,2,0 ; 3,2,0 ];% 2,1,0 ; 4,1,0 ];
%neumann   = [ 3,2,fscalar ];
%
%% First, import the mesh
%[ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/plate.msh' );
%nnodes = size(nodes,1); nelem = size(elements,1);
%
%% Then, build the stiffness matrix :
%[K,C,nbloq,node2c,c2node] = Krig2 (nodes,elements,[0, E, nu],order,boundary,dirichlet,1);
%% The right hand side :
%%f1 = loading(nbloq,nodes,boundary,neumann);
%f1 = volumicLoad( nbloq, nodes, elements, 2, 30 );
%udi = zeros( nnodes, 1 );
%ind = 2:2:2*nnodes;
%udi(ind,1) = 0.02;
%udi = keepField( udi, 3, boundary );
%f1 = [ zeros(2*nnodes,1) ; C'*udi ];
%
%%% Proportionnal loading
%%T = 1:1:10;
%%f = f1*T/T(end);
%
%% Cyclic loading
%T = 1:1:25;
%fa = f1*T/T(end); fb = f1*(1-T/T(end)); f = fa;
%T = 1:1:75; f = [fa,fb,-fa];
%%f = [fa,fb,-fa,-fb];%,fa];%,fb,-fa,-fb,fa];
%
%% Solve the problem :
%%uin = K\f;
%[uin,sigm,pe,lam] = Elastoplast( mat, K, f, T, nodes, elements, order, 1e-6, 10, 1, 0 );
%% Extract displacement :
%u = uin(1:2*nnodes,end);
%ui = reshape(u,2,[])';  ux = ui(:,1);  uy = ui(:,2);
%
%% pass sigma on the nodes
%Red = zeros(nnodes,nelem);
%for i=1:size(elements,1)
%   Xloc = nodes(elements(i,:),:);    % Extract coords
%   nno = size(Xloc,1);
%
%   ne = size(elements,2);
%   for j=1:ne
%      nod = elements(i,j);
%      Red( nod, i ) = Red( nod, i ) + 1/ntoelem(nod,1);
%   end
%end
%sigma(1:3:3*nnodes-2,:) = Red*sigm(1:3:3*nelem-2,:);
%sigma(2:3:3*nnodes-1,:) = Red*sigm(2:3:3*nelem-1,:);
%sigma(3:3:3*nnodes,:)   = Red*sigm(3:3:3*nelem,:);
%p                       = Red*pe;
%
%sigma = sigma(:,end); p = p(:,end);
%
%% Output :
%plotGMSH({ux,'U_x';uy,'U_y';u,'U_vect';sigma,'stress';p,'cumulated plasticity'},...
%          elements, nodes, 'output/elasto_plastic');
%          
%% Plot traction curve at node 3
%f = -[ lam(:,2:end) ; zeros( size(uin,1)-size(u,1) , size(T,2) ) ];
%up = uin(6,1:end); fp = [0,f(6,:)];
%plot(up, fp);
