%15/03/2016
%Code FEM 2D contraintes planes

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 200000; % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 250;    % N.mm-1 : Loading on the plate

% Boundary conditions
% first index  : index of the boundary
% second index : 1=x, 2=y
% third        : value
% [0,1,value] marks a dirichlet regularization therm on x
dirichlet = [0,1,0 ; 0,2,0 ; 0,3,0];
neumann1  = [3,2,1];
neumann2  = [1,2,-1];
neumann   = [3,2,1 ; 1,2,-1];

% First, import the mesh
[ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/plate.msh' );
nnodes = size(nodes,1);

[node2to, to2node] = mapBound(2, boundary, nnodes);

% Then, build the stiffness matrix :
[K,C,nbloq,node2c,c2node] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);

% The right hand side :
% f1 = loading(nbloq,nodes,boundary,neumann1);
% f2 = loading(nbloq,nodes,boundary,neumann2);
f  = loading(nbloq,nodes,boundary,neumann);

% % Regularize the problem :
% g1 = zeros(2*nnodes,1); g2 = g1; g3p = g1;
% ind = 2:2:2*nnodes;
% g1(ind,1) = 1; g2(ind-1,1) = 1;
% %g1 = keepField( g1, 1, boundary );
% %g2 = keepField( g2, 1, boundary );
% g3p(ind,1) = nodes(ind/2,1);
% g3p(ind-1,1) = -nodes(ind/2,2);
% %g3p = keepField( g3p, 1, boundary );
% % orthogonalize
% g3 = g3p - (g3p'*g1)/(g1'*g1)*g1 - (g3p'*g2)/(g2'*g2)*g2;
% 
% G = [g1,g2,g3];
% % G = null(full(K));
% Kr = [K,G;G',zeros(size(G,2))];

% Solve the problem :
% uin1 = Kr\[f1;zeros(size(G,2),1)];
% uin2 = Kr\[f2;zeros(size(G,2),1)];
uin = K\f;
% uin = uin1+uin2;

% Extract displacement :
u = uin(1:2*nnodes,1);
ui = reshape(u,2,[])';  ux = ui(:,1);  uy = ui(:,2);

% Compute stress :
sigma = stress(u,E,nu,nodes,elements,order,1,ntoelem);

% Output :
% plotNodes(u,elements,nodes); : TODO debug on Matlab r>2013
plotGMSH({ux,'U_x';uy,'U_y';u,'U_vect';sigma,'stress'}, elements, nodes, 'linear_field');

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