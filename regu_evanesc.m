% 24/03/2016
% Algo de régularisation evanescente

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 70000;  % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 1;      % N.mm-1 : Loading on the plate
niter   = 50;
mu      = .1;    % Regularization parameter
br      = .0;     % noise

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

% find the nodes in the corners and suppress the element :
xmax = max(nodes(:,1));
xmin = min(nodes(:,1));
ymax = max(nodes(:,2));
ymin = min(nodes(:,2));
no1  = findNode(xmin, ymin, nodes, 1e-5);
no2  = findNode(xmax, ymin, nodes, 1e-5);
no3  = findNode(xmax, ymax, nodes, 1e-5);
no4  = findNode(xmin, ymax, nodes, 1e-5);

% Suppress some nodes from the boundaries
boundaryp1 = suppressBound( boundary, [no3;no4], 3 );
boundaryp1 = suppressBound( boundaryp1, [no1;no2], 1 );

% Then, build the stiffness matrix :
[K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
Kinter = K(1:2*nnodes, 1:2*nnodes);
M      = mass_mat(nodes, elements);

% Some indices
[node2b4, b2node4] = mapBound(4, boundaryp1, nnodes);
[node2b3, b2node3] = mapBound(3, boundaryp1, nnodes);
[node2b1, b2node1] = mapBound(1, boundaryp1, nnodes);
[node2b2, b2node2] = mapBound(2, boundaryp1, nnodes);
b2node12  = [b2node1;b2node2];
b2node123 = [b2node1;b2node2;b2node3];    % nodes
bbound    = [2*b2node123-1; 2*b2node123]; % dof
nbound    = size(bbound,1);

% The right hand side :
f = loading(nbloq,nodes,boundary,neumann);

% Solve the problem :
uin = K\f;

% Extract displacement and Lagrange multiplicators :
uref = uin(1:2*nnodes,1);
lagr = uin(2*nnodes+1:end,1);

urefb = ( 1 + br*randn(2*nnodes,1) ) .* uref;
frefb  = f(1:2*nnodes);

sigma = stress(uref,E,nu,nodes,elements,order,1,ntoelem);
plotGMSH({uref,'Vect_U';sigma,'stress'}, elements, nodes, 'reference');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computation of the inner stiffness
dirichlet1 = [4,1,0;4,2,0];
[K1,C1,nbloq1,node2c1,c2node1] = Krig (nodes,elements,E,nu,order,boundary,dirichlet1);
neumann1   = []; % There is no alone Neumann
f1 = loading( nbloq, nodes, boundary, neumann1 );

%% Schur operator
[ S, b, map ] = schurComp2( Kinter, f(1:2*nnodes), b2node123 );

error    = zeros(niter,1);
residual = zeros(niter,1);
regulari = zeros(niter,1);

%% Mass matrices
Mr  = bMass_mat(nodes, boundary, [2;1]);
Mrt = 1/E*Mr;  % Ideally 1/EL
M   = bMass_mat(nodes, boundary, [3;2;1]);
Mt  = 1/E*M;
% Debug
%M   = eye(2*nnodes);  M([2*b2node3-1,2*b2node3],[2*b2node3-1,2*b2node3]) = 0;
%Mr  = eye(2*nnodes);
%Mt  = M;
%Mrt = Mr;
% Extract coords
Mr  = Mr(bbound, bbound);
Mrt = Mrt(bbound, bbound);
M   = M(bbound, bbound);
Mt  = Mt(bbound, bbound);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Evanescent regularization method
Itere  = zeros( 2*nnodes, 1 );
Iteref = zeros( 2*nnodes, 1 );

% Compute errors
error(1)    = (Itere(bbound)-uref(bbound))'*M*...
   (Itere(bbound)-uref(bbound)) / (uref(bbound)'*M*uref(bbound));
residual(1) = (Itere(bbound)-urefb(bbound))'*Mr*...
   (Itere(bbound)-urefb(bbound)) / (uref(bbound)'*Mr*uref(bbound));
regulari(1) = (Itere(bbound)-Itere(bbound))'*M*...
   (Itere(bbound)-Itere(bbound)) / (uref(bbound)'*M*uref(bbound));

% Build the fat matrix
Atot = [Mr+mu*M, zeros(size(M)), S'
        zeros(size(M)), Mrt+mu*Mt, -eye(size(M,1),size(S,1))
        S, -eye(size(S,1), size(M,1)), zeros(size(M))];
        
disp( [ 'Log of the cond of the problem :', num2str(log10(cond(Atot))) ] )

plot(log10(abs(eig(Atot))));
% plot(eig(Atot));
figure;

for i = 2:niter

   % Rhs
   btot = [Mr*urefb(bbound) + mu*M*Itere(bbound)
           Mrt*frefb(bbound) + mu*Mt*Iteref(bbound)
           f1(bbound)]; % Don't forget to add Kinter*uimp if needed

   % Solve and extract the relevant parts
   Iterep = Itere;
   xtot   = Atot\btot;
   Itere(bbound)  = xtot(1:nbound);
   Iteref(bbound) = xtot(nbound+1:2*nbound);
   
   % Compute errors
   error(i)    = (Itere(bbound)-uref(bbound))'*M*...
      (Itere(bbound)-uref(bbound)) / (uref(bbound)'*M*uref(bbound));
   residual(i) = (Itere(bbound)-urefb(bbound))'*Mr*...
      (Itere(bbound)-urefb(bbound)) / (uref(bbound)'*Mr*uref(bbound));
   regulari(i) = (Itere(bbound)-Iterep(bbound))'*M*...
      (Itere(bbound)-Iterep(bbound)) / (uref(bbound)'*M*uref(bbound));
   
end

hold on
plot(log10(error),'Color','blue')
plot(log10(residual),'Color','red')
legend('error (log)','residual (log)')
figure;
% L-curve :
loglog(residual,regulari);
figure
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
fdir3 = dirichletRhs(Itere, 3, C, boundary);
usoli = K \ assembleDirichlet( [fdir1+fdir3,fdir2] );

usol = usoli(1:2*nnodes,1);
fsol = Kinter*usol;

hold on;
plot(frefb(2*b2node3), 'Color', 'red')
plot(fsol(2*b2node3), 'Color', 'blue')
figure;
hold on;
plot(uref(2*b2node3), 'Color', 'red')
plot(usol(2*b2node3), 'Color', 'blue')

total_error = norm(usol-uref)/norm(uref);
% Compute stress :
sigma = stress(usol,E,nu,nodes,elements,order,1,ntoelem);
plotGMSH({usol,'U_vect';sigma,'stress'}, elements, nodes, 'solution');