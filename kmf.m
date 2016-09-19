%21/03/2016
%Algo KMF

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 70000;  % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 1;      % N.mm-1 : Loading on the plate
niter   = 20;
br      = 0.2;      % noise
relax   = 0;      % Use a relaxation paramter

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

% Extract the index of the boundary
[ node2b3, b2node3 ] = mapBound( 3, boundary, nnodes );
index    = 2*b2node3-1;
index    = index(size(index):-1:1);

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init :
u1    = uref-uref;
u2    = u1;
fri   = u1;
v     = u1;
theta = ones(niter+1,1); % First relaxation parameter

% DN problem
dirichlet1 = [4,1,0;4,2,0;
              3,1,0;3,2,0];
neumann1   = [1,2,fscalar;
              2,1,fscalar];
[K1,C1,nbloq1,node2c1,c2node1] = Krig (nodes,elements,E,nu,order,boundary,dirichlet1);
%[L1,U1] = lu(K1);
% ND problem
dirichlet2 = [4,1,0;4,2,0;
              1,1,0;1,2,0;
              2,1,0;2,2,0];
neumann2   = [];% [3,1,lagr1; 3,2,lagr1]
                % is managed by lagr2forces
[K2,C2,nbloq2,node2c2,c2node2] = Krig (nodes,elements,E,nu,order,boundary,dirichlet2);
%[L2,U2] = lu(K2);

error1   = zeros(niter,1);
error2   = zeros(niter,1);
residual = zeros(niter,1);
regulari = zeros(niter,1);
serror1   = zeros(niter,1); % Error for sigma
serror2   = zeros(niter,1);
sresidual = zeros(niter,1);
ferror1   = zeros(niter,1); % Error for reaction force
ferror2   = zeros(niter,1);
fresidual = zeros(niter,1);

for iter = 1:niter
    % Solve DN
    u1o = u1;
    f1in = loading(nbloq1,nodes,boundary,neumann1);
    %fdir = dirichletRhs(u2, 3, C1, boundary);
    fdir = dirichletRhs2(v, 3, c2node1, boundary, nnodes );
    f1 = f1in + fdir;

    uin1 = K1\f1;
    u1 = uin1(1:2*nnodes,1);
    lagr1 = uin1(2*nnodes+1:end,1);
    frio = fri; % Store fri for the residual computation
    fri = lagr2forces2( lagr1, c2node1, 3, boundary, nnodes );
    
%     sigma1 = stress(u1,E,nu,nodes,elements,order,1,ntoelem);
    
    error1(iter) = sqrt( myps(u1-uref,u1-uref,Kinter,boundary,M,nodes)/...
        myps(uref,uref,Kinter,boundary,M,nodes) );
%     ferror1(iter) = sqrt( myps(fri-fref,fri-fref,Kinter,boundary,M,nodes)/...
%         myps(fref,fref,Kinter,boundary,M,nodes) );
%     serror1(iter) = transpose(sigma-sigma1)*(sigma-sigma1) / (sigma'*sigma);
%     
%     fresidual(iter) = sqrt( myps(fri-frio,fri-frio,Kinter,boundary,M,nodes)/...
%                       sqrt(myps(fri,fri,Kinter,boundary,M,nodes)*...
%                       myps(frio,frio,Kinter,boundary,M,nodes)) );
    
    % Solve ND
    %fri = lagr2forces( lagr1, C1, 3, boundary );
    u2o = u2;
    fr = [ fri; zeros(size(C2,2),1) ]; % Give to fr the right size
    %fdir1 = dirichletRhs(uref, 1, C2, boundary);
    %fdir2 = dirichletRhs(uref, 2, C2, boundary);
    fdir1 = dirichletRhs2(urefb, 1, c2node2, boundary, nnodes);
    fdir2 = dirichletRhs2(urefb, 2, c2node2, boundary, nnodes);
    % Assembly the dirichlet RHS (ensure a redundant CL isn't imposed 2
    % times
    f2 = fr + assembleDirichlet( [fdir1,fdir2] );

    uin2 = K2\f2;
    u2 = uin2(1:2*nnodes,1);
    
    vo = v;
    v = theta(iter)*u2 + (1-theta(iter))*vo;
    
    if relax == 1 && iter > 1
        e1 = u1-u1o;
        e2 = u2-u2o;
        theta(iter+1) = myps(e1,e1-e2,Kinter,boundary,M,nodes) /...
            myps(e1-e2,e1-e2,Kinter,boundary,M,nodes);
    end

%     sigma2 = stress(u2,E,nu,nodes,elements,order,1,ntoelem);
    
    error2(iter) = sqrt( myps(u2-uref,u2-uref,Kinter,boundary,M,nodes)/...
        myps(uref,uref,Kinter,boundary,M,nodes) );
%     ferror2(iter) = sqrt( myps(frio-fref,frio-fref,Kinter,boundary,M,nodes)/...
%         myps(fref,fref,Kinter,boundary,M,nodes) );
%     serror2(iter) = transpose(sigma-sigma2)*(sigma-sigma2) / (sigma'*sigma);
    if iter > 1
       residual(iter) = sqrt( myps(u1-u2,u1-u2,Kinter,boundary,M,nodes)/...
                        sqrt(myps(u1,u1,Kinter,boundary,M,nodes)*...
                        myps(u2,u2,Kinter,boundary,M,nodes)) );
    end
%     sresidual(iter) = (sigma1-sigma2)'*(sigma1-sigma2) /...
%                        sqrt( (sigma1'*sigma1)*(sigma2'*sigma2) );
                 
    regulari(iter) = sqrt(u2'*regul(u2, nodes, boundary, 3));
%     Nu = regul(u2, nodes, boundary, 3);
end

residual(1) = 1; % tiny hack

hold on;
set(gca, 'fontsize', 15);
set(gca,'ylim',[-3 0])
plot(log10(error1),'Color','black')
plot(log10(error2),'Color','blue')
plot(log10(residual),'Color','red')
legend('error1 (log)','error2 (log)','residual (log)')
figure
%hold on;
% plot(log10(ferror1),'Color','black')
% plot(log10(ferror2),'Color','blue')
% plot(log10(fresidual),'Color','red')
% figure
% L-curve
loglog(residual,regulari);
% plot(regulari.*residual)

% Plot solution
% hold on;
% set(gca, 'fontsize', 15);
% set(gca,'ylim',[-3e-5 3e-5])
% plot(uref(index));
% plot(u2(index),'Color','red');

% Compute stress :
sigma = stress(u2,E,nu,nodes,elements,order,1,ntoelem);
% Output :
plotGMSH({u2,'U_vect';sigma,'stress'}, elements, nodes(:,[1,2]), 'field2');