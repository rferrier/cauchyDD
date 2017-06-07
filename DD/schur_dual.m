% 06/04/2016
% M�thode de Schur Duale

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 70000;  % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 1;      % N.mm-1 : Loading on the plate
niter   = 10;
precond = 1;      % Use of a primal precond

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
[Ka,Ca,nbloqa,node2ca,c2nodea] = ...
    Krig (nodea,elementa,E,nu,ordea,boundara,dirichlea);
[Kb,Cb,nbloqb,node2cb,c2nodeb] = ...
    Krig (nodeb,elementa,E,nu,ordea,boundara,dirichleb);
%Kab = penalEq( Ka, Kb, E*1e9, nodea, nodeb, 1e-5 );
[ Kab, naj] = lagrEq( Ka, Kb, nodea, nodeb, 1e-5,c2nodea,c2nodeb);
%Kab = [Ka,zeros(size(Ka,1), size(Kb,1)); zeros(size(Kb,1), size(Ka,1)),Kb];

% Some useful tables
[node2b1, b2node1] = mapBound(1, boundara, nnodea);
[node2b3, b2node3] = mapBound(3, boundara, nnodea);

% The right hand side :
fa = volumicLoad( nbloqa, nodea, elementa, 1, fscalar );
fb = volumicLoad( nbloqb, nodeb, elementa, 1, fscalar );

% Solve the problem :
uinab = Kab\[fa;fb;zeros(naj,1)];

% Extract displacement and reaction :
u    = uinab(1:2*nnodea,1);
ua   = uinab(1:2*nnodea+nbloqa,1);
fref = fa - Ka*ua;

% Compute stress :
sigma = stress(u,E,nu,nodea,elementa,ordea,1,ntoelea);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dual Schur algorithm
% Solve (D10+D20)(W) = -D1-D2 with a CG algo

% We work on a half-mesh
[ nodes1,elements,ntoelem,boundary,order ] = readmesh( 'meshes/plate.msh' );
nnodes = size(nodes1,1);

% find the nodes in the corners and suppress the elements :
xmax = max(nodes1(:,1));
xmin = min(nodes1(:,1));
ymax = max(nodes1(:,2));
ymin = min(nodes1(:,2));
no1  = findNode(xmin, ymin, nodes1, 1e-5);
no2  = findNode(xmax, ymin, nodes1, 1e-5);
no3  = findNode(xmax, ymax, nodes1, 1e-5);
no4  = findNode(xmin, ymax, nodes1, 1e-5);
boundaryp2 = suppressBound( boundary, [no1;no2], 1 );
boundaryp1 = suppressBound( boundary, [no3;no4], 3 );

%% Primal operators
% lower problem
dirichlet1 = [4,1,0;4,2,0;
              1,1,0;1,2,0;
              2,1,0;2,2,0;
              3,1,0;3,2,0];

[K1,C1,nbloq1,node2c1,c2node1] = Krig (nodes1,elements,E,nu,order,boundary,dirichlet1);
Kinter = K1(1:2*nnodes, 1:2*nnodes); % Extract the inner matrix

% upper problem
nodes2 = translateMesh(nodes1, 0, 10); % Move the mesh2
[ b1to2, b2to1 ] = superNodes( nodes1, nodes2, 1e-5 ); % build the connec table
dirichlet2 = [4,1,0;4,2,0;
              1,1,0;1,2,0;
              2,1,0;2,2,0;
              3,1,0;3,2,0];
[K2,C2,nbloq2,node2c2,c2node2] = Krig (nodes2,elements,E,nu,order,boundary,dirichlet2);

%% Dual problem
% lower problem
dirichlet1d = [4,1,0;4,2,0;
               1,1,0;1,2,0;
               2,1,0;2,2,0];

[K1d,C1d,nbloq1d,node2c1d,c2node1d] = Krig (nodes1,elements,E,nu,order,boundary,dirichlet1d);
f1in = volumicLoad( nbloq1d, nodes1, elements, 1, fscalar );

% upper problem
dirichlet2d = [4,1,0;4,2,0;
               2,1,0;2,2,0;
               3,1,0;3,2,0];

[K2d,C2d,nbloq2d,node2c2d,c2node2d] = Krig (nodes2,elements,E,nu,order,boundary,dirichlet2d);
f2in = volumicLoad( nbloq2d, nodes2, elements, 1, fscalar );

% Move reference solution on the second mesh
frni = projectField(fref, nodea, nodes2, 1e-5);
frn = dirichletRhs(frni, 1, C2, boundary);
%%%%
%%
error    = zeros(niter+1,1);
residual = zeros(niter+1,1);

Iter = zeros(2*nnodes,1);       % Iter is on the mesh 2
Res  = zeros(2*nnodes,niter+1);
Zed  = zeros(2*nnodes,niter+1);
d    = zeros(2*nnodes,niter+1);

% Computation of the error
f2n = dirichletRhs2(Iter, 1, c2node2, boundary, nnodes);
% frni = projectField(fref, nodea, nodes2, 1e-5);
% frn = dirichletRhs(frni, 1, C2, boundary);
error(1) = sqrt( (f2n-frn)'*(f2n-frn)/(frn'*frn) );

%% Compute RHS
% Solve lower
f1 = f1in;
uin1 = K1d\f1;
u1i = uin1(1:2*nnodes,1);
u1 = projectField( u1i, nodes1, nodes2, 1e-5 );
% Solve upper
f2 = f2in;
uin2 = K2d\f2;
u2 = keepField( uin2(1:2*nnodes,1), 1, boundary );
%
b = -u1+u2;
%%%%
%% Compute Ax0
% Solve lower
f1i = projectField( -Iter, nodes2, nodes1, 1e-5 );
f1 = [ f1i ; zeros(nbloq1d,1) ];
uin1 = K1d\f1;
u1i = uin1(1:2*nnodes,1);
u1 = projectField( u1i, nodes1, nodes2, 1e-5 );
% Solve upper
f2 = [ Iter ; zeros(nbloq2d,1) ];
uin2 = K2d\f2;
u2 = keepField( uin2(1:2*nnodes,1), 1, boundary );
%
Axz = u1-u2;
%%%%
%% Init
Res(:,1) = b - Axz;

if precond == 1
    % Solve lower primal
    ui2 = projectField( Res(:,1)/2, nodes2, nodes1, 1e-5 );
    f1 = dirichletRhs2(ui2, 3, c2node1, boundaryp1, nnodes);
    uin1 = K1\f1;
    lagr1 = uin1(2*nnodes+1:end,1);
    lambda1i = lagr2forces2( lagr1, c2node1, 3, boundaryp1, nnodes );
    lambda1 = projectField( lambda1i, nodes1, nodes2, 1e-5 );
    % Solve upper primal
    f2 = dirichletRhs2(Res(:,1)/2, 1, c2node2, boundaryp2, nnodes);
    uin2 = K2\f2;
    lagr2 = uin2(2*nnodes+1:end,1);
    lambda2 = lagr2forces2( lagr2, c2node2, 1, boundaryp2, nnodes );
    %
    Zed(:,1) = lambda1/2+lambda2/2;
else
    Zed(:,1) = Res(:,1);
end

%residual(1) = sqrt( Zed(:,1)'*Zed(:,1) );
residual(1) = sqrt( Res(:,1)'*Res(:,1) );

d(:,1)   = Zed(:,1);

for iter = 1:niter
    %% Optimal step
    num = Res(:,iter)'*d(:,iter);
    
    % Solve lower
    f1i = projectField( -d(:,iter), nodes2, nodes1, 1e-5 );
    f1 = [ f1i ; zeros(nbloq1d,1) ];
    uin1 = K1d\f1;
    u1i = uin1(1:2*nnodes,1);
    u1 = projectField( u1i, nodes1, nodes2, 1e-5 );
    % Solve upper
    f2 = [ d(:,iter) ; zeros(nbloq2d,1) ];
    uin2 = K2d\f2;
    u2 = keepField( uin2(1:2*nnodes,1), 1, boundary  );
    %
    Ad = u1-u2;
    %
    alpha = num/(d(:,iter)'*Ad);
    %%
    Iter          = Iter + d(:,iter)*alpha;
    Res(:,iter+1) = Res(:,iter) - Ad*alpha;
    
    if precond == 1
        % Solve lower primal
        ui2 = projectField( Res(:,iter+1)/2, nodes2, nodes1, 1e-5 );
        f1 = dirichletRhs2(ui2, 3, c2node1, boundaryp1, nnodes);
        uin1 = K1\f1;
        lagr1 = uin1(2*nnodes+1:end,1);
        lambda1i = lagr2forces2( lagr1, c2node1, 3, boundaryp1, nnodes );
        lambda1 = projectField( lambda1i, nodes1, nodes2, 1e-5 );
        % Solve upper primal
        f2 = dirichletRhs2(Res(:,iter+1)/2, 1, c2node2, boundaryp2, nnodes);
        uin2 = K2\f2;
        lagr2 = uin2(2*nnodes+1:end,1);
        lambda2 = lagr2forces2( lagr2, c2node2, 1, boundaryp2, nnodes );
        %
        Zed(:,iter+1) = lambda1/2 + lambda2/2;
    else
        Zed(:,iter+1) = Res(:,iter+1);
    end
    
    %residual(iter+1) = sqrt( Zed(:,iter+1)'*Zed(:,iter+1) );
    residual(iter+1) = sqrt( Res(:,iter+1)'*Res(:,iter+1) );
    
    beta = - ( Zed(:,iter+1)'*Ad )/( d(:,iter)'*Ad );
    d(:,iter+1) = Zed(:,iter+1) + d(:,iter)*beta;
    
    % Computation of the error
    f2n = dirichletRhs2(Iter, 1, c2node2, boundary, nnodes);
    error(iter+1) = sqrt( (f2n-frn)'*(f2n-frn)/(frn'*frn) );
end
%% Post-processing

% Computation of the error :
% Solve lower
f1i = projectField( -Iter, nodes2, nodes1, 1e-5 );
f1 = [ f1i ; zeros(nbloq1d,1)] + f1in;
uin1 = K1d\f1;
u1 = uin1(1:2*nnodes,1);

% Solve upper
f2 = [ Iter ; zeros(nbloq1d,1)] + f2in;
uin2 = K2d\f2;
u2 = uin2(1:2*nnodes,1);

sigma1 = stress(u1,E,nu,nodes1,elements,order,1,ntoelem);
sigma2 = stress(u2,E,nu,nodes2,elements,order,1,ntoelem);

plotGMSH({u1,'U_vect';sigma1,'stress'}, elements, nodes1(:,[1,2]), 'field1');
plotGMSH({u2,'U_vect';sigma2,'stress'}, elements, nodes1(:,[1,2]), 'field2');

hold on;
set(gca, 'fontsize', 15);
set(gca,'ylim',[-18 0]);
plot(log10(error(2:end)),'Color','blue')
plot(log10(residual(2:end)/residual(1)),'Color','red')
legend('error (log)','residual (log)')