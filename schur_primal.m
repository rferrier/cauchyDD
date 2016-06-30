% 22/03/2016
% Méthode de Schur Primale

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 70000;  % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 1;      % N.mm-1 : Loading on the plate
niter   = 10;
precond = 1;      % Use of a dual precond

dirichlea = [4,1,0;4,2,0;
             1,1,0;1,2,0;
             2,1,0;2,2,0];
dirichleb = [4,1,0;4,2,0;
             3,1,0;3,2,0;
             2,1,0;2,2,0];
neumann   = [];

% First, import the mesh
[ nodea,elementa,ntoelea,boundara,ordea ] = readmesh( 'meshes/plate.msh' );
[ no,el,nt,bo,or ] = createmesh( 5, 10, 0.5, 0.5, 'meshes/plate_reg.msh' );
nodeb = translateMesh( nodea, 0, 10 );
nnodea = size(nodea,1);

% Then, build the stiffness matrix :
[Ka,Ca,nbloqa,node2ca,c2nodea] =...
    Krig (nodea,elementa,E,nu,ordea,boundara,dirichlea);
[Kb,Cb,nbloqb,node2cb,c2nodeb] =...
    Krig (nodeb,elementa,E,nu,ordea,boundara,dirichleb);
%Kab = penalEq( Ka, Kb, E*1e9, nodea, nodeb, 1e-5);
[ Kab, naj] = lagrEq( Ka, Kb, nodea, nodeb, 1e-5,c2nodea,c2nodeb);
%Kab = [Ka,zeros(size(Ka,1), size(Kb,1)); zeros(size(Kb,1), size(Ka,1)),Kb];

% The right hand side :
fa = volumicLoad( nbloqa, nodea, elementa, 1, fscalar );
fb = volumicLoad( nbloqb, nodeb, elementa, 1, fscalar );

% Solve the problem :
uinab = Kab\[fa;fb;zeros(naj,1)];

% Extract displacement :
u = uinab(1:2*nnodea,1);

% Compute stress :
sigma = stress(u,E,nu,nodea,elementa,ordea,1,ntoelea);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% primal Schur algorithm
% Solve (S10+S20)(W) = -S1-S2 with a CG algo

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
%boundary = suppressBound( boundary, [no3;no4], 3 );

% Some useful tables
[node2b1, b2node1] = mapBound(1, boundaryp2, nnodea);
[node2b3, b2node3] = mapBound(3, boundaryp1, nnodea);

% lower problem
dirichlet1 = [4,1,0;4,2,0;
              1,1,0;1,2,0;
              2,1,0;2,2,0;
              3,1,0;3,2,0];

[K1,C1,nbloq1] = Krig (nodes1,elements,E,nu,order,boundaryp1,dirichlet1);
f1in = volumicLoad( nbloq1, nodes1, elements, 1, fscalar );
Kinter = K1(1:2*nnodes, 1:2*nnodes); % Extract the inner matrix
% upper problem
nodes2 = translateMesh(nodes1, 0, 10); % Move the mesh2
[ b1to2, b2to1 ] = superNodes( nodes1, nodes2, 1e-5 ); % build the connec table

dirichlet2 = [4,1,0;4,2,0;
              1,1,0;1,2,0;
              2,1,0;2,2,0;
              3,1,0;3,2,0];

[K2,C2,nbloq2] = Krig (nodes2,elements,E,nu,order,boundaryp2,dirichlet2);
f2in = volumicLoad( nbloq2, nodes2, elements, 1, fscalar );

%% Dual problem
% lower problem
dirichlet1d = [4,1,0;4,2,0;
               1,1,0;1,2,0;
               2,1,0;2,2,0];

[K1d,C1d,nbloq1d] = Krig (nodes1,elements,E,nu,order,boundary,dirichlet1d);

% upper problem
dirichlet2d = [4,1,0;4,2,0;
               2,1,0;2,2,0;
               3,1,0;3,2,0];

[K2d,C2d,nbloq2d] = Krig (nodes2,elements,E,nu,order,boundary,dirichlet2d);
%%%%
%% Assembly of the Schur operators
K1s = penalise(  Kinter, [1,2,4], boundaryp1, nnodes, E*1e9 );
K2s = penalise(  Kinter, [3,2,4], boundaryp2, nnodes, E*1e9 );

[ S1i, b1i, map1 ] = schurComp( K1s, f1in(1:2*nnodes,1), 3, boundaryp1, nnodes );
[ S2, b2, map2 ] = schurComp( K2s, f2in(1:2*nnodes,1), 1, boundaryp2, nnodes );
ndof = size(S1i,1); S1 = zeros(ndof); b1 = zeros(ndof,1);
% Re-index S2 in order to match S1
% S1 on b1 -> b2node1 -> S1 on n2 -> b2ot1 -> S1 on n1 -> node2b3 -> S1 on 3
i0 = 1:1:ndof/2;         % Indices on b1
i1 = b2node1(i0,1)';     % indices on n2
i2 = b2to1(i1,1)';       % indices on n1
i3 = node2b3(i2,1)';     % indices on b2
m0 = 1:1:ndof;           % The same with dofs (1node = 2dofs)
m3 = zeros(1,ndof);
m3(1, 2.*i0) = 2.*i3;
m3(1, 2.*i0-1) = 2.*i3-1;
% Actual re-indexation
S1(m3, m3) = S1i(m0, m0);
b1(m3, 1) = b1i(m0, 1);
%
Stot = S1+S2;
btot = b1+b2;
% u01 = S1\b1;
% u02 = S2\b2;
% u0 = Stot\btot;
%
% hold on;
% plot(u01,'Color','black');
% plot(u02,'Color','blue');
% plot(u0,'Color','yellow');
% plot(u(map1,1),'Color','red');
% figure;
% 
% % Compute the error with this method :
% error0 = norm(u(map1,1)-u0) / norm(u(map1,1));
% Dual Schur Complements
D1 = inv(S1);
D2 = inv(S2);
Dtot = D1+D2;
%%%%
%%
error   = zeros(niter+1,1);
residual = zeros(niter+1,1);

Iter = zeros(2*nnodes,1); % Iter is on the mesh 2
Res  = zeros(2*nnodes,niter+1);
Zed  = zeros(2*nnodes,niter+1);
d    = zeros(2*nnodes,niter+1);

% Computation of the error
u2n = dirichletRhs(Iter, 1, C2, boundary);
urni = projectField(u, nodea, nodes2, 1e-5);
urn = dirichletRhs(urni, 1, C2, boundary);
error(1) = sqrt( (u2n-urn)'*(u2n-urn)/(urn'*urn) );

%% Compute RHS
% Solve lower
f1 = f1in;
uin1 = K1\f1;
lagr1 = uin1(2*nnodes+1:end,1);
lambda1i = lagr2forces( lagr1, C1, 3, boundaryp1 );
lambda1 = projectField( lambda1i, nodes1, nodes2, 1e-5 );
% Solve upper
f2 = f2in;
uin2 = K2\f2;
lagr2 = uin2(2*nnodes+1:end,1);
lambda2 = lagr2forces( lagr2, C2, 1, boundaryp2 );

b = - lambda1 - lambda2;
bs = btot;   % with Schur Complement
%%%%
%% Compute Ax0
% Solve lower
ui2 = projectField( Iter, nodes2, nodes1, 1e-5 );
f1 = dirichletRhs(ui2, 3, C1, boundaryp1);
uin1 = K1\f1;
lagr1 = uin1(2*nnodes+1:end,1);
lambda1i = lagr2forces( lagr1, C1, 3, boundaryp1 );
lambda1 = projectField( lambda1i, nodes1, nodes2, 1e-5 );

% Solve upper
f2 = dirichletRhs(Iter, 1, C2, boundaryp2);
uin2 = K2\f2;
lagr2 = uin2(2*nnodes+1:end,1);
lambda2 = lagr2forces( lagr2, C2, 1, boundaryp2 );

Axz = lambda1+lambda2;
Axzs = Stot*Iter(map2,1);  % With Schur complement
%%%%
%% Init
Res(:,1) = b - Axz;

if precond == 1
    % Solve lower dual
    loadi = projectField( Res(:,1)/2, nodes2, nodes1, 1e-5 );
    load1 = [loadi ; zeros(nbloq1d,1)];
    uin1 = K1d\load1; % Res lives on the domain 2
    u1i = uin1(1:2*nnodes,1);
    u1 = projectField( u1i, nodes1, nodes2, 1e-5 );

    % Solve upper dual
    load2 = [Res(:,1)/2 ; zeros(nbloq2d,1)];
    uin2 = K2d\load2;
    u2 = uin2(1:2*nnodes,1);

    Zed(:,1) = u1/2 + u2/2;
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
    ui2 = projectField( d(:,iter), nodes2, nodes1, 1e-5 );
    f1 = dirichletRhs(ui2, 3, C1, boundaryp1);
    uin1 = K1\f1;
    lagr1 = uin1(2*nnodes+1:end,1);
    lambda1i = lagr2forces( lagr1, C1, 3, boundaryp1 );
    lambda1 = projectField( lambda1i, nodes1, nodes2, 1e-5 );
    % Solve upper
    f2 = dirichletRhs(d(:,iter), 1, C2, boundaryp2);
    uin2 = K2\f2;
    lagr2 = uin2(2*nnodes+1:end,1);
    lambda2 = lagr2forces( lagr2, C2, 1, boundaryp2 );
    %
    Ad = lambda1+lambda2;
%     Ad = lambda1-lambda1;
%     Ad(map2, 1) = Stot*d(map2,iter);  % With Schur complement
    Ads = Stot*d(map2,iter);  % With Schur complement
%     hold on;
%     plot(Ads,'Color','blue');
%     plot(Ad(map2,1),'Color','red');
%     figure;
    
    alpha = num/(d(:,iter)'*Ad);
    %%
    Iter          = Iter + d(:,iter)*alpha;
    Res(:,iter+1) = Res(:,iter) - Ad*alpha;
    
    if precond == 1
        % Solve lower dual
        loadi = projectField( Res(:,iter+1)/2, nodes2, nodes1, 1e-5 );
        load = [loadi ; zeros(nbloq1d,1)];
        uin1 = K1d\load; % Res lives on the domain 2
        u1i = uin1(1:2*nnodes,1);
        u1 = projectField( u1i, nodes1, nodes2, 1e-5 );

        % Solve upper dual
        load = [Res(:,iter+1)/2 ; zeros(nbloq2d,1)];
        uin2 = K2d\load;
        u2 = uin2(1:2*nnodes,1);

        Zed(:,iter+1) = u1/2 + u2/2;
    else
        Zed(:,iter+1) = Res(:,iter+1);
    end
    
    %residual(iter+1) = sqrt( Zed(:,iter+1)'*Zed(:,iter+1) );
    residual(iter+1) = sqrt( Res(:,iter+1)'*Res(:,iter+1) );
    
    beta = - ( Zed(:,iter+1)'*Ad )/( d(:,iter)'*Ad );
    d(:,iter+1) = Zed(:,iter+1) + d(:,iter)*beta;
    
    % Computation of the error
    u2n = dirichletRhs(Iter, 1, C2, boundary);
    error(iter+1) = sqrt( (u2n-urn)'*(u2n-urn)/(urn'*urn) );
end
%% Post-processing

% Computation of the error :
% Solve lower
ui2 = projectField( Iter, nodes2, nodes1, 1e-5 );
f1 = dirichletRhs(ui2, 3, C1, boundaryp1) + f1in;
uin1 = K1\f1;
u1 = uin1(1:2*nnodes,1);

% Solve upper
f2 = dirichletRhs(Iter, 1, C2, boundaryp2) + f2in;
uin2 = K2\f2;
u2 = uin2(1:2*nnodes,1);

%     u1n = dirichletRhs(u1, 3, C1, boundary);
%     urni = projectField(u, nodes, nodes1, 1e-5);
%     urn = dirichletRhs(urni, 3, C1, boundary);
%     error1(iter) = (u1n-urn)'*(u1n-urn)/(urn'*urn);
    % Computation of the error :
%     u2ni = projectField(u2, nodes2, nodes1, 1e-5);
%     u2n = dirichletRhs(u2ni, 3, C1, boundary);
%     %urni = projectField(u, nodes, nodes1, .5);
%     %urn = dirichletRhs(urni, 3, C1, boundary);
%     error2(iter) = (u2n-urn)'*(u2n-urn)/(urn'*urn);
%     residual(iter) = (u2n-u1n)'*(u2n-u1n)/sqrt((u2n'*u2n)*(u1n'*u1n));

sigma1 = stress(u1,E,nu,nodes1,elements,order,1,ntoelem);
sigma2 = stress(u2,E,nu,nodes2,elements,order,1,ntoelem);

plotGMSH({u1,'U_vect';sigma1,'stress'}, elements, nodes1(:,[1,2]), 'field1');
plotGMSH({u2,'U_vect';sigma2,'stress'}, elements, nodes1(:,[1,2]), 'field2');
% 
hold on;
set(gca, 'fontsize', 15);
set(gca,'ylim',[-18 0]);
plot(log10(error(2:end)),'Color','blue')
plot(log10(residual(2:end)/residual(1)),'Color','red')
legend('error (log)','residual (log)')