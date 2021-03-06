% 01/02/2017
% Algo Steklov-Poincaré primal bloc avec Gradient Conjugué

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 70000;  % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 1;      % N.mm-1 : Loading on the plate
niter   = 5;
precond = 1;      % 1 : Use a dual precond, 2 : use regul precond
mu      = 0.;    % Regularization parameter
br      = 0.;     % noise

mat = [0, E, nu];

% Boundary conditions
% first index  : index of the boundary
% second index : 1=x, 2=y
% third        : value
% [0,1,value] marks a dirichlet regularization therm on x
dirichlet = [1,1,0; 1,2,0 ;
             3,1,0; 3,2,0 ];
%neumann   = [2,1,fscalar,0,fscalar;
%             4,1,fscalar,0,-fscalar];
neumann   = [2,1,fscalar;
             4,1,fscalar];

% Import the mesh
[ nodes,elements,ntoelem,boundary,order ] = readmesh( 'meshes/platee.msh' );
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

boundaryp1 = suppressBound( boundary, no2, 2 );
boundaryp1 = suppressBound( boundaryp1, no3, 2 );
boundaryp1 = suppressBound( boundaryp1, no1, 4 );
boundaryp1 = suppressBound( boundaryp1, no4, 4 );
boundaryp2 = boundaryp1;

[~, b2node2] = mapBound(2, boundaryp1, nnodes);
indexa = [2*b2node2-1; 2*b2node2];

% Then, build the stiffness matrix :
[K,C,nbloq] = Krig2 (nodes,elements,mat,order,boundary,dirichlet);
Kinter = K(1:2*nnodes, 1:2*nnodes);
M      = mass_mat(nodes, elements);
[node2b4, b2node4] = mapBound(4, boundaryp1, nnodes);
[node2b3, b2node3] = mapBound(3, boundaryp1, nnodes);
[node2b1, b2node1] = mapBound(1, boundaryp1, nnodes);
[node2b2, b2node2] = mapBound(2, boundaryp1, nnodes);
index = [2*b2node2-1;2*b2node2];
% The right hand side :
f = loading(nbloq,nodes,boundary,neumann);

% Solve the problem :
uin = K\f;

% Extract displacement and Lagrange multiplicators :
uref = uin(1:2*nnodes,1);
lagr = uin(2*nnodes+1:end,1);

urefb = ( 1 + br*randn(2*nnodes,1) ) .* uref;

sigma = stress(uref,E,nu,nodes,elements,order,1,ntoelem);
plotGMSH({uref,'Vect_U';sigma,'stress'}, elements, nodes, 'reference');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition of the operators

% First problem
dirichlet1 = [4,1,0;4,2,0;
              3,1,0;3,2,0;
              2,1,0;2,2,0;
              1,1,0;1,2,0];

[K1,C1,nbloq1,node2c1,c2node1] = Krig (nodes,elements,E,nu,order,boundaryp1,dirichlet1);

% Second problem
dirichlet2 = [1,1,0;1,2,0;
              2,1,0;2,2,0;
              3,1,0;3,2,0];
%neumann2   = [4,1,fscalar,0,-fscalar];
neumann2   = [4,1,fscalar];
neumann0   = [];
[K2,C2,nbloq2,node2c2,c2node2] = Krig (nodes,elements,E,nu,order,boundaryp2,dirichlet2);

error    = zeros(niter+1,1);
residual = zeros(niter+1,1);
regulari = zeros(niter+1,1);

%% Dual problems
% First problem
dirichlet1d = [4,1,0;4,2,0;
               3,1,0;3,2,0;
               1,1,0;1,2,0];
[K1d,C1d,nbloq1d] = Krig (nodes,elements,E,nu,order,boundary,dirichlet1d);

% Second problem
dirichlet2d = [1,1,0;1,2,0;
               3,1,0;3,2,0];
[K2d,C2d,nbloq2d] = Krig (nodes,elements,E,nu,order,boundary,dirichlet2d);

%% Anti-cancellation trick
% K1(indexa,indexa) = 0;
% K2(indexa,indexa) = 0;
% K1d(indexa,indexa) = 0;
% K2d(indexa,indexa) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Conjugate Gradient for the problem : (S10-S20) x = S2-S1
Itere = zeros( 2*nnodes, 2 );
d     = zeros( 2*nnodes, 2*(niter+1) );  % 2 directions per step
Ad    = zeros( 2*nnodes, 2*(niter+1) );
Res   = zeros( 2*nnodes, 2*(niter+1) );
Zed   = zeros( 2*nnodes, 2*(niter+1) );
H12   = H1demi(size(Res,1), nodes, boundary, 2 );
Mgr   = Mgrad(size(Res,1), nodes, boundary, 2 );

%% Perform A x0 :
% Solve 1
f1 = dirichletRhs2( Itere(:,1), 2, c2node1, boundaryp1, nnodes );
uin1 = K1\f1;
lagr1 = uin1(2*nnodes+1:end,1);
lamb1 = lagr2forces( lagr1, C1, 2, boundaryp1 );
% Solve 2
f2 = dirichletRhs2( Itere(:,2), 2, c2node2, boundaryp2, nnodes );
uin2 = K2\f2;
lagr2 = uin2(2*nnodes+1:end,1);
lamb2 = lagr2forces( lagr2, C2, 2, boundaryp2 );
%
Axz = [lamb1,lamb2];
%%%%
%% Compute Rhs :
% Solve 1
f1 = dirichletRhs2( urefb, 4, c2node1, boundary, nnodes );
uin1 = K1\f1;
lagr1 = uin1(2*nnodes+1:end,1);
lamb1 = lagr2forces2( lagr1, c2node1, 2, boundaryp1, nnodes );
% Solve 2
f2 = loading(nbloq2,nodes,boundary,neumann2);
uin2 = K2\f2;
lagr2 = uin2(2*nnodes+1:end,1);
lamb2 = lagr2forces2( lagr2, c2node2, 2, boundaryp2, nnodes );
%
b = -[lamb1,lamb2];

%%
Res(:,[1,2]) = b - Axz;

if precond == 1
    % Solve 1 (Sd)
    f1 = [Res(:,[1,2]); zeros(nbloq1d,2)];
    uin1 = K1d\f1;
    u1i = uin1(1:2*nnodes,[1,2]);
    u1 = keepField( u1i, 2, boundaryp1 );
    %
    Zed(:,[1,2]) = u1;
elseif precond == 2
    Zed(index,1) = 1/E*H12(index,index)\Res(index,1);
elseif precond == 3
    Zed(index,1) = 1/E*Mgr(index,index)\Res(index,1);
else
    Zed(:,[1,2]) = Res(:,[1,2]);
end

d(:,[1,2]) = Zed(:,[1,2]);

residual(1) = norm( Res(indexa,1)-Res(indexa,2) );
error(1)    = norm(Itere(indexa,1) - Itere(indexa,2) - uref(indexa)) / ...
                                    norm(uref(indexa));
regulari(1) = sqrt( (Itere(:,1)'-Itere(:,2)')* ... 
                     regul( Itere(:,1)-Itere(:,2) , nodes, boundary, 2) );

%%
for iter = 1:niter
    %% Optimal step
    
    % Solve 1
    f11 = dirichletRhs2( d(:,2*iter-1), 2, c2node1, boundaryp1, nnodes );
    f12 = dirichletRhs2( d(:,2*iter), 2, c2node1, boundaryp1, nnodes );
    uin1 = K1\[f11,f12];
    lagr1 = uin1(2*nnodes+1:end,[1,2]);
    lamb11 = lagr2forces( lagr1(:,1), C1, 2, boundaryp1 );
    lamb12 = lagr2forces( lagr1(:,2), C1, 2, boundaryp1 );
    % Solve 2
    f21 = dirichletRhs2( d(:,2*iter-1), 2, c2node2, boundaryp2, nnodes );
    f22 = dirichletRhs2( d(:,2*iter), 2, c2node2, boundaryp2, nnodes );
    uin2 = K2\[f21,f22];
    lagr2 = uin2(2*nnodes+1:end,[1,2]);
    lamb21 = lagr2forces( lagr2(:,1), C2, 2, boundaryp2 );
    lamb22 = lagr2forces( lagr2(:,2), C2, 2, boundaryp2 );
    %
    Ad(:,[2*iter-1,2*iter]) = [lamb11-lamb21,lamb12-lamb22];
    
    den = d(indexa,[2*iter-1,2*iter])'*Ad(indexa,[2*iter-1,2*iter]);
    sqD = den^(1/2);
    d(:,[2*iter-1,2*iter]) = d(:,[2*iter-1,2*iter])*inv(sqD); %transpose( sqD\d(:,[2*iter-1,2*iter])' );
    Ad(:,[2*iter-1,2*iter]) = Ad(:,[2*iter-1,2*iter])*inv(sqD); %transpose( sqD\Ad(:,[2*iter-1,2*iter])' );
    num = Res(indexa,[2*iter-1,2*iter])'*Zed(indexa,[2*iter-1,2*iter]);
    num = inv(sqD)*num; % because of non-commutiativity
    alpha = den\num;
%    Itere = Itere + d(:,[2*iter-1,2*iter])*alpha;
%    Res(:,[2*iter+1,2*iter+2]) = Res(:,[2*iter-1,2*iter]) - ...
%                                    Ad(:,[2*iter-1,2*iter])*alpha;
    Itere = Itere + d(:,[2*iter-1,2*iter])*num;
    Res(:,[2*iter+1,2*iter+2]) = Res(:,[2*iter-1,2*iter]) - ...
                                    Ad(:,[2*iter-1,2*iter])*num;
    
    residual(iter+1) = norm( Res(indexa,2*iter+1)-Res(indexa,2*iter+2) );
    error(iter+1)    = norm(Itere(indexa,1) - Itere(indexa,2) - uref(indexa)) / ...
                                    norm(uref(indexa));
    regulari(iter+1) = sqrt( (Itere(:,1)'-Itere(:,2)')* ... 
                          regul( Itere(:,1)-Itere(:,2) , nodes, boundary, 2) );
    
    if precond == 1
       % Solve 1 (Sd)
       f1 = [Res(:,[2*iter+1,2*iter+2]); zeros(nbloq1d,2)];
       uin1 = K1d\f1;
       u1i = uin1(1:2*nnodes,[1,2]);
       u1 = keepField( u1i, 2, boundaryp1 );
       %
       Zed(:,[2*iter+1,2*iter+2]) = u1;
    elseif precond == 2
        Zed(index,iter+1) = 1/E*H12(index,index)\Res(index,iter+1);
    elseif precond == 3
        Zed(index,iter+1) = 1/E*Mgr(index,index)\Res(index,iter+1);
    else
        Zed(:,[2*iter+1,2*iter+2]) = Res(:,[2*iter+1,2*iter+2]);
    end
    
    %% Orthogonalization
    d(:,[2*iter+1,2*iter+2]) = Zed(:,[2*iter+1,2*iter+2]);
    
    for jter=1:iter
        betaij = ( d(indexa,[2*jter-1,2*jter])'*Ad(indexa,[2*jter-1,2*jter]) ) \ ...
            ( d(indexa,[2*iter+1,2*iter+2])'*Ad(indexa,[2*jter-1,2*jter]) );
        betaij = ( Ad(indexa,[2*jter-1,2*jter])'*d(indexa,[2*jter-1,2*jter]) ) \ ...
            ( Ad(indexa,[2*jter-1,2*jter])'*d(indexa,[2*iter+1,2*iter+2]) );
%        betaij = ( d(indexa,[2*iter+1,2*iter+2])'*Ad(indexa,[2*jter-1,2*jter]) );
        d(:,[2*iter+1,2*iter+2]) = d(:,[2*iter+1,2*iter+2]) - ...
                                   d(:,[2*jter-1,2*jter]) * betaij;
    end

%     betaij = ( Zed(indexa,iter+1)'*Ad(indexa,iter) )/...
%         ( d(indexa,iter)'*Ad(indexa,iter) );
%     d(:,iter+1) = d(:,iter+1) - d(:,iter) * betaij;
    
end

hold on
plot(log10(error),'Color','blue')
plot(log10(residual),'Color','red')
legend('error','residual')
figure;
%% L-curve :
%loglog(residual(2:end),regulari(2:end));
%figure
%%%%
%% Final problem : compute u
% DN problem
dirichlet = [4,1,0;4,2,0;
             3,1,0;3,2,0;
             1,1,0;1,2,0;
             2,1,0;2,2,0];
neumann   = [];
[K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
fdir2 = dirichletRhs(Itere(:,1)-Itere(:,2), 2, C, boundary);
fdir4 = dirichletRhs(urefb, 4, C, boundary);
usoli = K \ (fdir4 + fdir2);

usol = usoli(1:2*nnodes,1);
fsol = Kinter*usol;

hold on;
plot(usol(2*b2node2-1),'Color','red');
plot(uref(2*b2node2-1),'Color','green');
legend('solution','reference')

total_error = norm(usol-uref)/norm(uref);
% Compute stress :
sigma = stress(usol,E,nu,nodes,elements,order,1,ntoelem);
plotGMSH({usol,'U_vect';sigma,'stress'}, elements, nodes, 'solution');