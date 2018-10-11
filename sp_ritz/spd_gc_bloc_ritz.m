% 07/02/2017
% Algo Steklov-Poincaré dual bloc avec Gradient Conjugué

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 70000;  % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 1;      % N.mm-1 : Loading on the plate
niter   = 5;
precond = 1;      % 1/2 : Use a dual precond
mu      = 0.;     % Regularization parameter
br      = 0.;     % noise
ntrunc  = 5;
difmesh = 0;

mat = [0, E, nu];

% Boundary conditions
% first index  : index of the boundary
% second index : 1=x, 2=y
% third        : value
% [0,1,value] marks a dirichlet regularization therm on x
dirichlet = [1,1,0; 1,2,0 ;
             3,1,0; 3,2,0 ];
neumann   = [2,1,fscalar,0,fscalar;
             4,1,fscalar,0,-fscalar];

%% Import the unrefined mesh (that is in fact an overkill mesh)
if difmesh == 1
   [ nodesu,elementsu,ntoelemu,boundaryu,orderu ] = readmesh( 'meshes/plateu.msh' );
   nnodesu = size(nodesu,1);
   
   % Then, build the stiffness matrix :
   [K,C,nbloq] = Krig2 (nodesu,elementsu,mat,orderu,boundaryu,dirichlet);
   Kinter = K(1:2*nnodesu, 1:2*nnodesu);
   
   % Solve the problem
   f = loading(nbloq,nodesu,boundaryu,neumann);
   uin = K\f;
   urefu = uin(1:2*nnodesu,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             
% Import the mesh
[ nodes,elements,ntoelem,boundary,order ] = readmesh( 'meshes/plate.msh' );
nnodes = size(nodes,1);

noises = load('./noises/noise1.mat'); % Particular noise vector
noise  = noises.bruit1;
%noise = randn(2*nnodes,1);

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
if difmesh == 0
   uin = K\f;
   % Extract displacement and Lagrange multiplicators :
   uref = uin(1:2*nnodes,1);
   lagr = uin(2*nnodes+1:end,1);
   
else % Pass solution
   uref = passMesh2D (nodesu, elementsu, nodes, elements, urefu);
end

urefb = ( 1 + br*noise ) .* uref;
fref  = Kinter*uref;

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
neumann2   = [4,1,fscalar,0,-fscalar];
%neumann2   = [4,1,fscalar];
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
[K1d,C1d,nbloq1d,node2c1d,c2node1d] = Krig (nodes,elements,E,nu,order,boundary,dirichlet1d);

% Second problem
dirichlet2d = [1,1,0;1,2,0;
               3,1,0;3,2,0];
[K2d,C2d,nbloq2d,node2c2d,c2node2d] = Krig (nodes,elements,E,nu,order,boundary,dirichlet2d);

%% Anti-cancellation trick
% K1(indexa,indexa) = 0;
% K2(indexa,indexa) = 0;
% K1d(indexa,indexa) = 0;
% K2d(indexa,indexa) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Conjugate Gradient for the problem : (D10-D20) x = D2-D1
Itere = zeros( 2*nnodes, 2 );
d     = zeros( 2*nnodes, 2*(niter+1) );  % 2 directions per step
Ad    = zeros( 2*nnodes, 2*(niter+1) );
Res   = zeros( 2*nnodes, 2*(niter+1) );
Zed   = zeros( 2*nnodes, 2*(niter+1) );
AZed  = zeros( 2*nnodes, 2*(niter+1) );
H12   = H1demi(size(Res,1), nodes, boundary, 2 );
Mgr   = Mgrad(size(Res,1), nodes, boundary, 2 );
alpha = zeros( 2*(niter+1) );
beta  = zeros( 2*(niter+1) );

% Equivalent Saad matrices
unsuralpha   = zeros( 2*(niter+1) );
betasuralpha = zeros( 2*(niter+1) );
etaeta       = zeros( 2*(niter+1) );

%% Perform A x0 :
% Solve 1
f1 = [Itere(:,1); zeros(nbloq1d,1)];
uin1 = K1d\f1;
u1i = uin1(1:2*nnodes,1);
u1 = keepField( u1i, 2, boundaryp1 );
% Solve 2
f2 = [Itere(:,2); zeros(nbloq2d,1)];
uin2 = K2d\f2;
u2i = uin2(1:2*nnodes,1);
u2 = keepField( u2i, 2, boundaryp2 );
%
Axz = -[u1,u2];
%%%%
%% Compute Rhs :
% Solve 1
f1 = dirichletRhs2( urefb, 4, c2node1d, boundary, nnodes );
uin1 = K1d\f1;
u1i = uin1(1:2*nnodes,1);
u1 = keepField( u1i, 2, boundaryp1 );
% Solve 2
f2 = loading(nbloq2d,nodes,boundary,neumann2);
uin2 = K2d\f2;
u2i = uin2(1:2*nnodes,1);
u2 = keepField( u2i, 2, boundaryp2 );
%
b = [u1,u2];
%%
Res(:,[1,2]) = b - Axz;

if precond == 1
    % Solve 1 (Sd)
    f11 = dirichletRhs2( Res(:,1), 2, c2node1, boundaryp1, nnodes );
    f12 = dirichletRhs2( Res(:,2), 2, c2node1, boundaryp1, nnodes );
    uin1 = K1\[f11,f12];
    lagr1 = uin1(2*nnodes+1:end,:);
    lamb11 = lagr2forces2( lagr1(:,1), c2node1, 2, boundaryp1, nnodes );
    lamb12 = lagr2forces2( lagr1(:,2), c2node1, 2, boundaryp1, nnodes );
    %
    Zed(:,[1,2]) = [lamb11,lamb12];
elseif precond == 2
    % Solve 2 (Sn)
    f21 = dirichletRhs2( Res(:,1), 2, c2node2, boundaryp2, nnodes );
    f22 = dirichletRhs2( Res(:,2), 2, c2node2, boundaryp2, nnodes );
    uin2 = K2\[f21,f22];
    lagr2 = uin2(2*nnodes+1:end,:);
    lamb21 = lagr2forces2( lagr2(:,1), c2node2, 2, boundaryp2, nnodes );
    lamb22 = lagr2forces2( lagr2(:,2), c2node2, 2, boundaryp2, nnodes );
    %
    Zed(:,[1,2]) = [lamb21,lamb22];
else
    Zed(:,[1,2]) = Res(:,[1,2]);
end

d(:,[1,2]) = Zed(:,[1,2]);

residual(1) = norm( Res(indexa,1)-Res(indexa,2) );
error(1)    = norm(Itere(indexa,1) - Itere(indexa,2) - fref(indexa)) / ...
                                    norm(fref(indexa));
regulari(1) = sqrt( (Itere(:,1)'-Itere(:,2)')* ... 
                     regul( Itere(:,1)-Itere(:,2) , nodes, boundary, 2) );
%%
V  = zeros(2*nnodes, 2*niter);
AV = zeros(2*nnodes, 2*niter);
MV = zeros(2*nnodes, 2*niter);
H  = zeros(2*niter);

num = [0,0;0,0]; % useless, but eta needs initialization #lazy
den = [0,0;0,0];
%%
for iter = 1:niter
    %% Optimal step
    % Solve 1
    f1 = [d(:,[2*iter-1,2*iter]); zeros(nbloq1d,2)];
    uin1 = K1d\f1;
    u1i = uin1(1:2*nnodes,:);
    u1 = keepField( u1i, 2, boundaryp1 );
    % Solve 2
    f2 = [d(:,[2*iter-1,2*iter]); zeros(nbloq2d,2)];
    uin2 = K2d\f2;
    u2i = uin2(1:2*nnodes,:);
    u2 = keepField( u2i, 2, boundaryp2 );
    %
    Ad(:,[2*iter-1,2*iter]) = u2-u1;

    denprec = den; numprec = num; % Store those ones
    den = d(indexa,[2*iter-1,2*iter])'*Ad(indexa,[2*iter-1,2*iter]);
    sqD = den^(1/2);
    d(:,[2*iter-1,2*iter]) = d(:,[2*iter-1,2*iter]) * inv(sqD);
    Ad(:,[2*iter-1,2*iter]) = Ad(:,[2*iter-1,2*iter]) * inv(sqD);
    num = Res(indexa,[2*iter-1,2*iter])'*Zed(indexa,[2*iter-1,2*iter]);
    num = sqD\num; % because of Zed and not d
    
%    Itere = Itere + d(:,[2*iter-1,2*iter])*alpha;
%    Res(:,[2*iter+1,2*iter+2]) = Res(:,[2*iter-1,2*iter]) - ...
%                                    Ad(:,[2*iter-1,2*iter])*alpha;
    Itere = Itere + d(:,[2*iter-1,2*iter])*num;
    Res(:,[2*iter+1,2*iter+2]) = Res(:,[2*iter-1,2*iter]) - ...
                                    Ad(:,[2*iter-1,2*iter])*num;
    
    residual(iter+1) = norm( Res(indexa,2*iter+1)-Res(indexa,2*iter+2) );
    error(iter+1)    = norm(Itere(indexa,1) - Itere(indexa,2) - fref(indexa)) / ...
                                    norm(fref(indexa));
    regulari(iter+1) = sqrt( (Itere(:,1)'-Itere(:,2)')* ... 
                          regul( Itere(:,1)-Itere(:,2) , nodes, boundary, 2) );
    
    if precond == 1
       % Solve 1 (Sd)
       f11 = dirichletRhs2( Res(:,2*iter+1), 2, c2node1, boundaryp1, nnodes );
       f12 = dirichletRhs2( Res(:,2*iter+2), 2, c2node1, boundaryp1, nnodes );
       uin1 = K1\[f11,f12];
       lagr1 = uin1(2*nnodes+1:end,:);
       lamb11 = lagr2forces2( lagr1(:,1), c2node1, 2, boundaryp1, nnodes );
       lamb12 = lagr2forces2( lagr1(:,2), c2node1, 2, boundaryp1, nnodes );
       %
       Zed(:,[2*iter+1,2*iter+2]) = [lamb11,lamb12];
    elseif precond == 2
       % Solve 2 (Sn)
       f21 = dirichletRhs2( Res(:,2*iter+1), 2, c2node2, boundaryp2, nnodes );
       f22 = dirichletRhs2( Res(:,2*iter+2), 2, c2node2, boundaryp2, nnodes );
       uin2 = K2\[f21,f22];
       lagr2 = uin2(2*nnodes+1:end,:);
       lamb21 = lagr2forces2( lagr2(:,1), c2node2, 2, boundaryp2, nnodes );
       lamb22 = lagr2forces2( lagr2(:,2), c2node2, 2, boundaryp2, nnodes );
       %
       Zed(:,[2*iter+1,2*iter+2]) = [lamb21,lamb22];
    else
        Zed(:,[2*iter+1,2*iter+2]) = Res(:,[2*iter+1,2*iter+2]);
    end
    
    % Ritz variables (Saad++)
    unsuralpha( [2*iter-1,2*iter], [2*iter-1,2*iter] ) = ...
                                    (sqD*num)^(-1/2)*den*(sqD*num)^(-1/2);
    
    if iter > 1
       betasuralpha( [2*iter-3,2*iter-2], [2*iter-3,2*iter-2] ) = ...
                         (sqD*num)^(-1/2) * betaij'*betaij * (sqD*num)^(-1/2);
                         % use betaij from the previous iteration
                         
       etaeta( [2*iter-1,2*iter], [2*iter-1,2*iter] ) = ...
               (denprec^(1/2)*numprec)^(-1/2) * denprec * ...
               inv( denprec^(1/2)*numprec ) * (sqD*num)^(1/2);
    end
    
    % First Reorthogonalize the residual (as we use it next), in sense of M
    for jter=1:iter-1
        betac = (Res(indexa,[2*jter-1,2*jter])'*Zed(indexa,[2*jter-1,2*jter])) \...
                (Res(indexa,[2*jter-1,2*jter])'*Zed(indexa,[2*iter+1,2*iter+2])) ;

        Zed(:,[2*iter+1,2*iter+2]) = Zed(:,[2*iter+1,2*iter+2]) - ...
                                      Zed(:,[2*jter-1,2*jter]) * betac;
        Res(:,[2*iter+1,2*iter+2]) = Res(:,[2*iter+1,2*iter+2]) - ...
                                      Res(:,[2*jter-1,2*jter]) * betac;
    end

    %% Orthogonalization
    d(:,[2*iter+1,2*iter+2]) = Zed(:,[2*iter+1,2*iter+2]);
    for jter=iter:iter
        betaij = ( Ad(indexa,[2*jter-1,2*jter])'*d(indexa,[2*jter-1,2*jter]) ) \ ...
            ( Ad(indexa,[2*jter-1,2*jter])'*d(indexa,[2*iter+1,2*iter+2]) );

        d(:,[2*iter+1,2*iter+2]) = d(:,[2*iter+1,2*iter+2]) - ...
                                   d(:,[2*jter-1,2*jter]) * betaij;
    end

    %% The Ritz elements
    V(indexa,[2*iter-1,2*iter]) = (-1)^(iter-1) * Zed(indexa,[2*iter-1,2*iter]) * ...
      ((Res(indexa,[2*iter-1,2*iter])'*Zed(indexa,[2*iter-1,2*iter]))^(-1/2)) ;
                       
    delta  = unsuralpha( [2*iter-1,2*iter], [2*iter-1,2*iter] ) ;
    if iter > 1
       delta = delta + betasuralpha( [2*iter-3,2*iter-2], [2*iter-3,2*iter-2] );
    end

    eta = etaeta( [2*iter-1,2*iter], [2*iter-1,2*iter] ); % what a stupid variable name
    
    if iter > 1
       H( [2*iter-1,2*iter] , [2*iter-3,2*iter-2,2*iter-1,2*iter] ) = ... 
                                                              [eta', delta];
       H( [2*iter-3,2*iter-2] , [2*iter-1,2*iter] ) = eta;
    else
       H( [2*iter-1,2*iter] , [2*iter-1,2*iter] ) = delta;
    end
    
%    Hd = V'*AV;  % Debug H
    
    % Compute eigenelems of the Hessenberg :
    [Q,Theta1] = eig(H);
    theta = diag(Theta1);
    % Sort it
    [theta,Ind] = sort(theta,'descend');
    Q = Q(:,Ind);
    Theta1 = Theta1(Ind,Ind);
    Y = V*Q;
    
end

% Compute the solution
chi = inv(Theta1)*Y'*( b(:,1)-b(:,2) );
if ntrunc > 0
   chi(ntrunc:end) = 0;
end
ItereR = Y*chi;

% Build residual and such
regD = zeros(niter,1); resD = zeros(niter,1); bt = Y'*b;
for i = 1:2*iter+1
   chiD   = inv(Theta1)*Y'*(b(:,1)-b(:,2)); chiD(i:end) = 0;
   ItereD = Y*chiD;
   resD(i) = sqrt( sum( bt(i:end).^2) );  
   regD(i) = sqrt( ItereD'*regul(ItereD, nodes, boundary, 2) );
end

% Compute the number of pertinent Ritz modes
[indm2,pol] = findPicard2(log10(abs(chiD)), 3, 1);
n = size(pol,1);
t = 1:.05:2*niter; tt = zeros(n,20*(2*niter-1)+1);
for j=1:n
   tt(j,:) = t.^(n-j);
end
px = pol'*tt;

figure;
hold on
plot(error,'Color','blue')
plot(residual,'Color','red')
legend('error','residual')
%% L-curve :
%loglog(residual(2:end),regulari(2:end));
%figure

figure;
hold on;
plot(log10(theta),'Color','blue')
plot(log10(abs(Y'*( b(:,1)-b(:,2) ))),'Color','red')
plot(log10(abs(chiD)),'Color','black')
plot(t,px,'Color','cyan')
legend( 'Ritz Values','RHS values','solution coefficients', ...
        'polynomial approximation' )

figure;
hold on;
plot(Itere(2*b2node2-1,1)-Itere(2*b2node2-1,2),'Color','red');
plot(ItereR(2*b2node2-1),'Color','blue');
plot(fref(2*b2node2-1),'Color','green');
%legend('solution','reference')
legend('solution','Ritz solution','reference')

%%%%
%% Final problem : compute u
dirichlet = [4,1,0;4,2,0;
             3,1,0;3,2,0;
             1,1,0;1,2,0];
neumann   = [];
[K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
f2    = [keepField( ItereR, 2, boundary ); zeros(nbloq,1)];
fdir4 = dirichletRhs(urefb, 4, C, boundary);
usoli = K \ (fdir4 + f2);

usol = usoli(1:2*nnodes,1);
fsol = Kinter*usol;

figure;
hold on;
plot(usol(2*b2node2-1),'Color','red');
plot(uref(2*b2node2-1),'Color','green');
legend('solution','reference')

total_error = norm(usol-uref)/norm(uref);
% Compute stress :
sigma = stress(usol,E,nu,nodes,elements,order,1,ntoelem);
plotGMSH({usol,'U_vect';sigma,'stress'}, elements, nodes, 'solution');
