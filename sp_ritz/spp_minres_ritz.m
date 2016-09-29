% 19/09/2016
% Algo Steklov-Poincaré primal avec MINRES 
% (passke Sd-Sn n'est >0 que de justesse)

close all;
clear all;

addpath(genpath('./tools'))
addpath(genpath('./sp_ritz'))

% Parameters
E       = 70000;  % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 1;      % N.mm-1 : Loading on the plate
niter   = 5;
precond = 0;      % 1 : Use a dual precond
mu      = 0.;     % Regularization parameter
ratio   = 5e-200;    % Maximal ratio (for eigenfilter)
br      = 0.;      % noise
brt     = 0;      % "translation" noise
epsilon = 1e-1;   % Convergence criterion for ritz value

% Boundary conditions
% first index  : index of the boundary
% second index : 1=x, 2=y
% third        : value
% [0,1,value] marks a dirichlet regularization therm on x
dirichlet = [1,1,0; 1,2,0 ;
             3,1,0; 3,2,0 ];
neumann   = [2,1,fscalar,0,fscalar;
             4,1,fscalar,0,-fscalar];

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

boundaryp1 = suppressBound( boundary, no2, 2 );
boundaryp1 = suppressBound( boundaryp1, no3, 2 );
boundaryp2 = boundaryp1;

[~, b2node2] = mapBound(2, boundaryp1, nnodes);
indexa = [2*b2node2-1; 2*b2node2];

% Then, build the stiffness matrix :
[K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
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

urefb = ( 1 + brt + br*randn(2*nnodes,1) ) .* uref;

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
K1r = K1; K2r = K2; K1dr = K1d; K2dr = K2d;
%K1(indexa,indexa) = 0;
%K2(indexa,indexa) = 0;
%K1d(indexa,indexa) = 0;
%K2d(indexa,indexa) = 0;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Conjugate Gradient for the problem : (S10-S20) x = S2-S1
Itere  = zeros( 2*nnodes, 1 );
d      = zeros( 2*nnodes, niter+1 );
Ad     = zeros( 2*nnodes, niter+1 );
Res    = zeros( 2*nnodes, niter+1 );
MRes   = zeros( 2*nnodes, niter+1 );
AMRes  = zeros( 2*nnodes, niter+1 );  % AMRes = A*M*Res
Zed    = zeros( 2*nnodes, niter+1 );
ntrunc = 0;  % In case the algo finishes at niter

%% Perform A x0 :
% Solve 1
f1 = dirichletRhs2( Itere, 2, c2node1, boundaryp1, nnodes );
uin1 = K1\f1;
lagr1 = uin1(2*nnodes+1:end,1);
lamb1 = lagr2forces( lagr1, C1, 2, boundaryp1 );
% Solve 2
f2 = dirichletRhs2( Itere, 2, c2node2, boundaryp2, nnodes );
uin2 = K2\f2;
lagr2 = uin2(2*nnodes+1:end,1);
lamb2 = lagr2forces( lagr2, C2, 2, boundaryp2 );
% Regularization term
Nu = regul(Itere, nodes, boundary, 2);
%
Axz = mu*Nu+lamb1-lamb2;
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
b = lamb2-lamb1;
%plot(b(2*b2node2-1));
%figure;
%%
Res(:,1) = b - Axz;

if precond == 1
    % Solve 1
    f1 = [Res(:,1); zeros(nbloq1d,1)];
    uin1 = K1dr\f1;
    u1i = uin1(1:2*nnodes,1);
    u1 = keepField( u1i, 2, boundaryp1 );
    d(:,1) = u1;
else
    d(:,1) = Res(:,1);
end

MRes(:,1) = d(:,1);

% Solve 1
rhs1 = d(:,1);
f1 = dirichletRhs(rhs1, 2, C1, boundaryp1);
uin1 = K1\f1;
lagr1 = uin1(2*nnodes+1:end,1);
lamb1 = lagr2forces( lagr1, C1, 2, boundaryp1 );
% Solve 2
rhs2 = d(:,1);
f2 = dirichletRhs(rhs2, 2, C2, boundaryp2);
uin2 = K2\f2;
lagr2 = uin2(2*nnodes+1:end,1);
lamb2 = lagr2forces( lagr2, C2, 2, boundaryp2 );
%
Ad(:,1) = lamb1-lamb2;

AMRes(:,1) = Ad(:,1);
 
residual(1) = sqrt( norm(Res( indexa,1)));
error(1)    = sqrt( norm(Itere(indexa) - uref(indexa)) / norm(uref(indexa)));
regulari(1) = sqrt( Itere'*regul(Itere, nodes, boundary, 2) );

ritzval  = 0; % Last ritz value that converged
oldtheta = 0;
eta      = 0;
getmeout = 0; % utility

%%
for iter = 1:niter
    
    if precond == 1
        % Solve 1
        f1 = [Ad(:,iter); zeros(nbloq1d,1)];
        uin1 = K1dr\f1;
        u1i = uin1(1:2*nnodes,1);
        u1 = keepField( u1i, 2, boundaryp1 );
        Zed(:,iter) = u1;
    else
        Zed(:,iter) = Ad(:,iter);
    end
    
    den = (Ad(indexa,iter)'*Zed(indexa,iter));
    d(:,iter) = d(:,iter)/sqrt(den); Ad(:,iter) = Ad(:,iter)/sqrt(den); 
    Zed(:,iter) = Zed(:,iter)/sqrt(den);
    
    num = Res(indexa,iter)'*Zed(indexa,iter);
    Itere          = Itere + d(:,iter)*num;%/den;
    Res(:,iter+1)  = Res(:,iter) - Ad(:,iter)*num;%/den;
    MRes(:,iter+1) = MRes(:,iter) - Zed(:,iter)*num;%/den;
    
    residual(iter+1) = sqrt( norm(Res(indexa,iter+1)));
    error(iter+1)    = sqrt( norm(Itere(indexa) - uref(indexa)) / norm(uref(indexa)));
    regulari(iter+1) = sqrt( Itere'*regul(Itere, nodes, boundary, 2) );

    %% Orthogonalization
    d(:,iter+1) = Zed(:,iter); % Ad(:,iter+1);
    
    % Solve 1
    rhs1 = d(:,iter+1);
    f1 = dirichletRhs(rhs1, 2, C1, boundaryp1);
    uin1 = K1\f1;
    lagr1 = uin1(2*nnodes+1:end,1);
    lamb1 = lagr2forces( lagr1, C1, 2, boundaryp1 );
    % Solve 2
    rhs2 = d(:,iter+1);
    f2 = dirichletRhs(rhs2, 2, C2, boundaryp2);
    uin2 = K2\f2;
    lagr2 = uin2(2*nnodes+1:end,1);
    lamb2 = lagr2forces( lagr2, C2, 2, boundaryp2 );
    %
    Ad(:,iter+1) = lamb1-lamb2;
    
    AMRes(:,iter+1) = AMRes(:,iter) - Ad(:,iter+1)*num;%/den;  % For AV. /!\ in case of precond
    
%    for jter=1:iter   % Reorthogonalize the residual in sense of AM
%        betac = Res(indexa,iter+1)'*AMRes(indexa,jter) / (Res(indexa,jter)'*AMRes(indexa,jter));
%        MRes(:,iter+1) = MRes(:,iter+1) - MRes(:,jter) * betac;
%        Res(:,iter+1)  = Res(:,iter+1) - Res(:,jter) * betac;
%        AMRes(:,iter+1) = AMRes(:,iter+1) - AMRes(:,jter) * betac;
%    end
    
    for jter = 1:iter % max(1,iter-1):iter 
        betaij = ( Ad(indexa,iter+1)'*Zed(indexa,jter) ) / ( Ad(indexa,jter)'*Zed(indexa,jter) );
        d(:,iter+1)  = d(:,iter+1) - d(:,jter) * betaij;
        Ad(:,iter+1) = Ad(:,iter+1) - Ad(:,jter) * betaij;
    end
    
    %% Ritz algo : find the Ritz elements
    % Build the matrices
    V(:,iter) = zeros(2*nnodes,1);
    V(indexa,iter) = (-1)^(iter-1)*MRes(indexa,iter)/(sqrt(Res(indexa,iter)'*AMRes(indexa,iter)));
    Vo(:,iter) = zeros(2*nnodes,1);
    Vo(indexa,iter) = (-1)^(iter-1)*Res(indexa,iter)/(sqrt(Res(indexa,iter)'*AMRes(indexa,iter)));
    
    %norm(Zed(indexa,iter))

    if iter > 1
       Hb(iter,iter-1) = - AMRes(indexa,iter)'*AMRes(indexa,iter-1)/ ...
                             ( sqrt(Res(indexa,iter)'*AMRes(indexa,iter)) * ...
                             sqrt(Res(indexa,iter-1)'*AMRes(indexa,iter-1)) );
       Hb(iter-1,iter) = - AMRes(indexa,iter)'*AMRes(indexa,iter-1)/ ...
                             ( sqrt(Res(indexa,iter)'*AMRes(indexa,iter)) * ...
                             sqrt(Res(indexa,iter-1)'*AMRes(indexa,iter-1)) );
    end
    %H(iter,iter)  = AMRes(indexa,iter)'*Res(indexa,iter) / ...
                                %(Res(indexa,iter)'*AMRes(indexa,iter));  %= 1
    Hb(iter,iter) = AMRes(indexa,iter)'*AMRes(indexa,iter) / ...
                                (Res(indexa,iter)'*AMRes(indexa,iter));
    
    % Compute eigenelems of the Hessenberg :
    [Q,Theta1] = eig(Hb);
    theta = diag(Theta1);
    % Sort it
    [theta,Ind] = sort(theta,'descend');
    Q = Q(:,Ind);
    Theta1 = Theta1(Ind,Ind);
    Y = V*Q;
    
    % See if the current one converged
    if abs(theta(ritzval+1) - oldtheta) < epsilon*oldtheta
       % increment oldtheta
       ritzval = ritzval+1;
       if size(theta,1) > ritzval
          oldtheta = theta(ritzval+1);
       else
          oldtheta = 0;
       end
       
       % Check small value / hold value
       if ritzval > 1
          if theta(ritzval) < ratio*theta(1)
             ntrunc = ritzval
             getmeout = 1;
             break;
          end
       end
    else
       oldtheta = theta(ritzval+1);
    end
%    if getmeout == 1  % In case I use a while over there
%       break;
%    end
end

% Debug : compute He
AV = zeros( 2*nnodes, niter );
for i=1:niter
    rhs1 = V(:,i);
    f1 = dirichletRhs(rhs1, 2, C1, boundaryp1);
    uin1 = K1\f1;
    lagr1 = uin1(2*nnodes+1:end,1);
    lamb1 = lagr2forces( lagr1, C1, 2, boundaryp1 );
    % Solve 2
    rhs2 = V(:,i);
    f2 = dirichletRhs(rhs2, 2, C2, boundaryp2);
    uin2 = K2\f2;
    lagr2 = uin2(2*nnodes+1:end,1);
    lamb2 = lagr2forces( lagr2, C2, 2, boundaryp2 );
    %
    AV(:,i) = lamb1-lamb2;
end
He = V'*AV;
Heb = AV'*AV;

% Compute the solution
chi = Y'*b;%inv(Theta1)*Y'*b;
if ntrunc > 0
   chi(ntrunc:end) = 0;
end
ItereR = Y*chi;

%hold on;
%plot(log10(theta),'Color','blue')
%plot(log10(abs(Y'*b)),'Color','red')
%plot(log10(abs(chi)),'Color','black')
%legend('Ritz Values','RHS values','solution coefficients')
%figure;
%
hold on;
plot(Itere(2*b2node2-1),'Color','red')
plot(ItereR(2*b2node2-1),'Color','blue')
plot(uref(2*b2node2-1),'Color','green')
legend('brutal solution','filtred solution', 'reference')
figure;

hold on
plot(log10(error(2:end)),'Color','blue')
plot(log10(residual(2:end)),'Color','red')
legend('error (log)','residual (log)')
%figure;
% L-curve :
%loglog(residual(2:end),regulari(2:end));
%figure
%%%%%
%% Final problem : compute u
dirichlet = [4,1,0;4,2,0;
             3,1,0;3,2,0;
             1,1,0;1,2,0;
             2,1,0;2,2,0];
neumann   = [];
[K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
fdir2 = dirichletRhs(Itere, 2, C, boundary);
fdir4 = dirichletRhs(urefb, 4, C, boundary);
usoli = K \ (fdir4 + fdir2);

usol = usoli(1:2*nnodes,1);
fsol = Kinter*usol;

% With the reduced solution
dirichlet = [4,1,0;4,2,0;
             3,1,0;3,2,0;
             1,1,0;1,2,0;
             2,1,0;2,2,0];
neumann   = [];
[K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
fdir2 = dirichletRhs(ItereR, 2, C, boundary);
fdir4 = dirichletRhs(urefb, 4, C, boundary);
usoli = K \ (fdir4 + fdir2);

usolR = usoli(1:2*nnodes,1);
fsolR = Kinter*usolR;

%hold on;
%plot(fsol(2*b2node2-1),'Color','red')
%plot(fsolR(2*b2node2-1),'Color','blue')
%plot(f(2*b2node2-1),'Color','green')
%legend('brutal solution','filtred solution','reference')

total_error = norm(usol-uref)/norm(uref);
% Compute stress :
sigma = stress(usol,E,nu,nodes,elements,order,1,ntoelem);
plotGMSH({usol,'U_vect';sigma,'stress'}, elements, nodes, 'solution');