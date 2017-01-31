% 30/09/2016
% Algo Steklov-Poincaré primal avec RC 
% (passke Sd-Sn n'est >0 que de justesse)

close all;
clear all;

addpath(genpath('./tools'))
addpath(genpath('./sp_ritz'))

% Parameters
E       = 70000;  % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 1;      % N.mm-1 : Loading on the plate
niter   = 15;
precond = 1;      % 1 : Use a dual precond
mu      = 0.;     % Regularization parameter
ratio   = 5e-200; % Maximal ratio (for eigenfilter)
br      = 0.1;   % noise
brt     = 0;      % "translation" noise
epsilon = 1e-20;  % Convergence criterion for ritz value
ntrunc  = 7;      % In case the algo finishes at niter
inhomog = 0;      % inhomogeneous medium

if inhomog == 2  % load previously stored matrix
   mat = [0, E, nu]; % dummy material
   Kinter = load('./noises/stocrig1.mat'); Kinter = Kinter.Kinter;
elseif inhomog == 1  % /!\ IN THE INHOMOGENEOUS CASE, ALL THE SIGMAS ARE WRONG
   mat = [2, E, nu, .1, 1];
else
   mat = [0, E, nu];
end

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
boundaryp2 = boundaryp1;

[~, b2node2] = mapBound(2, boundaryp1, nnodes);
indexa = [2*b2node2-1; 2*b2node2];

% Then, build the stiffness matrix :
[K,C,nbloq] = Krig2 (nodes,elements,mat,order,boundary,dirichlet);
if inhomog == 2
   K(1:2*nnodes, 1:2*nnodes) = Kinter;
else
   Kinter = K(1:2*nnodes, 1:2*nnodes);
end
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

urefb = ( 1 + brt + br*noise ) .* uref;

sigma = stress(uref,E,nu,nodes,elements,order,1,ntoelem);
plotGMSH({uref,'Vect_U';sigma,'stress'}, elements, nodes, 'reference');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition of the operators

% First problem
dirichlet1 = [4,1,0;4,2,0;
              3,1,0;3,2,0;
              2,1,0;2,2,0;
              1,1,0;1,2,0];

[K1,C1,nbloq1,node2c1,c2node1] = Krig2 (nodes,elements,mat,order,boundaryp1,dirichlet1);
if inhomog >= 1  % Because of the random stuff
   K1(1:2*nnodes, 1:2*nnodes) = Kinter;
end
% Second problem
dirichlet2 = [1,1,0;1,2,0;
              2,1,0;2,2,0;
              3,1,0;3,2,0];
neumann2   = [4,1,fscalar,0,-fscalar];
neumann0   = [];
[K2,C2,nbloq2,node2c2,c2node2] = Krig2 (nodes,elements,mat,order,boundaryp2,dirichlet2);
if inhomog >= 1
   K2(1:2*nnodes, 1:2*nnodes) = Kinter;
end

error    = zeros(niter+1,1);
residual = zeros(niter+1,1);
regulari = zeros(niter+1,1);

%% Dual problems
% First problem
dirichlet1d = [4,1,0;4,2,0;
               3,1,0;3,2,0;
               1,1,0;1,2,0];
[K1d,C1d,nbloq1d] = Krig2 (nodes,elements,mat,order,boundary,dirichlet1d);
if inhomog >= 1
   K1d(1:2*nnodes, 1:2*nnodes) = Kinter;
end
% Second problem
dirichlet2d = [1,1,0;1,2,0;
               3,1,0;3,2,0];
[K2d,C2d,nbloq2d] = Krig2 (nodes,elements,mat,order,boundary,dirichlet2d);
if inhomog >= 1
   K2d(1:2*nnodes, 1:2*nnodes) = Kinter;
end

%% Anti-cancellation trick
K1r = K1; K2r = K2; K1dr = K1d; K2dr = K2d;
%K1(indexa,indexa) = 0;
%K2(indexa,indexa) = 0;
%K1d(indexa,indexa) = 0;
%K2d(indexa,indexa) = 0;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Conjugate Residual for the problem : (S10-S20) x = S2-S1
Itere  = zeros( 2*nnodes, 1 );
d      = zeros( 2*nnodes, niter+1 );
Ad     = zeros( 2*nnodes, niter+1 );
Res    = zeros( 2*nnodes, niter+1 );
Zed    = zeros( 2*nnodes, niter+1 );
AZed   = zeros( 2*nnodes, niter+1 );

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
    Zed(:,1) = u1;
else
    Zed(:,1) = Res(:,1);
end

d(:,1) = Zed(:,1);

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

AZed(:,1) = Ad(:,1);
 
residual(1) = sqrt( norm(Res( indexa,1)));
error(1)    = sqrt( norm(Itere(indexa) - uref(indexa)) / norm(uref(indexa)));
regulari(1) = sqrt( Itere'*regul(Itere, nodes, boundary, 2) );

ritzval  = 0; % Last ritz value that converged
oldtheta = 0;
eta      = 0;
getmeout = 0; % utility

V(:,1) = zeros(2*nnodes,1);
V(indexa,1) = Zed(indexa,1)/(sqrt(Zed(indexa,1)'*AZed(indexa,1)));

Hb(1,1) = AZed(indexa,1)'*AZed(indexa,1) / (Zed(indexa,1)'*AZed(indexa,1));

%%
for iter = 1:niter
    
    den = (Ad(indexa,iter)'*Ad(indexa,iter));
    d(:,iter) = d(:,iter)/sqrt(den); Ad(:,iter) = Ad(:,iter)/sqrt(den); 
    
    num            = Res(indexa,iter)'*Ad(indexa,iter);
    Itere          = Itere + d(:,iter)*num;%/den;
    Res(:,iter+1)  = Res(:,iter) - Ad(:,iter)*num;%/den;

    if precond == 1
        % Solve 1
        f1 = [Res(:,iter+1); zeros(nbloq1d,1)];
        uin1 = K1dr\f1;
        u1i = uin1(1:2*nnodes,1);
        u1 = keepField( u1i, 2, boundaryp1 );
        Zed(:,iter+1) = u1;
    else
        Zed(:,iter+1) = Res(:,iter+1);
    end
    
    residual(iter+1) = sqrt( norm(Res(indexa,iter+1)));
    error(iter+1)    = sqrt( norm(Itere(indexa) - uref(indexa)) / norm(uref(indexa)));
    regulari(iter+1) = sqrt( Itere'*regul(Itere, nodes, boundary, 2) );

    %% Orthogonalization
    d(:,iter+1) = Zed(:,iter+1);
    
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
    
    AZed(:,iter+1) = Ad(:,iter+1); % for AV
    
    for jter=1:iter   % Reorthogonalize the residual in sense of A
        betac = Zed(indexa,iter+1)'*AZed(indexa,jter) / (Zed(indexa,jter)'*AZed(indexa,jter));
        Zed(:,iter+1) = Zed(:,iter+1) - Zed(:,jter) * betac;
        Res(:,iter+1)  = Res(:,iter+1) - Res(:,jter) * betac;
        AZed(:,iter+1) = AZed(:,iter+1) - AZed(:,jter) * betac;
    end

    for jter = 1:iter % max(1,iter-1):iter 
        betaij = ( Ad(indexa,iter+1)'*Ad(indexa,jter) );% / ( Ad(indexa,jter)'*Ad(indexa,jter) ) = 1;
        d(:,iter+1)  = d(:,iter+1) - d(:,jter) * betaij;
        Ad(:,iter+1) = Ad(:,iter+1) - Ad(:,jter) * betaij;
    end
    
    %% Ritz algo : find the Ritz elements
    % Build the matrices
    V(:,iter+1) = zeros(2*nnodes,1);
    V(indexa,iter+1) = (-1)^(iter)*Zed(indexa,iter+1) / ...
                            (sqrt(Zed(indexa,iter+1)'*AZed(indexa,iter+1)));
    
    %norm(Zed(indexa,iter))

    Hb(iter+1,iter) = - AZed(indexa,iter+1)'*AZed(indexa,iter)/ ...
                          ( sqrt(Zed(indexa,iter+1)'*AZed(indexa,iter+1)) * ...
                          sqrt(Zed(indexa,iter)'*AZed(indexa,iter)) );
    Hb(iter,iter+1) = - AZed(indexa,iter+1)'*AZed(indexa,iter)/ ...
                          ( sqrt(Zed(indexa,iter+1)'*AZed(indexa,iter+1)) * ...
                          sqrt(Zed(indexa,iter)'*AZed(indexa,iter)) );
                          
    Hb(iter+1,iter+1) = AZed(indexa,iter+1)'*AZed(indexa,iter+1) / ...
                                   (Zed(indexa,iter+1)'*AZed(indexa,iter+1));
    
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

%% Debug : compute He
%AV = zeros( 2*nnodes, niter );
%for i=1:iter+1
%    rhs1 = V(:,i);
%    f1 = dirichletRhs(rhs1, 2, C1, boundaryp1);
%    uin1 = K1\f1;
%    lagr1 = uin1(2*nnodes+1:end,1);
%    lamb1 = lagr2forces( lagr1, C1, 2, boundaryp1 );
%    % Solve 2
%    rhs2 = V(:,i);
%    f2 = dirichletRhs(rhs2, 2, C2, boundaryp2);
%    uin2 = K2\f2;
%    lagr2 = uin2(2*nnodes+1:end,1);
%    lamb2 = lagr2forces( lagr2, C2, 2, boundaryp2 );
%    %
%    AV(:,i) = lamb1-lamb2;
%end
%He = V'*AV;
%Heb = AV'*AV;

% First, do A'b (=Ab) : x is the solution of A'Ax = A'b

% Solve 1
rhs1 = b;
f1 = dirichletRhs(rhs1, 2, C1, boundaryp1);
uin1 = K1\f1;
lagr1 = uin1(2*nnodes+1:end,1);
lamb1 = lagr2forces( lagr1, C1, 2, boundaryp1 );
% Solve 2
rhs2 = b;
f2 = dirichletRhs(rhs2, 2, C2, boundaryp2);
uin2 = K2\f2;
lagr2 = uin2(2*nnodes+1:end,1);
lamb2 = lagr2forces( lagr2, C2, 2, boundaryp2 );
%
Ab = lamb1-lamb2;

% Compute the solution
chi = Y'*b;%inv(Theta1)*Y'*b;
%chi = inv(Theta1)*Y'*Ab;
%chi = inv(Theta1)*Y'*Y*Y'*b;%inv(Theta1)*Y'*b;
if ntrunc > 0
   chi(ntrunc:end) = 0;
end
ItereR = Y*chi;

regS = zeros(niter,1);
resS = zeros(niter,1);
%% Build the L-curve regul, ntrunc
for i = 1:iter+1
%   chiS   = inv(Theta1)*Y'*Ab; chiS(i:end) = 0;
   chiS   = Y'*b; chiS(i:end) = 0;
   ItereS = Y*chiS;
   
   % Solve 1
   rhs1 = ItereS;
   f1 = dirichletRhs(rhs1, 2, C1, boundaryp1);
   uin1 = K1\f1;
   lagr1 = uin1(2*nnodes+1:end,1);
   lamb1 = lagr2forces( lagr1, C1, 2, boundaryp1 );
   % Solve 2
   rhs2 = ItereS;
   f2 = dirichletRhs(rhs2, 2, C2, boundaryp2);
   uin2 = K2\f2;
   lagr2 = uin2(2*nnodes+1:end,1);
   lamb2 = lagr2forces( lagr2, C2, 2, boundaryp2 );
   % Regularization term
   AI = lamb1-lamb2;
   
   ResS = AI-b;
   resS(i) = sqrt(norm(ResS));   
   regS(i) = sqrt( ItereS'*regul(ItereS, nodes, boundary, 2) );
end

hold on;
plot(log10(theta),'Color','blue')
plot(log10(abs(Y'*b)),'Color','red')
plot(log10(abs(chi)),'Color','black')
legend('Ritz Values','RHS values','solution coefficients')
figure;
%
hold on;
plot(Itere(2*b2node2-1),'Color','red')
plot(ItereR(2*b2node2-1),'Color','blue')
plot(uref(2*b2node2-1),'Color','green')
legend('brutal solution','filtred solution', 'reference')
figure;
%hold on;
%plot(Itere(2*b2node2),'Color','red')
%plot(ItereR(2*b2node2),'Color','blue')
%plot(uref(2*b2node2),'Color','green')
%legend('brutal solution','filtred solution', 'reference')
%figure;

hold on
plot(log10(error(2:end)),'Color','blue')
plot(log10(residual(2:end)),'Color','red')
legend('error (log)','residual (log)')
figure;

regD = zeros(niter,1); resD = zeros(niter,1); bt = Y'*b;
%errD = zeros(niter,1);
for i = 1:iter+1
%   chiD   = inv(Theta1)*Y'*Ab; chiD(i:end) = 0;
   chiD   = Y'*b; chiD(i:end) = 0;
   ItereD = Y*chiD;
   resD(i) = sqrt( sum( bt(i:end).^2) );  
   regD(i) = sqrt( ItereD'*regul(ItereD, nodes, boundary, 2) );
   %errD(i) = norm(ItereD(indexa) - fref(indexa)) / norm(fref(indexa));
end
% RL-curve
loglog(resD(2:iter),regD(2:iter),'-+');
legend('RL-curve')
figure

hold on;
set(gca, 'fontsize', 20);
loglog(residual(2:end),regulari(2:end),'-*','Color','red','linewidth',3);
loglog(resS(2:iter),regS(2:iter),'-+','linewidth',3);
legend('L-curve','RL-curve')
%legend('L-curve')
xlabel('residual')
ylabel('H1 norm')

%%%%%
%% Final problem : compute u
dirichlet = [4,1,0;4,2,0;
             3,1,0;3,2,0;
             1,1,0;1,2,0;
             2,1,0;2,2,0];
neumann   = [];
[K,C,nbloq] = Krig2 (nodes,elements,mat,order,boundary,dirichlet);
if inhomog >= 1
   K(1:2*nnodes, 1:2*nnodes) = Kinter;
end
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
[K,C,nbloq] = Krig2 (nodes,elements,mat,order,boundary,dirichlet);
if inhomog >= 1
   K(1:2*nnodes, 1:2*nnodes) = Kinter;
end
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