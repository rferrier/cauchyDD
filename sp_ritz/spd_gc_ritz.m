% 03/10/2016
% Algo Steklov-Poincaré dual avec Gradient Conjugué

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 70000;  % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 1;      % N.mm-1 : Loading on the plate
niter   = 10;
precond = 0;      % 1/2 : Use a dual precond, 3 : use H1/2 precond, 4 : use gradient precond
mu      = 0.;     % Regularization parameter
ratio   = 1e-300;    % Maximal ratio (for eigenfilter)
br      = .1;      % noise
brt     = 0;      % "translation" noise
epsilon = 1e-200;   % Convergence criterion for ritz value
ntrunc  = 0;      % In case the algo finishes at niter
inhomog = 0;      % inhomogeneous medium

if inhomog == 2  % load previously stored matrix
   mat = [2, E, nu, .1, 1];
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

noises = load('./noises/noise0.mat'); % Particular noise vector
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

urefb = ( 1 + brt + br*noise ) .* uref;
fref  = Kinter*uref;

sigma = stress(uref,E,nu,nodes,elements,order,1,ntoelem);
plotGMSH({uref,'Vect_U';sigma,'stress'}, elements, nodes, 'output/reference');

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
regulari2 = zeros(niter+1,1);

%% Dual problems
% First problem
dirichlet1d = [4,1,0;4,2,0;
               3,1,0;3,2,0;
               1,1,0;1,2,0];
[K1d,C1d,nbloq1d,node2c1d,c2node1d] = Krig2 (nodes,elements,mat,order,boundary,dirichlet1d);
if inhomog >= 1
   K1d(1:2*nnodes, 1:2*nnodes) = Kinter;
end

% Second problem
dirichlet2d = [1,1,0;1,2,0;
               3,1,0;3,2,0];
[K2d,C2d,nbloq2d,node2c2d,c2node2d] = Krig2 (nodes,elements,mat,order,boundary,dirichlet2d);
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
%% Conjugate Gradient for the problem : (D20-D10) x = D1-D2
Itere = zeros( 2*nnodes, 1 );
d     = zeros( 2*nnodes, niter+1 );
Ad    = zeros( 2*nnodes, niter+1 );
Res   = zeros( 2*nnodes, niter+1 );
Zed   = zeros( 2*nnodes, niter+1 );
alpha = zeros( niter+1, 1 );
beta  = zeros( niter+1, 1 );
alpha2 = zeros( niter+1, 1 );
H12   = H1demi(size(Res,1), nodes, boundary, 2 );
Mgr   = Mgrad(size(Res,1), nodes, boundary, 2 );

%% Perform A x0 :
% Solve 1
f1 = [Itere; zeros(nbloq1d,1)];
uin1 = K1d\f1;
u1i = uin1(1:2*nnodes,1);
u1 = keepField( u1i, 2, boundaryp1 );
% Solve 2
f2 = [Itere; zeros(nbloq2d,1)];
uin2 = K2d\f2;
u2i = uin2(1:2*nnodes,1);
u2 = keepField( u2i, 2, boundaryp2 );
% Regularization term
Nu = regul(Itere, nodes, boundary, 2);
%
Axz = mu*Nu+u2-u1;
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
b = u1-u2;
%plot(b(2*b2node2-1));
%figure;
%%
Res(:,1) = b - Axz;

if precond == 1
    % Solve 1
    f1 = dirichletRhs2( Res(:,1)/2, 2, c2node1, boundaryp1, nnodes );
    uin1 = K1\f1;
    lagr1 = uin1(2*nnodes+1:end,1);
    lamb1 = lagr2forces2( lagr1, c2node1, 2, boundaryp1, nnodes );
    %
    Zed(:,1) = lamb1/2;
elseif precond == 2
    % Solve 2
    f2 = dirichletRhs2( Res(:,1)/2, 2, c2node2, boundaryp2, nnodes );
    uin2 = K2\f2;
    lagr2 = uin2(2*nnodes+1:end,1);
    lamb2 = lagr2forces2( lagr2, c2node2, 2, boundaryp2, nnodes );
    %
    Zed(:,1) = lamb2/2;
elseif precond == 3
    Zed(index,1) = E*H12(index,index)*Res(index,1);
elseif precond == 4
    Zed(index,1) = 1/E*Mgr(index,index)\Res(index,1);
else
    Zed(:,1) = Res(:,1);
end

d(:,1) = Zed(:,1);

residual(1) = norm(Res( indexa,1));
error(1)    = norm(Itere(indexa) - fref(indexa)) / norm(fref(indexa));
regulari(1) = sqrt( Itere'*regul(Itere, nodes, boundary, 2) );
regulari2(1) = norm( Itere(indexa) );

ritzval  = 0; % Last ritz value that converged
oldtheta = 0;
eta      = 0;
getmeout = 0; % utility
%V = zeros(2*nnodes, iter);
%H = zeros(iter);
%%
EA2=[0]; XM2 = [0];
pMx = 0; xM2 = 0; eA2 = 0;
for iter = 1:niter

    if iter>1
       xM2 = xM2 + num^2/den*pM2 + 2*num/sqrt(den)*pMx;
       XM2 = [XM2,xM2];
       pMx = -num*betaij/den*pM2 - betaij/sqrt(den)*pMx; %-beta(iter-1)
    end

    %% Optimal step
    % Solve 1
    f1 = [d(:,iter); zeros(nbloq1d,1)];
    uin1 = K1d\f1;
    u1i = uin1(1:2*nnodes,1);
    u1 = keepField( u1i, 2, boundaryp1 );
    % Solve 2
    f2 = [d(:,iter); zeros(nbloq2d,1)];
    uin2 = K2d\f2;
    u2i = uin2(1:2*nnodes,1);
    u2 = keepField( u2i, 2, boundaryp2 );
    % Regularization term
    Nu = regul(d(:,iter), nodes, boundary, 2);
    %
    Ad(:,iter) = mu*Nu+u2-u1;
    
    den = (d(indexa,iter)'*Ad(indexa,iter));
    d(:,iter) = d(:,iter)/sqrt(den); Ad(:,iter) = Ad(:,iter)/sqrt(den);
    num = Res(indexa,iter)'*d(indexa,iter);
    
    if iter==1
       pM2 = num*sqrt(den);
    else
       pM2 = num*sqrt(den) + betaij^2/denp * pM2; betaij^2/denp; %-beta(iter-1)
    end
    eA2 = eA2 - num^2;
    EA2 = [EA2, eA2];

    Itere         = Itere + d(:,iter)*num;%/den;
    Res(:,iter+1) = Res(:,iter) - Ad(:,iter)*num;%/den;
    
    residual(iter+1) = norm(Res(indexa,iter+1));
    error(iter+1)    = norm(Itere(indexa) - fref(indexa)) / norm(fref(indexa));
    regulari(iter+1) = sqrt( Itere'*regul(Itere, nodes, boundary, 2) );
    regulari2(iter+1) = norm( Itere(indexa) );
    
    if precond == 1
        % Solve 1
        f1 = dirichletRhs2( Res(:,iter+1)/2, 2, c2node1, boundaryp1, nnodes );
        uin1 = K1\f1;
        lagr1 = uin1(2*nnodes+1:end,1);
        lamb1 = lagr2forces2( lagr1, c2node1, 2, boundaryp1, nnodes );
        %
        Zed(:,iter+1) = lamb1/2;
    elseif precond == 2
        % Solve 2
        f2 = dirichletRhs2( Res(:,iter+1)/2, 2, c2node2, boundaryp2, nnodes );
        uin2 = K2\f2;
        lagr2 = uin2(2*nnodes+1:end,1);
        lamb2 = lagr2forces2( lagr2, c2node2, 2, boundaryp2, nnodes );
        %
        Zed(:,iter+1) = lamb2/2;
    elseif precond == 3
        Zed(index,iter+1) = E*H12(index,index)*Res(index,iter+1);
    elseif precond == 4
        Zed(index,iter+1) = 1/E*Mgr(index,index)\Res(index,iter+1);
    else
        Zed(:,iter+1) = Res(:,iter+1);
    end
    
    % Needed values for the Ritz stuff
    alpha(iter) = num/sqrt(den);
    beta(iter)  = - Zed(indexa,iter+1)'*Ad(indexa,iter)/sqrt(den);
%    alpha(iter) = Res(indexa,iter)'*Res(indexa,iter) / den;
%    beta(iter)  = Zed(indexa,iter+1)'*Res(indexa,iter+1) /... 
%                                (Zed(indexa,iter)'*Res(indexa,iter));
    
    % First Reorthogonalize the residual (as we use it next), in sense of M
    for jter=1:iter
        betac = Zed(indexa,iter+1)'*Res(indexa,jter) / (Zed(indexa,jter)'*Res(indexa,jter));
        Zed(:,iter+1) = Zed(:,iter+1) - Zed(:,jter) * betac;
        Res(:,iter+1) = Res(:,iter+1) - Res(:,jter) * betac;
    end
    
    %% Orthogonalization
    d(:,iter+1) = Zed(:,iter+1);
    
    for jter=iter:iter % No need to reorthogonalize (see above). do it anyway
        betaij = ( Zed(indexa,iter+1)'*Ad(indexa,jter) );%/...
            %( d(indexa,jter)'*Ad(indexa,jter) );
        d(:,iter+1) = d(:,iter+1) - d(:,jter) * betaij;
    end
    
    %% Ritz algo : find the Ritz elements
    % Build the matrices
    V(:,iter) = zeros(2*nnodes,1);
    V(indexa,iter) = (-1)^(iter-1)*Zed(indexa,iter)/(sqrt(Res(indexa,iter)'*Zed(indexa,iter)));
    %norm(Zed(indexa,iter))
    etap   = eta;
    delta  = 1/alpha(iter);
    if iter > 1
       delta  = delta + beta(iter-1)/alpha(iter-1);
    end
    eta    = sqrt(beta(iter))/alpha(iter);
    
    if iter > 1
       H(iter,[iter-1,iter]) = [etap, delta];
       H(iter-1,iter)        = etap;
    else
       H(iter,iter) = delta;
    end
    
    % Compute eigenelems of the Hessenberg :
    [Q,Theta1] = eig(H);
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
   denp = den;
end

% Compute the solution
chi = inv(Theta1)*Y'*b;
if ntrunc > 0
   chi(ntrunc:end) = 0;
end
ItereR = Y*chi;

%%for i=1:niter
%   plot(V(2*b2node2-1,20:25))
%   figure;
%%end

regS = zeros(niter,1);
resS = zeros(niter,1);
%% Build the L-curve regul, ntrunc
for i = 1:iter+1
   chiS   = inv(Theta1)*Y'*b; chiS(i:end) = 0;
   ItereS = Y*chiS;
   
   % Solve 1
   f1 = [ItereS; zeros(nbloq1d,1)];
   uin1 = K1d\f1;
   u1i = uin1(1:2*nnodes,1);
   u1 = keepField( u1i, 2, boundaryp1 );
   % Solve 2
   f2 = [ItereS; zeros(nbloq2d,1)];
   uin2 = K2d\f2;
   u2i = uin2(1:2*nnodes,1);
   u2 = keepField( u2i, 2, boundaryp2 );
   %
   AI = u2-u1;
   
   ResS = AI-b;
   resS(i) = norm(ResS);   
   regS(i) = sqrt( ItereS'*regul(ItereS, nodes, boundary, 2) );
end

%% Build Pieter's RL-curve
EA2R = [0]; XM2R = [0]; ratio = [0];
for i=1:iter
   ui = (Y(:,i)'*Res(:,1)) / sqrt((Zed(:,1)'*Res(:,1)));
   EA2R = [ EA2R, EA2R(end) - ui^2/Theta1(i,i) ] ;
   XM2R = [ XM2R, XM2R(end) + ui^2/Theta1(i,i)^2 ] ;
   %ratio = [ratio, 
end
EA2R = (Zed(:,1)'*Res(:,1))*EA2R;
XM2R = (Zed(:,1)'*Res(:,1))*XM2R; % Coefficient (for offset)

traD = 1 - sum(theta(1:ntrunc-1))/sum(theta);
regD = zeros(niter,1); resD = zeros(niter,1); bt = Y'*b;
for i = 1:iter+1
   chiD   = inv(Theta1)*Y'*b; chiD(i:end) = 0;
   ItereD = Y*chiD;
   resD(i) = sqrt( sum( bt(i:end).^2) );  
   regD(i) = sqrt( ItereD'*regul(ItereD, nodes, boundary, 2) );
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
plot(f(2*b2node2-1),'Color','green')
legend('brutal solution','filtred solution', 'reference')
figure;

%hold on;
%plot(Itere(2*b2node2),'Color','red')
%plot(ItereR(2*b2node2),'Color','blue')
%plot(f(2*b2node2),'Color','green')
%legend('brutal solution','filtred solution', 'reference')
%figure;

hold on
plot(error(2:end),'Color','blue')
%plot(log10(residual(2:end)),'Color','red')
legend('error')%legend('error (log)','residual (log)')
figure;
%L-curve :
hold on;
loglog(residual(2:iter+1),regulari2(2:iter+1),'Color','red','-*');
%figure
loglog(resS(2:iter+1),regS(2:iter+1),'-+');
legend('L-curve','RL-curve')
%findCorner (residual(2:iter+1), regulari(2:iter+1), 3)
%findCorner (resS(2:iter+1), regS(2:iter+1), 3)

%% Pieter's L-curves
EA  = sqrt(EA2-min(EA2));   XM  = sqrt(XM2);
EAR = sqrt(EA2R-min(EA2R)); XMR = sqrt(XM2R);
figure;
hold on;
loglog( EA(2:end-1), XM(2:end), 'Color', 'red', '-+', 'linewidth', 3 );
loglog( EAR(2:end-1), XMR(2:end-1), '-+', 'linewidth', 3 );
%plot( EA2(2:end-1)-min(EA2), XM2(2:end), 'Color', 'red', '-+', 'linewidth', 3 );
%plot( EA2R(2:end-1)-min(EA2R), XM2R(2:end-1), '-+', 'linewidth', 3 );
legend('L-curve','RL-curve')
xlabel('residual (A error)')
ylabel('M norm')

%%hold on;
%%loglog(resS(2:iter+1),regS(2:iter+1),'Color','red','-*');
%loglog(resD(2:iter),regD(2:iter),'-+');
%%legend('RL-curve (natural basis)','RL-curve (diagonal basis)')
%figure
%findCorner (resD(2:iter), regD(2:iter), 3)
%%%%%
%% Final problem : compute u
dirichlet = [4,1,0;4,2,0;
             3,1,0;3,2,0;
             1,1,0;1,2,0];
neumann   = [];
[K,C,nbloq] = Krig2 (nodes,elements,mat,order,boundary,dirichlet);
if inhomog >= 1
   K(1:2*nnodes, 1:2*nnodes) = Kinter;
end
f2    = [keepField( Itere, 2, boundary ); zeros(nbloq,1)];
fdir4 = dirichletRhs(urefb, 4, C, boundary);
usoli = K \ (fdir4 + f2);

usol = usoli(1:2*nnodes,1);
fsol = Kinter*usol;

% With the reduced solution
dirichlet = [4,1,0;4,2,0;
             3,1,0;3,2,0;
             1,1,0;1,2,0];
neumann   = [];
[K,C,nbloq] = Krig2 (nodes,elements,mat,order,boundary,dirichlet);
if inhomog >= 1
   K(1:2*nnodes, 1:2*nnodes) = Kinter;
end
f2    = [keepField( ItereR, 2, boundary ); zeros(nbloq,1)];
fdir4 = dirichletRhs(urefb, 4, C, boundary);
usoli = K \ (fdir4 + f2);

usolR = usoli(1:2*nnodes,1);
fsolR = Kinter*usol;

figure
hold on;
%plot(usol(2*b2node2-1),'Color','red')
plot(usolR(2*b2node2-1),'Color','blue')
plot(uref(2*b2node2-1),'Color','green')
legend('filtred solution','reference')

total_error = norm(usol-uref)/norm(uref);
[a,b] = min(error);
erroru = norm(usol(indexa) - uref(indexa)) / norm(uref(indexa));
erroru2 = norm(usol(indexa) - uref(indexa));
% Compute the energy error
f2 = dirichletRhs2( usol-uref, 2, c2node2, boundaryp2, nnodes );
uin2 = K2\f2;
lagr2 = uin2(2*nnodes+1:end,1);
lamb2 = lagr2forces2( lagr2, c2node2, 2, boundaryp2, nnodes );
%
Su = lamb2;
errorE = sqrt( ( transpose(usol-uref)*Su ) / (fref(indexa)'*uref(indexa)) );

errorx = abs((usol(2*b2node2-1)-uref(2*b2node2-1)))/max(abs(uref(2*b2node2-1)));

% Compute stress :
sigma = stress(usol,E,nu,nodes,elements,order,1,ntoelem);
plotGMSH({usol,'U_vect';sigma,'stress'}, elements, nodes, 'solution');
