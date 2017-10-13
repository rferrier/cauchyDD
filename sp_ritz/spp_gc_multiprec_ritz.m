% 19/09/2016
% Algo Steklov-Poincaré primal avec Gradient Conjugué multipréconditionné : 
% critère d'arr\^et par valeurs de Ritz.

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 70000;  % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 1;      % N.mm-1 : Loading on the plate
niter   = 30;
br      = 0.;     % noise
epsilon = 1e-100;   % Convergence criterion for ritz value
ratio   = 1e-500;   % Max ratio between eigenvalues
ntrunc  = 0;
inhomog = 2;      % inhomogeneous medium

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
% The right hand side :
f = loading(nbloq,nodes,boundary,neumann);

% Solve the problem :
uin = K\f;

% Extract displacement and Lagrange multiplicators :
uref = uin(1:2*nnodes,1);
lagr = uin(2*nnodes+1:end,1);

urefb = ( 1 + br*noise ) .* uref;

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

%% Anti-cancellation trick (don't do that for the precond ! )
% K1(indexa,indexa) = 0;
% K2(indexa,indexa) = 0;
% K1d(indexa,indexa) = 0;
% K2d(indexa,indexa) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Conjugate Gradient for the problem : (S10-S20) x = S2-S1
Itere = zeros( 2*nnodes, 1 );
d1    = zeros( 2*nnodes, niter+1 );
d2    = zeros( 2*nnodes, niter+1 );
Ad1   = zeros( 2*nnodes, niter+1 );
Ad2   = zeros( 2*nnodes, niter+1 );
d1o   = zeros( 2*nnodes, niter+1 );
d2o   = zeros( 2*nnodes, niter+1 );
Ad1o  = zeros( 2*nnodes, niter+1 );
Ad2o  = zeros( 2*nnodes, niter+1 );
Res   = zeros( 2*nnodes, niter+1 );
Zed1  = zeros( 2*nnodes, niter+1 );
Zed2  = zeros( 2*nnodes, niter+1 );
alpha = zeros( 2, niter+1 );
alpha2= zeros( 2, niter+1 );
beta  = cell(niter+1,1);

ritzval  = 0;
oldtheta = 0;
getmeout = 0;

%% Perform A x0 :
% Solve 1
f1 = dirichletRhs2( Itere, 2, c2node1, boundaryp1, nnodes );
uin1 = K1\f1;
lagr1 = uin1(2*nnodes+1:end,1);
lamb1 = lagr2forces2( lagr1, c2node1, 2, boundaryp1, nnodes );
% Solve 2
f2 = dirichletRhs2( Itere, 2, c2node2, boundaryp2, nnodes );
uin2 = K2\f2;
lagr2 = uin2(2*nnodes+1:end,1);
lamb2 = lagr2forces2( lagr2, c2node2, 2, boundaryp2, nnodes );
%
Axz = lamb1-lamb2;
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
b = -lamb1+lamb2;

%%
Res(:,1) = b - Axz;

%% Precond
% Solve 1
f1 = [Res(:,1); zeros(nbloq1d,1)];
uin1 = K1d\f1;
u1i = uin1(1:2*nnodes,1);
u1 = keepField( u1i, 2, boundaryp1 );
% Solve 2
f2 = [Res(:,1); zeros(nbloq2d,1)];
uin2 = K2d\f2;
u2i = uin2(1:2*nnodes,1);
u2 = keepField( u2i, 2, boundaryp2 );
%
Zed1(:,1) = u1;
Zed2(:,1) = u2;

d1(:,1) = Zed1(:,1);
d2(:,1) = Zed2(:,1);

residual(1) = norm(Res( indexa,1));
error(1)    = norm(Itere(indexa) - uref(indexa)) / norm(uref(indexa));
regulari(1) = sqrt( Itere'*regul(Itere, nodes, boundary, 2) );

%%
for iter = 1:niter
    
    %% Ad1
    % Solve 1
    rhs1 = d1(:,iter);
    f1 = dirichletRhs2(rhs1, 2, c2node1, boundaryp1, nnodes);
    uin1 = K1\f1;
    lagr1 = uin1(2*nnodes+1:end,1);
    lamb1 = lagr2forces2( lagr1, c2node1, 2, boundaryp1, nnodes );
    % Solve 2
    rhs2 = d1(:,iter);
    f2 = dirichletRhs2(rhs2, 2, c2node2, boundaryp2, nnodes);
    uin2 = K2\f2;
    lagr2 = uin2(2*nnodes+1:end,1);
    lamb2 = lagr2forces2( lagr2, c2node2, 2, boundaryp2, nnodes );
    %
    Ad1(:,iter) = lamb1-lamb2;
    
    %% Ad2
    % Solve 1
    rhs1 = d2(:,iter);
    f1 = dirichletRhs2(rhs1, 2, c2node1, boundaryp1, nnodes);
    uin1 = K1\f1;
    lagr1 = uin1(2*nnodes+1:end,1);
    lamb1 = lagr2forces2( lagr1, c2node1, 2, boundaryp1, nnodes );
    % Solve 2
    rhs2 = d2(:,iter);
    f2 = dirichletRhs2(rhs2, 2, c2node2, boundaryp2, nnodes);
    uin2 = K2\f2;
    lagr2 = uin2(2*nnodes+1:end,1);
    lamb2 = lagr2forces2( lagr2, c2node2, 2, boundaryp2, nnodes );
    %
    Ad2(:,iter) = lamb1-lamb2;

    %%
    Delta           = ([Ad1(indexa,iter),Ad2(indexa,iter)]'*...
                                             [d1(indexa,iter),d2(indexa,iter)]);
%    sqDel           = Delta^(1/2);
%    dtot            = sqDel\[d1(:,iter) , d2(:,iter)]';
%    Adtot           = sqDel\[Ad1(:,iter) , Ad2(:,iter)]';
%    d1o(:,iter) = d1(:,iter);  d2o(:,iter) = d2(:,iter);  % Store old ones
%    Ad1o(:,iter) = Ad1o(:,iter);  Ad2o(:,iter) = Ad2o(:,iter);
%    d1(:,iter) = dtot(1,:)';  d2(:,iter) = dtot(2,:)';
%    Ad1(:,iter) = Adtot(1,:)';  Ad2(:,iter) = Adtot(2,:)';

%    gamma           = ([Zed1(indexa,iter),Zed2(indexa,iter)]'*Res(indexa,iter));
    gamma           = ([d1(indexa,iter),d2(indexa,iter)]'*Res(indexa,iter));
    alpha(:,iter)   = Delta\gamma;

    Itere         = Itere + [d1(:,iter),d2(:,iter)]*alpha(:,iter);
    Res(:,iter+1) = Res(:,iter) - [Ad1(:,iter),Ad2(:,iter)]*alpha(:,iter);
%    Itere         = Itere + [d1(:,iter),d2(:,iter)]*gamma;
%    Res(:,iter+1) = Res(:,iter) - [Ad1(:,iter),Ad2(:,iter)]*gamma;
    
    residual(iter+1) = norm(Res(indexa,iter+1));
    error(iter+1)    = norm(Itere(indexa) - uref(indexa)) / norm(uref(indexa));
    regulari(iter+1) = sqrt( Itere'*regul(Itere, nodes, boundary, 2) );
    
    %% Precond
    % Solve 1
    f1 = [Res(:,iter+1); zeros(nbloq1d,1)];
    uin1 = K1d\f1;
    u1i = uin1(1:2*nnodes,1);
    u1 = keepField( u1i, 2, boundaryp1 );
    % Solve 2
    f2 = [Res(:,iter+1); zeros(nbloq2d,1)];
    uin2 = K2d\f2;
    u2i = uin2(1:2*nnodes,1);
    u2 = keepField( u2i, 2, boundaryp2 );
    %
    Zed1(:,iter+1) = u1;
    Zed2(:,iter+1) = u2;
    
    % Needed values for Ritz stuff
%    alpha(:,iter)   = sqDel\gamma;% sqDel\gamma
%    beta(iter)      = sqDel \ ( [Ad1(indexa,iter),Ad2(indexa,iter)]'* ...
%                                  [Zed1(indexa,iter+1),Zed2(indexa,iter+1)] );
    
    % First Reorthogonalize the residual
    for jter=1:iter-1
        betac = Zed1(indexa,iter+1)'*Res(indexa,jter) /...
                   (Zed1(indexa,jter)'*Res(indexa,jter));
        Zed1(:,iter+1) = Zed1(:,iter+1) - Zed1(:,jter) * betac;
        betac = Zed2(indexa,iter+1)'*Res(indexa,jter) /...
                   (Zed2(indexa,jter)'*Res(indexa,jter));
        Zed2(:,iter+1) = Zed2(:,iter+1) - Zed2(:,jter) * betac;
    end
    
    %% Orthogonalization
    d1(:,iter+1) = Zed1(:,iter+1);
    d2(:,iter+1) = Zed2(:,iter+1);
    
    for jter=iter:iter
%        phiij      = ( [Ad1(indexa,jter),Ad2(indexa,jter)]'*[d1(indexa,iter+1),d2(indexa,iter+1)] );
        phiij      = ( [Ad1(indexa,jter),Ad2(indexa,jter)]'*[Zed1(indexa,iter+1),Zed2(indexa,iter+1)] );
       % matpr = inv([Ad1(indexa,jter),Ad2(indexa,jter)]'*[d1(indexa,jter),d2(indexa,jter)])
        betaij     = ([Ad1(indexa,jter),Ad2(indexa,jter)]'*[d1(indexa,jter),d2(indexa,jter)]) \ phiij;
        beta(iter) = betaij;
        Prov = [d1(:,iter+1),d2(:,iter+1)] - [d1(:,jter),d2(:,jter)] * betaij;%phiij;%* betaij;%phiij;
        
        d1(:,iter+1) = Prov(:,1);
        d2(:,iter+1) = Prov(:,2);
    end
    
    %% Build V
%    V(:,2*iter-1)  = zeros(4*nnodes,1);
%    V(:,2*iter)    = zeros(4*nnodes,1);
%    V(2*indexa-1,2*iter-1) = (-1)^(iter-1)*Zed1(indexa,iter)/...
%                             (sqrt(Res(indexa,iter)'*Zed1(indexa,iter)));
%    V(2*indexa,2*iter)     = (-1)^(iter-1)*Zed2(indexa,iter)/...
%                             (sqrt(Res(indexa,iter)'*Zed2(indexa,iter)));
    V(:,2*iter-1)  = zeros(2*nnodes,1);
    V(:,2*iter)    = zeros(2*nnodes,1);
%    V(indexa,2*iter-1) = (-1)^(iter-1)*Zed1(indexa,iter);%/...
%                             (sqrt(Res(indexa,iter)'*Zed1(indexa,iter)));
%    V(indexa,2*iter)   = (-1)^(iter-1)*Zed2(indexa,iter);%/...
%                             (sqrt(Res(indexa,iter)'*Zed2(indexa,iter)));
    V(indexa,2*iter-1) = Zed1(indexa,iter)/(sqrt(Res(indexa,iter)'*Zed1(indexa,iter)));
    V(indexa,2*iter)   = Zed2(indexa,iter)/(sqrt(Res(indexa,iter)'*Zed2(indexa,iter)));
%    V(indexa,2*iter-1) = Zed1(indexa,iter);
%    V(indexa,2*iter)   = Zed2(indexa,iter);
%                             
    if iter > 1
       betapa = cell2mat(beta(iter-1));
       zt1p = zt1; zt2p = zt2; Azt1p = Azt1; Azt2p = Azt2; 
       zt1  = d1(:,iter) + betapa(1,1)*d1(:,iter-1) + betapa(2,1)*d2(:,iter-1);
       Azt1 = Ad1(:,iter) + betapa(1,1)*Ad1(:,iter-1) + betapa(2,1)*Ad2(:,iter-1);
       zt2  = d2(:,iter) + betapa(1,2)*d1(:,iter-1) + betapa(2,2)*d2(:,iter-1);
       Azt2 = Ad2(:,iter) + betapa(1,2)*Ad1(:,iter-1) + betapa(2,2)*Ad2(:,iter-1);
       
       H(2*iter-1,2*iter-1) = d1(:,iter)'*Ad1(:,iter)...
                              + betapa(1,1)^2 * d1(:,iter-1)'*Ad1(:,iter-1)...
                              + betapa(1,2)^2 * d2(:,iter-1)'*Ad2(:,iter-1)...
                              + 2*betapa(1,1)*betapa(2,1) * d1(:,iter-1)'*Ad2(:,iter-1);

       H(2*iter,2*iter)     = d2(:,iter)'*Ad2(:,iter)...
                              + betapa(2,1)^2 * d1(:,iter-1)'*Ad1(:,iter-1)...
                              + betapa(2,2)^2 * d2(:,iter-1)'*Ad2(:,iter-1)...
                              + 2*betapa(1,2)*betapa(2,2) * d1(:,iter-1)'*Ad2(:,iter-1);
       ind = [2*iter-1;2*iter];
       H(ind,ind) = [zt1(indexa),zt2(indexa)]'*[Azt1(indexa),Azt2(indexa)];
       H(ind-2,ind) = [zt1p(indexa),zt2p(indexa)]'*[Azt1(indexa),Azt2(indexa)];
       H(ind,ind-2) = [zt1(indexa),zt2(indexa)]'*[Azt1p(indexa),Azt2p(indexa)];

%%       % Compute denominators
       de1  = sqrt(Res(indexa,iter)'*Zed1(indexa,iter));
       de2  = sqrt(Res(indexa,iter)'*Zed2(indexa,iter));
       de1p = sqrt(Res(indexa,iter-1)'*Zed1(indexa,iter-1));
       de2p = sqrt(Res(indexa,iter-1)'*Zed2(indexa,iter-1));
       
       % And apply it
       H(ind,ind) = H(ind,ind) ./ ([de1,de2]'*[de1,de2]);
       H(ind-2,ind) = H(ind-2,ind) ./ ([de1p,de2p]'*[de1,de2]);
       H(ind,ind-2) = H(ind,ind-2) ./ ([de1,de2]'*[de1p,de2p]);
       
    else
       zt1  = d1(:,iter);
       Azt1 = Ad1(:,iter);
       zt2  = d2(:,iter);
       Azt2 = Ad2(:,iter);
       ind = [2*iter-1;2*iter];
       H(ind,ind) = [d1(:,iter),d2(:,iter)]'*[Ad1(:,iter),Ad2(:,iter)];
       de1  = sqrt(Res(indexa,iter)'*Zed1(indexa,iter));
       de2  = sqrt(Res(indexa,iter)'*Zed2(indexa,iter));
       H(ind,ind) = H(ind,ind) ./ ([de1,de2]'*[de1,de2]);
    end
    
    % Compute eigenelems of the Hessenberg :
    [Q,Theta1] = eig(H);
    theta = diag(Theta1);
    % Sort it
    [theta,Ind] = sort(theta,'descend');
    Q = Q(:,Ind);
    Theta1 = Theta1(Ind,Ind);
    Y = V*Q;
    
    % See if the current two converged
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
    
    % Do it again ( I'm too lazy to debug my while )
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
 end
 
 % Compute the solution
chi = inv(Theta1)*Y'*b;
if ntrunc > 0
   chi(ntrunc:end) = 0;
end
ItereR = Y*chi;
 
%    %% Build H from V (DEBUG)
%for iter = 1:size(V,2)
%    %% Avi
%    % Solve 1
%    rhs1 = V(:,iter);
%    f1 = dirichletRhs(rhs1, 2, C1, boundaryp1);
%    uin1 = K1\f1;
%    lagr1 = uin1(2*nnodes+1:end,1);
%    lamb1 = lagr2forces( lagr1, C1, 2, boundaryp1 );
%    % Solve 2
%    rhs2 = V(:,iter);
%    f2 = dirichletRhs(rhs2, 2, C2, boundaryp2);
%    uin2 = K2\f2;
%    lagr2 = uin2(2*nnodes+1:end,1);
%    lamb2 = lagr2forces( lagr2, C2, 2, boundaryp2 );
%    %
%    AV(:,iter) = lamb1-lamb2;
%end
%
%He = V'*AV;

figure;
hold on;
plot(log10(abs(theta)),'Color','blue')  % Rem : I need abs on this one...
plot(log10(abs(Y'*b)),'Color','red')
plot(log10(abs(chi)),'Color','black')
legend('Ritz Values','RHS values','solution coefficients')

figure;
hold on;
plot(Itere(2*b2node2-1),'Color','red')
plot(ItereR(2*b2node2-1),'Color','blue')
plot(uref(2*b2node2-1),'Color','green')
legend('brutal solution','filtred solution', 'reference')

figure;
hold on
plot(log10(error),'Color','blue')
plot(log10(residual),'Color','red')
legend('error (log)','residual (log)')

% L-curve :
figure
loglog(residual,regulari,'-+');

%%%%
%% Final problem : compute u
% DN problem
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

%plot(fsol(2*b2node2-1))

total_error = norm(usol-uref)/norm(uref);
% Compute stress :
sigma = stress(usol,E,nu,nodes,elements,order,1,ntoelem);
plotGMSH({usol,'U_vect';sigma,'stress'}, elements, nodes, 'solution');