% 22/09/2017
% Algo Bayésien s'appuyant sur l'analyse de Ritz, u VS u, bloc primal

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E        = 70000;    % MPa : Young modulus
nu       = 0.3;      % Poisson ratio
fscalar  = 1;        % N.mm-1 : Loading on the plate
niter    = 10;
precond  = 1;        % 1 : Use a dual precond, 2 : use regul precond
mu       = 0.;       % Regularization parameter
br       = 0.1;      % noise
ntrunc   = 2;        % Ritz truncature level for the prior
nbayesp1 = 20;        % Ritz truncature level for Bayesian inversion
nMC       = 100000; % nb of MC samples
refined   = .5;      % Choose the mesh to use
chaosor   = 1;      % Order of the polynomial chaos max = 10
cornermes = 0;      % Take the encastred corners as measurement (or not)
reduction = 2;      % 0 : no reduction, 1 : replace G by its SVD, 2 : reduced stochastic parameters

szpri     = 0;   % Size of the prior distribution, 0 : automatic estimate

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
%neumann   = [2,1,fscalar;
%             4,1,fscalar];

% Import the mesh
if refined == 5
   [ nodes,elements,ntoelem,boundary,order ] = readmesh( 'meshes/plateUs/plate5.msh' );
    nnodes = size(nodes,1); noise = randn(2*nnodes,1);
elseif refined == 4
   [ nodes,elements,ntoelem,boundary,order ] = readmesh( 'meshes/plateUs/plate4.msh' );
    nnodes = size(nodes,1); noise = randn(2*nnodes,1);
elseif refined == 3
   [ nodes,elements,ntoelem,boundary,order ] = readmesh( 'meshes/plateUs/plate3.msh' );
    nnodes = size(nodes,1); noise = randn(2*nnodes,1);
elseif refined == 2
   [ nodes,elements,ntoelem,boundary,order ] = readmesh( 'meshes/plateUs/plate2.msh' );
    noises = load('./noises/noise0U.mat'); % Particular noise vector
    nnodes = size(nodes,1); noise  = noises.bruit1;
elseif refined == .5
   [ nodes,elements,ntoelem,boundary,order ] = readmesh( 'meshes/plateUs/plate05.msh' );
   noises = load('./noises/noise0.mat'); % Particular noise vector
   nnodes = size(nodes,1); noise  = noises.bruit1;
%   noise = randn(2*nnodes,1);
end

%noises = load('./noises/noise1.mat'); % Particular noise vector
%noise  = noises.bruit1;
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
[~, b2node4] = mapBound(4, boundaryp1, nnodes);
indexa = [2*b2node2-1; 2*b2node2]; nindexa = size(indexa,1);
indexb = [2*b2node4-1; 2*b2node4]; nindexb = size(indexb,1);

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
uref = uin(1:2*nnodes,1); fref  = Kinter*uref;
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
f1 = dirichletRhs2( Itere(:,1), 2, c2node1, boundaryp1, nnodes );
uin1 = K1\f1;
lagr1 = uin1(2*nnodes+1:end,1);
lamb1 = lagr2forces2( lagr1, c2node1, 2, boundaryp1, nnodes );
% Solve 2
f2 = dirichletRhs2( Itere(:,2), 2, c2node2, boundaryp2, nnodes );
uin2 = K2\f2;
lagr2 = uin2(2*nnodes+1:end,1);
lamb2 = lagr2forces2( lagr2, c2node2, 2, boundaryp2, nnodes );
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
V  = zeros(2*nnodes, 2*niter);
AV = zeros(2*nnodes, 2*niter);
MV = zeros(2*nnodes, 2*niter);
H  = zeros(2*niter);
%eta = [0,0;0,0];
num = [0,0;0,0]; % useless, but eta needs initialization #lazy
den = [0,0;0,0];
%%
for iter = 1:niter
    %% Optimal step
    
    % Solve 1
    f11 = dirichletRhs2( d(:,2*iter-1), 2, c2node1, boundaryp1, nnodes );
    f12 = dirichletRhs2( d(:,2*iter), 2, c2node1, boundaryp1, nnodes );
    uin1 = K1\[f11,f12];
    lagr1 = uin1(2*nnodes+1:end,[1,2]);
    lamb11 = lagr2forces2( lagr1(:,1), c2node1, 2, boundaryp1, nnodes );
    lamb12 = lagr2forces2( lagr1(:,2), c2node1, 2, boundaryp1, nnodes );
    % Solve 2
    f21 = dirichletRhs2( d(:,2*iter-1), 2, c2node2, boundaryp2, nnodes );
    f22 = dirichletRhs2( d(:,2*iter), 2, c2node2, boundaryp2, nnodes );
    uin2 = K2\[f21,f22];
    lagr2 = uin2(2*nnodes+1:end,[1,2]);
    lamb21 = lagr2forces2( lagr2(:,1), c2node2, 2, boundaryp2, nnodes );
    lamb22 = lagr2forces2( lagr2(:,2), c2node2, 2, boundaryp2, nnodes );
    %
    Ad(:,[2*iter-1,2*iter]) = [lamb11-lamb21,lamb12-lamb22];

    denprec = den; numprec = num; % Store those ones
    den = d(indexa,[2*iter-1,2*iter])'*Ad(indexa,[2*iter-1,2*iter]);
    sqD = den^(1/2);
    d(:,[2*iter-1,2*iter]) = d(:,[2*iter-1,2*iter]) * inv(sqD);
    Ad(:,[2*iter-1,2*iter]) = Ad(:,[2*iter-1,2*iter]) * inv(sqD);
    num = Res(indexa,[2*iter-1,2*iter])'*Zed(indexa,[2*iter-1,2*iter]);
    num = sqD\num; % because of Zed and not d
    
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
    
    % Compute eigenelems of the Hessenberg :
    [Q,Theta1] = eig(H);
    theta = diag(Theta1);
    % Sort it
    [theta,Ind] = sort(theta,'descend');
    Q = Q(:,Ind);
    Theta1 = Theta1(Ind,Ind);
    Y = V*Q;
    
end

residual = residual/residual(1); % Normalize the residual

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

%figure;
%hold on
%plot(error,'Color','blue')
%plot(residual,'Color','red')
%legend('error','residual')
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
%plot(Itere(2*b2node2-1,1)-Itere(2*b2node2-1,2),'Color','red');
plot(ItereR(2*b2node2-1),'Color','blue');
plot(uref(2*b2node2-1),'Color','green');
legend('solution','reference')
%legend('solution','Ritz solution','reference')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bayesian stuff
if szpri == 0
   szpri = max(ItereR(indexa))-min(ItereR(indexa));
end

% Approximate the rectangular prior with chaos expansion (1D)
Hermite = [ 1,0,0,0,0,0,0,0,0,0,0
            0,1,0,0,0,0,0,0,0,0,0
            -1,0,1,0,0,0,0,0,0,0,0
            0,-3,0,1,0,0,0,0,0,0,0
            3,0,-6,0,1,0,0,0,0,0,0
            0,15,0,-10,0,1,0,0,0,0,0
            -15,0,45,0,-15,0,1,0,0,0,0
            0,-105,0,105,0,-21,0,1,0,0,0
            105,0,-420,0,210,0,-28,0,1,0,0        % Source : Modeling multibody systems with uncertainties. Part I:Theoretica and computational aspects
            0,945,0,-1260,0,378,0,-36,0,1,0       %
            -945,0,4725,0,-3150,0,630,0,-45,0,1]; % Polynomial coefficients

theta  = rand(1,nMC);                        % random seed
prior  = szpri*theta-szpri/2;                        % rectangular law
thetaG = sqrt(2)*erfinv(2*theta-1);          % gaussian law N(0,1)

% Determinate the chaos coefficients (via Monte-Carlo)
coef = zeros(chaosor+1,1);
for i=0:chaosor
   esti = 0;
   for j=0:chaosor
      esti = esti + Hermite(i+1,j+1)*thetaG.^j;
   end
   coef(i+1) = mean(prior.*esti)/factorial(i);
end

%%%%%%%%%%%%%%%%%%%%%%
% Forward model
dirichletm = [1,1,0; 1,2,0 ;
              3,1,0; 3,2,0 ;
              2,1,0; 2,2,0];
neumannm   = [4,1,fscalar,0,-fscalar];
% (TODO : the inhomogeneous trick)
[Km,Cm,nbloqm,node2cm,c2nodem] = Krig2 (nodes,elements,mat,order,boundary,dirichletm);
%%%%%%%%%%%%%%%%%%%%%%

nbayes = nbayesp1 - 1;

% Compute the approximated TGSVD of the model
if reduction ~= 0 % G \simeq Ut*Sigmat*Vt'; Ut'*Zd*Ut = 1; Vt'*M*Vt = 1;
   Yt = Y(indexa,1:nbayes);
   Thetat = Theta1(1:nbayes,1:nbayes);
   
   Sigmat = sqrt(Thetat); %Vt = Yt;
   
   % Vt = MYt
   if precond == 1
      rhs = zeros(2*nnodes,nbayes); rhs(indexa,:) = Yt;
      f11 = dirichletRhs2( rhs, 2, c2node1, boundaryp1, nnodes );
      uin1 = K1\f11;
      lagr1 = uin1(2*nnodes+1:end,:);
      rhs = zeros(2*nnodes,nbayes);
      for i=1:nbayes
         rhs(:,i) = lagr2forces2( lagr1(:,i), c2node1, 2, boundaryp1, nnodes );
      end
   else
      rhs = zeros(2*nnodes,nbayes); rhs(indexa,:) = Yt;
   end
   Vt = rhs(indexa,:);
   
   % Ut = G*M^-1*Vt*Sigmat^-1
   if precond == 1
%      rhs = zeros(2*nnodes,size(Vt*Sigmat,2)); rhs(indexa,:) = Vt*inv(Sigmat);
%      f11 = dirichletRhs2( rhs, 2, c2node1, boundaryp1, nnodes );
%      uin1 = K1\f11;
%      lagr1 = uin1(2*nnodes+1:end,:);
%      rhs = zeros(2*nnodes,size(Vt*Sigmat,2));
%      for i=1:size(Vt*Sigmat,2)
%         rhs(:,i) = lagr2forces2( lagr1(:,i), c2node1, 2, boundaryp1, nnodes );
%      end
      rhs = zeros(2*nnodes,size(Vt*Sigmat,2)); rhs(indexa,:) = Vt*inv(Sigmat);
      rhs = K1d\[rhs;zeros(nbloq1d,size(Vt*Sigmat,2))];
   else
      rhs = zeros(2*nnodes,size(Vt*Sigmat,2)); rhs(indexa,:) = Vt*inv(Sigmat);
   end
   %Mdep = zeros(2*nnodes,size(Vt*Sigmat,2)); Mdep(indexa,:) = rhs;
   Mrhs = dirichletRhs2( rhs, 2, c2nodem, boundary, nnodes );
   Ut = Km\Mrhs; Ut = Ut(indexb,:);
   
%   % DEBUG : test the orthogonality of Ut
%   rhs = zeros(2*nnodes,size(Ut,2)); rhs(indexb,:) = Ut;
%   f11 = dirichletRhs2( rhs, 4, c2node1, boundaryp1, nnodes );
%   uin1 = K1\f11;
%   lagr1 = uin1(2*nnodes+1:end,:);
%   rhs = zeros(2*nnodes,size(Ut,2));
%   for i=1:size(Ut,2)
%      rhs(:,i) = lagr2forces2( lagr1(:,i), c2node1, 4, boundaryp1, nnodes );
%   end
%   rhs = rhs(indexb,:);
%   Ut'*rhs

   Gtilde = Ut*Sigmat*Vt';
end

if reduction == 2 % Linear Bayesian update (Kalman filter) 
   % Propagate the uncertainties into the reduced basis : multiply by Vt'
   prior = zeros(nindexa,1+(chaosor)*nindexa);
   for i=1:nindexa % See how optimized it is
      prior(i,(chaosor)*(i-1)+2:(chaosor)*i+1) = transpose(coef(2:end));  % Coefficients of the PCE for the prior
   end
   prior(:,1) = ItereR(indexa);
   pritilde = Vt'*prior; % Reduced prior
   
   % Reduced measurments
   f1 = dirichletRhs2( urefb, 4, c2node1, boundary, nnodes );
   uin1 = K1\f1; lagr1 = uin1(2*nnodes+1:end,1);
   Zu = lagr2forces2( lagr1, c2node1, 4, boundaryp1, nnodes );
   btilde = Ut'*Zu(indexb,:);
%   btilde = (Ut'*Ut)\Ut'*urefb;
   
   % Estimate additive noise level and project it on the basis
   nl = br*mean(abs(uref(indexb))); % indexb
   Ceps0 = nl*eye(nindexb);     % Correlation matrix of the noise
   
   Ceps1 = zeros(2*nnodes,size(Ceps0,2));
   % Multiply in order to get the noise on b from the noise on u
   for i = 1:size(Ceps0,1)
      rhs = zeros(2*nnodes); rhs(indexb) = Ceps0(:,i);
      f1 = dirichletRhs2( rhs, 4, c2node1, boundary, nnodes );
      uin1 = K1\f1;
      lagr1 = uin1(2*nnodes+1:end,1);
      Ceps1(:,i) = lagr2forces2( lagr1, c2node1, 4, boundaryp1, nnodes );
   end
   Ceps1 = Ceps1(indexb,:); Ceps = Ceps1*Ceps1';
   
   Ctilde = Ut'*Ceps*Ut;         % Reduced correlation matrix
   ntilde = Ut'*Ceps1;           % Reduced noise 

   % Solve the neumann part of the problem
   fm = loading(nbloqm,nodes,boundary,neumannm); fm = fm(indexb);
   ttilde = Ut'*fm;
%   Msol1 = Km\fm; Msol1 = Msol1(indexb,:);  % TOCHECK : I guess it is Zd-1
   
   Msol = Sigmat*pritilde;
   Msol(:,1) = Msol(:,1) + ttilde; % Generate the synthetic measurement PDF

   % Compute observation covariance
   Cy = zeros(nbayes);
   for i=1:chaosor % We're intentionnaly starting at 1 'cause uniform law has no covariance
      for j=1:nindexa
         Cy = Cy + factorial(i)*Msol(:,i+1+(j-1)*(chaosor))*Msol(:,i+1+(j-1)*(chaosor))';
      end
   end
   Cd = Ctilde+Cy; % Total covariance (from noise and prediction)
   
   % Estimate the Kalman factor
   ZmY = -Msol;        % Difference between measurement and prediction
   ZmY(:,1) = ZmY(:,1) + btilde;  % The mean
   for i=1:nindexa  % and probability measurement (again in PCE form)
      ZmY(:,(chaosor)*(i-1)+2) = ZmY(:,(chaosor)*(i-1)+2) + ...
              ntilde(:,i)/sqrt(nindexa);  % /!\ Rem : I guess there should be other RVs here
   end
   Geai = Cd\ZmY;
   
   % Covariance between a-priori and reconstructed observation
   Cqy = zeros(nbayes);
   for i=1:chaosor
      for j=1:nindexa
         Cqy = Cqy + factorial(i)*pritilde(:,i+1+(j-1)*(chaosor))*Msol(:,i+1+(j-1)*(chaosor))';
      end
   end

   % Build the PCE of the solution
   Qb = pritilde + Cqy*Geai;
   
   sigmav = zeros(nbayes,1);
   CD     = zeros(nbayes); % Posterior variance
   for j=1:nindexa 
      for i=1:chaosor
         CD = CD + factorial(i)*Qb(:,i+1+(j-1)*(chaosor))*Qb(:,i+1+(j-1)*(chaosor))';
      end
   end
   sigmav = sqrt(diag(CD));
   meanD = Qb(:,1);
   
   %% Back to the regular basis
   Qbu = Vt*((Vt'*Vt)\Qb); meanU = Qbu(:,1);
   CU    = zeros(nindexa); % Posterior variance
   for j=1:nindexa 
      for i=1:chaosor
         CU = CU + factorial(i)*Qbu(:,i+1+(j-1)*(chaosor))*Qbu(:,i+1+(j-1)*(chaosor))';
      end
   end
   sigmau = sqrt(diag(CU));
   
%   meanU = Yt*meanD; % Back to the regular basis
%   meanU = Vt*((Vt'*Vt)\meanD); % Back to the regular basis
%   sCD   = CD^(1/2); sCU = Vt*((Vt'*Vt)\sCD); CU = sCU*sCU';
%   CU    = Yt*CD*Yt';

%%   if precond == 1
%%      rhs = zeros(2*nnodes,1); rhs(indexa,:) = YmeanD;
%%      f11 = dirichletRhs2( rhs, 2, c2node1, boundaryp1, nnodes );
%%      uin1 = K1\f11; lagr1 = uin1(2*nnodes+1:end,:);
%%      meanU = lagr2forces2( lagr1, c2node1, 2, boundaryp1, nnodes );
%%      meanU = meanU(indexa,1);
%%      
%%      rhs = zeros(2*nnodes,size(sCD,2)); rhs(indexa,:) = YsCD;
%%      f11 = dirichletRhs2( rhs, 2, c2node1, boundaryp1, nnodes );
%%      uin1 = K1\f11;
%%      lagr1 = uin1(2*nnodes+1:end,:);
%%      sCU = zeros(2*nnodes,size(sCD,2));
%%      for i=1:size(sCD,2)
%%         sCU(:,i) = lagr2forces2( lagr1(:,i), c2node1, 2, boundaryp1, nnodes );
%%      end
%%      sCU = sCU(indexa,:); CU = sCU*sCU';
%%   else
%%      meanU = YmeanD; CU = YsCD*YsCD';
%%   end
%   sigmau = sqrt(diag(CU));
else % basic Kalman Filter
   %% Solve the stochastic forward problem to synthetize observations
   prior = zeros(nindexa,1+(chaosor)*nindexa);
   for i=1:nindexa % See how optimized it is
      prior(i,(chaosor)*(i-1)+2:(chaosor)*i+1) = transpose(coef(2:end));  % Coefficients of the PCE for the prior
      prior(i,1) = coef(1);
   end

   if reduction == 0
      % Solve the linear (and not affine) problems for each dof 
      % (in order to compute Msol0, matrix of the direct model)
      Mdepp = eye(nindexa); Mdep = zeros(2*nnodes,nindexa); Mdep(indexa,:) = Mdepp;
      Mrhs = dirichletRhs2( Mdep, 2, c2nodem, boundary, nnodes );
      Msol0 = Km\Mrhs;  Msol0 = Msol0(indexb,:);
   end

   % Solve the neumann part of the problem
   fm = loading(nbloqm,nodes,boundary,neumannm);
   Msol1 = Km\fm; Msol1 = Msol1(indexb,:);
   
   Msol = zeros(nindexb,1+(chaosor)*nindexa);
   if reduction == 0 Msol = Msol0*prior;
   else Msol = Gtilde*prior; end
   Msol(:,1) = Msol(:,1) + Msol1;

   % Estimate additive noise level
   nl = br*mean(abs(uref(indexb))); % indexb
   Ceps = nl^2*eye(nindexb);  % Correlation matrix of the noise
   % Compute observation covariance
   Cy = zeros(nindexb);
   for i=1:chaosor % We're intentionnaly starting at 1 'cause uniform law has no covariance
      for j=1:nindexa
         Cy = Cy + factorial(i)*Msol(:,i+1+(j-1)*(chaosor))*Msol(:,i+1+(j-1)*(chaosor))';
      end
   end
   Cd = Ceps+Cy; % Total covariance (from noise and prediction)
   
   % Estimate the Kalman factor
   ZmY = -Msol;        % Difference between measurement and prediction
   ZmY(:,1) = ZmY(:,1) + urefb(indexb);  % The mean
   for i=1:nindexa  % and probability measurement (again in PCE form)
      ZmY(:,(chaosor)*(i-1)+2) = ZmY(:,(chaosor)*(i-1)+2) + nl;%/sqrt(nindexa);  % /!\ Rem : I guess there should be other RVs here
   end
   Geai = Cd\ZmY;
   
   % Covariance between a-priori and reconstructed observation
   Cqy = zeros(nindexa,nindexb);
   for i=1:chaosor
      for j=1:nindexa
         Cqy = Cqy + factorial(i)*prior(:,i+1+(j-1)*(chaosor))*Msol(:,i+1+(j-1)*(chaosor))';
      end
   end
   
   % Build the PCE of the solution
   Qb = prior + Cqy*Geai; Qbu = Qb;
   
   sigmau = zeros(nindexa,1);
   CU     = zeros(nindexa); % Posterior variance
   for j=1:nindexa 
      for i=1:chaosor
         CU = CU + factorial(i)*Qb(:,i+1+(j-1)*(chaosor))*Qb(:,i+1+(j-1)*(chaosor))';
      end
   end
   sigmau = sqrt(diag(CU));
   meanU = Qb(:,1);
end

usoltot = zeros(2*nnodes,1); usoltot(indexa) = meanU;
usolmax = zeros(2*nnodes,1); usolmax(indexa) = meanU+sigmau;
usolmin = zeros(2*nnodes,1); usolmin(indexa) = meanU-sigmau;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%
%% Final problem : compute u and f : Rem : in the reduced case, it is possible to do it on M-1V (less dimensions)
% DN problem
dirichlet = [3,1,0;3,2,0;
             1,1,0;1,2,0;
             2,1,0;2,2,0];
neumann   = [4,1,fscalar,0,-fscalar];
[K,C,nbloq,node2c,c2node] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
rhs = zeros(2*nnodes,size(Qbu,2)); rhs(indexa,:) = Qbu;
fdir2 = dirichletRhs2(rhs, 2, c2node, boundary,nnodes);
flin4 = loading(nbloq,nodes,boundary,neumann);
usol0 = K \ fdir2; usol1 = K\fdir2;

usol = usol1(1:2*nnodes,1) + usol0(1:2*nnodes,1)*ones(1,size(Qbu,2));
fsol = Kinter*usol; Qf = fsol(indexa,:);
%
sigmaf = zeros(nindexa,1);
CF     = zeros(nindexa); % Posterior variance
for j=1:nindexa 
   for i=1:chaosor
      CF = CF + factorial(i)*Qf(:,i+1+(j-1)*(chaosor))*Qf(:,i+1+(j-1)*(chaosor))';
   end
end
sigmaf = sqrt(diag(CF));
meanF = Qf(:,1);
%
fsoltot = zeros(2*nnodes,1); fsoltot(indexa) = meanF;
fsolmax = zeros(2*nnodes,1); fsolmax(indexa) = meanF+sigmaf;
fsolmin = zeros(2*nnodes,1); fsolmin(indexa) = meanF-sigmaf;
%
figure;
hold on;
plot(usoltot(2*b2node2-1),'Color','red','linewidth',3);
plot(uref(2*b2node2-1),'Color','green','linewidth',3);
plot(usolmax(2*b2node2-1),'Color','red','linewidth',3);
plot(usolmin(2*b2node2-1),'Color','red','linewidth',3);
legend('solution','reference');
%
figure;
hold on;
plot(fsoltot(2*b2node2-1),'Color','red','linewidth',3);
plot(fref(2*b2node2-1),'Color','green','linewidth',3);
plot(fsolmax(2*b2node2-1),'Color','red','linewidth',3);
plot(fsolmin(2*b2node2-1),'Color','red','linewidth',3);
legend('solution (force)','reference');
%
error0 = norm(ItereR(indexa)-uref(indexa))/norm(uref(indexa));
error1 = norm(usoltot(indexa)-uref(indexa))/norm(uref(indexa));
%total_error = norm(usol-uref)/norm(uref);
%% Compute stress :
%sigma = stress(usol,E,nu,nodes,elements,order,1,ntoelem);
%plotGMSH({usol,'U_vect';sigma,'stress'}, elements, nodes, 'solution');