%13/09/2017
%Inversion bayesienne du pb de Cauchy avec les Polynomes du chaos, procédure iterative

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E         = 70000;  % MPa : Young modulus
nu        = 0.3;    % Poisson ratio
fscalar   = 1;      % N.mm-2 : Loading on the plate
br        = .1;      % noise
brt       = 0;      % "translation" noise
inhomog   = 0;      % inhomogeneous medium
nMC       = 100000; % nb of MC samples
refined   = .5;      % Choose the mesh to use
chaosor   = 1;      % Order of the polynomial chaos max = 10
cornermes = 0;      % Take the encastred corners as measurement (or not)
niter     = 1;      % Nb of iterations

upperpr   = 1e-3;   % Bounds of the prior distribution
lowerpr   = -1e-3;  %

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
             2,2,-10*fscalar,0,0;
             4,1,fscalar,0,-fscalar];

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
    noise  = noises.bruit1;
elseif refined == .5
   [ nodes,elements,ntoelem,boundary,order ] = readmesh( 'meshes/plateUs/plate05.msh' );
   noises = load('./noises/noise0.mat'); % Particular noise vector
   noise  = noises.bruit1;
end
nnodes = size(nodes,1);

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
if cornermes == 0 % Suppress the corners from the mesurement
   boundaryp1 = suppressBound( boundaryp1, no1, 4 );
   boundaryp1 = suppressBound( boundaryp1, no4, 4 );
end
boundaryp2 = boundaryp1;

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

% Ugly hack for 1dof case :
if refined == 5 b2node2 = 5;  end
if refined == 5 && cornermes == 0 b2node4 = 6; end

indexa = [2*b2node2-1; 2*b2node2]; nindexa = size(indexa,1);
indexb = [2*b2node4-1; 2*b2node4]; nindexb = size(indexb,1);
% The right hand side :
f = loading(nbloq,nodes,boundary,neumann);

% Solve the problem :
uin = K\f;

% Extract displacement and Lagrange multiplicators :
uref = uin(1:2*nnodes,1);
lagr = uin(2*nnodes+1:end,1);

%urefb = ( 1 + brt + br*noise ) .* uref;
urefb = uref + br*noise*mean(abs(uref(indexb)));
fref  = Kinter*uref;

sigma = stress(uref,E,nu,nodes,elements,order,1,ntoelem);
plotGMSH({uref,'Vect_U';sigma,'stress'}, elements, nodes, 'output/reference');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bayesian inversion
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
prior  = (upperpr-lowerpr)*theta + lowerpr;  % rectangular law
thetaG = sqrt(2)*erfinv(2*theta-1);                  % gaussian law N(0,1)

% Determinate the chaos coefficients (via Monte-Carlo)
coef = zeros(chaosor+1,1);
for i=0:chaosor
   esti = 0;
   for j=0:chaosor
      esti = esti + Hermite(i+1,j+1)*thetaG.^j;
   end
   coef(i+1) = mean(prior.*esti)/factorial(i);
end

% Note : there is no crossed coefficient and there will never be
% as the problem is linear !

% Determine the coeffs in the standard polynomial basis
coefS = Hermite(1:chaosor+1,1:chaosor+1)'*coef;

%% Plot the approximated theta
%thetae  = randn(1,nMC);
%toplo = 0;
%for j=0:chaosor
%   toplo = toplo + coefS(j+1)*thetae.^j;
%end
%figure;
%hist(toplo,nMC/1000);

%% Solve the stochastic forward problem to synthetize observations
% Build the LHS
dirichletm = [1,1,0; 1,2,0 ;
              3,1,0; 3,2,0 ;
              2,1,0; 2,2,0];
neumannm   = [4,1,fscalar,0,-fscalar];
% TODO : the inhomogeneous trick !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
[Km,Cm,nbloqm,node2cm,c2nodem] = Krig2 (nodes,elements,mat,order,boundary,dirichletm);

prior = zeros(nindexa,1+(chaosor)*nindexa);
for i=1:nindexa % See how optimized it is
   prior(i,(chaosor)*(i-1)+2:(chaosor)*i+1) = transpose(coef(2:end));  % Coefficients of the PCE for the prior
   prior(i,1) = coef(1);
end

% Solve the linear (and not affine) problems for each dof (in order to compute Msol0, matrix of the direct model)
Mdepp = eye(nindexa); Mdep = zeros(2*nnodes,nindexa); Mdep(indexa,:) = Mdepp;
Mrhs = dirichletRhs2( Mdep, 2, c2nodem, boundary, nnodes );
Msol0 = Km\Mrhs; Msol0 = Msol0(indexb,:);

% Solve the neumann part of the problem
fm = loading(nbloqm,nodes,boundary,neumannm);
Msol1 = Km\fm; Msol1 = Msol1(indexb,:);

for iter = 1:niter

   % New measurement
   noise = randn(2*nnodes,1);
   urefb = uref + br*noise*mean(abs(uref(indexb)));
   %
   Msol = zeros(nindexb,1+(chaosor)*nindexa);
   Msol = Msol0 * prior;
   Msol(:,1) = Msol(:,1) + Msol1;% + Msol0*prior(:,1); % The mean
%   for i=1:nindexa % Msol is the PCE of the measurements
%      Msol(:,(chaosor)*(i-1)+2:(chaosor)*i+1) = Msol0 *...
%                           prior(:,(chaosor)*(i-1)+2:(chaosor)*i+1);
%%      Msol(:,(chaosor)*(i-1)+2:(chaosor)*i+1) = Msol0(:,i) *...
%%                        prior(i,(chaosor)*(i-1)+2:(chaosor)*i+1);
%   end

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
   %prior(1,:)
   %Msol(1,:)
   norm(diag(Cy));
%   norm(diag(Ceps))
   Cd = Ceps+Cy; % Total covariance (from noise and prediction)
   %norm(diag(Cy))
   % Estimate the Kalman factor
   ZmY = -Msol;        % Difference between measurement and prediction
   ZmY(:,1) = ZmY(:,1) + urefb(indexb);  % The mean
   for i=1:nindexa  % and probability measurement (again in PCE form)
      ZmY(:,(chaosor)*(i-1)+2) = ZmY(:,(chaosor)*(i-1)+2) + nl;%/sqrt(nindexa);  % /!\ Rem : I guess there should be other RVs introduced here
   end
   Geai = Cd\ZmY;
   
   % Covariance between a-priori and reconstructed observation
   Cqy = zeros(nindexa,nindexb);
   for i=1:chaosor
      for j=1:nindexa
         Cqy = Cqy + factorial(i)*prior(:,i+1+(j-1)*(chaosor))*Msol(:,i+1+(j-1)*(chaosor))';
      end
   end
   
%   testF = Cqy*Geai; norm(testF(:,1))
   % Build the PCE of the solution
   Qb = prior + Cqy*Geai;
   
   sigmav = zeros(nindexa,1);
   CD     = zeros(nindexa); % Posterior variance
   for j=1:nindexa 
      for i=1:chaosor
         CD = CD + factorial(i)*Qb(:,i+1+(j-1)*(chaosor))*Qb(:,i+1+(j-1)*(chaosor))';
      end
   end
   sigmav = sqrt(diag(CD));
   meanD = Qb(:,1);
   
   % Diagonalize the covariance matrix
   %[~,CDt] = eig(CD); sigmav = sqrt(diag(CDt));
   
   figure;
   hold on;
   plot(meanD(1:nindexa/2),'Color','red');
   plot(uref(2*b2node2-1),'Color','green');
   plot(meanD(1:nindexa/2)-sigmav(1:nindexa/2),'Color','red');
   plot(meanD(1:nindexa/2)+sigmav(1:nindexa/2),'Color','red');
   legend('solution','reference')
   
%   figure;
%   hold on;
%   plot(meanD(nindexa/2+1:end),'Color','red');
%   plot(uref(2*b2node2),'Color','green');
%   plot(meanD(nindexa/2+1:end)-sigmav(nindexa/2+1:end),'Color','red');
%   plot(meanD(nindexa/2+1:end)+sigmav(nindexa/2+1:end),'Color','red');
%   legend('solution','reference')
   
   error = norm(uref(indexa)-meanD)/norm(uref(indexa));
   
   % Update the stuff
   prior2 = Qb;
   priorold = prior;
   prior = prior2-prior2;
   % Hack on the prior : make it uncoupled
   for j=1:nindexa 
      for i=1:nindexa
         prior(j,2+(j-1)*(chaosor):1+(j)*(chaosor)) =...
                          prior(j,2+(j-1)*(chaosor):1+(j)*(chaosor)) + ...
                          (prior2(j,2+(i-1)*(chaosor):1+(i)*(chaosor))).^2 ;
      end
   end
   prior = sqrt(prior);
   prior(:,1) = prior2(:,1);
%   prior = prior2;
end

% Plot the proba density on 1 dof
ddl = 1;
coefS1 = [];
for i=1:nindexa
   coefS1 = [ coefS1 , Hermite(2:chaosor+1,2:chaosor+1)' * transpose(Qb(ddl,2+(chaosor)*(i-1):1+(chaosor)*(i))) ];
end
% Plot the approximated theta
thetae  = randn(nindexa,nMC);
toplo = Qb(ddl,1);
for i=1:nindexa 
   for j=1:chaosor
      toplo = toplo + coefS1(j+chaosor*(i-1))*thetae(i,:).^j;
   end
end
figure;
hist(toplo,nMC/1000);
%% Plot an independant noise
%thetae  = randn(nindexa,nMC);
%toplo = 0;
%for i=1:nindexa 
%   toplo = toplo + 1/sqrt(nindexa)*thetae(i,:);
%end
%figure;
%hist(toplo,nMC/1000);