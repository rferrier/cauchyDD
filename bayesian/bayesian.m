%05/09/2017
%Inversion bayesienne basique du pb de Cauchy

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 70000; % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 1;    % N.mm-2 : Loading on the plate
br      = .1;      % noise
brt     = 0;      % "translation" noise
inhomog = 0;      % inhomogeneous medium
nMC     = 100000;     % nb of MC samples
refined = 2;      % Chooses the mesh to use

upperpr = 1e-3;   % Bounds of the prior distribution
lowerpr = -1e-3;  %

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
if refined == 5 b2node2 = 5; end

indexa = [2*b2node2-1; 2*b2node2]; nindexa = size(indexa,1);
indexb = [2*b2node4-1; 2*b2node4]; nindexb = size(indexb,1);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bayesian inversion
% Build the LHS
dirichletm = [1,1,0; 1,2,0 ;
              3,1,0; 3,2,0 ;
              2,1,0; 2,2,0];
neumannm   = [4,1,fscalar,0,-fscalar];
% TODO : the inhomogeneous trick !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
[Km,Cm,nbloqm,node2cm,c2nodem] = Krig2 (nodes,elements,mat,order,boundary,dirichletm);

prior = rand(nindexa,nMC);                  % Build the prior MRHS
prior = (upperpr-lowerpr)*prior + lowerpr;  %

Mdep  = zeros(2*nnodes,nMC);   %
Mdep(indexa,:) = prior;        % Prescripted displacement

%Mrhs = zeros(2*nnodes+nbloqm,nMC);
%for i=1:nMC
%   Mrhs(:,i) = dirichletRhs2( Mdep(:,i), 2, c2nodem, boundary, nnodes ); % Prior Mrhs
%end

fm = loading(nbloqm,nodes,boundary,neumannm);
Mrhs = dirichletRhs2( Mdep, 2, c2nodem, boundary, nnodes ); % Prior Mrhs

Mrhs = Mrhs + fm*ones(1,nMC); % Add the Neumann BC forces
% MC sampling model evaluation
tic
Msol = Km\Mrhs; Msol = Msol(indexb,:);
toc

% Estimate additive noise level  /!\/!\/!\ TOSEE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
nl = br*mean(uref(indexa)); % br
%nl = 1.; % Debug stuff
Co = nl^2*eye(nindexb);  % Correlation matrix
Co1 = inv(Co); factor = 1/((2*pi)^nindexb*det(Co));

% MC normalization factor evaluation
mesure = urefb(indexb); %mesure = mesure*ones(1,nMC);
pia = 0; meanD = zeros(nindexa,1); Vsq =  zeros(nindexa,1);
for i=1:nMC
   pia   = pia   + factor*exp(-1/2*(Msol(:,i)-mesure)'*Co1*(Msol(:,i)-mesure));
   meanD = meanD + Mdep(indexa,i).*factor*exp(-1/2*(Msol(:,i)-mesure)'*Co1*(Msol(:,i)-mesure));
end
meanD = meanD/pia; % normalize

for i=1:nMC % Variation
   Vsq = Vsq + ((Mdep(indexa,i)-meanD).^2).*factor*exp(-1/2*(Msol(:,i)-mesure)'*Co1*(Msol(:,i)-mesure));
end
Vsq = Vsq/pia;
sigmav = sqrt(Vsq);

figure;
hold on;
plot(meanD(1:nindexa/2),'Color','red');
plot(uref(2*b2node2-1),'Color','green');
plot(meanD(1:nindexa/2)-sigmav(1:nindexa/2),'Color','red');
plot(meanD(1:nindexa/2)+sigmav(1:nindexa/2),'Color','red');
legend('solution','reference')

figure;
hold on;
plot(meanD(nindexa/2+1:end),'Color','red');
plot(uref(2*b2node2),'Color','green');
plot(meanD(nindexa/2+1:end)-sigmav(nindexa/2+1:end),'Color','red');
plot(meanD(nindexa/2+1:end)+sigmav(nindexa/2+1:end),'Color','red');
legend('solution','reference')