% 18/09/2018 : Couplage entre Cauchy et décomposition de domaines,
% Surcharge de problèmes par pénalisation + Lagrange

clear all;
close all;

addpath(genpath('./tools'))

% Parameters
E       = 200000; % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 250;    % N.mm-2 : Loading on the plate
mat = [0, E, nu];
br      = 0.;      % noise

nsub    = 1; % nb subdomains = 2*nsub
niter   = 10;
precond = 0;
alpreac = 0;%1e-5; % Reaction penalization parameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Direct problem

%dirichlet   = [1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0];
dirichlet   = [2,1,0 ; 2,2,0];

% Import the meshes
[ nodes,elements,ntoelem,boundary,order ] =...
    readmesh( 'meshes/tube.msh' );
nnodes = size(nodes,1);

%noises = load('./noises/noisetube0.mat'); % Particular noise vector
%noise  = noises.bruit1;
noise  = randn(2*nnodes,1);

% patch('Faces',elements,'Vertices',nodes,'FaceAlpha',0);
% figure

% Then, build the stiffness matrix :
[K,C,nbloq,node2c,c2node] =...
    Krig2 (nodes,elements,mat,order,boundary,dirichlet,1);
Kinter = K(1:2*nnodes, 1:2*nnodes);
M      = mass_mat(nodes, elements);
[ node2b1, b2node1 ] = mapBound( 1, boundary, nnodes );
[ node2b2, b2node2 ] = mapBound( 2, boundary, nnodes );
[ node2b3, b2node3 ] = mapBound( 3, boundary, nnodes );
[ node2b9, b2node9 ] = mapBound( 9, boundary, nnodes );

% The right hand side :
f = pressureLoad( nbloq, nodes, boundary, fscalar*[10*sin(pi/6),0;-1,0], 8 );
plotGMSH({f(1:2*nnodes),'aife'}, elements, nodes(:,[1,2]), 'output/aife');

% Solve the problem :
uin = K\f;

% Extract displacement and Lagrange multiplicators :
uref = uin(1:2*nnodes,1);
fref = Kinter*uref;
lagr = uin(2*nnodes+1:end,1);
urefb = ( 1 + br*noise ) .* uref;
frefb = ( 1 + br*noise ) .* fref;

% Compute stress :
sigma = stress(uref,E,nu,nodes,elements,order,1,ntoelem,1);
% Output :
plotGMSH({uref,'U_vect';sigma,'stress'}, elements, nodes, 'output/reference');

%% Plot displacement on the interface :
%%index = 2*[b2node1;b2node2;b2node3];
index = 2*b2node3;
thetax = 0:2*pi/size(index,1):2*pi*(1-1/size(index,1));
try
figure
hold on
set(gca, 'fontsize', 15);
plot(thetax,uref(index,1));
plot(thetax,uref(index-1,1),'Color','red');
legend('uy','ux')
xlabel('angle(rad)')
end
try
figure
hold on
set(gca, 'fontsize', 15);
plot(thetax,fref(index,1));
plot(thetax,fref(index-1,1),'Color','red');
legend('fy','fx')
xlabel('angle(rad)')
end
%uref(index-1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Substructuring

% Find the dimensions of the plate
L1 = min(nodes(:,1)); L2 = max(nodes(:,1));
H1 = min(nodes(:,2)); H2 = max(nodes(:,2));

newelem  = {};
newnode  = {};
newbouns = {}; % global numerotation of the internal boundaries
bounsloc = {}; % local numerotation of the internal boundaries
map      = {};
boundar  = {}; % restriction of the external boundaries to each subdomain

newelem{1} = elements ;
newnode{1} = nodes ;
boundar{1} = boundary ;
map{1}     = 1:nnodes ;

%% First, split the domain

% separate in 2
elementsc = newelem{1}; nodesc = newnode{1}; boundaryc = boundar{1};
line = [ 0 ; H1 ; 0 ; H2];
[ nodes1, elements1, boundary1, map1,...
  nodes2, elements2, boundary2, map2, newbound ] =...
                                cutMesh (nodesc, elementsc, boundaryc, line);

b1to2 = (1:nnodes)';
                                
newelem{1}   = elements1; newelem{2} = elements2;
newnode{1}   = nodes1; newnode{2} = nodes2;
% Increment the map
oldmap = map{1}; map{1} = map1(oldmap); map{2} = map2(oldmap);
% get the local boundary
newbouns{1} = b1to2(newbound);
% The boundaries
boundar{1} = boundary1; boundar{2} = boundary2;
% store b1to2 (always useful)
b1to2s{1} = b1to2;

for i = 1:nsub-1
   elementsc1 = newelem{2*i-1}; elementsc2 = newelem{2*i};
   nodesc1    = newnode{2*i-1}; nodesc2    = newnode{2*i};
   boundaryc1 = boundar{2*i-1}; boundaryc2 = boundar{2*i};

   % Cut along a line
   theta1 = -pi/2 + i*pi/nsub; theta2 = -pi/2 - i*pi/nsub;
   line1 = [ 0 ; 0 ; cos(theta1) ; sin(theta1)];
   line2 = [ 0 ; 0 ; cos(theta2) ; sin(theta2)];
   [ nodes1, elements1, boundary1, map1,...
     nodes3, elements3, boundary3, map3, newbound1 ] =...
                               cutMesh (nodesc1, elementsc1, boundaryc1, line1);
   [ nodes4, elements4, boundary4, map4,...  % It's reversed on purpose
     nodes2, elements2, boundary2, map2, newbound2 ] =...
                               cutMesh (nodesc2, elementsc2, boundaryc2, line2);

   % Get the connectivity between the divided domain and the global domain
  [ b1to3, b3to1 ] = superNodes( nodesc1, nodes, 1e-6 );
  [ b2to4, b4to2 ] = superNodes( nodesc2, nodes, 1e-6 );
                                   
   newelem{2*i-1} = elements1; newelem{2*i} = elements2;
   newelem{2*i+1} = elements3; newelem{2*i+2} = elements4;
   newnode{2*i-1} = nodes1; newnode{2*i} = nodes2;
   newnode{2*i+1} = nodes3; newnode{2*i+2} = nodes4;
   
   % Increment the map
   oldmap1 = map{2*i-1};  oldmap2 = map{2*i};
   map{2*i-1} = map1(oldmap1); map{2*i} = map2(oldmap2);
   map{2*i+1} = map3(oldmap1); map{2*i+2} = map4(oldmap2);
   % get the local boundary
   newbouns{2*i} = b1to3(newbound1); newbouns{2*i+1} = b2to4(newbound2);
   % The boundaries
   boundar{2*i-1} = boundary1; boundar{2*i} = boundary2;
   boundar{2*i+1} = boundary3; boundar{2*i+2} = boundary4;
   % store b1to2 (always useful (or not))
%   b1to2s{i} = b1to2;
end

% Split the boundary 1
bo1 = newbouns{1}; nbo1 = nodes(bo1,:);
indup = find(nbo1(:,2)>0); inddo = find(nbo1(:,2)<0);
bon2 = bo1(indup); bon1 = bo1(inddo);
newbouns{1} = bon1; newbouns{2*nsub} = bon2;

% Link the boundaries to the subelemts
if nsub > 1
   el2bound = zeros(2*nsub, 2);
   el2bound([1,2],:) = [1,2;1,3];
   for i=3:2*nsub-1
      el2bound(i,:) = [i-1,i+1];
   end
   el2bound([2*nsub-1,2*nsub],:) = [2*nsub-2,2*nsub;2*nsub-1,2*nsub];
else
   el2bound = [1,2;1,2];
end

figure
hold on;
for i = 1:2*nsub
   elementsc = newelem{i};
   nodesc    = newnode{i};
   ret  = patch('Faces',elementsc(:,1:3),'Vertices',nodesc);
   col1 = rand(); col2 = rand(); col3 = rand(); coto = col1+col2+col3;
   set (ret, 'FaceColor', [col1,col2,col3]);
end
axis('equal');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin the serious stuff

% Define the Dirichlet's BCs
dirichlet1 = [2,1,0;2,2,0;
              1,1,0;1,2,0];
dirichlet2 = [2,1,0;2,2,0];
              
dirichlet1p = [2,1,0;2,2,0;
               1,1,0;1,2,0
               3,1,0;3,2,0];
dirichlet2p = [2,1,0;2,2,0;
               3,1,0;3,2,0];

% Stiffness, rigid modes & Co
K1 = {}; f1 = {}; nbloq1 = {}; 
Kp1 = {}; nbloqp1 = {}; Ct1 = {};
K2 = {}; f2 = {}; nbloq2 = {}; 
Kp2 = {}; nbloqp2 = {}; Ct2 = {};

K = {}; f = {}; nbloq = {}; 
Kp = {}; nbloqp = {}; Ct = {}; % trace operator on the boundary
G = {}; R = {}; Gglob = zeros( 6*nnodes, 3*nsub ); % Gglob's size will decrease
urb = {}; b1to2 = {}; cornoglob = {}; corno = {};

% Connectivity tables between the subdomains and the rigid modes
sdtorm = zeros(nsub,1); rmtosd = [];

Car = zeros(6*nnodes,0); % action/reaction matrix

j = 1; % index for the connectivity tables
for i = 1:2*nsub
   nodess    = newnode{i};
   elementss = newelem{i};
   boundarys = boundar{i};
   map1      = map{i};
   
   sdtorm(i) = 0;
   addboun = [];

   boun1       = map1(newbouns{el2bound(i,1)}); % Local boundaries
   bounsloc{2*i-1} = boun1;
   addboun = [addboun ; boun1];
   boun2       = map1(newbouns{el2bound(i,2)});
   bounsloc{2*i} = boun2;
   addboun = [addboun ; boun2];
   
   nno{i} = size(nodess,1);
   
   % Extract the local uref
   [ b1to2i, b2to1 ] = superNodes( nodess, nodes, 1e-6 ); % Yeah, I'm a pig
   b1to2{i} = b1to2i;
   urbi = zeros(2*nno{i},1);
   urbi( [1:2:2*nno{i}-1, 2:2:2*nno{i}] ) = urefb( [2*b1to2{i}-1 ; 2*b1to2{i}] );
   urb{i} = urbi;
   
   % Corners (ie nodes with DD and Cauchy)
   boun3 = boundarys( find( boundarys(:,1)==3 ) , 2:3 );
   boun3 = unique( [ boun3(:,1) ; boun3(:,2) ] );
   corno{i} = [ intersect(boun3,bounsloc{2*i-1}) , ...
                intersect(boun3,bounsloc{2*i}) ];
   cornoglob{i} = b1to2i(corno{i});

   % Stiffness matrices
   [K{i},C,nbloq{i},node2c,c2node] =...
                 Krig2 (nodess,elementss,mat,order,boundarys,dirichlet,1);
   [K1{i},C,nbloq1{i},node2c1,c2node1] =...
                 Krig2 (nodess,elementss,mat,order,boundarys,dirichlet1,1);
   [K2{i},C,nbloq2{i},node2c2,c2node2] =...
                 Krig2 (nodess,elementss,mat,order,boundarys,dirichlet2,1);

   f{i} = pressureLoad( nbloq{i}, nodess, boundarys, ...
                         fscalar*[10*sin(pi/6),0;-1,0], 8 );
   f1{i} = dirichletRhs2( urb{i}, 1, c2node1, boundarys, nno{i} );
   f2{i} = zeros(2*nno{i}+nbloq2{i},1);  % every known Neumann is 0
   
   R{i} = [];
   
   % Primal problem : 
   [Kpt,C,nbloqp{i},node2c,c2node] =...
                 Krig2 (nodess,elementss,mat,order,boundarys,dirichlet,1);
   [Kpt1,C1,nbloqp1{i},node2c1,c2node1] =...
                 Krig2 (nodess,elementss,mat,order,boundarys,dirichlet1p,1);
   [Kpt2,C2,nbloqp2{i},node2c2,c2node2] =...
                 Krig2 (nodess,elementss,mat,order,boundarys,dirichlet2p,1);
   Kpinter = Kpt(1:2*nno{i}, 1:2*nno{i});
   nbloqp{i} = nbloqp{i} + 2*size(addboun,1);
   
   % See witch dof has already Dirichlet conditions
   testC = C*ones(size(C,2),1);
   testC1 = C1*ones(size(C1,2),1);
   testC2 = C2*ones(size(C2,2),1);
   
   %% Lagrange multipliers for addboun
   %C1c = zeros( 2*nno{i}, 2*size(addboun) ); % intermediate matrix
   C1c = zeros( 2*nno{i}, 0 ); % intermediate matrix
   index = 0;
   for no=1:size(addboun,1)
      if testC(2*addboun(no)-1) == 0
         index = index+1;
         C1c(2*addboun(no)-1, index) = 1;
      end
      if testC(2*addboun(no)) == 0 
         index = index+1;
         C1c(2*addboun(no), index) = 1;
      end
   end
   Ct{i} = [ C, C1c ];
   Kp{i} = [Kpinter, Ct{i} ; Ct{i}', zeros(size(Ct{i},2))];
   
   C1c = zeros( 2*nno{i}, 0 ); % intermediate matrix
   index = 0;
   for no=1:size(addboun,1)
      if testC1(2*addboun(no)-1) == 0
         index = index+1;
         C1c(2*addboun(no)-1, index) = 1;
      end
      if testC1(2*addboun(no)) == 0 
         index = index+1;
         C1c(2*addboun(no), index) = 1;
      end
   end
   Ct1{i} = [ C1, C1c ];
   Kp1{i} = [Kpinter, Ct1{i} ; Ct1{i}', zeros(size(Ct1{i},2))];

   C1c = zeros( 2*nno{i}, 0 ); % intermediate matrix
   index = 0;
   for no=1:size(addboun,1)
      if testC2(2*addboun(no)-1) == 0
         index = index+1;
         C1c(2*addboun(no)-1, index) = 1;
      end
      if testC2(2*addboun(no)) == 0 
         index = index+1;
         C1c(2*addboun(no), index) = 1;
      end
   end
   Ct2{i} = [ C2, C1c ];
   Kp2{i} = [Kpinter, Ct2{i} ; Ct2{i}', zeros(size(Ct2{i},2))];
   
   sign = 1;
   if rem(i,4) == 0 || rem(i-1,4) == 0
      sign = -1;
   end

   % Management of the rigid modes
   if nbloq2{i} == 0 % nbloq < 3 is also problematic : ex if there is 1 encastred point /!\
   
      % Connectivity stuff (subdomains VS rigid modes)
      sdtorm(i) = j;
      rmtosd(j) = i;
   
      nnodess = nno{i};

      R{i} = null(K2{i}); Rloc = R{i};
      
      nbloq2{i} = size(R{i},2);
      [ node2b3s, b2node3s ] = mapBound( 3, boundarys, nnodess );
      b1to2s = b1to2{i};
      b2node3g = intersect(b1to2s, b2node3);

      Gloc = zeros( 2*nnodess, 3 ); % local rigid modes
      Gloc( [ 2*b2node3s-1, 2*b2node3s ], : ) = ...
                                   Rloc( [ 2*b2node3s-1, 2*b2node3s ], : );
      Gloc( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ], : ) = ...
                  Rloc( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ], : );
      Gloc( [ 2*bounsloc{2*i}-1, 2*bounsloc{2*i} ], : ) = ...
                  Rloc( [ 2*bounsloc{2*i}-1, 2*bounsloc{2*i} ], : );

      % Global rigid modes
      if sign == 1
         Gglob( [ 2*b2node3g-1, 2*b2node3g ], 3*j-2:3*j ) = ...
                                   Rloc( [ 2*b2node3s-1, 2*b2node3s ], : );
         Gglob( [ 2*newbouns{el2bound(i,1)}-1, 2*newbouns{el2bound(i,1)} ], 3*j-2:3*j ) = ...
                Rloc( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ], : );
         Gglob( [ 2*newbouns{el2bound(i,2)}-1, 2*newbouns{el2bound(i,2)} ], 3*j-2:3*j ) = ...
                Rloc( [ 2*bounsloc{2*i}-1, 2*bounsloc{2*i} ], : );
      else
         Gglob( 2*nnodes + [ 2*b2node3g-1, 2*b2node3g ], 3*j-2:3*j ) = ...
                                   Rloc( [ 2*b2node3s-1, 2*b2node3s ], : );
         Gglob( 2*nnodes + [ 2*newbouns{el2bound(i,1)}-1, 2*newbouns{el2bound(i,1)} ], 3*j-2:3*j ) = ...
                Rloc( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ], : );
         Gglob( 2*nnodes + [ 2*newbouns{el2bound(i,2)}-1, 2*newbouns{el2bound(i,2)} ], 3*j-2:3*j ) = ...
                Rloc( [ 2*bounsloc{2*i}-1, 2*bounsloc{2*i} ], : );
      end

      G{i} = Gloc;
      j = j+1;
   end
end
Gglob(:,3*j-2:end) = [];  % reshape Gglob

for i=1:2*nsub % Loop over the interfaces
   % Find the associated sub domains
   if i == 1
      sdom1 = 1; sdom2 = 2;
      ind1 = newbouns{el2bound(sdom1,1)}; % Find the corresponding dofs
      ind2 = newbouns{el2bound(sdom2,1)};
   elseif i == 2*nsub
      sdom1 = 2*nsub; sdom2 = 2*nsub-1;
      ind1 = newbouns{el2bound(sdom1,2)};
      ind2 = newbouns{el2bound(sdom2,2)};
   elseif rem(i,2) == 0
      sdom1 = i; sdom2 = i+2;
      ind1 = newbouns{el2bound(sdom1,2)};
      ind2 = newbouns{el2bound(sdom2,1)};
   else%if rem(i,2) == 0
      sdom1 = i-2; sdom2 = i;
      ind1 = newbouns{el2bound(sdom1,2)};
      ind2 = newbouns{el2bound(sdom2,1)};
   end
   ind1 = setdiff( ind1, cornoglob{sdom1} ); % Remove the Cauchy dofs (as load can be applied there)
   ind2 = setdiff( ind2, cornoglob{sdom2} );
   nnobo = size(ind1,1);

   sign1 = 1;
   if rem(sdom1,4) == 0 || rem(sdom1-1,4) == 0
      sign1 = -1;
   end
   sign2 = 1;
   if rem(sdom2,4) == 0 || rem(sdom2-1,4) == 0
      sign2 = -1;
   end

   % Management of the action/reaction principle
   Carloc = zeros( 6*nnodes, 2*nnobo );
   if sign1 == 1 % we know that sign1=-sign2
      Carloc( [2*ind1-1;2*ind1], : ) = eye(2*nnobo);
      Carloc( 2*nnodes + [2*ind2-1;2*ind2], : ) = eye(2*nnobo);
   else
      Carloc( [2*ind2-1;2*ind2], : ) = eye(2*nnobo);
      Carloc( 2*nnodes + [2*ind1-1;2*ind1], : ) = eye(2*nnobo);
   end
   Car = [ Car, Carloc ]; % Car'*Lambda = 0 <=> action/reaction principle
end

% rework the rigid modes
if j > 1
   eD = zeros( size(Gglob,2), 1 );  % 0 of the same size as Gglob
   
   for i = 1:2*nsub
      if isempty(R{i}) % We do not loop over empty matrices
         continue;
      end

      j = sdtorm(i);   % Recover the rigid mode no
      eD(3*j-2:3*j) = eD(3*j-2:3*j) + R{i}'*f2{i};

      K2{i} = [ K2{i}, G{i} ; G{i}', zeros(size(G{i},2)) ];
      f2{i} = [ f2{i} ; zeros( size(G{i},2), 1 ) ];

   end
   P = eye(6*nnodes) - Gglob*inv(Gglob'*Gglob)*Gglob';  % The projector
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ( constraint - ) CG algorithm
Lambda = zeros( 6*nnodes,1);       % Iterate (on the global mesh for convenience)
Res    = zeros( 6*nnodes,niter+1 );
Rez    = zeros( 6*nnodes,niter+1 ); % Residual without orthogonal projection
Zed    = zeros( 6*nnodes,niter+1 );
d      = zeros( 6*nnodes,niter+1 );
Ad     = zeros( 6*nnodes,niter+1 );
Az     = zeros( 6*nnodes,niter+1 ); % Ad without the projection
resid  = zeros( niter+1,1 );
error  = zeros( niter+1,1 );  % Useless for now (except for debug)

tic
% Compute Rhs
b = zeros( 6*nnodes, 1 );
for i = 1:2*nsub
   ui1 = K1{i}\f1{i};
   ui2 = K2{i}\f2{i};
   u  = zeros( 2*nnodes, 1 ); uc = zeros( 2*nnodes, 1 );
   ud = zeros( 2*nnodes, 1 );
   
   u1 = keepField( ui1, 3, boundar{i} ); % SP Fields
   u2 = keepField( ui2, 3, boundar{i} ); %
   
   sign = -1;
   if rem(i,4) == 0 || rem(i-1,4) == 0
      sign = 1; % Sign is the opposite as usually
   end
               
   % Cauchy operator (including DD bounds)
   uc( [2*b1to2{i}-1 ; 2*b1to2{i}] ) = ... %u( [2*b1to2{i}-1 ; 2*b1to2{i}] ) + ...
          u1( [1:2:2*nno{i}-1, 2:2:2*nno{i}] ) - u2( [1:2:2*nno{i}-1, 2:2:2*nno{i}] );

   uc( [ 2*newbouns{el2bound(i,1)}-1, 2*newbouns{el2bound(i,1)} ]  ) =...
          ui1( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ] ) - ui2( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ] );

   uc( [ 2*newbouns{el2bound(i,2)}-1, 2*newbouns{el2bound(i,2)} ]  ) =...
          ui1( [ 2*bounsloc{2*i}-1, 2*bounsloc{2*i} ] ) - ui2( [ 2*bounsloc{2*i}-1, 2*bounsloc{2*i} ] );
          
   % Domain decomposition part (on Neumann domains)
   ud( [ 2*newbouns{el2bound(i,1)}-1, 2*newbouns{el2bound(i,1)} ]  ) = ...
          - ui2( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ] );
   ud( [ 2*newbouns{el2bound(i,2)}-1, 2*newbouns{el2bound(i,2)} ]  ) = ...
          - ui2( [ 2*bounsloc{2*i}-1, 2*bounsloc{2*i} ] );

   b(4*nnodes+1:end) = b(4*nnodes+1:end) + sign*ud;
   if sign == -1
      b(1:2*nnodes) = b(1:2*nnodes) + uc;
   else
      b(2*nnodes+1:4*nnodes) = b(2*nnodes+1:4*nnodes) + uc;
   end
end
plotGMSH({b(1:2*nnodes),'be1';b(2*nnodes+1:4*nnodes),'be2'}, elements, nodes(:,[1,2]), 'output/be');
% initialize
if norm(Gglob) ~= 0
   Lambda = -Gglob*inv(Gglob'*Gglob)*eD ;
end

% Compute Ax0
Axz = zeros( 6*nnodes, 1 );
for i = 1:2*nsub
   floc   = zeros( 2*nno{i}, 1 );
   flocDD = zeros( 2*nno{i}+nbloq2{i}, 1 );

   sign = 1;
   if rem(i,4) == 0 || rem(i-1,4) == 0
      sign = -1;
   end
   
   % Cauchy part
   if sign == 1
      floc( [1:2:2*nno{i}-1, 2:2:2*nno{i}] ) = ...
                          Lambda( [2*b1to2{i}-1 ; 2*b1to2{i}] ) ;

      floc( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ] ) = ...
                  Lambda( [ 2*newbouns{el2bound(i,1)}-1, 2*newbouns{el2bound(i,1)} ] );
   else
      floc( [1:2:2*nno{i}-1, 2:2:2*nno{i}] ) = ...
                          Lambda( 2*nnodes + [2*b1to2{i}-1 ; 2*b1to2{i}] ) ;

      floc( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ] ) =...
                  Lambda( 2*nnodes + [ 2*newbouns{el2bound(i,1)}-1, 2*newbouns{el2bound(i,1)} ] );
   end

   flocDD( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ] ) = ...
               Lambda( 4*nnodes + [ 2*newbouns{el2bound(i,1)}-1, 2*newbouns{el2bound(i,1)} ] ); % RHS For the transpose problem
   
   floc1 = [ floc ; zeros(nbloq1{i},1) ]; floc2 = [ floc ; zeros(nbloq2{i},1) ];
   
   ui1  = K1{i}\floc1;
   ui2i = K2{i}\[floc2,flocDD];
   ui2 = ui2i(:,1); uiDD = ui2i(:,2);
   u = zeros( 2*nnodes, 1 ); uc = zeros( 2*nnodes, 1 );
   ud = zeros( 2*nnodes, 1 ); udt = zeros( 2*nnodes, 1 );
   
   u1  = keepField( ui1, 3, boundar{i} ); % SP Fields (keep only the bound)
   u2  = keepField( ui2, 3, boundar{i} ); %
   uDD = keepField( uiDD, 3, boundar{i} ); %

   % Cauchy operator (including DD bounds)
   uc( [2*b1to2{i}-1 ; 2*b1to2{i}] ) = ...
          u2( [1:2:2*nno{i}-1, 2:2:2*nno{i}] ) - u1( [1:2:2*nno{i}-1, 2:2:2*nno{i}] );
   uc( [ 2*newbouns{el2bound(i,1)}-1, 2*newbouns{el2bound(i,1)} ]  ) = ...
          ui2( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ] ) - ui1( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ] );
   uc( [ 2*newbouns{el2bound(i,2)}-1, 2*newbouns{el2bound(i,2)} ]  ) = ...
          ui2( [ 2*bounsloc{2*i}-1, 2*bounsloc{2*i} ] ) - ui1( [ 2*bounsloc{2*i}-1, 2*bounsloc{2*i} ] );

   % Transpose part
   udt( [2*b1to2{i}-1 ; 2*b1to2{i}] ) = uDD( [1:2:2*nno{i}-1, 2:2:2*nno{i}] );

   udt( [ 2*newbouns{el2bound(i,1)}-1, 2*newbouns{el2bound(i,1)} ]  ) = ...
          uiDD( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ] );
   udt( [ 2*newbouns{el2bound(i,2)}-1, 2*newbouns{el2bound(i,2)} ]  ) = ...
          uiDD( [ 2*bounsloc{2*i}-1, 2*bounsloc{2*i} ] );

   % Lagrange part
   ud( [ 2*newbouns{el2bound(i,1)}-1, 2*newbouns{el2bound(i,1)} ]  ) = ...
          ui2( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ] );
   ud( [ 2*newbouns{el2bound(i,2)}-1, 2*newbouns{el2bound(i,2)} ]  ) = ...
          ui2( [ 2*bounsloc{2*i}-1, 2*bounsloc{2*i} ] );

   Axz(4*nnodes+1:end) = Axz(4*nnodes+1:end) + sign*ud;
   if sign == 1
      Axz(1:2*nnodes) = Axz(1:2*nnodes) + uc + sign*udt;
   else
      Axz(2*nnodes+1:4*nnodes) = Axz(2*nnodes+1:4*nnodes) + uc + sign*udt;
   end
end

% Add the reaction condition
CarCar = Car*Car';
Axz = Axz + alpreac*CarCar*Lambda;

Rez(:,1) = b - Axz;
if norm(Gglob) ~= 0
   Res(:,1) = P'*(b - Axz) ;
else
   Res(:,1) = b - Axz;
end

if precond == 1
   Zed(:,1) = Res(:,1);
else
   Zed(:,1) = Res(:,1);
end

if norm(Gglob) ~= 0
   Zed(:,1) = P*Zed(:,1) ;
else
   Zed(:,1) = Zed(:,1);
end

resid(1) = norm( Res(:,1) );

d(:,1) = Zed(:,1);

%% Loop over the iterations
iter = 0;  % in case niter=0
for iter = 1:niter

   % Compute Ad
   Ad(:,iter) = zeros( 6*nnodes, 1 );
   for i = 1:2*nsub

      floc = zeros( 2*nno{i}, 1 );
      flocDD = zeros( 2*nno{i}+nbloq2{i}, 1 );

      % Now, we need a sign      
      sign = 1;
      if rem(i,4) == 0 || rem(i-1,4) == 0
         sign = -1;
      end
      
      % Cauchy part
      if sign == 1
         floc( [1:2:2*nno{i}-1, 2:2:2*nno{i}] ) = ...
                             d( [2*b1to2{i}-1 ; 2*b1to2{i}], iter ) ;
   
         floc( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ] ) =...
                     d( [ 2*newbouns{el2bound(i,1)}-1, 2*newbouns{el2bound(i,1)} ], iter );
      else
         floc( [1:2:2*nno{i}-1, 2:2:2*nno{i}] ) = ...
                             d( 2*nnodes+[2*b1to2{i}-1 ; 2*b1to2{i}], iter ) ;
   
         floc( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ] ) =...
                     d( 2*nnodes + [ 2*newbouns{el2bound(i,1)}-1, 2*newbouns{el2bound(i,1)} ], iter );
      end
   
      flocDD( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ] ) = ...  % RHS For the transpose problem
                  d( 4*nnodes + [ 2*newbouns{el2bound(i,1)}-1, 2*newbouns{el2bound(i,1)} ], iter );

      floc1 = [ floc ; zeros(nbloq1{i},1) ]; floc2 = [ floc ; zeros(nbloq2{i},1) ];
   
      ui1 = K1{i}\floc1;
      ui2i = K2{i}\[floc2,flocDD];
      ui2 = ui2i(:,1); uiDD = ui2i(:,2);
      u = zeros( 2*nnodes, 1 ); uc = zeros( 2*nnodes, 1 );
      ud = zeros( 2*nnodes, 1 ); udt = zeros( 2*nnodes, 1 );
   
      u1 = keepField( ui1, 3, boundar{i} ); % SP Fields
      u2 = keepField( ui2, 3, boundar{i} ); %
      uDD = keepField( uiDD, 3, boundar{i} );

      % Cauchy operator (including DD bounds)
      uc( [2*b1to2{i}-1 ; 2*b1to2{i}] ) = ...
             u2( [1:2:2*nno{i}-1, 2:2:2*nno{i}] ) - u1( [1:2:2*nno{i}-1, 2:2:2*nno{i}] );
      uc( [ 2*newbouns{el2bound(i,1)}-1, 2*newbouns{el2bound(i,1)} ]  ) = ...
             ui2( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ] ) - ui1( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ] );
      uc( [ 2*newbouns{el2bound(i,2)}-1, 2*newbouns{el2bound(i,2)} ]  ) = ...
             ui2( [ 2*bounsloc{2*i}-1, 2*bounsloc{2*i} ] ) - ui1( [ 2*bounsloc{2*i}-1, 2*bounsloc{2*i} ] );

      % Transpose part
      udt( [2*b1to2{i}-1 ; 2*b1to2{i}] ) = uDD( [1:2:2*nno{i}-1, 2:2:2*nno{i}] );

      udt( [ 2*newbouns{el2bound(i,1)}-1, 2*newbouns{el2bound(i,1)} ]  ) = ...
             uiDD( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ] );
      udt( [ 2*newbouns{el2bound(i,2)}-1, 2*newbouns{el2bound(i,2)} ]  ) = ...
             uiDD( [ 2*bounsloc{2*i}-1, 2*bounsloc{2*i} ] );

      % Lagrange part (only on the DD bounds)
      ud( [ 2*newbouns{el2bound(i,1)}-1, 2*newbouns{el2bound(i,1)} ]  ) = ...
             ui2( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ] );
      ud( [ 2*newbouns{el2bound(i,2)}-1, 2*newbouns{el2bound(i,2)} ]  ) = ...
             ui2( [ 2*bounsloc{2*i}-1, 2*bounsloc{2*i} ] );

      Ad(4*nnodes+1:end,iter) = Ad(4*nnodes+1:end,iter) + sign*ud;
      if sign == 1
         Ad(1:2*nnodes,iter) = Ad(1:2*nnodes,iter) + uc + sign*udt;
      else
         Ad(2*nnodes+1:4*nnodes,iter) = Ad(2*nnodes+1:4*nnodes,iter) + uc + sign*udt;
      end
   end

   % Add the reaction condition
   Ad(:,iter) = Ad(:,iter) + alpreac*CarCar*d(:,iter);

   Az(:,iter) = Ad(:,iter);
   if norm(Gglob) ~= 0
      Ad(:,iter) = P'*Ad(:,iter);
   end

   den = (d(:,iter)'*Ad(:,iter));
   d(:,iter) = d(:,iter)/sqrt(den); Ad(:,iter) = Ad(:,iter)/sqrt(den);
   Az(:,iter) = Az(:,iter)/sqrt(den);
   num = Res(:,iter)'*d(:,iter);
   Lambda        = Lambda + d(:,iter)*num;
   Res(:,iter+1) = Res(:,iter) - Ad(:,iter)*num;
   Rez(:,iter+1) = Rez(:,iter) - Az(:,iter)*num;
   
   if precond == 1
      Zed(:,iter+1) = Res(:,iter+1);
   else
      Zed(:,iter+1) = Res(:,iter+1);
   end

   if norm(Gglob) ~= 0
      Zed(:,iter+1) = P*Zed(:,iter+1);
      resid(iter+1) = norm( P*Res(:,iter+1) );
   else
      Zed(:,iter+1) = Zed(:,iter+1);
      resid(iter+1) = norm(Res(:,iter+1));
   end
   
   % Orthogonalization
   d(:,iter+1) = Zed(:,iter+1);
   for jter=1:iter
       betaij = ( Zed(:,iter+1)'*Ad(:,jter) );
       d(:,iter+1) = d(:,iter+1) - d(:,jter) * betaij;
   end
   
end
disp([ 'End of the CG ', num2str(toc) ]);

try
figure; % plot the residual
plot(log10(resid/resid(1)));
end

plotGMSH({Lambda(1:2*nnodes),'Lambda1';Lambda(2*nnodes+1:4*nnodes),'Lambda2';Lambda(1:2*nnodes)+Lambda(2*nnodes+1:4*nnodes),'SumLambda'}, elements, nodes(:,[1,2]), 'output/lambda');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the final solution and assemble it
usol = zeros( 2*nnodes, 1 );

if norm(Gglob) ~= 0
   alphaD = inv(Gglob'*Gglob)*Gglob'*Rez(:,niter+1);
end

for i = 1:2*nsub
   nodess    = newnode{i};
%   elementss = newelem{i};
   nnodess   = nno{i};
   
   sign = 1;
   if rem(i,4) == 0 || rem(i-1,4) == 0
      sign = -1;
   end
%   sign=1;
   floc = zeros( 2*size(newnode{i},1), 1 );

   if sign == 1
      floc( [1:2:2*nno{i}-1, 2:2:2*nno{i}] ) = ...
                                Lambda( [2*b1to2{i}-1 ; 2*b1to2{i}] ) ;
      floc( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ] ) =...
                  Lambda( [ 2*newbouns{el2bound(i,1)}-1, 2*newbouns{el2bound(i,1)} ] );
      floc( [ 2*bounsloc{2*i}-1, 2*bounsloc{2*i} ] ) =...
                  Lambda( [ 2*newbouns{el2bound(i,2)}-1, 2*newbouns{el2bound(i,2)} ] );
   else             
      floc( [1:2:2*nno{i}-1, 2:2:2*nno{i}] ) = ...
                                Lambda( 2*nnodes+[2*b1to2{i}-1 ; 2*b1to2{i}] ) ;
      floc( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ] ) =...
                  Lambda( 2*nnodes + [ 2*newbouns{el2bound(i,1)}-1, 2*newbouns{el2bound(i,1)} ] );
      floc( [ 2*bounsloc{2*i}-1, 2*bounsloc{2*i} ] ) =...
                  Lambda( 2*nnodes + [ 2*newbouns{el2bound(i,2)}-1, 2*newbouns{el2bound(i,2)} ] );
   end

   floc1 = [ floc ; zeros(nbloq1{i},1) ];
   floc2 = [ floc ; zeros(nbloq2{i},1) ];

   % In case Cauchy and DD contribute to these points
%   floc( [ 2*corno{i}-1, 2*corno{i} ] ) = floc( [ 2*corno{i}-1, 2*corno{i} ] ) / 2;

   ui1 = K1{i}\(floc1+f1{i});
   ui2 = K2{i}\(floc2+f2{i});
   if ~isempty(R{i}) % Rigid modes
      j = sdtorm(i);
      ui2 = ui2 + [ R{i}*alphaD(3*j-2:3*j) ; zeros(nbloq2{i},1) ];
   end
   
   u = zeros( 2*nnodes, 1 );
%   [ b1to2, b2to1 ] = superNodes( nodess, nodes, 1e-6 );
   u( [2*b1to2{i}-1 ; 2*b1to2{i}] ) = ui2( [1:2:2*nnodess-1, 2:2:2*nnodess] );
   
   % Special treatment for the interface nodes

   u( [ 2*newbouns{el2bound(i,1)}-1, 2*newbouns{el2bound(i,1)} ]  ) =...
               .5 * ui2( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ] );
   u( [ 2*newbouns{el2bound(i,2)}-1, 2*newbouns{el2bound(i,2)} ] ) =...
               .5 * ui2( [ 2*bounsloc{2*i}-1, 2*bounsloc{2*i} ] );
   
   usol = usol + u;
end

uti = reshape(usol,2,[])';  ux = uti(:,1);  uy = uti(:,2);
plotGMSH({ux,'U_x';uy,'U_y';usol,'U_vect'}, elements, nodes(:,[1,2]), 'output/solution');
total_error = norm(usol-uref)/norm(uref);
plotGMSH({usol-uref,'Error'}, elements, nodes(:,[1,2]), 'output/error');

% Plot the solution
fsol = Kinter*usol;

index = 2*b2node3;
try
figure
hold on
set(gca, 'fontsize', 15);
plot(thetax,usol(index,1));
plot(thetax,usol(index-1,1),'Color','red');
legend('uy','ux')
xlabel('angle(rad)')
end
%try
%figure
%hold on
%set(gca, 'fontsize', 15);
%plot(thetax,fsol(index,1));
%plot(thetax,fsol(index-1,1),'Color','red');
%legend('fy','fx')
%xlabel('angle(rad)')
%end
lambd = Lambda(1:2*nnodes)+Lambda(2*nnodes+1:4*nnodes);
try
figure
hold on
set(gca, 'fontsize', 15);
plot(thetax,lambd(index,1));
plot(thetax,lambd(index-1,1),'Color','red');
legend('fy','fx')
xlabel('angle(rad)')
end
