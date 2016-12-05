% 01/12/2016 : Algo de Schur Dual avec plusieurs sous-domaines

clear all;
close all;

addpath(genpath('./tools'))

% Parameters
E       = 200000; % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 250;    % N.mm-2 : Loading on the plate
mat = [0, E, nu];

nsub  = 5;
niter = 20;

[ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/plate.msh' );
nnodes = size(nodes,1);

%dirichlet = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ; 4,1,0 ; 4,2,0 ];
dirichlet = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ];
%dirichlet = [ 1,1,0 ; 1,2,0 ];
neumann   = [ 3,2,fscalar ];

%figure
%patch('Faces',elements(:,1:3),'Vertices',nodes,'FaceAlpha',0);
%axis('equal');

%% Submeshes

% Find the dimensions of the plate
L1 = min(nodes(:,1)); L2 = max(nodes(:,1));
H1 = min(nodes(:,2)); H2 = max(nodes(:,2));

newelem  = {};
newnode  = {};
newbouns = {};
bounsloc = {}; % local numerotation of the internal boundaries
map      = {};
boundar  = {};

newelem{1} = elements ;
newnode{1} = nodes ;
boundar{1} = boundary ;
map{1}     = 1:nnodes ;

for i = 1:nsub-1
   elementsc = newelem{i};
   nodesc    = newnode{i};
   boundaryc = boundar{i};

   % Cut along a line
   line = [ L1 ; i*(H1+H2)/nsub ; L2 ; i*(H1+H2)/nsub];
   [ nodes1, elements1, boundary1, map1,...
     nodes2, elements2, boundary2, map2, newbound ] =...
                                   cutMesh (nodesc, elementsc, boundaryc, line);

   % Get the connectivity between the divided domain and the global domain
   [ b1to2, b2to1 ] = superNodes( nodesc, nodes, 1e-6 );
                                   
   newelem{i}   = elements1;
   newelem{i+1} = elements2;
   newnode{i}   = nodes1;
   newnode{i+1} = nodes2;
   
   % Increment the map
   oldmap = map{i};
   map{i} = map1(oldmap);
   map{i+1} = map2(oldmap);

   % get the local boundary
   %map1 = map{i}; map2 = map{i+1};
%   newbound1   = map1(newbound);
%   newbound2   = map2(newbound);
%   newbouns{2*i-1} = newbound1; newbouns{2*i} = newbound2;
   newbouns{i} = b1to2(newbound);
   
   % The boundaries
   boundar{i} = boundary1; boundar{i+1} = boundary2;
end

figure
hold on;
for i = 1:nsub
   elementsc = newelem{i};
   nodesc    = newnode{i};
   ret  = patch('Faces',elementsc(:,1:3),'Vertices',nodesc);
   col1 = rand(); col2 = rand(); col3 = rand(); coto = col1+col2+col3;
   set (ret, 'FaceColor', [col1/coto,col2/coto,col3/coto]);
end
axis('equal');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reference solution
[K,C,nbloq,node2c,c2node] = Krig2 (nodes,elements,mat,order,boundary,dirichlet);
f = loading(nbloq,nodes,boundary,neumann);
u = K\f; u = u(1:2*nnodes,1);
sigma = stress(u,E,nu,nodes,elements,order,1,ntoelem);
plotGMSH({u,'U_vect';sigma,'stress'}, elements, nodes(:,[1,2]), 'reference');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin the serious stuff

% Stiffness, rigid modes & Co
K = {}; f = {}; nbloq = {};
for i = 1:nsub
   nodess    = newnode{i};
   elementss = newelem{i};
   boundarys = boundar{i};
   map1      = map{i};
   %map2      = map{i+1};
   
   if i-1 > 0
      boun1       = map1(newbouns{i-1}); % Local boundaries
      bounsloc{2*i-1} = boun1;
   end
   if i < nsub
      boun2       = map1(newbouns{i});
      bounsloc{2*i} = boun2;
   end
   
   nno{i} = size(nodess,1);
   [K{i},C,nbloq{i},node2c,c2node] =...
                 Krig2 (nodess,elementss,mat,order,boundarys,dirichlet);
   f{i} = loading(nbloq{i},nodess,boundarys,neumann);
   
   % Management of the rigid modes (TODO)
%   if nbloq{i} == 0 % nbloq == 1 is also problematic, but ...
%      nnodess = nno{i};
%      r1 = zeros(2*nnodess,1); r2 = r1; r3p = r1; 
%      g1 = r1; g2 = r1; g3p = r1;
%      ind = 2:2:2*nnodess;
%      r1(ind-1,1) = 1; r2(ind,1) = 1;
%      
%
%      r3p(ind-1,1) = -nodess(ind/2,2);
%      r3p(ind,1) = nodess(ind/2,1);
%      g1(boun1) = r1(boun1); g2(boun1) = r2(boun1); g3p(boun1) = r3p(boun1);
%      % orthonormalize G and ensure G = trace(R)
%      g3 = g3p - (g3p'*g1)/(g1'*g1)*g1 - (g3p'*g2)/(g2'*g2)*g2;
%      r3 = r3p - (r3p'*g1)/(g1'*g1)*r1 - (r3p'*g2)/(g2'*g2)*r2;
%      
%      ng1 = norm(g1); ng2 = norm(g2); ng3 = norm(g3);
%      g1 = g1/ng1; g2 = g2/ng2; g3 = g3/ng3;
%      r1 = r1/ng1; r2 = r2/ng2; r3 = r3/ng3;
%   
%      G = [g1,g2,g3];
%      R = [r1,r2,r3];
%      K{i} = [ K{i}, R ; R', zeros(size(R,2)) ];
%      nbloq{i} = size(R,2);
%   end
   
end

%% CG algorithm
Lambda = zeros(2*nnodes,1);       % Iterate (on the global mesh for convenience)
Res    = zeros(2*nnodes,niter+1);
Zed    = zeros(2*nnodes,niter+1);
d      = zeros(2*nnodes,niter+1);
resid  = zeros( niter+1,1 );
error  = zeros( niter+1,1 );

% Compute Rhs
b = zeros( 2*nnodes, 1 );
for i = 1:nsub
   ui = K{i}\f{i};
   u = zeros( 2*nnodes, 1 );
   if i>1 % Lower boundary
      u(newbouns{i-1}) = ui(bounsloc{2*i-1});
   end
   if i<nsub
      u(newbouns{i}) = ui(bounsloc{2*i});
   end
   b = b+u;
end

% Compute Ax0
Axz = zeros( 2*nnodes, 1 );
for i = 1:nsub
   floc = zeros( 2*size(newnode{i},1) + nbloq{i}, 1 );
   if i>1 % Lower boundary
      floc(bounsloc{2*i-1}) = Lambda(newbouns{i-1});
   end
   if i<nsub
      floc(bounsloc{2*i}) = Lambda(newbouns{i});
   end
   
   ui = K{i}\floc;
   u = zeros( 2*nnodes, 1 );
   if i>1 % Lower boundary
      u(newbouns{i-1}) = ui(bounsloc{2*i-1});
   end
   if i<nsub
      u(newbouns{i}) = ui(bounsloc{2*i});
   end
   Axz = Axz+u;
end

Res(:,1) = b-Axz;
resid(1) = norm( Res(:,1) );

d(:,1) = Res(:,1);

%% Loop over the iterations
for iter = 1:niter

   % Compute Ad
   Ad(:,iter) = zeros( 2*nnodes, 1 );
   for i = 1:nsub
      floc = zeros( 2*size(newnode{i},1) + nbloq{i}, 1 );
      if i>1 % Lower boundary
         floc(bounsloc{2*i-1}) = d(newbouns{i-1},iter);
      end
      if i<nsub
         floc(bounsloc{2*i}) = d(newbouns{i},iter);
      end
      
      ui = K{i}\floc;
      u = zeros( 2*nnodes, 1 );
      if i>1 % Lower boundary
         u(newbouns{i-1}) = ui(bounsloc{2*i-1});
      end
      if i<nsub
         u(newbouns{i}) = ui(bounsloc{2*i});
      end
      Ad(:,iter) = Ad(:,iter) + u;
   end

   den = (d(:,iter)'*Ad(:,iter));
   d(:,iter) = d(:,iter)/sqrt(den); Ad(:,iter) = Ad(:,iter)/sqrt(den);
   num = Res(:,iter)'*d(:,iter);
   Lambda        = Lambda + d(:,iter)*num;
   Res(:,iter+1) = Res(:,iter) - Ad(:,iter)*num;
   
   resid(iter+1) = norm(Res(:,iter+1));
   
   % Orthogonalization
   d(:,iter+1) = Res(:,iter+1);
   for jter=1:iter
       betaij = ( Res(:,iter+1)'*Ad(:,jter) );
       d(:,iter+1) = d(:,iter+1) - d(:,jter) * betaij;
   end
   
end

figure;
plot(log10(resid));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the final solution and assemble it
usol = zeros( 2*nnodes, 1 );
for i = 1:nsub
   floc = zeros( 2*size(newnode{i},1) + nbloq{i}, 1 );
   if i>1 % Lower boundary
      floc(bounsloc{2*i-1}) = Lambda(newbouns{i-1});
   end
   if i<nsub
      floc(bounsloc{2*i}) = Lambda(newbouns{i});
   end
   
   ui = K{i}\floc;
   u = zeros( 2*nnodes, 1 );
   %u(map{i}) = ui;
   usol = usol + u;
end

plotGMSH({usol,'U_vect'}, elements, nodes(:,[1,2]), 'solution');