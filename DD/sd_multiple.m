% 01/12/2016 : Algo de Schur Dual avec plusieurs sous-domaines

clear all;
close all;

addpath(genpath('./tools'))

% Parameters
E       = 200000; % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 250;    % N.mm-2 : Loading on the plate
mat = [0, E, nu];

nsub = 4;

[ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/plate.msh' );
nnodes = size(nodes,1);

dirichlet = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ; 4,1,0 ; 4,2,0 ];
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
   newelem{i}   = elements1;
   newelem{i+1} = elements2;
   newnode{i}   = nodes1;
   newnode{i+1} = nodes2;
   
   % Increment the map
   oldmap = map{i};
   map{i} = map1(oldmap);
   map{i+1} = map2;
   
   % get the local boundary
   map1 = map{i}; map2 = map{i+1};
   newbound1   = map1(newbound);
   newbound2   = map2(newbound);
   newbouns{2*i-1} = newbound1; newbouns{2*i} = newbound2;
   
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

%% Begin the serious stuff

% Stiffness, rigid modes & Co
for i = 1:nsub
   nodess    = newnode{i};
   elementss = newelem{i};
   boundarys = boundar{i};
   
   if 2*i-2 > 0
      boun1 = newbouns{2*i-2};
   else
      boun1 = [];
   end
   
   if 2*i-1 < 2*nsub-1
      boun2 = newbouns{2*i-1};
   else
      boun2 = [];
   end
   
   nno{i} = size(nodess,1);
   [K{i},C,nbloq,node2c,c2node] =...
                 Krig2 (nodess,elementss,mat,order,boundarys,dirichlet);
   f{i}  = loading(nbloq,nodess,boundarys,neumann);
   
   % Management of the rigid modes (TODO)
%   if nbloq == 0 % nbloq == 1 is also problematic, but ...
%      nnodess = nno{i};
%      r1 = zeros(2*nnodess,1); r2 = r1; r3p = r1;
%      ind = 2:2:2*nnodess;
%      r1(ind-1,1) = 1; r2(ind,1) = 1;
%      
%      g1 = r1(boun1);
%      g2 = keepField( r2, 3, boundary );
%      r3p(ind-1,1) = -nodes(ind/2,2);
%      r3p(ind,1) = nodes(ind/2,1);
%      g3p = keepField( r3p, 3, boundary );
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
%      K2d = [ K2d, R ; R', zeros(size(R,2)) ];
%      nbloq2d = size(R,2);
%   end
   
end

% global indices of the boundary nodes
for i = 1:2*nsub-2
   nodess    = newnode{i};
   [ b1to2, b2to1 ] =...
             superNodes( nodess(newbouns{floor(i/2)+1},:), nodes, 1e-6 );
   globbo{i} = b1to2;
end

%% CG algorithm
Lambda = zeros(2*nnodes,1);       % Iterate (on the global mesh for convenience)
Res    = zeros(2*nnodes,niter+1);
Zed    = zeros(2*nnodes,niter+1);
d      = zeros(2*nnodes,niter+1);