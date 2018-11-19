% 02/11/2018
% ID fissure qcq par RG Petrov-Galerkin

tic
close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E           = 210000; % MPa : Young modulus
nu          = 0.3;    % Poisson ratio
fscalar     = 250;    % N.mm-1 : Loading on the plate
mat         = [0, E, nu];
br          = .0;      % Noise level
mur         = 1e1;%2e3;    % Regularization parameter
regular     = 1;      % Use the derivative regularization matrix (0 : Id, 1 : derivative)
froreg      = 1;      % Frobenius preconditioner
Npg         = 2;      % Nb Gauss points
ordertest   = 20;     % Order of test fonctions
zerobound   = 1;      % Put the boundaries of the crack to 0
nuzawa      = 100;     % (nuzawa = 1 means no Uzawa)
kuzawa      = 0;%1e2;     % Parameters of the Uzawa algorithm (kuzawa = 0 means best parameter)
teskase     = 20;       % Choice of the test case 20 => 2 cracks
resmooth    = 0;      % Resmooth coefficient (for the partitioning)
myceil      = .75;     % Ceil for the longth of the crack
newtonNor   = 0;       % Use Newton to find the normal

nbDirichlet = [];
%nbDirichlet = [ 1,10 ; 2,11 ; 3,11 ; 4,11 ];
%nbDirichlet = [ 1,5 ; 2,5 ; 3,5 ; 4,5 ]; % Nb of displacement measure points on the boundaries (0=all, /!\ NEVER go above the nb of dofs)
%nbDirichlet = [ 1,6 ; 2,6 ; 3,6 ; 4,6 ];
%nbDirichlet = [ 1,11 ; 2,0 ; 3,11 ; 4,0 ];

% Boundary conditions
dirichlet  = [0,1,0 ; 0,2,0 ; 0,3,0];
neumann1   = [3,2,fscalar ; 1,2,-fscalar];
neumann1p  = [3,1,-5*fscalar,10*fscalar,0];
neumann2   = [2,1,fscalar ; 4,1,-fscalar];
neumann3   = [3,1,fscalar ; 3,2,fscalar ; 2,1,fscalar ; 2,2,fscalar ; ...
              1,1,-fscalar ; 1,2,-fscalar ; 4,1,-fscalar ; 4,2,-fscalar];
neumann4   = [3,1,-fscalar ; 3,2,fscalar ; 2,1,fscalar ; 2,2,-fscalar ; ...
              1,1,fscalar ; 1,2,-fscalar ; 4,1,-fscalar ; 4,2,fscalar];
%neumann4   = [3,1,-fscalar ; 2,2,-fscalar ; ...
%              1,1,fscalar ; 4,2,fscalar];

dirichlet0 = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ; 3,1,0 ; 3,2,0 ; 4,1,0 ; 4,2,0];
neumann0   = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ; 3,1,0 ; 3,2,0 ; 4,1,0 ; 4,2,0];              
%dirichlet0 = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ; 3,1,0 ; 3,2,0 ; 4,1,0 ; 4,2,0];
%neumann0   = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ; 4,1,0 ; 4,2,0];
%dirichlet0 = [ 1,1,0 ; 1,2,0 ; 3,1,0 ; 3,2,0 ; 4,1,0 ; 4,2,0];
%neumann0   = [ 1,1,0 ; 1,2,0 ; 3,1,0 ; 3,2,0 ; 4,1,0 ; 4,2,0];
%dirichlet0 = [ 1,1,0 ; 1,2,0 ; 3,1,0 ; 3,2,0];
%neumann0   = [ 1,1,0 ; 1,2,0 ; 3,1,0 ; 3,2,0];
%dirichlet0 = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ; 3,1,0 ; 3,2,0 ; 4,1,0 ; 4,2,0];
%neumann0   = [ 1,1,0 ; 1,2,0 ; 3,1,0 ; 3,2,0];
%dirichlet0 = [ 1,1,0 ; 1,2,0 ; 3,1,0 ; 3,2,0 ; 4,1,0 ; 4,2,0];
%neumann0   = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ; 3,1,0 ; 3,2,0 ; 4,1,0 ; 4,2,0];
%dirichlet0 = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ];
%neumann0   = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ; 3,1,0 ; 3,2,0 ; 4,1,0 ; 4,2,0];
%dirichlet0 = [ 1,1,0 ; 1,2,0 ];
%neumann0   = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ; 3,1,0 ; 3,2,0 ; 4,1,0 ; 4,2,0];

% First, import the mesh
if teskase == 3
   [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_c_squared3.msh' );
elseif teskase == 4
   [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_c_squared4.msh' );
elseif teskase == 2
   [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_c_squared2.msh' );
elseif teskase == 5
   [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_c_squared5.msh' );
elseif teskase == 20
   [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc_r_2.msh' );
end

nnodes = size(nodes,1);

% mapBounds
[ node2b1, b2node1 ]   = mapBound( 1, boundary, nnodes );
[ node2b2, b2node2 ]   = mapBound( 2, boundary, nnodes );
[ node2b3, b2node3 ]   = mapBound( 3, boundary, nnodes );
[ node2b4, b2node4 ]   = mapBound( 4, boundary, nnodes );
[ node2b5, b2node5 ]   = mapBound( 5, boundary, nnodes );
[ node2b6, b2node6 ]   = mapBound( 6, boundary, nnodes );
               
% Then, build the stiffness matrix :
[K,C,nbloq,node2c,c2node] = Krig2 (nodes,elements,mat,order,boundary,dirichlet);
Kinter = K( 1:2*nnodes, 1:2*nnodes );

Xmax = max(nodes(:,1)); Xmin = min(nodes(:,1)); Xmoy = (Xmax+Xmin)/2;
Ymax = max(nodes(:,2)); Ymin = min(nodes(:,2)); Ymoy = (Ymax+Ymin)/2;
Lx = Xmax-Xmin; Ly = Ymax-Ymin;

f1  = loading(nbloq,nodes,boundary,neumann1);
f2  = loading(nbloq,nodes,boundary,neumann2);
f3  = loading(nbloq,nodes,boundary,neumann3);
f4  = loading(nbloq,nodes,boundary,neumann4);

uin = K\[f1,f2,f3,f4];
u1 = uin(1:2*nnodes,1); u2 = uin(1:2*nnodes,2);
u3 = uin(1:2*nnodes,3); u4 = uin(1:2*nnodes,4);

%% Remove rigid modes (useless)
%R  = rigidModes(nodes);
%UU = [u1,u2,u3,u4] - R*((R'*R)\(R'*[u1,u2,u3,u4]));
%u1 = UU(:,1); u2 = UU(:,2); u3 = UU(:,3); u4 = UU(:,4);

f1 = Kinter*u1; f2 = Kinter*u2; f3 = Kinter*u3; f4 = Kinter*u4;

ui = reshape(u1,2,[])';  ux = ui(:,1);  uy = ui(:,2);

u1ref = u1; u2ref = u2; u3ref = u3; u4ref = u4;
f1ref = f1; f2ref = f2; f3ref = f3; f4ref = f4;

% Compute stress :
sigma = stress(u1,E,nu,nodes,elements,order,1,ntoelem);

% Output :
plotGMSH({u1,'U1';u2,'U2';u3,'U3';u4,'U4'}, elements, nodes, 'output/reference');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ nodes2,elements2,ntoelem2,boundary2,order] = readmesh( 'meshes/rg_refined/plate_nu_r.msh' );
nnodes2 = size(nodes2,1);
[K2,C2,nbloq2,node2c2,c2node2] = Krig2 (nodes2,elements2,mat,order,boundary2,dirichlet);
Kinter2 = K2( 1:2*nnodes2, 1:2*nnodes2 );

%% Mapbounds
%[ node2b12, b2node12 ] = mapBound( 1, boundary2, nnodes2 );
%[ node2b22, b2node22 ] = mapBound( 2, boundary2, nnodes2 );
%[ node2b32, b2node32 ] = mapBound( 3, boundary2, nnodes2 );
%[ node2b42, b2node42 ] = mapBound( 4, boundary2, nnodes2 );
%
%indexbound  = [2*b2node1-1 ; 2*b2node1 ; 2*b2node2-1 ; 2*b2node2 ;...
%               2*b2node3-1 ; 2*b2node3 ; 2*b2node4-1 ; 2*b2node4];
%
% % Sort the nodes (in case order > 1)
%[~,order5]  = sort( nodes( b2node5, 1 ) );
%[~,order6]  = sort( nodes( b2node6, 1 ) );
%
%icrack5x    = [ 2*b2node5(order5)-1 ];
%icrack6x    = [ 2*b2node6(order6)-1 ];
%icrack5y    = [ 2*b2node5(order5) ];
%icrack6y    = [ 2*b2node6(order6) ];
%
%indexbound2 = [2*b2node12-1 ; 2*b2node12 ; 2*b2node22-1 ; 2*b2node22 ;...
%               2*b2node32-1 ; 2*b2node32 ; 2*b2node42-1 ; 2*b2node42];
%% Pass f and u on the uncracked mesh
%ur1 = zeros( 2*nnodes2, 1 ); fr1 = zeros( 2*nnodes2, 1 );
%ur1(indexbound2) = u1(indexbound);
%fr1(indexbound2) = f1(indexbound);
%% Same for the second one
%ur2 = zeros( 2*nnodes2, 1 ); fr2 = zeros( 2*nnodes2, 1 );
%ur2(indexbound2) = u2(indexbound);
%fr2(indexbound2) = f2(indexbound);
%% 3
%ur3 = zeros( 2*nnodes2, 1 ); fr3 = zeros( 2*nnodes2, 1 );
%ur3(indexbound2) = u3(indexbound);
%fr3(indexbound2) = f3(indexbound);
%% 4
%ur4 = zeros( 2*nnodes2, 1 ); fr4 = zeros( 2*nnodes2, 1 );
%ur4(indexbound2) = u4(indexbound);
%fr4(indexbound2) = f4(indexbound);

%% Pass meshes
UF = [u1,u2,u3,u4,f1,f2,f3,f4];
UR = passMesh2D(nodes, elements, nodes2, elements2, UF, 0);

ur1 = UR(:,1); fr1 = UR(:,5)*25; % 25 is the ratio between mesh sizes (hackissimo)
ur2 = UR(:,2); fr2 = UR(:,6)*25;
ur3 = UR(:,3); fr3 = UR(:,7)*25;
ur4 = UR(:,4); fr4 = UR(:,8)*25;

% Add the noise
u1n = ur1; u2n = ur2; u3n = ur3; u4n = ur4;
am1 = sqrt(mean(ur1.^2)); am2 = sqrt(mean(ur2.^2));
am3 = sqrt(mean(ur3.^2)); am4 = sqrt(mean(ur4.^2));
br1 = randn(2*nnodes2,1); br2 = randn(2*nnodes2,1);
br3 = randn(2*nnodes2,1); br4 = randn(2*nnodes2,1);
%noise = load('noises/cauchyRG_dirichlet.mat');
%br1 = noise.br1; br2 = noise.br2; br3 = noise.br3; br4 = noise.br4;
%ur1 = ( 1 + br*br1 ) .* ur1; ur2 = ( 1 + br*br2 ) .* ur2;
%ur3 = ( 1 + br*br3 ) .* ur3; ur4 = ( 1 + br*br4 ) .* ur4;
ur1 = ur1 + am1*br*br1; ur2 = ur2 + am2*br*br2;
ur3 = ur3 + am3*br*br3; ur4 = ur4 + am4*br*br4;

plotGMSH({ur1,'U1';ur2,'U2';ur3,'U3';ur4,'U4'},elements2, nodes2, 'output/bound');

%% Test the discrete measure point stuff
if min(size(nbDirichlet))>0
   nmin = [];
   for i=1:4
      nbdiscrete = find(nbDirichlet(:,1)==i);
      if min(size(nbdiscrete)) > 0
         nbi = nbDirichlet(nbdiscrete,2);
         if nbi>0
            [~, b2node] = mapBound(i, boundary2, nnodes2);
            xmin = min(nodes2(b2node,1)); xmax = max(nodes2(b2node,1)); Lx1 = xmax-xmin; sx = Lx1/(nbi-1);
            ymin = min(nodes2(b2node,2)); ymax = max(nodes2(b2node,2)); Ly1 = ymax-ymin; sy = Ly1/(nbi-1);
            if sx==0 xy = [ xmin*ones(1,nbi) ; ymin:sy:ymax ];
            elseif sy==0 xy = [ xmin:sx:xmax ; ymin*ones(1,nbi) ];
            else xy = [ xmin:sx:xmax ; ymin:sy:ymax ];
            end
         
            nm = zeros(nbi,1);
            for j=1:nbi % Find the closest nodes
               xyj = xy(:,j)';
               nodiff = nodes2(b2node,:)-xyj;
               normdiff = nodiff(:,1).^2+nodiff(:,2).^2;
               [~,nm(j)] = min(normdiff);
            end
            nm = unique(nm); nmin = [nmin;b2node(nm)];
         end
      end
   end
   nmin = unique(nmin);
else
   nmin = 1:nnodes2;
end

%% Build the maps of the Dirichlet and Neumann bounds
%b2nodesD = cell(4,2); % nb of boundaries * nb of dimensions
%b2nodesN = cell(4,2);
b2nodesD = []; b2nodesN = []; b2nodesTOT = [];
for i=1:4
   [~, b2node] = mapBound(i, boundary2, nnodes2);
   dofD = dirichlet0(find(dirichlet0(:,1)==i),2); % dof of the given dirichlet
   dofN = neumann0(find(neumann0(:,1)==i),2); % dof of the given neumann
   b2nodesDo = []; % Stores the added nodes
   for j=1:size(dofD,1)
      b2nodesD = [b2nodesD ; 2*b2node-2+dofD(j)];
   end
   b2nodesTOT = [b2nodesTOT ; 2*b2node-1 ; 2*b2node];
%   if max(size(dofN)) > 0
   for j=1:size(dofN,1)
      b2nodesN = [b2nodesN ; 2*b2node-2+dofN(j)];
   end
end

% Discrete boundary stuff
b2nodesD = intersect(b2nodesD,[2*nmin-1;2*nmin]);

b2nodesnoD = setdiff(b2nodesTOT,b2nodesD);  % Missing Dirichlet conditions
b2nodesnoN = setdiff(b2nodesTOT,b2nodesN);  % Missing Neumann conditions (unused except for visualization)
% Put u to 0 at those missing BCs
% u1(b2nodesnoD) = 0; u2(b2nodesnoD) = 0; u3(b2nodesnoD) = 0; u4(b2nodesnoD) = 0;

%% Preliminary stuff : find the volumic elements corresponding to the boundaries
nboun2 = size(boundary2,1); nelem2 = size(elements2,1);
boun2vol2 = zeros( nboun2, 1 ); extnorm2 = zeros( nboun2, 2 );
for i=1:nboun2
   % Volumic element
   no1 = boundary2(i,2); no2 = boundary2(i,3); % only with 2 nodes even if order > 1
   cand1 = rem( find(elements2==no1)-1,nelem2 )+1; % find gives line + column*size
   cand2 = rem( find(elements2==no2)-1,nelem2 )+1;
   boun2vol2(i) = intersect(cand1, cand2); % If everything went well, there is only one
   
   % Exterior normal
   elt = boun2vol2(i); no3 = setdiff( elements2( elt, 1:3 ), [no1,no2]);
   x1 = nodes2(no1,1); y1 = nodes2(no1,2);
   x2 = nodes2(no2,1); y2 = nodes2(no2,2);
   x3 = nodes2(no3,1); y3 = nodes2(no3,2);
   extnorm2(i,:) = [ y1-y2 , -(x1-x2) ];
   extnorm2(i,:) = extnorm2(i,:)/norm(extnorm2(i,:));
   if (x3-x2)*(y1-y2) - (y3-y2)*(x1-x2) > 0 % Check that the normal is exterior
      extnorm2(i,:) = -extnorm2(i,:);
   end
end

nboun1 = size(boundary2,1); nelem1 = size(elements2,1);
boun2vol1 = zeros( nboun1, 1 ); extnorm1 = zeros( nboun1, 2 );
sigr1 = zeros( nboun1, 3*order ); sigr2 = zeros( nboun1, 3*order );
urr1  = zeros( nboun1, 2+2*order ); urr2 = zeros( nboun1, 2+2*order );
urr3  = zeros( nboun1, 2+2*order ); urr4 = zeros( nboun1, 2+2*order );
for i=1:nboun1 % TODO : rationalize the stuff (don't do 2 times the same work)
   % Volumic element
   no1 = boundary2(i,2); no2 = boundary2(i,3); % only with 2 nodes even if order > 1
   cand1 = rem( find(elements2==no1)-1,nelem2 )+1; % find gives line + column*size
   cand2 = rem( find(elements2==no2)-1,nelem2 )+1;
   boun2vol1(i) = intersect(cand1, cand2); % If everything went well, there is only one
   
   % Exterior normal
   elt = boun2vol1(i); no3 = setdiff( elements2( elt, 1:3 ), [no1,no2]);
   x1 = nodes2(no1,1); y1 = nodes2(no1,2);
   x2 = nodes2(no2,1); y2 = nodes2(no2,2);
   x3 = nodes2(no3,1); y3 = nodes2(no3,2);
   extnorm1(i,:) = [ y1-y2 , -(x1-x2) ];
   extnorm1(i,:) = extnorm1(i,:)/norm(extnorm1(i,:));
   if (x3-x2)*(y1-y2) - (y3-y2)*(x1-x2) > 0 % Check that the normal is exterior
      extnorm1(i,:) = -extnorm1(i,:);
   end
   
   % ur
   urr1(i,1:4) = u1( [2*no1-1,2*no1,2*no2-1,2*no2] );
   urr2(i,1:4) = u2( [2*no1-1,2*no1,2*no2-1,2*no2] );
   urr3(i,1:4) = u3( [2*no1-1,2*no1,2*no2-1,2*no2] );
   urr4(i,1:4) = u4( [2*no1-1,2*no1,2*no2-1,2*no2] );
   if order == 2
      no4 = boundary(i,4);
         urr1(i,5:6) = u1( [2*no4-1,2*no4] );
         urr2(i,5:6) = u2( [2*no4-1,2*no4] );
         urr3(i,5:6) = u3( [2*no4-1,2*no4] );
         urr4(i,5:6) = u4( [2*no4-1,2*no4] );
   end
end

% Give a dof to each node...
nodes2dof = zeros(2*nnodes2,1); dof2nodes = []; szD = 0;
for i=1:2*nnodes2
   if min(size(find(b2nodesD==i)))==0
      nodes2dof(i) = szD+1;
      dof2nodes(szD+1) = i;
      szD = szD+1;
   end
end

% ...and then give a dof to each element
boun2dof = zeros(nboun2,1);
for i=1:nboun2
   boun2dof(i,1) = nodes2dof(2*boundary2(i,2)-1);
   boun2dof(i,2) = nodes2dof(2*boundary2(i,2));
   boun2dof(i,3) = nodes2dof(2*boundary2(i,3)-1);
   boun2dof(i,4) = nodes2dof(2*boundary2(i,3));
end

%% Order the dofs : force and displacement on the boundary
nodesbound2 = unique(boundary2(:)); % Nodes on the boundary
nnbound2    = size(nodesbound2,1);
boun2doftot = zeros(nboun2,2);

for i=1:nboun2 % Stores the nodes associated to each element of bound2 (in boundary numbering)
   boun2doftot(i,1) = find(nodesbound2==boundary2(i,2)); % Rem : boundary2(i,1) is the no of boundary
   boun2doftot(i,2) = find(nodesbound2==boundary2(i,3));
end

knownD = [];
knownN = [];
for i=1:nboun2 % Chooses the known displacement dofs
   index = boundary2(i,1);
   isimposed = find( dirichlet0(:,1) == index );
   dofs = dirichlet0(isimposed,2);
   knownD = [ knownD ; 2*boun2doftot(i,1)-2+dofs ];
   knownD = [ knownD ; 2*boun2doftot(i,2)-2+dofs ];
   
   isimposed = find( neumann0(:,1) == index );
   dofs = neumann0(isimposed,2);
   knownN = [ knownN ; 2*i-2+dofs ];
end

knownD = unique(knownD); knownN = unique(knownN); % Remove redondnacy
knownD = intersect( knownD, [2*nmin-1;2*nmin] ); % Discrete measurement    /!\ There is a global/local confusion /!\

tofindD = setdiff( 1:2*nnbound2 , knownD );
tofindN = setdiff( 1:2*nnbound2 , knownN );
%%

[ nodesu,elementsu,ntoelemu,boundaryu,orderu] = readmesh( 'meshes/rg_refined/plate_nu_r.msh' );
nnodesu = size(nodesu,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test-functions shit
% Build the polynomial test functions.
if ordertest == 10
   load('conditions10_2d.mat','-ascii');
   M = spconvert(conditions10_2d); clear('conditions10_2d');
   nmax = 10;
elseif ordertest == 20
   load('conditions20_2d.mat','-ascii');
   M = spconvert(conditions20_2d); clear('conditions20_2d');
   nmax = 20;
end
ncoef =2*(nmax+1)^2; neq = ncoef;

% Suppress the superior terms
for i=0:nmax
   for j=0:nmax
      if i+j>nmax
         M((nmax+1)*i+j+1,:) = 0;
         M((nmax+1)*i+j+1,(nmax+1)*i+j+1) = 1;
      end
   end
end

% Suppress zeros rows
toremove = [];
for k=1:size(M,1)
  if norm(M(k,:)) == 0
     toremove(end+1) = k;
  end
end
M(toremove,:) = [];

coef   = null(full(M));
nftest = size(coef,2); % Nb of avaliable test functions

%% Build the (Petrov-Galerkin) Left Hand Side
nelemu = size(elementsu,1); nn = size(elementsu,2); %nn=3

% Build the list of segment elements
nseg = nn*nelemu; segel = zeros(nseg,2); nbu = size(boundaryu,1);
for j=1:nelemu
   segel(nn*(j-1)+1,1) = elementsu(j,2);
   segel(nn*(j-1)+1,2) = elementsu(j,3);
   segel(nn*(j-1)+2,1) = elementsu(j,1);
   segel(nn*(j-1)+2,2) = elementsu(j,3);
   segel(nn*(j-1)+3,1) = elementsu(j,2);
   segel(nn*(j-1)+3,2) = elementsu(j,1);
end
j = 1;
while j <= nseg % Remove redundancy (I don't use unique because of inverted values)
   lis1 = find( segel(1:j-1,:) == segel(j,1) );
   lis1( find(lis1>j-1) ) = lis1( find(lis1>j-1) ) - j + 1;
   lis2 = find( segel(1:j-1,:) == segel(j,2) );
   lis2( find(lis2>j-1) ) = lis2( find(lis2>j-1) ) - j + 1;
   
   % also remove boundary elements
   lis3 = find( boundaryu(:,2:3) == segel(j,1) );
   lis3( find(lis3>nbu) ) = lis3( find(lis3>nbu) ) - nbu;
   lis4 = find( boundaryu(:,2:3) == segel(j,2) );
   lis4( find(lis4>nbu) ) = lis4( find(lis4>nbu) ) - nbu;
   
   if min(size( intersect(lis1,lis2) )) > 0 || min(size( intersect(lis3,lis4) )) > 0% || (size(lis3,1) > 0 || size(lis4,1) > 0)%
      segel(j,:) = [];
      j = j-1;
      nseg = nseg-1;
   end
   j = j+1;
end

%%
% Transform the given data in the proper format
f_known1 = zeros(2*nboun2,1); u_known1 = zeros(2*nnbound2,1);
f_known2 = zeros(2*nboun2,1); u_known2 = zeros(2*nnbound2,1);
f_known3 = zeros(2*nboun2,1); u_known3 = zeros(2*nnbound2,1);
f_known4 = zeros(2*nboun2,1); u_known4 = zeros(2*nnbound2,1);
   
for i=1:nboun2
   bonod = boundary2(i,:); exno = extnorm2(i,:)';

   if exno(1) == 1 % Bound 2
      fer1 = [0;0]; fer2 = [fscalar;0];
      fer3 = [fscalar;fscalar]; fer4 = [fscalar;-fscalar];
   elseif exno(1) == -1 % Bound 4
      fer1 = [0;0]; fer2 = -[fscalar;0];
      fer3 = [-fscalar;-fscalar]; fer4 = [-fscalar;fscalar];
   elseif exno(2) == 1 % Bound 3
      fer1 = [0;fscalar]; fer2 = [0;0];
      fer3 = [fscalar;fscalar]; fer4 = [-fscalar;fscalar];
   elseif exno(2) == -1 % Bound 1
      fer1 = -[0;fscalar]; fer2 = [0;0];
      fer3 = [-fscalar;-fscalar]; fer4 = [fscalar;-fscalar];
   end
     
   indicesLoc = [2*boun2doftot(i,:)-1,2*boun2doftot(i,:)]; % Local Displacement dofs
   indicesGlo = [2*boundary2(i,[2,3])-1,2*boundary2(i,[2,3])]; % Global Displacement dofs
      
   u_known1(indicesLoc) = ur1(indicesGlo);
   u_known2(indicesLoc) = ur2(indicesGlo);
   u_known3(indicesLoc) = ur3(indicesGlo);
   u_known4(indicesLoc) = ur4(indicesGlo);
      
   f_known1([2*i-1,2*i]) = fer1;
   f_known2([2*i-1,2*i]) = fer2;
   f_known3([2*i-1,2*i]) = fer3;
   f_known4([2*i-1,2*i]) = fer4;
end

disp([ 'Direct problem solved and data management ', num2str(toc) ]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the operators
tic
Ng = Npg;
[ Xg, Wg ] = gaussPt1d( Ng );
Ndots = size(Wg,1); % (= Ng by the way)

%% Build the list of Gauss points, and of construction-functions
Xxg = zeros( nboun2*Ndots,1 ); Yyg = zeros( nboun2*Ndots,1 );
Wwg = zeros( nboun2*Ndots,1 );
Phi = sparse( 2*nboun2*Ndots, 2*nnbound2 );
Psi = sparse( 2*nboun2*Ndots, 2*nboun2 );
exnor = zeros( nboun2*Ndots,2); % Normal

for i=1:nboun2
   bonod = boundary2(i,:); exno = extnorm2(i,:)';
   
   no1 = bonod(2); no2 = bonod(3); boname = bonod(1);
   x1 = nodes2(no1,1); y1 = nodes2(no1,2);
   x2 = nodes2(no2,1); y2 = nodes2(no2,2);
      
   len = sqrt( (x1-x2)^2 + (y1-y2)^2 );

   indDtot = [2*boun2doftot(i,1)-1,2*boun2doftot(i,1),...
              2*boun2doftot(i,2)-1,2*boun2doftot(i,2)]; % U dofs the element is associated to
   indNtot = [2*i-1,2*i]; % F dofs the element is associated to
                
   for j=1:Ndots
      xg = Xg(j,:); wg = Wg(j);
      xgr  = ( (1-xg)*[x1;y1] + xg*[x2;y2] ); % abscissae

      X = xgr(1)/Lx; Y = xgr(2)/Lx; % Normalize
      Swg   = len * wg;

      indexg = Ndots*(i-1) + j;
      Xxg( indexg ) = X; Yyg( indexg ) = Y; Wwg( indexg ) = Swg;
      exnor( indexg, :) = exno';

      Phi( 2*indexg-1, indDtot)  = [ 1-xg, 0, xg, 0 ];
      Phi( 2*indexg  , indDtot)  = [ 0, 1-xg, 0, xg ];
      Psi( 2*indexg-1, indNtot)  = [1,0];
      Psi( 2*indexg  , indNtot)  = [0,1];
   end
end

%% Build the matrix of test-functions
Vv = zeros( size(coef,1), 2*nboun2*Ndots );
Sv = zeros( size(coef,1), 2*nboun2*Ndots );
exnoX = exnor(:,1); exnoY = exnor(:,2);

for ii=0:nmax
   nd1 = nmax-ii;
   for jj=0:nd1
                  
         sloc11a = zeros(nboun2*Ndots,1); sloc12a = zeros(nboun2*Ndots,1);
         sloc22a = zeros(nboun2*Ndots,1); sloc11b = zeros(nboun2*Ndots,1);
         sloc12b = zeros(nboun2*Ndots,1); sloc22b = zeros(nboun2*Ndots,1);
         if ii>0
            XY = ii.*Xxg.^(ii-1).*Yyg.^jj;
            sloc11a = E/(1-nu^2)*XY;
            sloc22a = nu*E/(1-nu^2)*XY;
            sloc12b = E/(2*(1+nu))*XY;
         end
         if jj>0
            XY = jj.*Xxg.^ii.*Yyg.^(jj-1);
            sloc11b = nu*E/(1-nu^2)*XY;
            sloc22b = E/(1-nu^2)*XY;
            sloc12a = E/(2*(1+nu))*XY;
         end

         fpaax = sloc11a.*exnoX + sloc12a.*exnoY;
         fpaay = sloc12a.*exnoX + sloc22a.*exnoY;

         fpabx = sloc11b.*exnoX + sloc12b.*exnoY;
         fpaby = sloc12b.*exnoX + sloc22b.*exnoY;

         XY = Xxg.^ii.*Yyg.^jj;
         vpaax = XY; vpaby = XY;

         index = (nmax+1)*ii + jj + 1;

         Sv( index, 1:2:2*nboun2*Ndots-1 )              = Wwg' .* fpaax';
         Sv( (nmax+1)^2 + index, 1:2:2*nboun2*Ndots-1 ) = Wwg' .* fpabx';

         Sv( index, 2:2:2*nboun2*Ndots )                = Wwg' .* fpaay';
         Sv( (nmax+1)^2 + index, 2:2:2*nboun2*Ndots )   = Wwg' .* fpaby';

         Vv( index, 1:2:2*nboun2*Ndots-1 )              = Wwg' .* vpaax';
         Vv( (nmax+1)^2 + index, 2:2:2*nboun2*Ndots )   = Wwg' .* vpaby';
   end
end

Rupij = Sv*Phi; clear Sv; clear Phi;
Rfpij = Vv*Psi; clear Vv; clear Psi;
Ru = coef'*Rupij; Rf = coef'*Rfpij;

%% Restrict the data on the given part (not necessary, but cleaner)
f_known1(tofindN) = 0; f_known2(tofindN) = 0;
f_known3(tofindN) = 0; f_known4(tofindN) = 0;
u_known1(tofindD) = 0; u_known2(tofindD) = 0;
u_known3(tofindD) = 0; u_known4(tofindD) = 0;

ndofs = size(f_known1(tofindN),1) + size(u_known1(tofindD),1) + 3*nnodes2;

%% Restrict the operators
Rfm = Rf(:,tofindN); Rfr = Rf(:,knownN);
Rum = Ru(:,tofindD); Rur = Ru(:,knownD);

LhsA = -Rum;
LhsB = Rfm;
Rhs1 = Rur*u_known1(knownD) - Rfr*f_known1(knownN);
Rhs2 = Rur*u_known2(knownD) - Rfr*f_known2(knownN);
Rhs3 = Rur*u_known3(knownD) - Rfr*f_known3(knownN);
Rhs4 = Rur*u_known4(knownD) - Rfr*f_known4(knownN);
   
% Build the matrix that passes f on the nodes from the bound
Fntob = zeros( 2*nnbound2, 2*nboun2 );
nodeMass = zeros(2*nnbound2); elemMass = zeros(2*nboun2);
for i=1:nboun2
   coef1 = [ boun2doftot(i,1), boun2doftot(i,2) ];
   len = norm(nodes2(boundary2(i,2),:) - nodes2(boundary2(i,3),:));
   Fntob( 2*i-1, 2*coef1-1 ) = len/2 * [1,1];
   Fntob( 2*i, 2*coef1 )     = len/2 * [1,1];

   ico = [ 2*i-1, 2*i ];
   cco = [ 2*coef1-1, 2*coef1 ];
   elemMass(ico,ico) = len*eye(2);
   nodeMass(cco,cco) = nodeMass(cco,cco) + len/2 * [ eye(2), zeros(2) ; ...
                                                     zeros(2), eye(2), ];
end
Mfm = elemMass(tofindN,tofindN); Mfr = elemMass(knownN,knownN);
Mum = nodeMass(tofindD,tofindD); Mur = nodeMass(knownD,knownD);

% If needed, use the differential regularization matrix
if regular == 1
   Du = zeros(2*nnbound2);
   Df = eye(2*nboun2); % No derivative for f
   for i=1:nboun2
      coefU = [ boun2doftot(i,1) , boun2doftot(i,2) ];
      len = norm(nodes2(boundary2(i,2),:) - nodes2(boundary2(i,3),:));
      Du( 2*coefU-1, 2*coefU-1 ) = Du( 2*coefU-1, 2*coefU-1 ) + 1/len*[1,-1;-1,1];
      Du( 2*coefU, 2*coefU )     = Du( 2*coefU, 2*coefU ) + 1/len*[1,-1;-1,1];
   end

   Duu  = Du(tofindD,tofindD); Dfu  = Df(tofindN,tofindN);
   Duk  = Du(knownD,knownD);   Dfk  = Df(knownN,knownN);
   Duuk = Du(tofindD,knownD);  Dfuk = Df(tofindN,knownN);
end
  
Rucij = zeros(size(coef,1),2*nnodes2);
   
if order>1 warning('Mesh order is too high') end
Ng = Npg;
[ Xg, Wg ] = gaussPt( Ng ); Ndots = size(Wg,1);

%%%% Build the operator
%% Build the list of Gauss points, and of construction-functions
Xxg = zeros( nelem2*Ndots,1 ); Yyg = zeros( nelem2*Ndots,1 );
Wwg = zeros( nelem2*Ndots,1 );
Phi = sparse( 3*nelem2*Ndots, 3*nnodes2 );
     
for i=1:nelem2
   bonod = elements2(i,:);
        
   no1 = bonod(1); no2 = bonod(2); no3 = bonod(3);
   x1 = nodes2(no1,1); y1 = nodes2(no1,2);
   x2 = nodes2(no2,1); y2 = nodes2(no2,2);
   x3 = nodes2(no3,1); y3 = nodes2(no3,2);
           
   nsurf = (x1-x2)*(y1-y3)-(x1-x3)*(y1-y2);  %
   surfa = .5*abs(nsurf);
          
   % Dofs in the numerotation of the boundary nodes
   N = elements2(i,:);
   indDtot = [ 3*N(1)-2, 3*N(1)-1, 3*N(1),...
               3*N(2)-2, 3*N(2)-1, 3*N(2),...
               3*N(3)-2, 3*N(3)-1, 3*N(3) ];
                       
   % Dofs in the numerotation of the boundary elements
   indNtot = [ 3*i-2, 3*i-1, 3*i ];
                     
   for j=1:Ndots
      xg = Xg(j,:); wg = Wg(j); umxgx = 1-xg(1)-xg(2);
      xgr  = umxgx*[x1;y1] + xg(1)*[x2;y2] + xg(2)*[x3;y3] ; % abscissae
     
      X = xgr(1)/Lx; Y = xgr(2)/Lx; % Normalize
      %Swg   = surfa * wg;
     
      indexg = Ndots*(i-1) + j;
      Xxg( indexg ) = X; Yyg( indexg ) = Y; Wwg( indexg ) = surfa * wg;
     
      Phi( 3*indexg-2, indDtot)  = [ umxgx, 0, 0, xg(1), 0, 0, xg(2), 0, 0 ];
      Phi( 3*indexg-1, indDtot)  = [ 0, umxgx, 0, 0, xg(1), 0, 0, xg(2), 0 ];
      Phi( 3*indexg  , indDtot)  = [ 0, 0, umxgx, 0, 0, xg(1), 0, 0, xg(2) ];
   end
end
     
%% Build the matrix of test-functions
Sv = zeros( size(coef,1), 3*nelem2*Ndots );
for ii=0:nmax
   nd1 = nmax-ii;
   for jj=0:nd1
                  
      sloc11a = zeros(nelem2*Ndots,1); sloc12a = zeros(nelem2*Ndots,1);
      sloc22a = zeros(nelem2*Ndots,1); sloc11b = zeros(nelem2*Ndots,1);
      sloc12b = zeros(nelem2*Ndots,1); sloc22b = zeros(nelem2*Ndots,1);
      if ii>0
         XY = ii.*Xxg.^(ii-1).*Yyg.^jj;
         sloc11a = E/(1-nu^2)*XY;
         sloc22a = nu*E/(1-nu^2)*XY;
         sloc12b = E/(2*(1+nu))*XY;
      end
      if jj>0
         XY = jj.*Xxg.^ii.*Yyg.^(jj-1);
         sloc11b = nu*E/(1-nu^2)*XY;
         sloc22b = E/(1-nu^2)*XY;
         sloc12a = E/(2*(1+nu))*XY;
      end
     
      index = (nmax+1)*ii + jj + 1;
                       
      Sv( index, 1:3:3*nelem2*Ndots-2 )              = Wwg' .* sloc11a';
      Sv( index, 2:3:3*nelem2*Ndots-1 )              = Wwg' .* sloc22a';
      Sv( index, 3:3:3*nelem2*Ndots )                = Wwg' .* sloc12a';

      Sv( (nmax+1)^2 + index, 1:3:3*nelem2*Ndots-2 ) = Wwg' .* sloc11b';
      Sv( (nmax+1)^2 + index, 2:3:3*nelem2*Ndots-1 ) = Wwg' .* sloc22b';
      Sv( (nmax+1)^2 + index, 3:3:3*nelem2*Ndots )   = Wwg' .* sloc12b';
   end
end
     
Rucij = Sv*Phi; clear Sv; clear Phi;
Ruc = coef'*Rucij; clear Rucij; % don't clear coef as we reuse it

disp([ 'Building of the operators ', num2str(toc) ]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%      nodeMass3 = zeros(3*nelem2+2);
%      for i=1:nelem2
%         coef1 = [ i, i+1 ];
%         len = norm(nodes2(elements2(i,2),:) - nodes2(elements2(i,3),:));

%         ico = [ 2*i-1, 2*i ];
%         cco = [ 2*coef1-1, 2*coef1 ];
%         nodeMass3(cco,cco) = nodeMass3(cco,cco) + len/2 * [ eye(2), zeros(2) ; ...
%                                                          zeros(2), eye(2), ];
%      end

tic;

kB = norm(Ru,'fro')/norm(Rf,'fro');
Lhs = [LhsA,kB*LhsB,Ruc];

Rhs = Rhs1;
A = Lhs'*Lhs; b = Lhs'*Rhs; sA = size(A,1);
   
%Z12    = zeros(size(Mum,1),size(Mfm,2)); Z13 = zeros(size(Mum,1),3*nnodes2);
%Z23    = zeros(size(Mfm,1),3*nnodes2);
%      Mtot   = [ Mum, Z12, Z13  ; Z12', Mfm, Z23 ; Z13', Z23', nodeMass3 ];
%      Msmall = [ 0*Mum, Z12, Z13  ; Z12', 0*Mfm, Z23 ; Z13', Z23', nodeMass3 ]; % Only on the crack
%      Mtop   = Mtot-Msmall; % Debug asset

%% And the derivative operator
D3  = sparse(6*nelem2,3*nnodes2);
D3m = sparse(2*nelem2,nnodes2);
for i=1:nelem2
   no1 = elements2(i,1); no2 = elements2(i,2); no3 = elements2(i,3);
   x1 = nodes2(no1,1); y1 = nodes2(no1,2);
   x2 = nodes2(no2,1); y2 = nodes2(no2,2);
   x3 = nodes2(no3,1); y3 = nodes2(no3,2);
      
   n = (x1-x2)*(y1-y3)-(x1-x3)*(y1-y2);
   S = .5*abs(n);

   xg = (x1+x2+x3)/3; yg = (y1+y2+y3)/3; % Gauss point
            
   Be = 1/(2*S)*[ x3-x2, 0, 0, x1-x3, 0, 0, x2-x1, 0, 0 ;...
                  y2-y3, 0, 0, y3-y1, 0, 0, y1-y2, 0, 0 ;...
                  0, x3-x2, 0, 0, x1-x3, 0, 0, x2-x1, 0 ;...
                  0, y2-y3, 0, 0, y3-y1, 0, 0, y1-y2, 0 ;...
                  0, 0, x3-x2, 0, 0, x1-x3, 0, 0, x2-x1 ;...
                  0, 0, y2-y3, 0, 0, y3-y1, 0, 0, y1-y2 ];

   Bem = 1/(2*S)*[ x3-x2, x1-x3, x2-x1 ;...
                   y2-y3, y3-y1, y1-y2 ];
      
   indelem = [ 6*i-5, 6*i-4, 6*i-3, 6*i-2, 6*i-1, 6*i ];
   coefU   = [ 3*no1-2, 3*no1-1, 3*no1, 3*no2-2, 3*no2-1, 3*no2, 3*no3-2, 3*no3-1, 3*no3 ];

   indelemm = [ 2*i-1, 2*i ];
   coefUm   = [ no1, no2, no3, ];
                        
   D3( indelem, coefU )    = D3( indelem, coefU )    + Be*sqrt(S);
   D3m( indelemm, coefUm ) = D3m( indelemm, coefUm ) + Bem*sqrt(S);
end

if regular == 1
   D3u = D3'*D3;
   Z12 = zeros( size(Duu,1) , size(Dfu,2) ); Z13 = zeros( size(Duu,1), size(D3,2) );
   Z23 = zeros( size(Dfu,1), size(D3,2) );
   Dtot = [ Duu ,Z12, Z13 ; Z12', Dfu, Z23 ; Z13', Z23', D3u ];
   L = Dtot;

   L12 = [ Du(tofindD,knownD) ; zeros(size(Dfu,2),size(Duk,2)) ; zeros(3*nnodes2,size(Duk,2)) ];
   L121 = L12*u_known1(knownD); L122 = L12*u_known2(knownD);
   L123 = L12*u_known3(knownD); L124 = L12*u_known4(knownD);
else
   L = Mtot;
   %L = eye(size(A))*nodeMass(1,1);
   L121 = zeros(ndofs,1); L122 = zeros(ndofs,1); L123 = zeros(ndofs,1); L124 = zeros(ndofs,1);
end

%sL = real(L^(1/2));   % Square root of L (as usual, real should not be) /!\ takes too much time

MAT = Lhs'*Lhs + mur*L;
VE1 = Lhs'*Rhs1-mur*L121; VE2 = Lhs'*Rhs2-mur*L122;
VE3 = Lhs'*Rhs3-mur*L123; VE4 = Lhs'*Rhs4-mur*L124;

      % Build the contact inequation condition
C = zeros(nnodes2, 3*nnodes2);
for i=1:nelem2 %ux*nx + uy*ny > 0
   no1 = elements2(i,1); no2 = elements2(i,2); no3 = elements2(i,3);
   C(no1,3*no1-2) = C(no1,3*no1-2) + 1;
   C(no2,3*no2-2) = C(no2,3*no2-2) + 1;
   C(no3,3*no3-2) = C(no3,3*no3-2) + 1;

   C(no1,3*no1-1) = C(no1,3*no1-1) + 1;
   C(no2,3*no2-1) = C(no2,3*no2-1) + 1;
   C(no3,3*no3-1) = C(no3,3*no3-1) + 1;
end
for i=1:nnodes2
   C(i,:) = C(i,:)/norm(C(i,:));
end
C = [ zeros(nnodes2,ndofs-3*nnodes2), C ];
f = zeros(nnodes2,4); Ctf = C'*f;
respos = zeros(nuzawa,1); df = zeros(nuzawa,1);

if zerobound == 1
   toremove = [ 3*boundary2(:,2)-2, 3*boundary2(:,2)-1, 3*boundary2(:,2) ];
   tokeep = setdiff(1:size(L,1),toremove);
   [Ll, Uu, Pp] = lu (MAT(tokeep,tokeep));
   kuzawa1 = .999*min(eig(MAT(tokeep,tokeep)));
else
   [Ll, Uu, Pp] = lu (MAT);  % Pp*MAT = Ll*Uu
   kuzawa1 = .999*min(eig(MAT));
end

for i=1:nuzawa % Uzawa for the contact problems
   if zerobound == 1
      SoluRG = zeros(size(L,1),4);
      SoluRG(tokeep,:) = Uu \ ( Ll \ ( Pp * ( [VE1(tokeep),VE2(tokeep),...
                                              VE3(tokeep),VE4(tokeep)] + Ctf(tokeep,:) ) ) );
   else
      SoluRG = Uu \ ( Ll \ ( Pp * ( [VE1,VE2,VE3,VE4] + Ctf ) ) );
   end
   respos(i) = norm(C*SoluRG - abs(C*SoluRG),'fro');
   fp = f;
   if kuzawa == 0
      f = f - kuzawa1*C*SoluRG;
   else
      f = f - kuzawa*C*SoluRG;
   end
   f = .5*(f + abs(f)); Ctf = C'*f;
   df(i) = norm(f-fp);
end

%if zerobound == 1
%   SoluRG = zeros(size(L,1),4);
%   SoluRG(tokeep,:) = Uu \ ( Ll \ ( Pp * ( [VE1(tokeep),VE2(tokeep),...
%                                            VE3(tokeep),VE4(tokeep)] ) ) );
%else
%   SoluRG = Uu \ ( Ll \ ( Pp * ( [VE1,VE2,VE3,VE4] ) ) );
%end

Solu1 = SoluRG(:,1); Solu2 = SoluRG(:,2); Solu3 = SoluRG(:,3); Solu4 = SoluRG(:,4);

disp(['Resolution terminated ', num2str(toc) ]);

%% Find the normal to the crack at each point
tic;
ucrsol1 = Solu1(end-3*nnodes2+1:end);
ucrsol2 = Solu2(end-3*nnodes2+1:end);
ucrsol3 = Solu3(end-3*nnodes2+1:end);
ucrsol4 = Solu4(end-3*nnodes2+1:end);

normal = zeros(2*nnodes2,1); % Normal
ucrnor = zeros(nnodes2,1);   % Norm of the gap

if newtonNor == 1
   resceil = 1e-7; %duu = 0;
   nhma = 5; uuu = zeros(10,nhma+1);
   for i=1:nnodes2
      Nu1 = ucrsol1( [ 3*i-2 ; 3*i-1 ; 3*i ] );
      Nu2 = ucrsol2( [ 3*i-2 ; 3*i-1 ; 3*i ] );
      Nu3 = ucrsol3( [ 3*i-2 ; 3*i-1 ; 3*i ] );
      Nu4 = ucrsol4( [ 3*i-2 ; 3*i-1 ; 3*i ] );

      nn  = [1;0];
      uu1 = [1;1]; uu2 = [1;1]; uu3 = [1;1]; uu4 = [1;1]; % Initialize
      uuu(:,1) = [uu1;uu2;uu3;uu4;nn];

      stopme = 0; %uuu = zeros(nhma+1,1);
      for j=1:nhma
         if norm([Nu1;Nu2;Nu3;Nu4]) == 0 % Means Nu=0
            %nn = [1,0]; % Arbitrary value
            %uu1 = [0,0]; uu2 = [0,0]; uu3 = [0,0]; uu4 = [0,0];
            uuu(:,2) = [0;0;0;0;0;0;0;0;1;0]; ress = 0;
            break;
         end

         Nn1 = .5*(nn*uu1'+uu1*nn'); Nun1 = [Nn1(1,1);Nn1(2,2);2*Nn1(1,2)];
         Nn2 = .5*(nn*uu2'+uu2*nn'); Nun2 = [Nn2(1,1);Nn2(2,2);2*Nn2(1,2)];
         Nn3 = .5*(nn*uu3'+uu3*nn'); Nun3 = [Nn3(1,1);Nn3(2,2);2*Nn3(1,2)];
         Nn4 = .5*(nn*uu4'+uu4*nn'); Nun4 = [Nn4(1,1);Nn4(2,2);2*Nn4(1,2)];

         RN1 = Nu1-Nun1; RN2 = Nu2-Nun2; RN3 = Nu3-Nun3; RN4 = Nu4-Nun4;
         ress(j) = norm([RN1;RN2;RN3;RN4]);

%         if j>1
%            if ress(j) - ress(j-1) > 0 % Cancel iteration that failed
%               uuu = uuup;
%               uu1 = uuu([1,2]); uu2 = uuu([3,4]); uu3 = uuu([5,6]); uu4 = uuu([7,8]); nn = uuu([9,10]);
%               stopme = 1;
%            end
%         end

%         if stopme==1
%            break;
%         end

         Au1 = [ uu1(1),0 ; 0,uu1(2) ; uu1(2),uu1(1) ];
         Au2 = [ uu2(1),0 ; 0,uu2(2) ; uu2(2),uu2(1) ];
         Au3 = [ uu3(1),0 ; 0,uu3(2) ; uu3(2),uu3(1) ];
         Au4 = [ uu4(1),0 ; 0,uu4(2) ; uu4(2),uu4(1) ];

         Z2 = [ 0,0;0,0;0,0 ];
         An = [ nn(1),0 ; 0,nn(2) ; nn(2), nn(1) ];
         At = [ [An;Z2;Z2;Z2] , [Z2;An;Z2;Z2] , [Z2;Z2;An;Z2] , [Z2;Z2;Z2;An] , [Au1;Au2;Au3;Au4] ];

         %C = [0,0,0,0,0,0,0,0,1,0]; % Arbitrarly impose nn(1) Cte (because singularity because of ||n||=1)
         C = [0,0,0,0,0,0,0,0,nn(1),nn(2)]; % correction is orthogonal to nn, because of ||n||=1
         Kt = [At'*At,C';C,0]; ft = [At'*[RN1;RN2;RN3;RN4];0];
         %Kt = At'*At; ft = At'*[RN1;RN2;RN3;RN4];

         %duu0 = duu;
         duu  = Kt \ ft; %uuup = uuu;
         uuu(:,j+1) = uuu(:,j) + duu([1:10]);

         uu1 = uuu([1,2],j+1); uu2 = uuu([3,4],j+1);
         uu3 = uuu([5,6],j+1); uu4 = uuu([7,8],j+1); nn = uuu([9,10],j+1);
      end

      %[~,ze] = min(ress); uuu = uuu(:,ze+1); % Select the best one
      uuu = uuu(:,end);
      uu1 = uuu([1,2]); uu2 = uuu([3,4]); uu3 = uuu([5,6]); uu4 = uuu([7,8]); nn = uuu([9,10]);

      uu1 = uu1*norm(nn); uu2 = uu2*norm(nn); uu3 = uu3*norm(nn); uu4 = uu4*norm(nn);
      nn = nn/norm(nn);
      if nn(1)<0, nn = -nn; end % Arbitrarly impose positive nn(1)
      normal([2*i-1,2*i]) = nn;

%      if i==527
%         ze
%         figure; plot(log10(ress));
%      end

      ucrnor(i) = ( uu1(1)^2 + uu1(2)^2 + uu2(1)^2 + uu2(2)^2 + uu3(1)^2 + uu3(2)^2 + uu4(1)^2 + uu4(2)^2 )^(1/2);
   end
else % newtonNor == 0
   % Test all the values for n and choose the best one
   nn = zeros(2,1); nst = 50;
   for i=1:nnodes2
      Nu1 = ucrsol1( [ 3*i-2 ; 3*i-1 ; 3*i ] ); %Nu1(3) = 2*Nu1(3); % No 2 because Voight
      Nu2 = ucrsol2( [ 3*i-2 ; 3*i-1 ; 3*i ] ); %Nu2(3) = 2*Nu2(3);
      Nu3 = ucrsol3( [ 3*i-2 ; 3*i-1 ; 3*i ] ); %Nu3(3) = 2*Nu3(3);
      Nu4 = ucrsol4( [ 3*i-2 ; 3*i-1 ; 3*i ] ); %Nu4(3) = 2*Nu4(3);

      theta = 0:pi/nst:pi;

      for j=1:nst+1
         nn(1) = cos(theta(j)); nn(2) = sin(theta(j));
         An = [ nn(1),0 ; 0,nn(2) ; nn(2),nn(1) ];
         uu = An \ [Nu1,Nu2,Nu3,Nu4];
         uu1 = uu(:,1); uu2 = uu(:,2); uu3 = uu(:,3); uu4 = uu(:,4);

         Nn1 = .5*(nn*uu1'+uu1*nn'); Nun1 = [Nn1(1,1);Nn1(2,2);2*Nn1(1,2)];
         Nn2 = .5*(nn*uu2'+uu2*nn'); Nun2 = [Nn2(1,1);Nn2(2,2);2*Nn2(1,2)];
         Nn3 = .5*(nn*uu3'+uu3*nn'); Nun3 = [Nn3(1,1);Nn3(2,2);2*Nn3(1,2)];
         Nn4 = .5*(nn*uu4'+uu4*nn'); Nun4 = [Nn4(1,1);Nn4(2,2);2*Nn4(1,2)];
         RN1 = Nu1-Nun1; RN2 = Nu2-Nun2; RN3 = Nu3-Nun3; RN4 = Nu4-Nun4;
         res1(j) = norm(RN1); res2(j) = norm(RN2); res3(j) = norm(RN3); res4(j) = norm(RN4);
         ress(j) = norm([RN1;RN2;RN3;RN4]);
      end
      [~,jj] = min(ress);

      nn(1) = cos(theta(jj)); nn(2) = sin(theta(jj));
      An = [ nn(1),0 ; 0,nn(2) ; nn(2),nn(1) ];
      uu = An \ [Nu1,Nu2,Nu3,Nu4];
      uu1 = uu(:,1); uu2 = uu(:,2); uu3 = uu(:,3); uu4 = uu(:,4);

%      if i==113
%         figure; hold on;
%         plot(theta/pi,ress,'Color','black');
%         plot(theta/pi,res1,'Color','blue');
%         plot(theta/pi,res2,'Color','red');
%         plot(theta/pi,res3,'Color','green');
%         plot(theta/pi,res4,'Color','magenta');
%      end

      normal([2*i-1,2*i]) = nn;
      ucrnor(i) = ( uu1(1)^2 + uu1(2)^2 + uu2(1)^2 + uu2(2)^2 + uu3(1)^2 + uu3(2)^2 + uu4(1)^2 + uu4(2)^2 )^(1/2);
   end
end

% % Compute the smoother ucrnor field
%M = D3m'*D3m; k = norm(M,'fro')/sqrt(nnodes2);
%OP = speye(nnodes2,nnodes2) + (resmooth/k)*M;
%toremove = boundary2(:,2); tokeep = setdiff(1:nnodes2,toremove);
%usmoot = zeros(nnodes2,1);
%usmoot(tokeep) = OP(tokeep,tokeep) \ ucrnor(tokeep);

usmoot = zeros(nnodes2,1);
lcor = .1;
for i=1:nnodes2
   x = nodes2(i,1); y = nodes2(i,2);
   dist = (x-nodes2(:,1)).^2 + (y-nodes2(:,2)).^2;
   zenodes = find(dist<lcor^2); nnod = max(size(zenodes));
   usmoot(i) = sum(ucrnor(zenodes))/nnod;
end
%usmoot = ucrnor;

%gradu = D3m*usmoot; % Gradient vector (on the elements)

% Find the local maxima
eps = 1e-12; maxima = [];
for i=1:nnodes2
   [elt,~] = find(elements2==i); nel = max(size(elt));
   zegreatest = 1; valu = usmoot(i);
   for j=1:nel
      myelem = elements2(elt(j),:);
      for k=1:3
         if usmoot(myelem(k)) > valu*(1+eps) % Tere is someone greater
            zegreatest = 0;
            break;
         end
      end
   end
   if zegreatest==1
      maxima = [maxima,i];
   end
end

nbmax = max(size(maxima));

% Filter the maxima
lcor = .2;
maxima2 = [];
for i=1:nbmax
   val = usmoot(maxima(i));
   xi = nodes2(maxima(i),1); yi = nodes2(maxima(i),2);
   keephim = 1;
   for j=1:nbmax
      if i==j, continue; end
      xj = nodes2(maxima(j),1); yj = nodes2(maxima(j),2);
      dist = (xi-xj)^2 + (yi-yj)^2;
      if dist <= lcor^2 && usmoot(maxima(j)) > usmoot(maxima(i))
         keephim = 0;
      end
   end
   if keephim == 1 % Report it
      maxima2 = [maxima2,maxima(i)];
   end
end
maxima = maxima2; nbmax = max(size(maxima));
%%% Ceil : L-curve
%nste = 20; epsi = 1e-12;
%ceils = 0:max(usmoot)/nste:max(usmoot);
%resT = zeros(nste+1,1); nZero = zeros(nste+1,1);

%for i=1:nste+1
%   cce       = ceils(i);
%   tozeroify = find(usmoot<cce); nZero(i) = nnodes2-max(size(tozeroify));
%   utrans1 = ucrsol1; utrans1([3*tozeroify-2,3*tozeroify-1,3*tozeroify]) = 0;% = Solu1(end-3*nnodes2+1:end);
%   utrans2 = ucrsol2; utrans2([3*tozeroify-2,3*tozeroify-1,3*tozeroify]) = 0;
%   utrans3 = ucrsol3; utrans3([3*tozeroify-2,3*tozeroify-1,3*tozeroify]) = 0;
%   utrans4 = ucrsol4; utrans4([3*tozeroify-2,3*tozeroify-1,3*tozeroify]) = 0;

%   SolT1 = Solu1; SolT1(end-3*nnodes2+1:end) = utrans1;
%   SolT2 = Solu2; SolT2(end-3*nnodes2+1:end) = utrans2;
%   SolT3 = Solu3; SolT3(end-3*nnodes2+1:end) = utrans3;
%   SolT4 = Solu4; SolT4(end-3*nnodes2+1:end) = utrans4;

%   resT(i) = sqrt( norm(Lhs*SolT1 - Rhs1)^2 + norm(Lhs*SolT2 - Rhs2)^2 + ...
%                   norm(Lhs*SolT3 - Rhs3)^2 + norm(Lhs*SolT4 - Rhs4)^2 );
%end

%try
%figure;
%%loglog(resT,nZero);
%plot(resT,nZero,'-*');
%legend('L-curve');
%end
%usmoot( find(usmoot<max(usmoot)/1.5) ) = 0;

disp(['Postpro terminated ', num2str(toc) ]);

%% Reconstruct the line
[~,mmm] = max(usmoot(maxima)); ind = maxima(mmm); x = nodes2(ind,1); y = nodes2(ind,2); %nn = normal([2*ind-1,2*ind]);
again = 1; i=0; step = .02; vmax = usmoot(ind);
ucrnorx = zeros(2*nnodes2); ucrnorx(1:2:2*nnodes2-1) = usmoot; % Small hack to have the right size
xx = [x,y];
while again
   UR = passMesh2D(nodes2, elements2, [x,y], [], [normal,ucrnorx]); % TODO : ponderation
   nn = UR(:,1);
   x  = x - step*nn(2); y = y + step*nn(1);
   xx = [ xx ; [x,y] ];
   if i>=100
      again = 0;
  % elseif UR(1,2) < myceil*vmax
      %again = 0;
   elseif x>=1 || x<=0 || y>=1 || y<=0
      again = 0;
   end
   i=i+1;
end
again = 1; i=0; x = nodes2(ind,1); y = nodes2(ind,2); % The other direction
xx2 = [x,y];
while again
   UR = passMesh2D(nodes2, elements2, [x,y], [], [normal,ucrnorx]); % TODO : ponderation
   nn = UR(:,1);
   x  = x + step*nn(2); y = y - step*nn(1);
   xx2 = [ xx2 ; [x,y] ];
   if i>=100
      again = 0;
  % elseif UR(1,2) < myceil*vmax
      %again = 0;
   elseif x>=1 || x<=0 || y>=1 || y<=0
      again = 0;
   end
   i=i+1;
end

if teskase == 20
   % Second line
   maxima(mmm) = [];
   [~,mmm] = max(usmoot(maxima)); ind = maxima(mmm); x = nodes2(ind,1); y = nodes2(ind,2); %nn = normal([2*ind-1,2*ind]);
   again = 1; i=0; step = .02; vmax = usmoot(ind);
   ucrnorx = zeros(2*nnodes2); ucrnorx(1:2:2*nnodes2-1) = usmoot; % Small hack to have the right size
   xxp = [x,y];
   while again
      UR = passMesh2D(nodes2, elements2, [x,y], [], [normal,ucrnorx]); % TODO : ponderation
      nn = UR(:,1);
      x  = x - step*nn(2); y = y + step*nn(1);
      xxp = [ xxp ; [x,y] ];
      if i>=100
         again = 0;
    %  elseif UR(1,2) < myceil*vmax
         %again = 0;
      elseif x>=1 || x<=0 || y>=1 || y<=0
         again = 0;
      end
      i=i+1;
   end
   again = 1; i=0; x = nodes2(ind,1); y = nodes2(ind,2); % The other direction
   xx2p = [x,y];
   while again
      UR = passMesh2D(nodes2, elements2, [x,y], [], [normal,ucrnorx]); % TODO : ponderation
      nn = UR(:,1);
      x  = x + step*nn(2); y = y - step*nn(1);
      xx2p = [ xx2p ; [x,y] ];
      if i>=100
         again = 0;
    %  elseif UR(1,2) < myceil*vmax
         %again = 0;
      elseif x>=1 || x<=0 || y>=1 || y<=0
         again = 0;
      end
      i=i+1;
   end
end

% Post-process : reconstruct u and f from Solu
fsolu1 = zeros(2*nboun2,1); usolu1 = zeros(2*nnodes2,1);
fsolu2 = zeros(2*nboun2,1); usolu2 = zeros(2*nnodes2,1);
fsolu3 = zeros(2*nboun2,1); usolu3 = zeros(2*nnodes2,1);
fsolu4 = zeros(2*nboun2,1); usolu4 = zeros(2*nnodes2,1);
szB = 0;
szD = size(LhsA,2);
for i=1:szD
   usolu1(dof2nodes(i)) = Solu1(i);
   usolu2(dof2nodes(i)) = Solu2(i);
   usolu3(dof2nodes(i)) = Solu3(i);
   usolu4(dof2nodes(i)) = Solu4(i);
end
for i=1:nboun2
   boname = boundary2(i,1);
   
   if max(size(neumann0))>0
      frdof = neumann0(find(neumann0(:,1)==boname),2) ;  % Dof of the given f
   else
      frdof = [];
   end
   fmdof = setdiff([1;2],frdof);  % Dof of the unknown f

   indicesB = szD + (szB+1:szB+size(fmdof,1));
   fsolu1(2*i-2+fmdof) = kB*Solu1(indicesB);
   fsolu2(2*i-2+fmdof) = kB*Solu2(indicesB);
   fsolu3(2*i-2+fmdof) = kB*Solu3(indicesB);
   fsolu4(2*i-2+fmdof) = kB*Solu4(indicesB);
   szB = szB + size(fmdof,1);
end

% Compute the values of f at the dofs
fsolu1no = zeros(2*nnodes2);
for i=1:size(b2nodesnoN)
   b2nodesnoNi = b2nodesnoN(i);
   noode  = floor((b2nodesnoNi+1)/2); % node
   boound = rem( find(boundary2(:,2:3)==noode)-1, nboun2 )+1; % boundaries
   
   % Compute lengths
   n1 = boundary2(boound(1),2); n2 = boundary2(boound(1),3);
   len1 = sqrt( (nodes2(n1,1)-nodes2(n2,1))^2 + (nodes2(n1,2)-nodes2(n2,2))^2 );
   n1 = boundary2(boound(2),2); n2 = boundary2(boound(2),3);
   len2 = sqrt( (nodes2(n1,1)-nodes2(n2,1))^2 + (nodes2(n1,2)-nodes2(n2,2))^2 );
   
   if rem(b2nodesnoNi,2)==0
      fsolu1no(b2nodesnoNi) = .5*( fsolu1(2*boound(1))*len1 + fsolu1(2*boound(2))*len2 );
      fsolu2no(b2nodesnoNi) = .5*( fsolu2(2*boound(1))*len1 + fsolu2(2*boound(2))*len2 );
      fsolu3no(b2nodesnoNi) = .5*( fsolu3(2*boound(1))*len1 + fsolu3(2*boound(2))*len2 );
      fsolu4no(b2nodesnoNi) = .5*( fsolu4(2*boound(1))*len1 + fsolu4(2*boound(2))*len2 );
   else
      fsolu1no(b2nodesnoNi) = .5*( fsolu1(2*boound(1)-1)*len1 + fsolu1(2*boound(2)-1)*len2 );
      fsolu2no(b2nodesnoNi) = .5*( fsolu2(2*boound(1)-1)*len1 + fsolu2(2*boound(2)-1)*len2 );
      fsolu3no(b2nodesnoNi) = .5*( fsolu3(2*boound(1)-1)*len1 + fsolu3(2*boound(2)-1)*len2 );
      fsolu4no(b2nodesnoNi) = .5*( fsolu4(2*boound(1)-1)*len1 + fsolu4(2*boound(2)-1)*len2 );
   end
end

%% Graph for u
%toplot = usolu1(b2nodesnoD); toplot2 = ur1(b2nodesnoD);
%figure;
%hold on;
%plot(toplot(2:2:end),'Color','red');
%plot(toplot2(2:2:end),'Color','blue');
%legend('Uy identified', 'Uy reference');

%try
%figure;
%hold on;
%plot(ind1+ind2+ind3+ind4,'Color','blue');
%plot(indp1+indp2+indp3+indp4,'Color','black');
%plot(indd1+indd2+indd3+indd4,'Color','green');
%plot(inds1+inds2+inds3+inds4,'Color','magenta');
%legend('stop Picard','stop L-curve','stop DL-curve','stop SL-curve');
%end

%try
%figure;
%hold on;
%plot(log10(phiP),'Color','blue');
%plot(log10(phiL),'Color','black');
%plot(log10(phiD),'Color','green');
%plot(log10(phiS),'Color','magenta');
%plot(log10(phi),'Color','red');
%legend('residual Picard','residual L-curve','residual DL-curve','residual SL-curve','residual');
%end

% Graph for f
if min(size(b2nodesnoN))>0
   toplot = fsolu4no(b2nodesnoN); toplot2 = fr4(b2nodesnoN);
   try
   figure;
   hold on;
   plot(toplot(1:2:end-1),'Color','red');
   plot(toplot2(1:2:end-1),'Color','blue');
   legend('Fy identified (4)', 'Fy reference (4)');
   end
end
% Graph for u
if min(size(b2nodesnoD))>0
   toplot = usolu4(b2nodesnoD); toplot2 = ur4(b2nodesnoD);
   try
   figure;
   hold on;
   plot(toplot(1:2:end-1),'Color','red');
   plot(toplot2(1:2:end-1),'Color','blue');
   legend('Uy identified (4)', 'Uy reference (4)');
   end
end

% Graph for [[u]]
x6 = nodes(6,1); y6 = nodes(6,2); x5 = nodes(5,1); y5 = nodes(5,2);
if teskase == 20
   x8 = nodes(8,1); y8 = nodes(8,2); x7 = nodes(7,1); y7 = nodes(7,2);
end
%try
%figure; hold on;
%set(gca, 'fontsize', 20);
%patch('Faces',elements2,'Vertices',nodes2,'FaceVertexCData',ucrsol1(1:3:3*nnodes2-2),'FaceColor','interp');
%plot( [x5,x6], [y5,y6] ,'Color', 'magenta', 'LineWidth',5);
%set(colorbar, 'fontsize', 20);
%end
%try
%figure; hold on;
%set(gca, 'fontsize', 20);
%patch('Faces',elements2,'Vertices',nodes2,'FaceVertexCData',ucrsol1(2:3:3*nnodes2-1),'FaceColor','interp');
%plot( [x5,x6], [y5,y6] ,'Color', 'magenta', 'LineWidth',5);
%set(colorbar, 'fontsize', 20);
%end
%try
%figure; hold on;
%set(gca, 'fontsize', 20);
%patch('Faces',elements2,'Vertices',nodes2,'FaceVertexCData',ucrsol1(3:3:3*nnodes2),'FaceColor','interp');
%plot( [x5,x6], [y5,y6] ,'Color', 'magenta', 'LineWidth',5);
%set(colorbar, 'fontsize', 20);
%end

try
figure; hold on;
set(gca, 'fontsize', 20);
%patch('Faces',elements2,'Vertices',nodes2,'FaceAlpha',0);
patch('Faces',elements2,'Vertices',nodes2,'FaceVertexCData',ucrnor,'FaceColor','interp');
plot( [x5,x6], [y5,y6] ,'Color', 'magenta', 'LineWidth',5);
if teskase == 20
   plot( [x7,x8], [y7,y8] ,'Color', 'magenta', 'LineWidth',5);
end
set(colorbar, 'fontsize', 20);
for i=1:nnodes2
   nn = normal([2*i-1,2*i]);
   plot( [nodes2(i,1),nodes2(i,1)+nn(1)/20] , [nodes2(i,2),nodes2(i,2)+nn(2)/20] , 'Color' , 'white', 'LineWidth',3 );
end
axis equal;
end

try
figure; hold on;
set(gca, 'fontsize', 20);
%patch('Faces',elements2,'Vertices',nodes2,'FaceAlpha',0);
patch('Faces',elements2,'Vertices',[nodes2,usmoot],'FaceVertexCData',usmoot,'FaceColor','interp');
%plot( [x5,x6], [y5,y6] ,'Color', 'magenta', 'LineWidth',5);
set(colorbar, 'fontsize', 20);
%axis equal;
end

try
figure; hold on;
set(gca, 'fontsize', 20);
%patch('Faces',elements2,'Vertices',nodes2,'FaceAlpha',0);
patch('Faces',elements2,'Vertices',nodes2,'FaceVertexCData',ucrnor,'FaceColor','interp');
plot( [x5,x6], [y5,y6] ,'Color', 'magenta', 'LineWidth',5);
if teskase == 20
   plot( [x7,x8], [y7,y8] ,'Color', 'magenta', 'LineWidth',5);
end
plot( xx(:,1) , xx(:,2)  ,'Color', 'green', 'LineWidth',3);
plot( xx2(:,1), xx2(:,2) ,'Color', 'green', 'LineWidth',3);
if teskase == 20
   plot( xxp(:,1) , xxp(:,2)  ,'Color', 'green', 'LineWidth',3);
   plot( xx2p(:,1), xx2p(:,2) ,'Color', 'green', 'LineWidth',3);
end
set(colorbar, 'fontsize', 20);
axis equal;
end

%try
%figure; hold on;
%set(gca, 'fontsize', 20);
%patch('Faces',elements2,'Vertices',nodes2,'FaceVertexCData',ucrsol2(1:3:3*nnodes2-2),'FaceColor','interp');
%plot( [x5,x6], [y5,y6] ,'Color', 'magenta', 'LineWidth',5);
%set(colorbar, 'fontsize', 20);
%end
%try
%figure; hold on;
%set(gca, 'fontsize', 20);
%patch('Faces',elements2,'Vertices',nodes2,'FaceVertexCData',ucrsol2(2:3:3*nnodes2-1),'FaceColor','interp');
%plot( [x5,x6], [y5,y6] ,'Color', 'magenta', 'LineWidth',5);
%set(colorbar, 'fontsize', 20);
%end
%try
%figure; hold on;
%set(gca, 'fontsize', 20);
%patch('Faces',elements2,'Vertices',nodes2,'FaceVertexCData',ucrsol2(3:3:3*nnodes2),'FaceColor','interp');
%plot( [x5,x6], [y5,y6] ,'Color', 'magenta', 'LineWidth',5);
%set(colorbar, 'fontsize', 20);
%end

%erroru = norm(usolu1(b2nodesnoD)-ur1(b2nodesnoD))   / norm(ur1(b2nodesnoD));
%errorf = norm(fsolu1no(b2nodesnoN)-fr1(b2nodesnoN)) / norm(fr1(b2nodesnoN));
