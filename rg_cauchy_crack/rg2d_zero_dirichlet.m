% 22/05/2018
% ID fissure par RG Petrov-Galerkin et fct-test nulle sur le bord de Dirichlet
% Il s'agit d'une particularisation de l'algo général : pas mal d'opérateurs sont construits inutilement
% The top boundary is Dirichlet only and the other are redondant

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
mur         = 1e1;%1e1;%2e3;    % Regularization parameter
regular     = 1;      % Use the derivative regularization matrix (0 : Id, 1 : derivative)
froreg      = 1;      % Frobenius preconditioner
recompute   = 0;      % Recompute the operators
theta1      = pi;     %3.7296;%pi;%pi/2; 3.8273
theta2      = 0;      %0.58800;%0%3*pi/2;5.7608  % Initial angles of the crack
anglestep   = 0;%pi/1000;  % Step in angle for Finite Differences anglestep = 0 means auto-adaptation
kauto       = 20;     % Coefficient for the auto-adaptation
nbstep      = 20;      % Nb of Newton Iterations
Npg         = 2;      % Nb Gauss points
ordertest   = 20;     % Order of test fonctions
zerobound   = 1;      % Put the boundaries of the crack to 0
nuzawa      = 100;     % (nuzawa = 1 means no Uzawa)
kuzawa      = 0;%1e2;     % Parameters of the Uzawa algorithm (kuzawa = 0 means best parameter)
ndofcrack   = 20;      % Nb of elements on the crack
teskase     = 4;       % Choice of the test case

nbDirichlet = [];

% Boundary conditions
dirichlet  = [3,1,0 ; 3,2,0];

%neumann1   = [1,2,-fscalar];
%neumann2   = [2,1,fscalar ; 4,1,-fscalar];
%neumann3   = [1,2,-fscalar;2,1,fscalar;4,1,-fscalar];
%neumann4   = [1,1,fscalar];

neumann1   = [1,2,-fscalar];
neumann2   = [2,1,fscalar ; 4,1,-fscalar];
neumann3   = [1,1,-fscalar ; 1,2,-fscalar ; 4,1,-fscalar ; 4,2,-fscalar];
neumann4   = [2,1,fscalar ; 2,2,-fscalar ; 1,1,fscalar ; 1,2,-fscalar];

%neumann1   = [3,2,fscalar ; 1,2,-fscalar];
%neumann2   = [2,1,fscalar ; 4,1,-fscalar];
%neumann3   = [3,1,fscalar ; 3,2,fscalar ; 2,1,fscalar ; 2,2,fscalar ; ...
%              1,1,-fscalar ; 1,2,-fscalar ; 4,1,-fscalar ; 4,2,-fscalar];
%neumann4   = [3,1,-fscalar ; 3,2,fscalar ; 2,1,fscalar ; 2,2,-fscalar ; ...
%              1,1,fscalar ; 1,2,-fscalar ; 4,1,-fscalar ; 4,2,fscalar];
           
dirichlet0 = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ; 3,1,0 ; 3,2,0 ; 4,1,0 ; 4,2,0 ];
neumann0   = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ; 4,1,0 ; 4,2,0 ];

% First, import the mesh
if teskase == 3
   [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_c_squared3.msh' );
elseif teskase == 4
   [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_c_squared4.msh' );
elseif teskase == 2
   [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_c_squared2.msh' );
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
f1 = Kinter*u1; f2 = Kinter*u2; f3 = Kinter*u3; f4 = Kinter*u4;

ui = reshape(u1,2,[])';  ux = ui(:,1);  uy = ui(:,2);

u1ref = u1; u2ref = u2; u3ref = u3; u4ref = u4;
f1ref = f1; f2ref = f2; f3ref = f3; f4ref = f4;

% Compute stress :
sigma = stress(u1,E,nu,nodes,elements,order,1,ntoelem);

% Output :
plotGMSH({u1,'U1';u2,'U2';u3,'U3';u4,'U4'}, elements, nodes, 'output/reference');
%plotGMSH({f1,'F1';f2,'F2';f3,'F3';f4,'F4'}, elements, nodes, 'output/reference');

%% Find the reference angles of the crack
x1 = nodes(5,1); y1 = nodes(5,2);
x2 = nodes(6,1); y2 = nodes(6,2);

xmin = min(nodes(:,1)); xmax = max(nodes(:,1));
ymin = min(nodes(:,2)); ymax = max(nodes(:,2));
xb = .5*(xmin+xmax); yb = .5*(ymin+ymax);

% 4 intersections of the line with the boundaries
t1 = (xmin-x1)/(x2-x1); t2 = (xmax-x1)/(x2-x1);
t3 = (ymin-y1)/(y2-y1); t4 = (ymax-y1)/(y2-y1);
xy11 = [xmin;y1+(y2-y1)*t1];
xy22 = [xmax;y1+(y2-y1)*t2];
xy33 = [x1+(x2-x1)*t3;ymin];
xy44 = [x1+(x2-x1)*t4;ymax];
xy1234 = [xy11,xy22,xy33,xy44];
% limit to those that are inside the square
elim1 = find(xy1234(1,:)>1); elim2 = find(xy1234(2,:)>1);
elim3 = find(xy1234(1,:)<0); elim4 = find(xy1234(2,:)<0);
total = setdiff( [1,2,3,4] , union(union(elim1,elim2),union(elim3,elim4)) );
xyt = xy1234(:,total); % Normally, this one should be of size 2
xy1 = xyt(:,1); xy2 = xyt(:,2); xy1r = xy1; xy2r = xy2;
theta1ref = atan2( xy1(2)-yb , xy1(1)-xb );
theta2ref = atan2( xy2(2)-yb , xy2(1)-xb );

theta1ref = mod(theta1ref,2*pi);
theta2ref = mod(theta2ref,2*pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ nodes2,elements2,ntoelem2,boundary2,order] = readmesh( 'meshes/rg_refined/plate_nu_r3.msh' );
nnodes2 = size(nodes2,1);
[K2,C2,nbloq2,node2c2,c2node2] = Krig2 (nodes2,elements2,mat,order,boundary2,dirichlet);
Kinter2 = K2( 1:2*nnodes2, 1:2*nnodes2 );

%% Pass meshes
UF = [u1,u2,u3,u4,f1,f2,f3,f4];
UR = passMesh2D(nodes, elements, nodes2, elements2, UF, 0);

ur1 = UR(:,1); fr1 = UR(:,5)*50; % 50 is the ratio between mesh sizes (hack)
ur2 = UR(:,2); fr2 = UR(:,6)*50;
ur3 = UR(:,3); fr3 = UR(:,7)*50;
ur4 = UR(:,4); fr4 = UR(:,8)*50;

% Recover the Neumann on 3
fr3p = fr3 - loading(0,nodes2,boundary2,neumann3);
fe3  = elem_forces(fr3p,nodes2,boundary2,3,order);

% Add the noise
u1n = ur1; u2n = ur2; u3n = ur3; u4n = ur4;
am1 = sqrt(mean(ur1.^2)); am2 = sqrt(mean(ur2.^2));
am3 = sqrt(mean(ur3.^2)); am4 = sqrt(mean(ur4.^2));
br1 = randn(2*nnodes2,1); br2 = randn(2*nnodes2,1);
br3 = randn(2*nnodes2,1); br4 = randn(2*nnodes2,1);
%noise = load('noises/cauchyRG.mat');
%br1 = noise.br1; br2 = noise.br2; br3 = noise.br3; br4 = noise.br4;
%ur1 = ( 1 + br*br1 ) .* ur1; ur2 = ( 1 + br*br2 ) .* ur2;
%ur3 = ( 1 + br*br3 ) .* ur3; ur4 = ( 1 + br*br4 ) .* ur4;
ur1 = ur1 + am1*br*br1; ur2 = ur2 + am1*br*br2;
ur3 = ur3 + am3*br*br3; ur4 = ur4 + am4*br*br4;

plotGMSH({ur1,'U1';ur2,'U2';ur3,'U3';ur4,'U4'},elements2, nodes2, 'output/bound');

%% Remove rigid modes (useless)
%R  = rigidModes(nodes2);
%UU = [ur1,ur2,ur3,ur4] - R*((R'*R)\(R'*[ur1,ur2,ur3,ur4]));
%ur1 = UU(:,1); ur2 = UU(:,2); ur3 = UU(:,3); ur4 = UU(:,4);

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

[ nodesu,elementsu,ntoelemu,boundaryu,orderu] = readmesh( 'meshes/rg_refined/plate_nu.msh' );
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
ncoef = 2*(nmax+1)^2; neq = ncoef;

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

% Add the condition : v(y=1) = 0
Ma = zeros( nmax+1, (nmax+1)^2 );
Mb = zeros( nmax+1, (nmax+1)^2 );
Z  = zeros( nmax+1, (nmax+1)^2 );
for i=0:nmax % \forall i, \sum a_{ij} = 0
   for j=0:nmax
      Ma(i+1,(nmax+1)*i+j+1) = 1; % For a
      Mb(i+1,(nmax+1)*i+j+1) = 1; % And for b
   end
end

M = [M;Ma,Z;Z,Mb];

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
   
   if min(size( intersect(lis1,lis2) )) > 0 || min(size( intersect(lis3,lis4) )) > 0  % || (size(lis3,1) > 0 || size(lis4,1) > 0)%
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
   % Reference force from the BC's
%   if exno(1) == 1 % Bound 2
%      fer1 = [0;0]; fer2 = [fscalar;0];
%      fer3 = [fscalar;0]; fer4 = [0;0];
%   elseif exno(1) == -1 % Bound 4
%      fer1 = [0;0]; fer2 = -[fscalar;0];
%      fer3 = [-fscalar;0]; fer4 = [0;0];
%   elseif exno(2) == 1 % Bound 3
%      fer1 = [0;0]; fer2 = [0;0];
%      fer3 = [0;0]; fer4 = [0;0];
%   elseif exno(2) == -1 % Bound 1
%      fer1 = -[0;fscalar]; fer2 = [0;0];
%      fer3 = -[0;fscalar]; fer4 = [fscalar;0];
%   end

   if exno(1) == 1 % Bound 2
      fer1 = [0;0]; fer2 = [fscalar;0];
      fer3 = [0;0]; fer4 = [fscalar;-fscalar];
   elseif exno(1) == -1 % Bound 4
      fer1 = [0;0]; fer2 = -[fscalar;0];
      fer3 = [-fscalar;-fscalar]; fer4 = [0;0];
   elseif exno(2) == 1 % Bound 3
      fer1 = [0;0]; fer2 = [0;0];
      fer3 = [0;0]; fer4 = [0;0];
   elseif exno(2) == -1 % Bound 1
      fer1 = -[0;fscalar]; fer2 = [0;0];
      fer3 = [-fscalar;-fscalar]; fer4 = [fscalar;-fscalar];
   end

%   if exno(1) == 1 % Bound 2
%      fer1 = [0;0]; fer2 = [fscalar;0];
%      fer3 = [fscalar;fscalar]; fer4 = [fscalar;-fscalar];
%   elseif exno(1) == -1 % Bound 4
%      fer1 = [0;0]; fer2 = -[fscalar;0];
%      fer3 = [-fscalar;-fscalar]; fer4 = [-fscalar;fscalar];
%   elseif exno(2) == 1 % Bound 3
%      fer1 = [0;0]; fer2 = [0;0];
%      fer3 = [0;0]; fer4 = [0;0];
%   elseif exno(2) == -1 % Bound 1
%      fer1 = -[0;fscalar]; fer2 = [0;0];
%      fer3 = [-fscalar;-fscalar]; fer4 = [fscalar;-fscalar];
%   end
     
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

%% Restrict the data on the given part
f_known1(tofindN) = 0; f_known2(tofindN) = 0;
f_known3(tofindN) = 0; f_known4(tofindN) = 0;
u_known1(tofindD) = 0; u_known2(tofindD) = 0;
u_known3(tofindD) = 0; u_known4(tofindD) = 0;

ndofs = size(f_known1(tofindN),1) + size(u_known1(tofindD),1) + 2*(ndofcrack+1);

%% Restrict the operators
Rfm = Rf(:,tofindN); Rfr = Rf(:,knownN);
Rum = Ru(:,tofindD); Rur = Ru(:,knownD);

LhsA = -Rum;
LhsB = Rfm;
Rhs1 = Ru*u_known1 - Rf*f_known1; % Full size despite all the zeros
Rhs2 = Ru*u_known2 - Rf*f_known2;
Rhs3 = Ru*u_known3 - Rf*f_known3;
Rhs4 = Ru*u_known4 - Rf*f_known4;
   
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

theta1rec = theta1;
theta2rec = theta2;
theta1b   = theta1;
theta2b   = theta2;
nbnotwork = 0; % Records the number of failures
previous = []; % Stores the index of the previous solutions
tic
if anglestep == 0 % Initialize the step
   anglestep1 = pi/kauto;
   anglestep2 = pi/kauto;
else
   anglestep1 = anglestep;
   anglestep2 = anglestep;
end

for iter = 1:nbstep % Newton loop

   pbmax = 2;
%   if iter == nbstep
%      pbmax == 0;
%   end
   if anglestep == 0 && iter > 1
      anglestep1 = dtheta(1)/kauto; anglestep2 = dtheta(2)/kauto;
   end

   for pb = 0:pbmax % Construct all the problems

      if pb == 0
         theta1c = theta1; theta2c = theta2;
      elseif pb == 1
         theta1c = theta1+anglestep1; theta2c = theta2;
      else
         theta1c = theta1; theta2c = theta2+anglestep2;
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% Compute the crack RHS (from angle)

      % Intersections
      %xb = mean(nodes2(:,1)); yb = mean(nodes2(:,2)); % Barycenter
      xmin = min(nodes2(:,1)); xmax = max(nodes2(:,1));
      ymin = min(nodes2(:,2)); ymax = max(nodes2(:,2));
      xb = .5*(xmin+xmax); yb = .5*(ymin+ymax);
      ma11 = [ 0, -cos(theta1c) ; ymin-ymax, -sin(theta1c) ]; bh1 = [ xb-xmax ; yb-ymax];
      ma12 = [ 0, -cos(theta2c) ; ymin-ymax, -sin(theta2c) ];
      ma21 = [ xmax-xmin, -cos(theta1c) ; 0, -sin(theta1c) ]; bh2 = [ xb-xmin ; yb-ymax];
      ma22 = [ xmax-xmin, -cos(theta2c) ; 0, -sin(theta2c) ];
      ma31 = [ 0, -cos(theta1c) ; ymax-ymin, -sin(theta1c) ]; bh3 = [ xb-xmin ; yb-ymin];
      ma32 = [ 0, -cos(theta2c) ; ymax-ymin, -sin(theta2c) ];
      ma41 = [ xmin-xmax, -cos(theta1c) ; 0, -sin(theta1c) ]; bh4 = [ xb-xmax ; yb-ymin];
      ma42 = [ xmin-xmax, -cos(theta2c) ; 0, -sin(theta2c) ];

      tr1s = zeros(2,4); % All the intersections
      tr2s = zeros(2,4);
      ran1 = [ rank(ma11) , rank(ma21) , rank(ma31) , rank(ma41) ]; % Ranks
      ran2 = [ rank(ma12) , rank(ma22) , rank(ma32) , rank(ma42) ]; %

      % Warning free implementation
      if ran1(1) == 2, tr1s(:,1) = ma11\bh1; end
      if ran1(2) == 2, tr1s(:,2) = ma21\bh2; end
      if ran1(3) == 2, tr1s(:,3) = ma31\bh3; end
      if ran1(4) == 2, tr1s(:,4) = ma41\bh4; end

      if ran2(1) == 2, tr2s(:,1) = ma12\bh1; end
      if ran2(2) == 2, tr2s(:,2) = ma22\bh2; end
      if ran2(3) == 2, tr2s(:,3) = ma32\bh3; end
      if ran2(4) == 2, tr2s(:,4) = ma42\bh4; end

      tr1s = tr1s( :, find(tr1s(2,:)>0) ); % Radius must be > 0
      [~,num] = min(tr1s(2,:));
      tr1 = tr1s(:,num); r = tr1(2);
      xy1 = [ xb + r*cos(theta1c) ; yb + r*sin(theta1c) ]; % coords of the right intersection
   
      tr2s = tr2s( :, find(tr2s(2,:)>0), : ); % Radius must be > 0
      [~,num] = min(tr2s(2,:));
      tr2 = tr2s(:,num); r = tr2(2);
      xy2 = [ xb + r*cos(theta2c) ; yb + r*sin(theta2c) ]; % coords of the right intersection

      %% Build the elements on the crack's line
      step = (xy2-xy1)/ndofcrack;

      if step(1)*step(2) ~= 0
         nodes3 = [ xy1(1):step(1):xy2(1) ; xy1(2):step(2):xy2(2) ];
      elseif step(1) == 0
         nodes3 = [ xy1(1)*ones(1,ndofcrack+1) ; xy1(2):step(2):xy2(2) ];
      else %if step(2) == 0
         nodes3 = [ xy1(1):step(1):xy2(1) ; xy1(2)*ones(1,ndofcrack+1) ];
      end

      nodes3 = nodes3'; nnodes3 = size(nodes3,1);
      boundary3 = zeros(nnodes3-1,3);
      boundary3(:,2) = 1:nnodes3-1; boundary3(:,3) = 2:nnodes3;
      nboun3 = size(boundary3,1);
   
      extnorm3 = zeros(nboun3,2);
      for i=1:nboun3
         x1 = nodes3( boundary3(i,2) , 1 ); y1 = nodes3( boundary3(i,2) , 2 );
         x2 = nodes3( boundary3(i,3) , 1 ); y2 = nodes3( boundary3(i,3) , 2 );
         extnorm3(i,:) = [ y1-y2 , -(x1-x2) ];
         extnorm3(i,:) = extnorm3(i,:)/norm(extnorm3(i,:));
      end
   
      Rucij = zeros(size(coef,1),2*nnodes3);
   
      if order>1 warning('Mesh order is too high') end
      Ng = Npg;
      [ Xg, Wg ] = gaussPt1d( Ng );
      Ndots = size(Wg,1); % (= Ng by the way)

      %% Build the list of Gauss points, and of construction-functions
      Xxg = zeros( nboun3*Ndots,1 ); Yyg = zeros( nboun3*Ndots,1 );
      Wwg = zeros( nboun3*Ndots,1 );
      Phi = sparse( 2*nboun3*Ndots, 2*nnodes3 );
      exnor = zeros( nboun3*Ndots,2); % Normal

      for i=1:nboun3
         bonod = boundary3(i,:); exno = extnorm3(i,:)';
   
         no1 = bonod(2); no2 = bonod(3); boname = bonod(1);
         x1 = nodes3(no1,1); y1 = nodes3(no1,2);
         x2 = nodes3(no2,1); y2 = nodes3(no2,2);
      
         len = sqrt( (x1-x2)^2 + (y1-y2)^2 );
     
         % Dofs in the numerotation of the boundary nodes
         indDtot = [ 2*i-1, 2*i, 2*(i+1)-1, 2*(i+1) ];   
                
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
         end
      end

      %% Build the matrix of test-functions
      Sv = zeros( size(coef,1), 2*nboun3*Ndots );
      exnoX = exnor(:,1); exnoY = exnor(:,2);

      for ii=0:nmax
         nd1 = nmax-ii;
         for jj=0:nd1
                  
               sloc11a = zeros(nboun3*Ndots,1); sloc12a = zeros(nboun3*Ndots,1);
               sloc22a = zeros(nboun3*Ndots,1); sloc11b = zeros(nboun3*Ndots,1);
               sloc12b = zeros(nboun3*Ndots,1); sloc22b = zeros(nboun3*Ndots,1);
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

               index = (nmax+1)*ii + jj + 1;

               Sv( index, 1:2:2*nboun3*Ndots-1 )              = Wwg' .* fpaax';
               Sv( (nmax+1)^2 + index, 1:2:2*nboun3*Ndots-1 ) = Wwg' .* fpabx';

               Sv( index, 2:2:2*nboun3*Ndots )                = Wwg' .* fpaay';
               Sv( (nmax+1)^2 + index, 2:2:2*nboun3*Ndots )   = Wwg' .* fpaby';
         end
      end

      Rucij = Sv*Phi; clear Sv; clear Phi;
      Ruc = coef'*Rucij; noRuc = norm(Ruc,'fro');
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      nodeMass3 = zeros(2*nboun3+2);
      for i=1:nboun3
         coef1 = [ i, i+1 ];
         len = norm(nodes3(boundary3(i,2),:) - nodes3(boundary3(i,3),:));

         ico = [ 2*i-1, 2*i ];
         cco = [ 2*coef1-1, 2*coef1 ];
         nodeMass3(cco,cco) = nodeMass3(cco,cco) + len/2 * [ eye(2), zeros(2) ; ...
                                                             zeros(2), eye(2), ];
      end
   
      %% Pure RG : Solve the linear system and recover the unknowns
%      if iter == 1 && pb == 0 % kB must be defined only once
         if froreg == 1 && min(size(LhsB))>0
            kB = sqrt(norm(LhsA,'fro')^2+norm(Ruc,'fro')^2)/norm(LhsB,'fro'); % In order to regularize the stuff
         %   kB = norm(LhsA,'fro')/norm(LhsB,'fro');
         else
            kB = 1;
         end
%      end
      Lhs = Ruc;%[LhsA,kB*LhsB,Ruc];

      A = Lhs'*Lhs; sA = size(A,1);
   
      Z12    = zeros(size(Mum,1),size(Mfm,2)); Z13 = zeros(size(Mum,1),size(nodeMass3,2));
      Z23    = zeros(size(Mfm,1),size(nodeMass3,2));
      Mtot   = nodeMass3;%[ Mum, Z12, Z13  ; Z12', Mfm, Z23 ; Z13', Z23', nodeMass3 ];
      Msmall = [ 0*Mum, Z12, Z13  ; Z12', 0*Mfm, Z23 ; Z13', Z23', nodeMass3 ]; % Only on the crack

      % Differential regularization matrix
      D3 = zeros(2*nboun3+2);
      for i=1:nboun3
         coefU = [ i , i+1 ];
         len = norm(nodes3(i,:) - nodes3(i+1,:));
         D3( 2*coefU-1, 2*coefU-1 ) = D3( 2*coefU-1, 2*coefU-1 ) + 1/len*[1,-1;-1,1];
         D3( 2*coefU, 2*coefU )     = D3( 2*coefU, 2*coefU )     + 1/len*[1,-1;-1,1];
      end
      Dsmall = D3;%[ 0*Mum, Z12, Z13  ; Z12', 0*Mfm, Z23 ; Z13', Z23', D3 ];

      if regular == 1
         Z12 = zeros( size(Duu,1) , size(Dfu,2) ); Z13 = zeros( size(Duu,1), size(D3,2) );
         Z23 = zeros( size(Dfu,1), size(D3,2) );
         Dtot = D3;%[ Duu ,Z12, Z13 ; Z12', Dfu, Z23 ; Z13', Z23', D3 ];
         L = Dtot;% + Mtot * norm(Dtot,'fro')/(10*norm(Mtot,'fro'));

         L12 = [ Du(tofindD,knownD) ; zeros(size(Dfu,2),size(Duk,2)) ; zeros(2*nboun3+2,size(Duk,2)) ];
         L2 = Du(knownD,knownD);
         L121 = L12*u_known1(knownD); L21 = u_known1(knownD)'*Du(knownD,knownD)*u_known1(knownD);
         L122 = L12*u_known2(knownD); L22 = u_known2(knownD)'*Du(knownD,knownD)*u_known2(knownD);
         L123 = L12*u_known3(knownD); L23 = u_known3(knownD)'*Du(knownD,knownD)*u_known3(knownD);
         L124 = L12*u_known4(knownD); L24 = u_known4(knownD)'*Du(knownD,knownD)*u_known4(knownD);
      else
%      L = eye(size(A));
         L = Mtot;
         L21 = 0; L22 = 0; L23 = 0; L24 = 0;
         L121 = zeros(ndofs,1); L122 = zeros(ndofs,1);
         L123 = zeros(ndofs,1); L124 = zeros(ndofs,1);
      end
      ninfty = 0;
   %L = eye(size(A));

      sL = real(L^(1/2));   % Square root of L (as usual, real should not be)

      MAT = Lhs'*Lhs + mur*L;
      VE1 = Lhs'*Rhs1; VE2 = Lhs'*Rhs2;
      VE3 = Lhs'*Rhs3; VE4 = Lhs'*Rhs4;

      % Build the contact inequation condition
      C = zeros(nnodes3, 2*nnodes3);
      C(1,2*1-1) = extnorm3(1,1); C(1,2*1) = extnorm3(1,2);
      C(nnodes3,2*nnodes3-1) = extnorm3(end,1); C(nnodes3,2*nnodes3) = extnorm3(end,2);
      for i=2:nnodes3-1
         C(i,2*i-1) = .5*(extnorm3(i-1,1) + extnorm3(i,1));
         C(i,2*i)   = .5*(extnorm3(i-1,2) + extnorm3(i,2));
      end
      %C = [ zeros(nnodes3,ndofs-2*nnodes3), C ];
      f = zeros(nnodes3,4); Ctf = C'*f;
      respos = zeros(nuzawa,1); df = zeros(nuzawa,1);

      if zerobound == 1
         toremove = [ 1, 2, 2*ndofcrack+1, 2*ndofcrack+2 ];
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
         respos(i) = norm(C*SoluRG - abs(C*SoluRG));
         fp = f;
         if kuzawa == 0
            f = f - kuzawa1*C*SoluRG;
         else
            f = f - kuzawa*C*SoluRG;
         end
         f = .5*(f + abs(f)); Ctf = C'*f;
         df(i) = norm(f-fp);
      end

      Solu1 = SoluRG(:,1); Solu2 = SoluRG(:,2); Solu3 = SoluRG(:,3); Solu4 = SoluRG(:,4);

      if pb == 0
         nor1 = Solu1'*Lhs'*Lhs*Solu1 - 2*Solu1'*Lhs'*Rhs1 + Rhs1'*Rhs1 + mur*Solu1'*L*Solu1;% + 2*mur*Solu1'*L121 + mur*L21;
         nor2 = Solu2'*Lhs'*Lhs*Solu2 - 2*Solu2'*Lhs'*Rhs2 + Rhs2'*Rhs2 + mur*Solu2'*L*Solu2;% + 2*mur*Solu2'*L122 + mur*L22;
         nor3 = Solu3'*Lhs'*Lhs*Solu3 - 2*Solu3'*Lhs'*Rhs3 + Rhs3'*Rhs3 + mur*Solu3'*L*Solu3;% + 2*mur*Solu3'*L123 + mur*L23;
         nor4 = Solu4'*Lhs'*Lhs*Solu4 - 2*Solu4'*Lhs'*Rhs4 + Rhs4'*Rhs4 + mur*Solu4'*L*Solu4;% + 2*mur*Solu4'*L124 + mur*L24;
         phi( iter )  = nor1 + nor2 + nor3 + nor4;

         res1 = Lhs*Solu1 - Rhs1; % Residuals
         res2 = Lhs*Solu2 - Rhs2; %
         res3 = Lhs*Solu3 - Rhs3; %
         res4 = Lhs*Solu4 - Rhs4; %

         rel1 = sL*Solu1; rel2 = sL*Solu2; rel3 = sL*Solu3; rel4 = sL*Solu4; 

         res = [ res1 ; res2 ; res3 ; res4 ];
         rel = [ rel1 ; rel2 ; rel3 ; rel4 ];
         Ax  = [ Lhs*Solu1 ; Lhs*Solu2 ; Lhs*Solu3 ; Lhs*Solu4 ];
         xt  = [ Solu1 ; Solu2 ; Solu3 ; Solu4 ];
         sLx = [ sL*Solu1 ; sL*Solu2 ; sL*Solu3 ; sL*Solu4 ];
         Lh0 = Lhs; sL0 = sL; L1210 = L121; L1220 = L122; L1230 = L123; L1240 = L124;

         Solu10 = Solu1; Solu20 = Solu2; Solu30 = Solu3; Solu40 = Solu4; % Store the values

         if phi(iter) == min(phi)
            nodes3b = nodes3;
         end

      elseif pb == 1
         D1    = ([ Lhs*Solu1 ; Lhs*Solu2; Lhs*Solu3; Lhs*Solu4 ] - Ax)/anglestep1;
         DL1   = ([ sL*Solu1 ; sL*Solu2 ; sL*Solu3 ; sL*Solu4 ] - sLx)/anglestep1;

%         D1    = ((Lhs-Lh0) * [Solu10,Solu20,Solu30,Solu40] )/anglestep1; D1 = D1(:);
%         DL1   = ((sL -sL0) * [Solu10,Solu20,Solu30,Solu40] )/anglestep1; DL1 = DL1(:);
      elseif pb == 2
         D2    = ([ Lhs*Solu1 ; Lhs*Solu2; Lhs*Solu3; Lhs*Solu4 ] - Ax)/anglestep2;
         DL2   = ([ sL*Solu1 ; sL*Solu2 ; sL*Solu3 ; sL*Solu4 ] - sLx)/anglestep2;

%         D2    = ((Lhs-Lh0) * [Solu10,Solu20,Solu30,Solu40] )/anglestep2; D2 = D2(:);
%         DL2   = ((sL -sL0) * [Solu10,Solu20,Solu30,Solu40] )/anglestep2; DL2 = DL2(:);
      end
   end
   D = [D1,D2]; DL = [DL1,DL2];% DL12 = [DL121,DL122];
   dtheta = - ( D'*D + mur*DL'*DL ) \ ( D'*res + mur*DL'*rel );

   theta1 = theta1 + dtheta(1); theta1 = mod(theta1,2*pi);
   theta2 = theta2 + dtheta(2); theta2 = mod(theta2,2*pi);

   theta1rec(iter+1) = theta1;
   theta2rec(iter+1) = theta2;
end
disp(['Iterative method terminated ', num2str(toc) ]);

Solu1 = Solu10; Solu2 = Solu20; Solu3 = Solu30; Solu4 = Solu40;

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

ucrsol1 = Solu1(end-2*nnodes3+1:end);
ucrsol2 = Solu2(end-2*nnodes3+1:end);
ucrsol3 = Solu3(end-2*nnodes3+1:end);
ucrsol4 = Solu4(end-2*nnodes3+1:end);

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

theta1pi = theta1rec/pi;
theta2pi = theta2rec/pi;

[~,thisone] = min(phi);
try
figure;
hold on;
plot(theta1pi,'Color','red');
plot(theta2pi,'Color','blue');
plot(theta1ref/pi*theta1pi./theta1pi,'Color','black'); % I love disgusting hacks
plot(theta2ref/pi*theta1pi./theta1pi,'Color','black');
plot( [thisone,thisone], [0,2], 'Color', 'black' );
legend('theta1/pi', 'theta2/pi','theta1ref/pi', 'theta2ref/pi', 'best iterate');
end

%% Recover the reference
step = (xy2r-xy1r)/30;
n = [-step(2);step(1)]; n = n/norm(n); % Normal
nodes3r = [ xy1r(1):step(1):xy2r(1) ; xy1r(2):step(2):xy2r(2) ];
nodes3r = nodes3r'; nnodes3r = size(nodes3r,1);
% Check if need to reverse the nodes (for visu purposes)
if norm(nodes3r(end,:)-nodes3b(1,:)) < norm(nodes3r(1,:)-nodes3b(1,:))
   nodes3r = nodes3r(end:-1:1,:);
end
nodes3s = nodes3r + 1e-3*ones(size(nodes3r,1),1)*n';
nodes3i = nodes3r - 1e-3*ones(size(nodes3r,1),1)*n';
urs = passMesh2D( nodes, elements, nodes3s, [], [u1,u2,u3,u4] );
uri = passMesh2D( nodes, elements, nodes3i, [], [u1,u2,u3,u4] );
urg = -(uri-urs)*(extnorm3(1,:)*n);  % Vectorial gap
urg([1,2,end-1,end],:) = 0; % Overwrite the strange stuff that comes from the fact that we get out of the domain
curvr = sqrt( (nodes3r(:,1)-nodes3r(1,1)).^2 + (nodes3r(:,2)-nodes3r(1,2)).^2 );

% Vizualize the crack's line
try
figure; hold on;
x1 = nodes3b(1,1);   y1 = nodes3b(1,2);
x2 = nodes3b(end,1); y2 = nodes3b(end,2);
x1r = xy1r(1); y1r = xy1r(2);
x2r = xy2r(1); y2r = xy2r(2);
plot( [xmin,xmax], [ymin,ymin], 'Color', 'black');
plot( [xmax,xmax], [ymin,ymax], 'Color', 'black');
plot( [xmax,xmin], [ymax,ymax], 'Color', 'black');
plot( [xmin,xmin], [ymax,ymin], 'Color', 'black');
plot( [x1r,x2r], [y1r,y2r], 'Color', 'black', 'LineWidth', 3 );
plot( [x1,x2], [y1,y2], 'Color', 'red', 'LineWidth', 3 );
axis equal;
end

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

% Graph for [[u]]
curv = sqrt( (nodes3b(:,1)-nodes3b(1,1)).^2 + (nodes3b(:,2)-nodes3b(1,2)).^2 );
%toplot1  = sqrt( ucrsol1(1:2:end-1).^2 + ucrsol1(2:2:end).^2 );
%toplot2  = sqrt( ucrsol2(1:2:end-1).^2 + ucrsol2(2:2:end).^2 );
%toplot3  = sqrt( ucrsol3(1:2:end-1).^2 + ucrsol3(2:2:end).^2 );
%toplot4  = sqrt( ucrsol4(1:2:end-1).^2 + ucrsol4(2:2:end).^2 );
%toplotr1 = sqrt( urg(1:2:end-1,1).^2 + urg(2:2:end,1).^2 );
%toplotr2 = sqrt( urg(1:2:end-1,2).^2 + urg(2:2:end,2).^2 );
%toplotr3 = sqrt( urg(1:2:end-1,3).^2 + urg(2:2:end,3).^2 );
%toplotr4 = sqrt( urg(1:2:end-1,4).^2 + urg(2:2:end,4).^2 );

n = extnorm3(1,:)';
toplot1  = ucrsol1(1:2:end-1)*n(1) + ucrsol1(2:2:end)*n(2);
toplot2  = ucrsol2(1:2:end-1)*n(1) + ucrsol2(2:2:end)*n(2);
toplot3  = ucrsol3(1:2:end-1)*n(1) + ucrsol3(2:2:end)*n(2);
toplot4  = ucrsol4(1:2:end-1)*n(1) + ucrsol4(2:2:end)*n(2);
%toplot4  = ucrsol4(1:2:end-1);
toplotr1 = urg(1:2:end-1,1)*n(1) + urg(2:2:end,1)*n(2);
toplotr2 = urg(1:2:end-1,2)*n(1) + urg(2:2:end,2)*n(2);
toplotr3 = urg(1:2:end-1,3)*n(1) + urg(2:2:end,3)*n(2);
toplotr4 = urg(1:2:end-1,4)*n(1) + urg(2:2:end,4)*n(2);
%toplotr4 = urg(1:2:end-1,4);


%try
%figure;
%hold on;
%set(gca, 'fontsize', 20);
%plot( curv, toplot1,'Color','green' );
%plot( curv, toplot2,'Color','black' );
%plot( curv, toplot3,'Color','blue' );
%plot( curv, toplot4,'Color','red' );
%plot( curvr, toplotr1,'Color','green' );
%plot( curvr, toplotr2,'Color','black' );
%plot( curvr, toplotr3,'Color','blue' );
%plot( curvr, toplotr4,'Color','red' );
%legend('[[u]] identified (1)', '[[u]] identified (2)', '[[u]] identified (3)', '[[u]] identified (4)',...
%       '[[u]] reference (1)', '[[u]] reference (2)', '[[u]] reference (3)', '[[u]] reference (4)');
%end

try
figure;
hold on;
set(gca, 'fontsize', 20);
plot( curv, toplot4,'Color','red','LineWidth',3 );
plot( curvr, toplotr4,'Color','blue','LineWidth',3 );
legend('[[u]] identified (4)', '[[u]] reference (4)');
end

try
figure;
plot(log10(phi));
legend('cost function');
end

%erroru = norm(usolu1(b2nodesnoD)-ur1(b2nodesnoD))   / norm(ur1(b2nodesnoD));
%errorf = norm(fsolu1no(b2nodesnoN)-fr1(b2nodesnoN)) / norm(fr1(b2nodesnoN));
