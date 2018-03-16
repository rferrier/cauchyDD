% 21/02/2018
% Recalage de RdC par écart à la réciprocité Newton

tic

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E1         = 70000;   %
E2         = 210000;  % MPa : Young modulus
E          = E1;      % Reference young modulus
nu         = 0.3;     % Poisson ratio
fscalar    = 250;     % N.mm-1 : Loading on the plate
br         = .0;      % Noise level
jmax       = 20;      % Eigenvalues truncation number
upper_term = 0;
recompute  = 1;
ordertest  = 20;
Npg        = 2;
niter      = 5;       % Nb of Newton iterations
dofull     = 1;       % Use the full matrix inversion
Cc         = [.4,.3]; % Center of the inclusion (for reference)
Rad        = .2;      % Radius of the inclusion (for reference)

mat = {}; mat{1} = [0, E1, nu]; mat{2} = [0, E2, nu];

nbDirichlet = []; % Deactivate discrete measure points

neumann1   = [3,2,fscalar ; 1,2,-fscalar];
neumann2   = [2,1,fscalar ; 4,1,-fscalar];
neumann3   = [3,1,fscalar ; 3,2,fscalar ; 2,1,fscalar ; 2,2,fscalar ; ...
              1,1,-fscalar ; 1,2,-fscalar ; 4,1,-fscalar ; 4,2,-fscalar];
neumann4   = [3,1,-fscalar ; 3,2,fscalar ; 2,1,fscalar ; 2,2,-fscalar ; ...
              1,1,fscalar ; 1,2,-fscalar ; 4,1,-fscalar ; 4,2,fscalar];

dirichlet  = [ 0,1,0 ; 0,2,0 ; 0,3,0 ];

dirichlet0 = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ; 3,1,0 ; 3,2,0 ; 4,1,0 ; 4,2,0];
neumann0   = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ; 3,1,0 ; 3,2,0 ; 4,1,0 ; 4,2,0];              
%dirichlet0 = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ; 3,1,0 ; 3,2,0 ; 4,1,0 ; 4,2,0];
%neumann0 = [];
%dirichlet0 = [ 1,1,0 ; 1,2,0 ; 3,1,0 ; 3,2,0];
%neumann0   = [ 2,1,0 ; 2,2,0 ; 4,1,0 ; 4,2,0];
%dirichlet0 = [ 1,1,0 ; 1,2,0 ; 3,1,0 ; 3,2,0 ; 4,1,0 ; 4,2,0];
%neumann0   = [ 4,1,0 ; 4,2,0];
%dirichlet0 = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ; 3,1,0 ; 3,2,0 ; 4,1,0 ; 4,2,0];
%neumann0   = [ 1,1,0 ; 1,2,0 ; 3,1,0 ; 3,2,0 ; 4,1,0 ; 4,2,0];
%dirichlet0 = [ 1,1,0 ; 1,2,0 ; 3,1,0 ; 3,2,0 ; 4,1,0 ; 4,2,0];
%neumann0   = [ 1,1,0 ; 1,2,0 ; 3,1,0 ; 3,2,0 ; 4,1,0 ; 4,2,0];
%dirichlet0 = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ; 3,1,0 ; 3,2,0 ; 4,1,0 ; 4,2,0];
%neumann0   = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0];
%dirichlet0 = [ 1,1,0 ; 1,2,0 ; 3,1,0 ; 3,2,0 ; 4,1,0 ; 4,2,0];
%neumann0   = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ; 3,1,0 ; 3,2,0 ; 4,1,0 ; 4,2,0];
%dirichlet0 = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ];
%neumann0   = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ; 3,1,0 ; 3,2,0 ; 4,1,0 ; 4,2,0];
%dirichlet0 = [ 1,1,0 ; 1,2,0 ];
%neumann0   = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ; 3,1,0 ; 3,2,0 ; 4,1,0 ; 4,2,0];

[ nodes,elements,ntoelem,boundary,order,physical] = readmesh( 'meshes/inclusion/plate_incl.msh' );
nnodes = size(nodes,1);
nelems = size(elements,1);

% mapBounds
[ node2b1, b2node1 ]   = mapBound( 1, boundary, nnodes );
[ node2b2, b2node2 ]   = mapBound( 2, boundary, nnodes );
[ node2b3, b2node3 ]   = mapBound( 3, boundary, nnodes );
[ node2b4, b2node4 ]   = mapBound( 4, boundary, nnodes );
[ node2b5, b2node5 ]   = mapBound( 5, boundary, nnodes );
[ node2b6, b2node6 ]   = mapBound( 6, boundary, nnodes );
indexbound  = [2*b2node1-1 ; 2*b2node1 ; 2*b2node2-1 ; 2*b2node2 ;...
               2*b2node3-1 ; 2*b2node3 ; 2*b2node4-1 ; 2*b2node4];
               
% Then, build the stiffness matrix :
[K,C,nbloq,node2c,c2node] = Krig2 (nodes,[physical,elements],mat,order,boundary,dirichlet);
Kinter = K( 1:2*nnodes, 1:2*nnodes );

Xmax = max(nodes(:,1)); Xmin = min(nodes(:,1)); Xmoy = (Xmax+Xmin)/2;
Ymax = max(nodes(:,2)); Ymin = min(nodes(:,2)); Ymoy = (Ymax+Ymin)/2;
Lx = Xmax-Xmin; Ly = Ymax-Ymin;

% Rigid body modes
R = rigidModes(nodes);

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

% Compute the reference sigma
mat1 = {}; mat1{1} = [0, 1e-6*E1, nu]; mat1{2} = [0, E2-E1, nu]; % Differential properties
%sigma1n = stress2(u1,nodes,[physical,elements],mat1,order,1,ntoelem);
%sigma2n = stress2(u2,nodes,[physical,elements],mat1,order,1,ntoelem);
%sigma3n = stress2(u3,nodes,[physical,elements],mat1,order,1,ntoelem);
%sigma4n = stress2(u4,nodes,[physical,elements],mat1,order,1,ntoelem);
sigma1e = stress2(u1,nodes,[physical,elements],mat1,order,0,ntoelem);
sigma2e = stress2(u2,nodes,[physical,elements],mat1,order,0,ntoelem);
sigma3e = stress2(u3,nodes,[physical,elements],mat1,order,0,ntoelem);
sigma4e = stress2(u4,nodes,[physical,elements],mat1,order,0,ntoelem);
%plotGMSH({u1,'U1';sigma1n,'S1';u2,'U2';sigma2n,'S2'; ...
%          u3,'U3';sigma3n,'S3';u4,'U4';sigma4n,'S4'}, elements, nodes, 'output/reference');
plotGMSH({u1,'U1';u2,'U2';u3,'U3';u4,'U4'}, elements, nodes, 'output/reference');
%plotGMSH({sigma1n(1:3:end-2),'S11';sigma1n(2:3:end-1),'S22';...
%          sigma1n(3:3:end),'S12'}, elements, nodes, 'output/sigref');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ nodes2,elements2,ntoelem2,boundary2,order] = readmesh( 'meshes/inclusion/plate_plain.msh' );
nnodes2 = size(nodes2,1);
%[K2,C2,nbloq2,node2c2,c2node2] = Krig2 (nodes2,elements2,mat,order,boundary2,dirichlet);
%Kinter2 = K2( 1:2*nnodes2, 1:2*nnodes2 );

%% Pass meshes
UF = [u1,u2,u3,u4,f1,f2,f3,f4];
UR = passMesh2D(nodes, elements, nodes2, elements2, UF, 0);

ur1 = UR(:,1); fr1 = UR(:,5)*2; % 2 is the ratio between mesh sizes (hack)
ur2 = UR(:,2); fr2 = UR(:,6)*2;
ur3 = UR(:,3); fr3 = UR(:,7)*2;
ur4 = UR(:,4); fr4 = UR(:,8)*2;

% Add the noise
u1n = ur1; u2n = ur2; u3n = ur3; u4n = ur4;
am1 = sqrt(mean(ur1.^2)); am2 = sqrt(mean(ur2.^2));
am3 = sqrt(mean(ur3.^2)); am4 = sqrt(mean(ur4.^2));
br1 = randn(2*nnodes2,1); br2 = randn(2*nnodes2,1);
br3 = randn(2*nnodes2,1); br4 = randn(2*nnodes2,1);
%noise = load('noises/cauchyRG.mat');
%br1 = noise.br1; br2 = noise.br2; br3 = noise.br3; br4 = noise.br4;
%u1 = ( 1 + br*br1 ) .* u1; u2 = ( 1 + br*br2 ) .* u2;
%u3 = ( 1 + br*br3 ) .* u3; u4 = ( 1 + br*br4 ) .* u4;
ur1 = ur1 + am1*br*br1; ur2 = ur2 + am1*br*br2;
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
   cand1 = rem( find(elements2==no1)-1,nelem1 )+1; % find gives line + column*size
   cand2 = rem( find(elements2==no2)-1,nelem1 )+1;
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

[ nodesu,elementsu,ntoelemu,boundaryu,orderu] = readmesh( 'meshes/inclusion/plate_plain.msh' );
nnodesu = size(nodesu,1);
[Ku,Cu,nbloqu,node2cu,c2nodeu] = Krig2 (nodesu,elementsu,[0,E,nu],orderu,boundaryu,[0,1,0;0,2,0;0,3,0]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
if upper_term == 0 % Rem : there will be lots of 0 in coef
   for i=0:nmax
      for j=0:nmax
         if i+j>nmax
            M((nmax+1)*i+j+1,:) = 0;
            M((nmax+1)*i+j+1,(nmax+1)*i+j+1) = 1;
         end
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

disp([ 'Direct problem solved and data management ', num2str(toc) ]);

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the Lhs and the Rhs
Rhs1  = zeros(nftest,1); Rhs2  = zeros(nftest,1);
Rhs3  = zeros(nftest,1); Rhs4  = zeros(nftest,1);
   
%Ru = zeros(nftest,2*nnbound2); Rf = zeros(nftest,2*nboun2);
Ruij = zeros(size(coef,1),2*nnbound2);
Rfij = zeros(size(coef,1),2*nboun2);

for i=1:nboun2
   bonod = boundary2(i,:); exno = extnorm2(i,:)';
  
   no1 = bonod(2); no2 = bonod(3); boname = bonod(1);
   x1 = nodes2(no1,1); y1 = nodes2(no1,2);
   x2 = nodes2(no2,1); y2 = nodes2(no2,2);
   len = sqrt( (x1-x2)^2 + (y1-y2)^2 );
        
   indDtot = [2*boun2doftot(i,1)-1,2*boun2doftot(i,1),...
              2*boun2doftot(i,2)-1,2*boun2doftot(i,2)]; % U dofs the element is associated to
   indNtot = [2*i-1,2*i]; % F dofs the element is associated to
      
   if order==1
      Ng = Npg; % (anyway, it's inexact)
   elseif order==2
      Ng  = Npg; no3 = bonod(4);
      x3  = nodes2(no3,1); y3 = nodes2(no3,2);
   end
   [ Xg, Wg ] = gaussPt1d( Ng );
                
   for j=1:Ng
      xg = Xg(j); wg = Wg(j);
      xgr  = ( (1-xg)*[x1;y1] + xg*[x2;y2] ); % abscissae of the Gauss point in the physical space
      X = xgr(1)/Lx; Y = xgr(2)/Lx; % Use the standard basis
            
      for ii=0:nmax
         nimax = nmax-ii;
         for jj=0:nimax
            sloc11a = 0; sloc12a = 0; sloc22a = 0;
            sloc11b = 0; sloc12b = 0; sloc22b = 0;
            if ii>0
               XY = ii*X^(ii-1)*Y^jj;
               sloc11a = E/(1-nu^2)*XY;
               sloc22a = nu*E/(1-nu^2)*XY;
               sloc12b = E/(2*(1+nu))*XY;
            end
            if jj>0
               XY = jj*X^ii*Y^(jj-1);
               sloc11b = nu*E/(1-nu^2)*XY;
               sloc22b = E/(1-nu^2)*XY;
               sloc12a = E/(2*(1+nu))*XY;
            end

            sloca = 1/Lx*[sloc11a,sloc12a;sloc12a,sloc22a];
            slocb = 1/Lx*[sloc11b,sloc12b;sloc12b,sloc22b];
            spaa = sloca; spab = slocb;
            fpaa = spaa*exno; fpab = spab*exno;

            vpaa = [ X^ii*Y^jj ; 0 ];
            vpab = [ 0 ; X^ii*Y^jj ];
               
            fpaTimesPhitest0a = [fpaa(1)*(1-xg);fpaa(2)*(1-xg);fpaa(1)*xg;fpaa(2)*xg];
            fpaTimesPhitest0b = [fpab(1)*(1-xg);fpab(2)*(1-xg);fpab(1)*xg;fpab(2)*xg]; 
              
            Ruij((nmax+1)*ii + jj+1,indDtot) = Ruij((nmax+1)*ii + jj+1,indDtot) ...
                                                + len*wg * fpaTimesPhitest0a';
            Ruij((nmax+1)^2 + (nmax+1)*ii + jj+1,indDtot) =...
                                   Ruij((nmax+1)^2 + (nmax+1)*ii + jj+1,indDtot) ...
                                                + len*wg * fpaTimesPhitest0b';
            Rfij((nmax+1)*ii + jj+1,indNtot) = Rfij((nmax+1)*ii + jj+1,indNtot) ...
                                               + len*wg * vpaa';
            Rfij((nmax+1)^2 + (nmax+1)*ii + jj+1,indNtot) = ...
                          Rfij((nmax+1)^2 + (nmax+1)*ii + jj+1,indNtot) + len*wg * vpab';
         end
      end
   end
end
Ru = coef'*Ruij; Rf = coef'*Rfij;

%%
% Transform the given data in the proper format
f_known1 = zeros(2*nboun2,1); u_known1 = zeros(2*nnbound2,1);
f_known2 = zeros(2*nboun2,1); u_known2 = zeros(2*nnbound2,1);
f_known3 = zeros(2*nboun2,1); u_known3 = zeros(2*nnbound2,1);
f_known4 = zeros(2*nboun2,1); u_known4 = zeros(2*nnbound2,1);
   
for i=1:nboun2
   bonod = boundary2(i,:); exno = extnorm2(i,:)';
   % Reference force from the BC's
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

%% Restrict the data on the given part (not necessary, but cleaner)
f_known1(tofindN) = 0; f_known2(tofindN) = 0;
f_known3(tofindN) = 0; f_known4(tofindN) = 0;
u_known1(tofindD) = 0; u_known2(tofindD) = 0;
u_known3(tofindD) = 0; u_known4(tofindD) = 0;

%% Restrict the operators
Rfm = Rf(:,tofindN); Rfr = Rf(:,knownN);
Rum = Ru(:,tofindD); Rur = Ru(:,knownD);

LhsA = -Rum;
LhsB = Rfm;
Rhs1 = - ( Rur*u_known1(knownD) - Rfr*f_known1(knownN) );
Rhs2 = - ( Rur*u_known2(knownD) - Rfr*f_known2(knownN) );
Rhs3 = - ( Rur*u_known3(knownD) - Rfr*f_known3(knownN) );
Rhs4 = - ( Rur*u_known4(knownD) - Rfr*f_known4(knownN) );

disp([ 'Right hand side generated ', num2str(toc) ]);

% Start algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the differential values
k12 = zeros(nelemu,1);
k33 = zeros(nelemu,1);

% Initialize Uu by direct resolutions (TOSEE : is it a good idea ?)
Ur1 = Ku\[fr1;zeros(3,1)]; Ur2 = Ku\[fr2;zeros(3,1)];
Ur3 = Ku\[fr3;zeros(3,1)]; Ur4 = Ku\[fr4;zeros(3,1)];
Uu1 = Ur1(1:2*nnodesu); Uu2 = Ur2(1:2*nnodesu);
Uu3 = Ur3(1:2*nnodesu); Uu4 = Ur4(1:2*nnodesu);

nm2 = (nmax+1)^2;

Ke12 = [1,1,0;1,1,0;0,0,0]; % Elemental derivative stiffnesses
Ke33 = [2,0,0;0,2,0;0,0,1];

% Ponderated preconditionner
Lu  = zeros( 2*nnodesu );
L12 = zeros( nelemu );
L33 = zeros( nelemu );
Zue = zeros( 2*nnodesu, nelemu );
Zuu = zeros( 2*nnodesu, 2*nnodesu );
Zee = zeros(nelemu);
for i=1:nelemu
   bonod = elementsu(i,:);

   no1 = bonod(1); no2 = bonod(2); no3 = bonod(3);
   x1 = nodesu(no1,1); y1 = nodesu(no1,2);
   x2 = nodesu(no2,1); y2 = nodesu(no2,2);
   x3 = nodesu(no3,1); y3 = nodesu(no3,2);

   S = .5 * ((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
   ind = [ 2*no1-1, 2*no1, 2*no2-1, 2*no2, 2*no3-1, 2*no3 ];
   Lu(ind,ind) = S/3*eye(6);
   L12(i,i) = S;
   L33(i,i) = S;
end

Res = zeros(niter,1); Deltau = zeros(niter,1);
Deltak12 = zeros(niter,1); Deltak33 = zeros(niter,1);

for iter = 1:niter

   tic
   Lhsij    = zeros( size(coef,1), 2*nnodesu ); % Decomposition functions are constant per element
   Lhsij121 = sparse( size(coef,1), nelemu );   % Derivative of the operator
   Lhsij331 = sparse( size(coef,1), nelemu );   %
   Lhsij122 = sparse( size(coef,1), nelemu );   %
   Lhsij332 = sparse( size(coef,1), nelemu );   %
   Lhsij123 = sparse( size(coef,1), nelemu );   %
   Lhsij333 = sparse( size(coef,1), nelemu );   %
   Lhsij124 = sparse( size(coef,1), nelemu );   %
   Lhsij334 = sparse( size(coef,1), nelemu );   %

   Ng = Npg; [ Xg, Wg ] = gaussPt( Ng ); sWg = size(Wg,1); % Could even be outside the iter loop

   %% Build the list of Gauss points, and of construction-functions
   Xxg = zeros( nelemu*sWg,1 ); Yyg = zeros( nelemu*sWg,1 ); Wwg = zeros( nelemu*sWg,1 );
   Phi = sparse( 3*nelemu*sWg, 2*nnodesu );

   Phi121 = sparse( 3*nelemu*sWg, nelemu );
   Phi122 = sparse( 3*nelemu*sWg, nelemu );
   Phi123 = sparse( 3*nelemu*sWg, nelemu );
   Phi124 = sparse( 3*nelemu*sWg, nelemu );

   Phi331 = sparse( 3*nelemu*sWg, nelemu );
   Phi332 = sparse( 3*nelemu*sWg, nelemu );
   Phi333 = sparse( 3*nelemu*sWg, nelemu );
   Phi334 = sparse( 3*nelemu*sWg, nelemu );

   for j=1:nelemu % Compute the integrals
      bonod = elementsu(j,:);

      no1 = bonod(1); no2 = bonod(2); no3 = bonod(3);
      x1 = nodesu(no1,1); y1 = nodesu(no1,2);
      x2 = nodesu(no2,1); y2 = nodesu(no2,2);
      x3 = nodesu(no3,1); y3 = nodesu(no3,2);

      S = abs(.5 * ((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)) );
   
      Be = [y2-y3,0,y3-y1,0,y1-y2,0; % Be*u = epsilon
            0,x3-x2,0,x1-x3,0,x2-x1;
            x3-x2,y2-y3,x1-x3,y3-y1,x2-x1,y1-y2]/(2*S);

      Ke   = [ 2*k33(j)+k12(j), k12(j), 0 ;... % Elemental stiffness matrix
               k12(j), 2*k33(j)+k12(j), 0 ;...
               0, 0, k33(j) ];

      KeBe   = Ke*Be; % Gives sigma from u
      KeBe12 = Ke12*Be;
      KeBe33 = Ke33*Be;
      
      ind  = [ 2*no1-1, 2*no1, 2*no2-1, 2*no2, 2*no3-1, 2*no3 ];
      uloc1 = Uu1(ind); uloc2 = Uu2(ind);
      uloc3 = Uu3(ind); uloc4 = Uu4(ind);

      K12u1 = KeBe12 * uloc1; K33u1 = KeBe33 * uloc1;
      K12u2 = KeBe12 * uloc2; K33u2 = KeBe33 * uloc2;
      K12u3 = KeBe12 * uloc3; K33u3 = KeBe33 * uloc3;
      K12u4 = KeBe12 * uloc4; K33u4 = KeBe33 * uloc4;

      for k=1:sWg
         xg = Xg(k,:); wg = Wg(k);
         swg = S*wg;
         xgr  = (1-xg(1)-xg(2))*[x1;y1] + xg(1)*[x2;y2] + xg(2)*[x3;y3] ; % abscissae
         X = xgr(1)/Lx; Y = xgr(2)/Lx;

         indexg = sWg*(j-1) + k;
         Xxg( indexg ) = X; Yyg( indexg ) = Y; Wwg( indexg ) = swg;
         Phi( [3*indexg-2,3*indexg-1,3*indexg], ind )    = KeBe;

         Phi121( [3*indexg-2,3*indexg-1,3*indexg], j ) = K12u1;
         Phi122( [3*indexg-2,3*indexg-1,3*indexg], j ) = K12u2;
         Phi123( [3*indexg-2,3*indexg-1,3*indexg], j ) = K12u3;
         Phi124( [3*indexg-2,3*indexg-1,3*indexg], j ) = K12u4;

         Phi331( [3*indexg-2,3*indexg-1,3*indexg], j ) = K33u1;
         Phi332( [3*indexg-2,3*indexg-1,3*indexg], j ) = K33u2;
         Phi333( [3*indexg-2,3*indexg-1,3*indexg], j ) = K33u3;
         Phi334( [3*indexg-2,3*indexg-1,3*indexg], j ) = K33u4;
      end
   end

   %% Build the matrix of test-functions
   Vv = zeros( size(coef,1), 3*nelemu*sWg );
   for ii=0:nmax
      nimax = nmax-ii;
      for jj=0:nimax
         Exxa = zeros( nelemu*sWg, 1 );
         Exya = zeros( nelemu*sWg, 1 );
         Eyyb = zeros( nelemu*sWg, 1 );
         Exyb = zeros( nelemu*sWg, 1 );

         if ii>0 % epsilon_ij on the surfacic element
            Exxa = ii*Xxg.^(ii-1).*Yyg.^jj;
            Exyb = .5 * ii*Xxg.^(ii-1).*Yyg.^jj;
         end
         if jj>0
            Eyyb = jj*Xxg.^ii.*Yyg.^(jj-1);
            Exya = .5 * jj*Xxg.^ii.*Yyg.^(jj-1);
         end

%         Ea = [ Exxa ; 0 ; Exya ];
%         Eb = [ 0 ; Eyyb ; Exyb ];

         inda = (nmax+1)*ii + jj+1;
         indb = nm2 + inda;

%         sEa = swg * Ea';
%         sEb = swg * Eb';

         Vv(inda,1:3:3*nelemu*sWg-2) = Wwg .* Exxa;
         %Vv(inda,2:3:3*nelemu*sWg-1) = Wwg .* Eyya;
         Vv(inda,3:3:3*nelemu*sWg)   = Wwg .* Exya;
         %Vv(indb,1:3:3*nelemu*sWg-2) = Wwg .* Exxb;
         Vv(indb,2:3:3*nelemu*sWg-1) = Wwg .* Eyyb;
         Vv(indb,3:3:3*nelemu*sWg)   = Wwg .* Exyb;
      end
   end
   %%

   Lhsijp    = Vv*Phi;
   Lhsij121p = Vv*Phi121; Lhsij331p = Vv*Phi331;
   Lhsij122p = Vv*Phi122; Lhsij332p = Vv*Phi332;
   Lhsij123p = Vv*Phi123; Lhsij333p = Vv*Phi333;
   Lhsij124p = Vv*Phi124; Lhsij334p = Vv*Phi334;

%   toc
%   tic

%   for j=1:nelemu % Compute the integrals
%      bonod = elementsu(j,:);

%      no1 = bonod(1); no2 = bonod(2); no3 = bonod(3);
%      x1 = nodesu(no1,1); y1 = nodesu(no1,2);
%      x2 = nodesu(no2,1); y2 = nodesu(no2,2);
%      x3 = nodesu(no3,1); y3 = nodesu(no3,2);

%      S = abs(.5 * ((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)) );
%   
%      Be = [y2-y3,0,y3-y1,0,y1-y2,0; % Be*u = epsilon
%            0,x3-x2,0,x1-x3,0,x2-x1;
%            x3-x2,y2-y3,x1-x3,y3-y1,x2-x1,y1-y2]/(2*S);

%      Ke   = [ 2*k33(j)+k12(j), k12(j), 0 ;... % Elemental stiffness matrix
%               k12(j), 2*k33(j)+k12(j), 0 ;...
%               0, 0, k33(j) ];

%      KeBe   = Ke*Be; % Gives sigma from u
%      KeBe12 = Ke12*Be;
%      KeBe33 = Ke33*Be;
%      
%      ind  = [ 2*no1-1, 2*no1, 2*no2-1, 2*no2, 2*no3-1, 2*no3 ];
%      uloc1 = Uu1(ind); uloc2 = Uu2(ind);
%      uloc3 = Uu3(ind); uloc4 = Uu4(ind);

%      K12u1 = KeBe12 * uloc1; K33u1 = KeBe33 * uloc1;
%      K12u2 = KeBe12 * uloc2; K33u2 = KeBe33 * uloc2;
%      K12u3 = KeBe12 * uloc3; K33u3 = KeBe33 * uloc3;
%      K12u4 = KeBe12 * uloc4; K33u4 = KeBe33 * uloc4;

%      for k=1:sWg
%         xg = Xg(k,:); wg = Wg(k);
%         swg = S*wg;
%         xgr  = (1-xg(1)-xg(2))*[x1;y1] + xg(1)*[x2;y2] + xg(2)*[x3;y3] ; % abscissae
%         X = xgr(1)/Lx; Y = xgr(2)/Lx;
%      
%         for ii=0:nmax
%            nimax = nmax-ii;
%            for jj=0:nimax
%               Exxa = 0; Exya = 0;
%               Eyyb = 0; Exyb = 0;
%               if ii>0 % integral of epsilon_ij on the surfacic element
%                  Exxa = ii*X^(ii-1)*Y^jj;
%                  Exyb = .5 * ii*X^(ii-1)*Y^jj;
%               end
%               if jj>0
%                  Eyyb = jj*X^ii*Y^(jj-1);
%                  Exya = .5 * jj*X^ii*Y^(jj-1);
%               end

%               Ea = [ Exxa ; 0 ; Exya ];
%               Eb = [ 0 ; Eyyb ; Exyb ];

%               inda = (nmax+1)*ii + jj+1;
%               indb = nm2 + inda;

%               sEa = swg * Ea';
%               sEb = swg * Eb';

%               Lhsij(inda,ind) = Lhsij(inda,ind) + sEa * KeBe;
%               Lhsij(indb,ind) = Lhsij(indb,ind) + sEb * KeBe;

%               Lhsij121(inda,j) = Lhsij121(inda,j) + sEa * K12u1;
%               Lhsij121(indb,j) = Lhsij121(indb,j) + sEb * K12u1;
%               Lhsij331(inda,j) = Lhsij331(inda,j) + sEa * K33u1;
%               Lhsij331(indb,j) = Lhsij331(indb,j) + sEb * K33u1;

%               Lhsij122(inda,j) = Lhsij122(inda,j) + sEa * K12u2;
%               Lhsij122(indb,j) = Lhsij122(indb,j) + sEb * K12u2;
%               Lhsij332(inda,j) = Lhsij332(inda,j) + sEa * K33u2;
%               Lhsij332(indb,j) = Lhsij332(indb,j) + sEb * K33u2;

%               Lhsij123(inda,j) = Lhsij123(inda,j) + sEa * K12u3;
%               Lhsij123(indb,j) = Lhsij123(indb,j) + sEb * K12u3;
%               Lhsij333(inda,j) = Lhsij333(inda,j) + sEa * K33u3;
%               Lhsij333(indb,j) = Lhsij333(indb,j) + sEb * K33u3;

%               Lhsij124(inda,j) = Lhsij124(inda,j) + sEa * K12u4;
%               Lhsij124(indb,j) = Lhsij124(indb,j) + sEb * K12u4;
%               Lhsij334(inda,j) = Lhsij334(inda,j) + sEa * K33u4;
%               Lhsij334(indb,j) = Lhsij334(indb,j) + sEb * K33u4;
%            end
%         end
%      end
%   end
%   Lhs0   = coef'*Lhsij;
%   Lhs121 = coef'*Lhsij121; Lhs331 = coef'*Lhsij331;
%   Lhs122 = coef'*Lhsij122; Lhs332 = coef'*Lhsij332;
%   Lhs123 = coef'*Lhsij123; Lhs333 = coef'*Lhsij333;
%   Lhs124 = coef'*Lhsij124; Lhs334 = coef'*Lhsij334;

   Lhs0   = coef'*Lhsijp;
   Lhs121 = coef'*Lhsij121p; Lhs331 = coef'*Lhsij331p;
   Lhs122 = coef'*Lhsij122p; Lhs332 = coef'*Lhsij332p;
   Lhs123 = coef'*Lhsij123p; Lhs333 = coef'*Lhsij333p;
   Lhs124 = coef'*Lhsij124p; Lhs334 = coef'*Lhsij334p;

   disp([ 'Left hand side generated ', num2str(toc) ]);

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% Solve the linear system
   if norm(Lhs0,'fro') > 0
      K121 = norm(Lhs0,'fro')/norm(Lhs121,'fro'); % Weight operators
      K331 = norm(Lhs0,'fro')/norm(Lhs331,'fro'); %
      K122 = norm(Lhs0,'fro')/norm(Lhs122,'fro'); %
      K332 = norm(Lhs0,'fro')/norm(Lhs332,'fro'); %
      K123 = norm(Lhs0,'fro')/norm(Lhs123,'fro'); %
      K333 = norm(Lhs0,'fro')/norm(Lhs333,'fro'); %
      K124 = norm(Lhs0,'fro')/norm(Lhs124,'fro'); %
      K334 = norm(Lhs0,'fro')/norm(Lhs334,'fro'); %
   else
      K121 = 1; K122 = 1; K123 = 1; K124 = 1;
      K331 = 1; K332 = 1; K333 = 1; K334 = 1;
   end

   %Lhs = [ Lhs0, Lhs0, Lhs0, Lhs0, K12*Lhs12, K33*Lhs33 ]; % Order is [u1;u2;u3;u4;K12*k12;K33*k33]
   tic;
   Z0 = zeros(size(Lhs0));

   if dofull == 1
      Lhs = [ Lhs0, Z0, Z0, Z0, K121*Lhs121, K331*Lhs331 ;...
              Z0, Lhs0, Z0, Z0, K121*Lhs122, K331*Lhs332 ;...
              Z0, Z0, Lhs0, Z0, K121*Lhs123, K331*Lhs333 ;... %/!\ There is the same K on purpose, on order to reconstruct k12 and k33
            Z0, Z0, Z0, Lhs0, K121*Lhs124, K331*Lhs334 ];
      Rhs = [ Rhs1 - Lhs0*Uu1 ; Rhs2 - Lhs0*Uu2 ; Rhs3 - Lhs0*Uu3 ; Rhs4 - Lhs0*Uu4 ];
      Res(iter) = norm(Rhs);
   
      A = Lhs'*Lhs; sA = size(A,1);
      L = [ Lu, Zuu, Zuu, Zuu, Zue, Zue ; ...
            Zuu, Lu, Zuu, Zuu, Zue, Zue ; ...
            Zuu, Zuu, Lu, Zuu, Zue, Zue ; ...
            Zuu, Zuu, Zuu, Lu, Zue, Zue ; ...
            Zue',Zue',Zue',Zue', L12, Zee ; ...
            Zue',Zue',Zue',Zue', Zee, L33 ];
      Msmall = [ Zuu, Zuu, Zuu, Zuu, Zue, Zue ; ...
                 Zuu, Zuu, Zuu, Zuu, Zue, Zue ; ...
                 Zuu, Zuu, Zuu, Zuu, Zue, Zue ; ...
                 Zuu, Zuu, Zuu, Zuu, Zue, Zue ; ...
                 Zue',Zue',Zue',Zue', L12, Zee ; ...
                 Zue',Zue',Zue',Zue', Zee, L33 ];

      [Q,Theta] = eig(A,L); Q = Q*real((Q'*L*Q)^(-1/2)); % real should not be, but you know, numerical shit...
      thetas = diag(Theta);
      [thetas,Ind] = sort( thetas,'descend' );
      Q = Q(:,Ind);
      Thetas = diag(thetas); 
      Theta = Theta(Ind,Ind);
   
      % Plot the Picard stuff
      imax = min( find(thetas/thetas(1)<1e-16) );
      if size(imax,1) == 0
         imax = size(thetas,1);
      end
   
      tplo = thetas(1:imax);
      bplo = Q'*Lhs'*Rhs; bplo = bplo(1:imax);
      rplp = (Q'*Lhs'*Rhs)./thetas; rplo = rplp(1:imax);

      % Remove Zeros in rploi (why on shit are there zeros in the first place ?)
      me = mean(abs(rplp))/1e5; arplo = max(me,abs(rplp));
      ind1 = findPicard2 (log10(arplo(1:imax)), ceil(imax/7), 1, 3);

%   try
%   figure
%   hold on;
%   plot(log10(abs(tplo)),'Color','green');
%   plot(log10(abs(bplo)),'Color','red');
%   plot(log10(abs(rplo)),'Color','black');
%   legend('Singular values','Rhs','sol');
%   end

%   % Filter eigenvalues
%   if jmax == 0
%      jmax = size(Thetas,1);
%    jmax1 = ind1;
%   else
%      jmax0 = min( size(Thetas,1) , jmax );
%      jmax1 = jmax;
%   end

      if iter == 1
         jmax1 = jmax;
      else
         jmax1 = 100;
      end

      nor = zeros(sA,1); res = zeros(sA,1);
      for i=1:sA
         xnor = Q(:,1:i)*rplp(1:i); nor(i) = xnor'*Msmall*xnor;
         res(i)  = norm(Lhs*xnor - Rhs);
      end
      nares = find(isnan(log10(res))); nanor = find(isnan(log10(nor))); % Manage the nans
      natot = union(nares,nanor); res(natot) = []; nor(nanor) = []; % We assume the nans are the last ones, so no impact on the numeroataion
      indd = findCorner (res, nor, ceil(sA/4), 1, 7);

      try
      figure;
      hold on;
      loglog(res,nor,'Color','red','-*');
      loglog(res(indd),nor(indd),'Color','red','o','markersize',15);
      legend('L-curve','');
      end

      ThetaT = Thetas( 1:jmax1 , 1:jmax1 );
      bU = Q'*Lhs'*Rhs; bT = bU(1:jmax1);
      Solu = Q(:,1:jmax1) * (ThetaT\bT);

      Uu1 = Uu1 + Solu(1:2*nnodesu);
      Uu2 = Uu2 + Solu(2*nnodesu+1:4*nnodesu);
      Uu3 = Uu3 + Solu(4*nnodesu+1:6*nnodesu);
      Uu4 = Uu4 + Solu(6*nnodesu+1:8*nnodesu);

      nnu8 = 8*nnodesu;
      k12 = k12 + K121*Solu(nnu8+1:nnu8+nelemu);
      k33 = k33 + K331*Solu(nnu8+nelemu+1:end);
      disp([ 'System pinverted ', num2str(toc) ]);

      Deltau(iter)   = norm(Solu(1:nnu8));
      Deltak12(iter) = norm(K121*Solu(nnu8+1:nnu8+nelemu));
      Deltak33(iter) = norm(K331*Solu(nnu8+nelemu+1:end));

   else
      %% Separate resolutions : first the Uu (if the norm of Lhs0 is not zero)
      if norm(Lhs0,'fro') ~= 0
         Lhsu  = Lhs0;
         Rhs11 = [ Rhs1 ];% - Lhs0*Uu1 ];% - Lhs121*k12 - Lhs331*k33 ];
         Rhs22 = [ Rhs2 ];% - Lhs0*Uu2 ];% - Lhs122*k12 - Lhs332*k33 ];
         Rhs33 = [ Rhs3 ];% - Lhs0*Uu3 ];% - Lhs123*k12 - Lhs333*k33 ];
         Rhs44 = [ Rhs4 ];% - Lhs0*Uu4 ];% - Lhs124*k12 - Lhs334*k33 ];

         A = Lhsu'*Lhsu; L = Lu;

         [Q,Theta] = eig(A,L); Q = Q*real((Q'*L*Q)^(-1/2));
         thetas = diag(Theta);
         [thetas,Ind] = sort( thetas,'descend' );
         Q = Q(:,Ind);
         Thetas = diag(thetas); 
         Theta = Theta(Ind,Ind);

         % Plot the Picard stuff
         imax = min( find(thetas/thetas(1)<1e-16) );
         if size(imax,1) == 0
            imax = size(thetas,1);
         end
   
         tplo = thetas(1:imax);
         bplo1 = Q'*Lhsu'*Rhs11; bplo1 = bplo1(1:imax);
         rplp1 = (Q'*Lhsu'*Rhs11)./thetas; rplo1 = rplp1(1:imax);
         bplo2 = Q'*Lhsu'*Rhs22; bplo2 = bplo2(1:imax);
         rplp2 = (Q'*Lhsu'*Rhs22)./thetas; rplo2 = rplp2(1:imax);
         bplo3 = Q'*Lhsu'*Rhs33; bplo3 = bplo3(1:imax);
         rplp3 = (Q'*Lhsu'*Rhs33)./thetas; rplo3 = rplp3(1:imax);
         bplo4 = Q'*Lhsu'*Rhs44; bplo4 = bplo4(1:imax);
         rplp4 = (Q'*Lhsu'*Rhs44)./thetas; rplo4 = rplp4(1:imax);

         try
         figure
         hold on;
         plot(log10(abs(rplo1)),'Color','red');
         plot(log10(abs(rplo2)),'Color','black');
         plot(log10(abs(rplo3)),'Color','magenta');
         plot(log10(abs(rplo4)),'Color','blue');
         legend('sol1','sol2','sol3','sol4');
         end

         % Remove Zeros in rploi (why on shit are there zeros in the first place ?)
         me1 = mean(abs(rplp1))/1e5; arplo1 = max(me1,abs(rplp1));
         me2 = mean(abs(rplp2))/1e5; arplo2 = max(me2,abs(rplp2));
         me3 = mean(abs(rplp3))/1e5; arplo3 = max(me3,abs(rplp3));
         me4 = mean(abs(rplp4))/1e5; arplo4 = max(me4,abs(rplp4));
         ind1 = findPicard2 (log10(arplo1(1:imax)), ceil(imax/7), 1, 3);
         ind2 = findPicard2 (log10(arplo2(1:imax)), ceil(imax/7), 1, 3);
         ind3 = findPicard2 (log10(arplo3(1:imax)), ceil(imax/7), 1, 3);
         ind4 = findPicard2 (log10(arplo4(1:imax)), ceil(imax/7), 1, 3);

         % Filter eigenvalues
         if jmax == 0
            jmax = size(Thetas,1);
            jmax1 = ind1; jmax2 = ind2; jmax3 = ind3; jmax4 = ind4;
         else
            jmax0 = min( size(Thetas,1) , jmax );
            jmax1 = jmax; jmax2 = jmax; jmax3 = jmax; jmax4 = jmax;
         end

         %jmax1 = 45; jmax2 = 45; jmax3 = 45; jmax4 = 45;

         ThetaT1 = Thetas( 1:jmax1 , 1:jmax1 );
         ThetaT2 = Thetas( 1:jmax2 , 1:jmax2 );
         ThetaT3 = Thetas( 1:jmax3 , 1:jmax3 );
         ThetaT4 = Thetas( 1:jmax4 , 1:jmax4 );
   
         bU1 = Q'*Lhsu'*Rhs11; bT1 = bU1(1:jmax1);
         bU2 = Q'*Lhsu'*Rhs22; bT2 = bU2(1:jmax2);
         bU3 = Q'*Lhsu'*Rhs33; bT3 = bU3(1:jmax3);
         bU4 = Q'*Lhsu'*Rhs44; bT4 = bU4(1:jmax4);
  
         Solu1 = Q(:,1:jmax1) * (ThetaT1\bT1);
         Solu2 = Q(:,1:jmax2) * (ThetaT2\bT2);
         Solu3 = Q(:,1:jmax3) * (ThetaT3\bT3);
         Solu4 = Q(:,1:jmax4) * (ThetaT4\bT4);

%         Uu1 = Uu1 + Solu1;
%         Uu2 = Uu2 + Solu2;
%         Uu3 = Uu3 + Solu3;
%         Uu4 = Uu4 + Solu4;

         Deltau(iter) = norm( [Solu1-Uu1;Solu2-Uu2;Solu3-Uu3;Solu4-Uu4] );

         Uu1 = Solu1; Uu2 = Solu2; Uu3 = Solu3; Uu4 = Solu4;

      end

      %% Then, solve for k12 and k33
      Lhs = [ Lhs121, Lhs331 ; Lhs122, Lhs332 ; Lhs123, Lhs333 ; Lhs124, Lhs334 ];
      Rhs = [ Rhs1 - Lhs0*Uu1; Rhs2 - Lhs0*Uu2; Rhs3 - Lhs0*Uu3 ; Rhs4 - Lhs0*Uu4 ];
      Res(iter) = norm(Rhs);
   
      A = Lhs'*Lhs; sA = size(A,1);
      L = [ L12, Zee ; Zee, L33 ];

      [Q,Theta] = eig(A,L); Q = Q*real((Q'*L*Q)^(-1/2)); % real should not be, but you know, numerical shit...
      thetas = diag(Theta);
      [thetas,Ind] = sort( thetas,'descend' );
      Q = Q(:,Ind);
      Thetas = diag(thetas); 
      Theta = Theta(Ind,Ind);
   
      % Plot the Picard stuff
      imax = min( find(thetas/thetas(1)<1e-16) );
      if size(imax,1) == 0
         imax = size(thetas,1);
      end
   
      tplo = thetas(1:imax);
      bplo = Q'*Lhs'*Rhs; bplo = bplo(1:imax);
      rplp = (Q'*Lhs'*Rhs)./thetas; rplo = rplp(1:imax);

      % Remove Zeros in rploi (why on shit are there zeros in the first place ?)
      me = mean(abs(rplp))/1e5; arplo = max(me,abs(rplp));
      ind1 = findPicard2 (log10(arplo(1:imax)), ceil(imax/7), 1, 3);

      try
      figure
      hold on;
      plot(log10(abs(tplo)),'Color','green');
      plot(log10(abs(bplo)),'Color','red');
      plot(log10(abs(rplo)),'Color','black');
      legend('Singular values','Rhs','sol');
      end

      % Filter eigenvalues
      if jmax == 0
         jmax = size(Thetas,1);
         jmax1 = ind1;
      else
         jmax0 = min( size(Thetas,1) , jmax );
         jmax1 = jmax;
      end

%      if iter == 1
%         jmax1 = jmax;
%      else
%         jmax1 = 100;
%      end

      ThetaT = Thetas( 1:jmax1 , 1:jmax1 );
      bU = Q'*Lhs'*Rhs; bT = bU(1:jmax1);
      Solu = Q(:,1:jmax1) * (ThetaT\bT);

      k12 = k12 + Solu(1:nelemu);
      k33 = k33 + Solu(nelemu+1:end);
      disp([ 'System pinverted ', num2str(toc) ]);

      Deltak12(iter) = norm(Solu(1:nelemu));
      Deltak33(iter) = norm(Solu(nelemu+1:end));
   end

end

lambda = 2*(k33.*k12)./(2*k33-k12);
mu     = k33;
Ed     = mu.*(3*lambda+2*mu)./(lambda+mu);
nud    = lambda./(2*lambda+mu);

try
figure;
hold on;
plot(log10(Res/Res(1)),'Color','black');
plot(2:niter,log10(Deltau(2:end)/Deltau(2)),'Color','blue');
plot(log10(Deltak12/Deltak12(1)),'Color','red');
plot(log10(Deltak33/Deltak33(1)),'Color','green');
legend('residual', '\Delta u', '\Delta k_{12}', '\Delta k_{33}');
end

try
figure;
hold on
patch('Faces',elementsu(:,1:3),'Vertices',nodesu,'FaceVertexCData',Ed+E,'FaceColor','flat');
teta = 0:.1:2*pi+.1; ixe = Rad*cos(teta)+Cc(1); igrec = Rad*sin(teta)+Cc(2);
plot( ixe, igrec, 'Color', 'black',  'LineWidth', 3 );
colorbar; axis equal;
caxis([70000,210000]); legend('E');
end

try
figure;
hold on
patch('Faces',elementsu(:,1:3),'Vertices',nodesu,'FaceVertexCData',nud+nu,'FaceColor','flat');
teta = 0:.1:2*pi+.1; ixe = Rad*cos(teta)+Cc(1); igrec = Rad*sin(teta)+Cc(2);
plot( ixe, igrec, 'Color', 'black',  'LineWidth', 3 );
colorbar; axis equal;
caxis([0,1]); legend('\nu');
end

try
figure;
hold on
patch('Faces',elementsu(:,1:3),'Vertices',nodesu,'FaceVertexCData',Uu1(1:2:end-1),'FaceColor','interp');
teta = 0:.1:2*pi+.1; ixe = Rad*cos(teta)+Cc(1); igrec = Rad*sin(teta)+Cc(2);
plot( ixe, igrec, 'Color', 'black',  'LineWidth', 3 );
colorbar; axis equal;
legend('u_1^x');
end
