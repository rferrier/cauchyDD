% 08/01/2018
% Problèmes directs et de Cauchy par écart à la réciprocité, ajout d'une fissure

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
jmaxRG      = 0;     % Eigenvalues truncation number
regular     = 0;      % Use the derivative regularization matrix (0 : Id, 1 : derivative)
upper_term  = 0;      % 1 : use i=0:10, j=0:10, 0 : use i>=0,j>=0,i+j<=10
froreg      = 1;      % frobenius preconditioner
recompute   = 0;      % Recompute the operators
RGorSP      = 1;      % Deactivated
theta1      = pi;     %3.7296;%pi;%pi/2;
theta2      = 0;      %0.58800;%0%3*pi/2; % Initial angles of the crack
anglestep   = pi/15;  % Step in angle for Markov search
limnotwork  = 1001;   % Nb of failures to start adapting the step
nbstep      = 100;
Npg         = 2;      % Nb Gauss points
ordertest   = 20;     % Order of test fonctions
paramchoice = 3;      % Choice of the regularization parameter : 1=L-curve, 2=Picard, 3=DL-curve, 4=SL-curve

weights = [0,1,0,0]; % Weights of the different regularization methods (order is the same as in paramchoice)

nbDirichlet = [];
%nbDirichlet = [ 1,10 ; 2,11 ; 3,11 ; 4,11 ];
%nbDirichlet = [ 1,5 ; 2,5 ; 3,5 ; 4,5 ]; % Nb of displacement measure points on the boundaries (0=all, /!\ never go above the nb of dofs)
%nbDirichlet = [ 1,6 ; 2,6 ; 3,6 ; 4,6 ];

% Boundary conditions
dirichlet  = [0,1,0 ; 0,2,0 ; 0,3,0];
neumann1   = [3,2,fscalar ; 1,2,-fscalar];
neumann2   = [2,1,fscalar ; 4,1,-fscalar];
neumann3   = [3,1,fscalar ; 3,2,fscalar ; 2,1,fscalar ; 2,2,fscalar ; ...
              1,1,-fscalar ; 1,2,-fscalar ; 4,1,-fscalar ; 4,2,-fscalar];
neumann4   = [3,1,-fscalar ; 3,2,fscalar ; 2,1,fscalar ; 2,2,-fscalar ; ...
              1,1,fscalar ; 1,2,-fscalar ; 4,1,-fscalar ; 4,2,fscalar];

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

% First, import the mesh
[ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_c_squared2.msh' );
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

%% Find the reference angles of the crack
x1 = nodes(5,1); y1 = nodes(5,2);
x2 = nodes(6,1); y2 = nodes(6,2);

xmin = min(nodes(:,1)); xmax = max(nodes(:,1));
ymin = min(nodes(:,2)); ymax = max(nodes(:,2));
xb = .5*(xmin+xmax); yb = .5*(ymin+ymax);
ma11 = [ 0, -(xb-x1) ; ymin-ymax, -(yb-y1) ]; bh11 = [ x1-xmax ; y1-ymax];
ma12 = [ 0, -(xb-x2) ; ymin-ymax, -(yb-y2) ]; bh12 = [ x2-xmax ; y2-ymax];
ma21 = [ xmax-xmin, -(xb-x1) ; 0, -(yb-y1) ]; bh21 = [ x1-xmin ; y1-ymax];
ma22 = [ xmax-xmin, -(xb-x2) ; 0, -(yb-y2) ]; bh22 = [ x2-xmin ; y2-ymax];
ma31 = [ 0, -(xb-x1) ; ymax-ymin, -(yb-y1) ]; bh31 = [ x1-xmin ; y1-ymin];
ma32 = [ 0, -(xb-x2) ; ymax-ymin, -(yb-y2) ]; bh32 = [ x2-xmin ; y2-ymin];
ma41 = [ xmin-xmax, -(xb-x1) ; 0, -(yb-y1) ]; bh41 = [ x1-xmax ; y1-ymin];
ma42 = [ xmin-xmax, -(xb-x2) ; 0, -(yb-y2) ]; bh42 = [ x2-xmax ; y2-ymin];

warning('off','Octave:singular-matrix-div'); % Yep, there will be singularity
tr1s = [ ma11\bh11, ma21\bh21, ma31\bh31, ma41\bh41 ]; % All the intersections
tr2s = [ ma12\bh12, ma22\bh22, ma32\bh32, ma42\bh42 ];
warning('on','Octave:singular-matrix-div');

% <Ugly_hack>
ran1 = [ rank(ma11)-1 , rank(ma21)-1 , rank(ma31)-1 , rank(ma41)-1 ]; % Ranks-1
ran2 = [ rank(ma12)-1 , rank(ma22)-1 , rank(ma32)-1 , rank(ma42)-1 ]; %
tr1s = tr1s.*ran1; % De-activate singular ones
tr2s = tr2s.*ran2; %
% </Ugly_hack>

if max(ran1) == 0 % It's the center !
   xy1 = [ xb ; yb ]; xy1r = xy1;% TODO : something to get theta1ref = -theta2ref;
else
   tr1s = tr1s( :, find(tr1s(2,:)>1) ); % Abscissa must be > 1
   [~,num] = min(tr1s(2,:));
   tr1 = tr1s(:,num); u = tr1(2);
   xy1 = [ (1-u)*x1 + u*xb ; (1-u)*y1 + u*yb ]; % coords of the right intersection
   xy1r = xy1; % Store it
   theta1ref = atan2( xy1(2)-yb , xy1(1)-xb );
end

if max(ran2) == 0
   xy2 = [ xb ; yb ]; xy2r = 2*xy2-xy1;
   theta2ref = theta1ref + pi;
else
   tr2s = tr2s( :, find(tr2s(2,:)>1), : ); % Radius must be > 0
   [~,num] = min(tr2s(2,:));
   tr2 = tr2s(:,num); u = tr2(2);
   xy2 = [ (1-u)*x2 + u*xb ; (1-u)*y2 + u*yb ]; % coords of the right intersection
   xy2r = xy2; % Store it
   theta1ref = atan2( xy2(2)-yb , xy2(1)-xb );
end

theta1ref = mod(theta1ref,2*pi);
theta2ref = mod(theta2ref,2*pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ nodes2,elements2,ntoelem2,boundary2,order] = readmesh( 'meshes/rg_refined/plate_nu.msh' );
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

ur1 = UR(:,1); fr1 = UR(:,5)*50; % 50 is the ratio between mesh sizes (hack)
ur2 = UR(:,2); fr2 = UR(:,6)*50;
ur3 = UR(:,3); fr3 = UR(:,7)*50;
ur4 = UR(:,4); fr4 = UR(:,8)*50;

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
   cand1 = rem( find(elements2==no1),nelem2 ); % find gives line + column*size
   cand2 = rem( find(elements2==no2),nelem2 );
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
   cand1 = rem( find(elements2==no1),nelem1 ); % find gives line + column*size
   cand2 = rem( find(elements2==no2),nelem1 );
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
if recompute == 1
   tic
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% Compute the RG and the Lhs
   Rhs1  = zeros(nftest,1); Rhs2  = zeros(nftest,1);
   Rhs3  = zeros(nftest,1); Rhs4  = zeros(nftest,1);
   LhsA  = zeros(nftest,szD); LhsB  = zeros(nftest,0);
   
   Ru = zeros(nftest,2*nnbound2); Rf = zeros(nftest,2*nboun2);
   
   for k=1:nftest
      coefa = coef(1:(nmax+1)^2,k);
      coefb = coef((nmax+1)^2+1:2*(nmax+1)^2,k);
      
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
            
            % Interpolations
            if order == 1
               xgr  = ( (1-xg)*[x1;y1] + xg*[x2;y2] ); % abscissae of the Gauss point in the physical space
            elseif order == 2
               xgr  = ( [x1;y1] + xg*(4*[x3;y3]-3*[x1;y1]-[x2;y2]) + ...
                      xg^2*(2*[x2;y2]+2*[x1;y1]-4*[x3;y3]) );
            end
            
            X = xgr(1)/Lx; Y = xgr(2)/Lx; % Use the standard basis
            
            % Build the test field
            vloc1 = 0; vloc2 = 0;
            sloc11 = 0; sloc12 = 0; sloc22 = 0;
   
            % It is vital to separate the cases ii=0 and jj=0 to avoid 1/0
            ii = 0; jj = 0;
            aij = coefa( (nmax+1)*ii + jj+1 );
            bij = coefb( (nmax+1)*ii + jj+1 );
            vloc1 = vloc1 + aij*X^ii*Y^jj;
            vloc2 = vloc2 + bij*X^ii*Y^jj;
         
            ii = 0;
            for jj=1:nmax
               aij = coefa( (nmax+1)*ii + jj+1 );
               bij = coefb( (nmax+1)*ii + jj+1 );
               vloc1 = vloc1 + aij*X^ii*Y^jj;
               vloc2 = vloc2 + bij*X^ii*Y^jj;
               sloc11 = sloc11 + nu*E/(1-nu^2)*jj*X^ii*Y^(jj-1)*bij;
               sloc22 = sloc22 + E/(1-nu^2)*jj*X^ii*Y^(jj-1)*bij;
               sloc12 = sloc12 + E/(2*(1+nu)) * jj*X^ii*Y^(jj-1)*aij;
            end
            
            jj = 0;
            for ii=1:nmax
               aij = coefa( (nmax+1)*ii + jj+1 );
               bij = coefb( (nmax+1)*ii + jj+1 );
               vloc1 = vloc1 + aij*X^ii*Y^jj;
               vloc2 = vloc2 + bij*X^ii*Y^jj;
               sloc11 = sloc11 + E/(1-nu^2)*ii*X^(ii-1)*Y^jj*aij;
               sloc22 = sloc22 + nu*E/(1-nu^2)*ii*X^(ii-1)*Y^jj*aij;
               sloc12 = sloc12 + E/(2*(1+nu)) * ii*X^(ii-1)*Y^jj*bij;
            end
            
            for ii=1:nmax
               for jj=1:nmax
                  if ii+jj > nmax continue; end % No need to add an other 0
                  aij = coefa( (nmax+1)*ii + jj+1 );
                  bij = coefb( (nmax+1)*ii + jj+1 );
                  vloc1 = vloc1 + aij*X^ii*Y^jj;
                  vloc2 = vloc2 + bij*X^ii*Y^jj;
                  sloc11 = sloc11 + E/(1-nu^2)*ii*X^(ii-1)*Y^jj*aij ...
                                  + nu*E/(1-nu^2)*jj*X^ii*Y^(jj-1)*bij;
                  sloc22 = sloc22 + nu*E/(1-nu^2)*ii*X^(ii-1)*Y^jj*aij ...
                                  + E/(1-nu^2)*jj*X^ii*Y^(jj-1)*bij;
                  sloc12 = sloc12 + E/(2*(1+nu)) * ...
                                    ( jj*X^ii*Y^(jj-1)*aij + ii*X^(ii-1)*Y^jj*bij );
               end
            end
            
            vloc = [ vloc1 ; vloc2 ];
            vpa = vloc;
   
            sloc = 1/Lx*[sloc11,sloc12;sloc12,sloc22];
            spa = sloc;
            fpa = spa*exno;
            
            fpaTimesPhitest0 = [fpa(1)*(1-xg);fpa(2)*(1-xg);fpa(1)*xg;fpa(2)*xg]; 
            
            Ru(k,indDtot) = Ru(k,indDtot) + len * wg * fpaTimesPhitest0';
            Rf(k,indNtot) = Rf(k,indNtot) + len * wg * vpa';
         end
      end
   end

   % Cut LhsA
   LhsA = LhsA( :, find(sum(LhsA.^2)) );
   
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
   
   disp([ 'Right hand side generated ', num2str(toc) ]);
else
   Anb  = load('fields/rg_cauchy_crack/reciprocity_crack_hat22.mat');
   Rf  = Anb.Rf; Ru  = Anb.Ru;
   u_known1 = Anb.u_known1; u_known2 = Anb.u_known2; u_known3 = Anb.u_known3; u_known4 = Anb.u_known4;
   f_known1 = Anb.f_known1; f_known2 = Anb.f_known2; f_known3 = Anb.f_known3; f_known4 = Anb.f_known4;
end
tic

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

ind1 = zeros(nbstep,1); ind2 = zeros(nbstep,1); % Store the stopping indices
ind3 = zeros(nbstep,1); ind4 = zeros(nbstep,1); %
indp1 = zeros(nbstep,1); indp2 = zeros(nbstep,1); % Store the stopping indices (L-curve)
indp3 = zeros(nbstep,1); indp4 = zeros(nbstep,1); %
indd1 = zeros(nbstep,1); indd2 = zeros(nbstep,1); % Store the stopping indices (L-curve)
indd3 = zeros(nbstep,1); indd4 = zeros(nbstep,1); %


theta1rec = theta1;
theta2rec = theta2;
theta1b   = theta1;
theta2b   = theta2;
nbnotwork = 0; % Records the number of failures
tic
for iter = 1:nbstep % Markov chain loop
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% Compute the crack RHS (from angle)
   
   % Intersections
   %xb = mean(nodes2(:,1)); yb = mean(nodes2(:,2)); % Barycenter
   xmin = min(nodes2(:,1)); xmax = max(nodes2(:,1));
   ymin = min(nodes2(:,2)); ymax = max(nodes2(:,2));
   xb = .5*(xmin+xmax); yb = .5*(ymin+ymax);
   ma11 = [ 0, -cos(theta1) ; ymin-ymax, -sin(theta1) ]; bh1 = [ xb-xmax ; yb-ymax];
   ma12 = [ 0, -cos(theta2) ; ymin-ymax, -sin(theta2) ];
   ma21 = [ xmax-xmin, -cos(theta1) ; 0, -sin(theta1) ]; bh2 = [ xb-xmin ; yb-ymax];
   ma22 = [ xmax-xmin, -cos(theta2) ; 0, -sin(theta2) ];
   ma31 = [ 0, -cos(theta1) ; ymax-ymin, -sin(theta1) ]; bh3 = [ xb-xmin ; yb-ymin];
   ma32 = [ 0, -cos(theta2) ; ymax-ymin, -sin(theta2) ];
   ma41 = [ xmin-xmax, -cos(theta1) ; 0, -sin(theta1) ]; bh4 = [ xb-xmax ; yb-ymin];
   ma42 = [ xmin-xmax, -cos(theta2) ; 0, -sin(theta2) ];
   
   %cond(ma11), cond(ma21), cond(ma31), cond(ma41), cond(ma12), cond(ma22), cond(ma32), cond(ma42)
   warning('off','Octave:singular-matrix-div'); % Yep, there will be singularity
   tr1s = [ ma11\bh1, ma21\bh2, ma31\bh3, ma41\bh4 ]; % All the intersections
   tr2s = [ ma12\bh1, ma22\bh2, ma32\bh3, ma42\bh4 ];
   warning('on','Octave:singular-matrix-div');
   
   % <Ugly_hack>
   ran1 = [ rank(ma11)-1 , rank(ma21)-1 , rank(ma31)-1 , rank(ma41)-1 ]; % Ranks-1
   ran2 = [ rank(ma12)-1 , rank(ma22)-1 , rank(ma32)-1 , rank(ma42)-1 ]; %
   tr1s = tr1s.*ran1; % De-activate singular ones
   tr2s = tr2s.*ran2; %
   % </Ugly_hack>
   
   tr1s = tr1s( :, find(tr1s(2,:)>0) ); % Radius must be > 0
   [~,num] = min(tr1s(2,:));
   tr1 = tr1s(:,num); r = tr1(2);
   xy1 = [ xb + r*cos(theta1) ; yb + r*sin(theta1) ]; % coords of the right intersection
   
   tr2s = tr2s( :, find(tr2s(2,:)>0), : ); % Radius must be > 0
   [~,num] = min(tr2s(2,:));
   tr2 = tr2s(:,num); r = tr2(2);
   xy2 = [ xb + r*cos(theta2) ; yb + r*sin(theta2) ]; % coords of the right intersection

%   % Detect if both points are on the same bound
%   if (xy1(1) == xy2(1) && (xy1(1) == xmax || xy1(1) == xmin)) || ...
%      (xy1(2) == xy2(2) && (xy1(1) == ymax || xy1(1) == ymin))
%      
%      theta1 = theta1b + (2*rand()-1)*anglestep; % Cancel this choice
%      theta2 = theta2b + (2*rand()-1)*anglestep;
%      continue;
%   end
   
   %% Build the elements on the crack's line
   step = (xy2-xy1)/10;
   nodes3 = [ xy1(1):step(1):xy2(1) ; xy1(2):step(2):xy2(2) ];
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

   for i=1:nboun3
      bonod = boundary3(i,:); exno = extnorm3(i,:)';
   
      no1 = bonod(2); no2 = bonod(3); boname = bonod(1);
      x1 = nodes3(no1,1); y1 = nodes3(no1,2);
      x2 = nodes3(no2,1); y2 = nodes3(no2,2);
      len = sqrt( (x1-x2)^2 + (y1-y2)^2 );
      
      indDtot = [ 2*i-1, 2*i, 2*(i+1)-1, 2*(i+1) ];      
                
      for j=1:Ng
         xg = Xg(j); wg = Wg(j);
         xgr  = ( (1-xg)*[x1;y1] + xg*[x2;y2] ); % abscissae of the Gauss point in the physical space
         X = xgr(1)/Lx; Y = xgr(2)/Lx; % Use the standard basis
         lw = len*wg;

         for ii=0:nmax
            nimax = nmax-ii;
            for jj=0:nimax
               %if ii+jj > nmax continue; end % No need to add an other 0
               
               % Build the test field
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
               
               fpaTimesPhitest0a = [fpaa(1)*(1-xg);fpaa(2)*(1-xg);fpaa(1)*xg;fpaa(2)*xg];
               fpaTimesPhitest0b = [fpab(1)*(1-xg);fpab(2)*(1-xg);fpab(1)*xg;fpab(2)*xg]; 
               
               Rucij((nmax+1)*ii + jj+1,indDtot) = Rucij((nmax+1)*ii + jj+1,indDtot) ...
                                                   + lw * fpaTimesPhitest0a';
               Rucij((nmax+1)^2 + (nmax+1)*ii + jj+1,indDtot) =...
                                      Rucij((nmax+1)^2 + (nmax+1)*ii + jj+1,indDtot) ...
                                                   + lw * fpaTimesPhitest0b';
            end
         end
      end
   end
   
   Ruc = coef'*Rucij;

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
   jmax = jmaxRG;
   if froreg == 1 && min(size(LhsB))>0
      kB = sqrt(norm(LhsA,'fro')^2+norm(Ruc,'fro')^2)/norm(LhsB,'fro'); % In order to regularize the stuff
   %   kB = norm(LhsA,'fro')/norm(LhsB,'fro');
   else
      kB = 1;
   end
   %Lhs = [LhsA,kB*LhsB];
   Lhs = [LhsA,kB*LhsB,Ruc];
   
   Rhs = Rhs1;
   A = Lhs'*Lhs; b = Lhs'*Rhs; sA = size(A,1);
   
   Z12    = zeros(size(Mum,1),size(Mfm,2)); Z13 = zeros(size(Mum,1),size(nodeMass3,2));
   Z23    = zeros(size(Mfm,1),size(nodeMass3,2));
   Mtot   = [ Mum, Z12, Z13  ; Z12', Mfm, Z23 ; Z13', Z23', nodeMass3 ];
   Msmall = [ 0*Mum, Z12, Z13  ; Z12', 0*Mfm, Z23 ; Z13', Z23', nodeMass3 ]; % Only on the crack
   Mtop   = Mtot-Msmall; % Debug asset

   % Differential regularization matrix
   D3 = zeros(2*nboun3+2);
   for i=1:nboun3
      coefU = [ i , i+1 ];
      len = norm(nodes3(i,:) - nodes3(i+1,:));
      D3( 2*coefU-1, 2*coefU-1 ) = D3( 2*coefU-1, 2*coefU-1 ) + 1/len*[1,-1;-1,1];
      D3( 2*coefU, 2*coefU )     = D3( 2*coefU, 2*coefU )     + 1/len*[1,-1;-1,1];
   end
   Dsmall = [ 0*Mum, Z12, Z13  ; Z12', 0*Mfm, Z23 ; Z13', Z23', D3 ];

   if regular == 1
      Z12 = zeros( size(Duu,1) , size(Dfu,2) ); Z13 = zeros( size(Duu,1), size(D3,2) );
      Z23 = zeros( size(Dfu,1), size(D3,2) );
      Dtot = [ Duu ,Z12, Z13 ; Z12', Dfu, Z23 ; Z13', Z23', D3 ];
      L = Dtot + Mtot * norm(Dtot,'fro')/(10*norm(Mtot,'fro'));
   else
%      L = eye(size(A));
      L = Mtot;
   end
   ninfty = 0;
   %L = eye(size(A));
   
   %[~,Theta,Q] = svd(Lhs,L); Q = Q*real((Q'*L*Q)^(-1/2)); Theta = Theta.^2;
   %[Q,Theta] = eig(Lhs'*L*Lhs); Q = Q*real((Q'*L*Q)^(-1/2)); Theta = Theta.^2;
   [Q,Theta] = eig(A,L); Q = Q*real((Q'*L*Q)^(-1/2)); % real should not be, but you know, numerical shit...
   thetas = diag(Theta);
   [thetas,Ind] = sort( thetas,'descend' );
   Q = Q(:,Ind);
   Thetas = diag(thetas); 
   Theta = Theta(Ind,Ind);
   
   %disp([ num2str(size(A,1)), ' - dofs rectangular system pinversed ', num2str(toc) ]);
   % Plot the Picard stuff
   imax = min( find(thetas/thetas(ninfty+1)<1e-16) );
   if size(imax,1) == 0
       imax = size(thetas,1);
   end

   % Eigenmode on the crack
   Qsmall = Q; Qsmall( 1:end-2*nnodes3+1, : ) = 0;
   
   tplo = thetas(ninfty+1:imax);
   bplo1 = Q'*Lhs'*Rhs1; bplo1 = bplo1(ninfty+1:imax);
   rplp1 = (Q'*Lhs'*Rhs1)./thetas; rplo1 = rplp1(1:imax);
   bplo2 = Q'*Lhs'*Rhs2; bplo2 = bplo2(ninfty+1:imax);
   rplp2 = (Q'*Lhs'*Rhs2)./thetas; rplo2 = rplp2(1:imax);
   bplo3 = Q'*Lhs'*Rhs3; bplo3 = bplo3(ninfty+1:imax);
   rplp3 = (Q'*Lhs'*Rhs3)./thetas; rplo3 = rplp3(1:imax);
   bplo4 = Q'*Lhs'*Rhs4; bplo4 = bplo4(ninfty+1:imax);
   rplp4 = (Q'*Lhs'*Rhs4)./thetas; rplo4 = rplp4(1:imax);

   nor1 = zeros(sA,1); nor2 = zeros(sA,1);
   nor3 = zeros(sA,1); nor4 = zeros(sA,1);
   nod1 = zeros(sA,1); nod2 = zeros(sA,1);
   nod3 = zeros(sA,1); nod4 = zeros(sA,1);
   nop1 = zeros(sA,1); nop2 = zeros(sA,1);
   nop3 = zeros(sA,1); nop4 = zeros(sA,1);
   res1 = zeros(sA,1); res2 = zeros(sA,1);
   res3 = zeros(sA,1); res4 = zeros(sA,1);
   for i=1:sA
      xnor1 = Q(:,1:i)*rplp1(1:i); nor1(i) = xnor1'*Msmall*xnor1; nod1(i) = xnor1'*Dsmall*xnor1;
      xnor2 = Q(:,1:i)*rplp2(1:i); nor2(i) = xnor2'*Msmall*xnor2; nod2(i) = xnor2'*Dsmall*xnor2;
      xnor3 = Q(:,1:i)*rplp3(1:i); nor3(i) = xnor3'*Msmall*xnor3; nod3(i) = xnor3'*Dsmall*xnor3;
      xnor4 = Q(:,1:i)*rplp4(1:i); nor4(i) = xnor4'*Msmall*xnor4; nod4(i) = xnor4'*Dsmall*xnor4;
      xnop1 = Q(:,1:i)*rplp1(1:i); nop1(i) = xnor1'*Mtop*xnor1;
      xnop2 = Q(:,1:i)*rplp2(1:i); nop2(i) = xnor2'*Mtop*xnor2;
      xnop3 = Q(:,1:i)*rplp3(1:i); nop3(i) = xnor3'*Mtop*xnor3;
      xnop4 = Q(:,1:i)*rplp4(1:i); nop4(i) = xnor4'*Mtop*xnor4;
%      xnop1 = Q(:,i)*rplo1(i); nop1(i) = sqrt(xnop1'*Msmall*xnop1); noq1(i) = sqrt(xnop1'*Mtop*xnop1);
%      xnop2 = Q(:,i)*rplo2(i); nop2(i) = sqrt(xnop2'*Msmall*xnop2); noq2(i) = sqrt(xnop2'*Mtop*xnop2);
%      xnop3 = Q(:,i)*rplo3(i); nop3(i) = sqrt(xnop3'*Msmall*xnop3); noq3(i) = sqrt(xnop3'*Mtop*xnop3);
%      xnop4 = Q(:,i)*rplo4(i); nop4(i) = sqrt(xnop4'*Msmall*xnop4); noq4(i) = sqrt(xnop4'*Mtop*xnop4);
      res1(i)  = norm(Lhs*xnor1 - Rhs1);
      res2(i)  = norm(Lhs*xnor2 - Rhs2);
      res3(i)  = norm(Lhs*xnor3 - Rhs3);
      res4(i)  = norm(Lhs*xnor4 - Rhs4);
   end
   not1 = nor1+nop1; not2 = nor2+nop2; not3 = nor3+nop3; not4 = nor4+nop4; % Total norm
   
   % Remove Zeros in rploi (why on hell are there zeros in the first place ?)
   me1 = mean(abs(rplp1))/1e5; arplo1 = max(me1,abs(rplp1));
   me2 = mean(abs(rplp2))/1e5; arplo2 = max(me2,abs(rplp2));
   me3 = mean(abs(rplp3))/1e5; arplo3 = max(me3,abs(rplp3));
   me4 = mean(abs(rplp4))/1e5; arplo4 = max(me4,abs(rplp4));
   ind1(iter) = findPicard2 (log10(arplo1), ceil(sA/7), 1, 3);
   ind2(iter) = findPicard2 (log10(arplo2), ceil(sA/7), 1, 3);
   ind3(iter) = findPicard2 (log10(arplo3), ceil(sA/7), 1, 3);
   ind4(iter) = findPicard2 (log10(arplo4), ceil(sA/7), 1, 3);

   indp1(iter) = findCorner (res1, nor1, ceil(sA/4), 1, 4);
   indp2(iter) = findCorner (res2, nor2, ceil(sA/4), 1, 4);
   indp3(iter) = findCorner (res3, nor3, ceil(sA/4), 1, 4);
   indp4(iter) = findCorner (res4, nor4, ceil(sA/4), 1, 4);

   indd1(iter) = findCorner (res1, nod1, ceil(sA/4), 1, 4);
   indd2(iter) = findCorner (res2, nod2, ceil(sA/4), 1, 4);
   indd3(iter) = findCorner (res3, nod3, ceil(sA/4), 1, 4);
   indd4(iter) = findCorner (res4, nod4, ceil(sA/4), 1, 4);

   if norm(nop1)~=0 % In case Mtop = []
      inds1(iter) = findCorner (res1, nop1, ceil(sA/4), 1, 4);
      inds2(iter) = findCorner (res2, nop2, ceil(sA/4), 1, 4);
      inds3(iter) = findCorner (res3, nop3, ceil(sA/4), 1, 4);
      inds4(iter) = findCorner (res4, nop4, ceil(sA/4), 1, 4);
   else
      inds1(iter) = 1; inds2(iter) = 1; inds3(iter) = 1; inds4(iter) = 1;
   end
%   indt1(iter) = findCorner (res1, not1, ceil(sA/4), 1, 4);
%   indt2(iter) = findCorner (res2, not2, ceil(sA/4), 1, 4);
%   indt3(iter) = findCorner (res3, not3, ceil(sA/4), 1, 4);
%   indt4(iter) = findCorner (res4, not4, ceil(sA/4), 1, 4);

   %if ind1(iter)==1
%      [~,p] = findPicard2 (log10(arplo1), 7, 1, 3);
%      n = size(p,1);
%      t = 1:.05:imax; tt = zeros(n,20*(imax-1)+1);
%      for j=1:n
%         tt(j,:) = t.^(n-j);
%      end
%      px = p'*tt;
%      
%      figure
%      hold on;
%      plot(log10(arplo1),'Color','black');
%      plot(t,px,'Color','red');
%      legend('sol1','poly1');

%      bug;
   %end

%   figure
%   hold on;
%   plot(log10(abs(tplo)),'Color','green');
%   plot(log10(abs(bplo1)),'Color','red');
%   plot(log10(abs(rplo1)),'Color','black');
%   plot(log10(abs(bplo2)),'Color','magenta');
%   plot(log10(abs(rplo2)),'Color','blue');
%   legend('Singular values','Rhs1','sol1','Rhs2','sol2');
   
%   try
%   figure
%   hold on;
%   plot(log10(abs(rplo1)),'Color','red');
%   plot(log10(abs(rplo2)),'Color','black');
%   plot(log10(abs(rplo3)),'Color','magenta');
%   plot(log10(abs(rplo4)),'Color','blue');
%   legend('sol1','sol2','sol3','sol4');
%   end

%   figure
%   hold on;
%   plot(log10(nop1),'Color','red');
%   plot(log10(nop2),'Color','black');
%   plot(log10(nop3),'Color','magenta');
%   plot(log10(nop4),'Color','blue');
%   legend('sol1','sol2','sol3','sol4');

%   figure
%   hold on;
%   plot(log10(noq1),'Color','red');
%   plot(log10(noq2),'Color','black');
%   plot(log10(noq3),'Color','magenta');
%   plot(log10(noq4),'Color','blue');
%   legend('sol1','sol2','sol3','sol4');

%   try
%   figure;
%   hold on;
%   loglog(res1,nor1,'Color','red','-*');
%   loglog(res1(indp1),nor1(indp1),'Color','red','o','markersize',15);
%   loglog(res2,nor2,'Color','black','-*');
%   loglog(res2(indp2),nor2(indp2),'Color','black','o','markersize',15);
%   loglog(res3,nor3,'Color','magenta','-*');
%   loglog(res3(indp3),nor3(indp3),'Color','magenta','o','markersize',15);
%   loglog(res4,nor4,'Color','blue','-*');
%   loglog(res4(indp4),nor4(indp4),'Color','blue','o','markersize',15);
%   legend('L1','','L2','','L3','','L4','');
%   end

%   try
%   figure;
%   hold on;
%   loglog(res1,nod1,'Color','red','-*');
%   loglog(res1(indd1),nod1(indd1),'Color','red','o','markersize',15);
%   loglog(res2,nod2,'Color','black','-*');
%   loglog(res2(indd2),nod2(indd2),'Color','black','o','markersize',15);
%   loglog(res3,nod3,'Color','magenta','-*');
%   loglog(res3(indd3),nod3(indd3),'Color','magenta','o','markersize',15);
%   loglog(res4,nod4,'Color','blue','-*');
%   loglog(res4(indd4),nod4(indd4),'Color','blue','o','markersize',15);
%   legend('L1','','L2','','L3','','L4','');
%   end

%   try
%   figure;
%   hold on;
%   loglog(res1,nop1,'Color','red','-*');
%   loglog(res1(inds1),nop1(inds1),'Color','red','o','markersize',15);
%   loglog(res2,nop2,'Color','black','-*');
%   loglog(res2(inds2),nop2(inds2),'Color','black','o','markersize',15);
%   loglog(res3,nop3,'Color','magenta','-*');
%   loglog(res3(inds3),nop3(inds3),'Color','magenta','o','markersize',15);
%   loglog(res4,nop4,'Color','blue','-*');
%   loglog(res4(inds4),nop4(inds4),'Color','blue','o','markersize',15);
%   legend('L1','','L2','','L3','','L4','');
%   end

%   try
%   figure;
%   hold on;
%   loglog(res1,not1,'Color','red','-*');
%   loglog(res1(indt1),not1(indt1),'Color','red','o','markersize',15);
%   loglog(res2,not2,'Color','black','-*');
%   loglog(res2(indt2),not2(indt2),'Color','black','o','markersize',15);
%   loglog(res3,not3,'Color','magenta','-*');
%   loglog(res3(indt3),not3(indt3),'Color','magenta','o','markersize',15);
%   loglog(res4,not4,'Color','blue','-*');
%   loglog(res4(indt4),not4(indt4),'Color','blue','o','markersize',15);
%   legend('L1','','L2','','L3','','L4','');
%   end
   
   % Filter eigenvalues
   if jmax == 0
      jmax = size(Thetas,1);
      if paramchoice == 1
         jmax1 = indp1(iter); jmax2 = indp2(iter); jmax3 = indp3(iter); jmax4 = indp4(iter);
      elseif paramchoice == 2
         jmax1 = ind1(iter); jmax2 = ind2(iter); jmax3 = ind3(iter); jmax4 = ind4(iter);
      elseif paramchoice == 3
         jmax1 = indd1(iter); jmax2 = indd2(iter); jmax3 = indd3(iter); jmax4 = indd4(iter);
      else % paramchoice == 4
         jmax1 = inds1(iter); jmax2 = inds2(iter); jmax3 = inds3(iter); jmax4 = inds4(iter);
      end
   else
      jmax0 = min( size(Thetas,1) , jmax );
      jmax1 = jmax; jmax2 = jmax; jmax3 = jmax; jmax4 = jmax;
   end
   ThetaT1 = Thetas( ninfty+1:jmax1 , ninfty+1:jmax1 );
   ThetaT2 = Thetas( ninfty+1:jmax2 , ninfty+1:jmax2 );
   ThetaT3 = Thetas( ninfty+1:jmax3 , ninfty+1:jmax3 );
   ThetaT4 = Thetas( ninfty+1:jmax4 , ninfty+1:jmax4 );
   
   bU1 = Q'*Lhs'*Rhs1; bT1 = bU1(ninfty+1:jmax1);
   bU2 = Q'*Lhs'*Rhs2; bT2 = bU2(ninfty+1:jmax2);
   bU3 = Q'*Lhs'*Rhs3; bT3 = bU3(ninfty+1:jmax3);
   bU4 = Q'*Lhs'*Rhs4; bT4 = bU4(ninfty+1:jmax4);
   
   Solu1RG = Q(:,ninfty+1:jmax1) * (ThetaT1\bT1);
   Solu2RG = Q(:,ninfty+1:jmax2) * (ThetaT2\bT2);
   Solu3RG = Q(:,ninfty+1:jmax3) * (ThetaT3\bT3);
   Solu4RG = Q(:,ninfty+1:jmax4) * (ThetaT4\bT4);
   
   Solu1 = Solu1RG; Solu2 = Solu2RG; Solu3 = Solu3RG; Solu4 = Solu4RG;

   %% Compute phi with 4 different ways
   jmax1 = indp1(iter); jmax2 = indp2(iter); jmax3 = indp3(iter); jmax4 = indp4(iter);
   ThetaT1 = Thetas( ninfty+1:jmax1 , ninfty+1:jmax1 ); bT1 = bU1(ninfty+1:jmax1);
   ThetaT2 = Thetas( ninfty+1:jmax2 , ninfty+1:jmax2 ); bT2 = bU2(ninfty+1:jmax2);
   ThetaT3 = Thetas( ninfty+1:jmax3 , ninfty+1:jmax3 ); bT3 = bU3(ninfty+1:jmax3);
   ThetaT4 = Thetas( ninfty+1:jmax4 , ninfty+1:jmax4 ); bT4 = bU4(ninfty+1:jmax4);
   Solu1L = Q(:,ninfty+1:jmax1) * (ThetaT1\bT1);
   Solu2L = Q(:,ninfty+1:jmax2) * (ThetaT2\bT2);
   Solu3L = Q(:,ninfty+1:jmax3) * (ThetaT3\bT3);
   Solu4L = Q(:,ninfty+1:jmax4) * (ThetaT4\bT4);

   jmax1 = ind1(iter); jmax2 = ind2(iter); jmax3 = ind3(iter); jmax4 = ind4(iter);
   ThetaT1 = Thetas( ninfty+1:jmax1 , ninfty+1:jmax1 ); bT1 = bU1(ninfty+1:jmax1);
   ThetaT2 = Thetas( ninfty+1:jmax2 , ninfty+1:jmax2 ); bT2 = bU2(ninfty+1:jmax2);
   ThetaT3 = Thetas( ninfty+1:jmax3 , ninfty+1:jmax3 ); bT3 = bU3(ninfty+1:jmax3);
   ThetaT4 = Thetas( ninfty+1:jmax4 , ninfty+1:jmax4 ); bT4 = bU4(ninfty+1:jmax4);
   Solu1P = Q(:,ninfty+1:jmax1) * (ThetaT1\bT1);
   Solu2P = Q(:,ninfty+1:jmax2) * (ThetaT2\bT2);
   Solu3P = Q(:,ninfty+1:jmax3) * (ThetaT3\bT3);
   Solu4P = Q(:,ninfty+1:jmax4) * (ThetaT4\bT4);

   jmax1 = indd1(iter); jmax2 = indd2(iter); jmax3 = indd3(iter); jmax4 = indd4(iter);
   ThetaT1 = Thetas( ninfty+1:jmax1 , ninfty+1:jmax1 ); bT1 = bU1(ninfty+1:jmax1);
   ThetaT2 = Thetas( ninfty+1:jmax2 , ninfty+1:jmax2 ); bT2 = bU2(ninfty+1:jmax2);
   ThetaT3 = Thetas( ninfty+1:jmax3 , ninfty+1:jmax3 ); bT3 = bU3(ninfty+1:jmax3);
   ThetaT4 = Thetas( ninfty+1:jmax4 , ninfty+1:jmax4 ); bT4 = bU4(ninfty+1:jmax4);
   Solu1D = Q(:,ninfty+1:jmax1) * (ThetaT1\bT1);
   Solu2D = Q(:,ninfty+1:jmax2) * (ThetaT2\bT2);
   Solu3D = Q(:,ninfty+1:jmax3) * (ThetaT3\bT3);
   Solu4D = Q(:,ninfty+1:jmax4) * (ThetaT4\bT4);

   jmax1 = inds1(iter); jmax2 = inds2(iter); jmax3 = inds3(iter); jmax4 = inds4(iter);
   ThetaT1 = Thetas( ninfty+1:jmax1 , ninfty+1:jmax1 ); bT1 = bU1(ninfty+1:jmax1);
   ThetaT2 = Thetas( ninfty+1:jmax2 , ninfty+1:jmax2 ); bT2 = bU2(ninfty+1:jmax2);
   ThetaT3 = Thetas( ninfty+1:jmax3 , ninfty+1:jmax3 ); bT3 = bU3(ninfty+1:jmax3);
   ThetaT4 = Thetas( ninfty+1:jmax4 , ninfty+1:jmax4 ); bT4 = bU4(ninfty+1:jmax4);
   Solu1S = Q(:,ninfty+1:jmax1) * (ThetaT1\bT1);
   Solu2S = Q(:,ninfty+1:jmax2) * (ThetaT2\bT2);
   Solu3S = Q(:,ninfty+1:jmax3) * (ThetaT3\bT3);
   Solu4S = Q(:,ninfty+1:jmax4) * (ThetaT4\bT4);

   phiL(iter) = norm(Lhs*Solu1L-Rhs1)^2 + norm(Lhs*Solu2L-Rhs2)^2 +...
               norm(Lhs*Solu3L-Rhs3)^2 + norm(Lhs*Solu4L-Rhs4)^2;
   phiP(iter) = norm(Lhs*Solu1P-Rhs1)^2 + norm(Lhs*Solu2P-Rhs2)^2 +...
               norm(Lhs*Solu3P-Rhs3)^2 + norm(Lhs*Solu4P-Rhs4)^2;
   phiD(iter) = norm(Lhs*Solu1D-Rhs1)^2 + norm(Lhs*Solu2D-Rhs2)^2 +...
               norm(Lhs*Solu3D-Rhs3)^2 + norm(Lhs*Solu4D-Rhs4)^2;
   phiS(iter) = norm(Lhs*Solu1S-Rhs1)^2 + norm(Lhs*Solu2S-Rhs2)^2 +...
               norm(Lhs*Solu3S-Rhs3)^2 + norm(Lhs*Solu4S-Rhs4)^2;

   phi(iter) = phiL(iter)^weights(1)*phiP(iter)^weights(2)*phiD(iter)^weights(3)*phiS(iter)^weights(4);
   
%   phi(iter) = sqrt( norm(Lhs*Solu1RG-Rhs1)^2 + norm(Lhs*Solu2RG-Rhs2)^2 +...
%                     norm(Lhs*Solu3RG-Rhs3)^2 + norm(Lhs*Solu4RG-Rhs4)^2 ); % residual

%   phi(iter) = norm(Lhs*Solu2RG-Rhs2)^2 +...
%               norm(Lhs*Solu3RG-Rhs3)^2 + norm(Lhs*Solu4RG-Rhs4)^2; % residual
               
   if phi(iter) == min(phi) % It's the minimum : continue from this one
%      iter
      Solu1B    = Solu1; Solu2B = Solu2; Solu3B = Solu3; Solu4B = Solu4;
      theta1b   = theta1;
      theta2b   = theta2;
      nodes3b   = nodes3;
      nbnotwork = 0; % Yay ! Got optimized !
      reloop = 1;
      while reloop % Check if we aren't in the same quarter
         theta1 = theta1 + (2*rand()-1)*anglestep; theta1 = mod(theta1,2*pi);
         theta2 = theta2 + (2*rand()-1)*anglestep; theta2 = mod(theta2,2*pi);
         
         in11 = theta1 <= pi/4   || theta1 >= 7*pi/4;
         in12 = theta1 >= pi/4   && theta1 <= 3*pi/4;
         in13 = theta1 >= 3*pi/4 && theta1 <= 5*pi/4;
         in14 = theta1 >= 5*pi/4 && theta1 <= 7*pi/4;;
         
         in21 = theta2 <= pi/4   || theta2 >= 7*pi/4;
         in22 = theta2 >= pi/4   && theta2 <= 3*pi/4;
         in23 = theta2 >= 3*pi/4 && theta2 <= 5*pi/4;
         in24 = theta2 >= 5*pi/4 && theta2 <= 7*pi/4;
         
         reloop = in11&&in21 || in12&&in22 || in13&&in23 || in14&&in24;
      end
   else % otherwise : continue from the previous one
      nbnotwork = nbnotwork+1; % Records the nb of consequent failures (to decrease step)
      if nbnotwork > limnotwork
         anglestep = anglestep/2;
         nbnotwork = 0;%ceil(nbnotwork/2); % Let a chance to this new value
      end

      reloop = 1;
      while reloop % Check if we aren't in the same quarter
         theta1 = theta1b + (2*rand()-1)*anglestep; theta1 = mod(theta1,2*pi);
         theta2 = theta2b + (2*rand()-1)*anglestep; theta2 = mod(theta2,2*pi);
         
         in11 = theta1 <= pi/4   || theta1 >= 7*pi/4;
         in12 = theta1 >= pi/4   && theta1 <= 3*pi/4;
         in13 = theta1 >= 3*pi/4 && theta1 <= 5*pi/4;
         in14 = theta1 >= 5*pi/4 && theta1 <= 7*pi/4;;
         
         in21 = theta2 <= pi/4   || theta2 >= 7*pi/4;
         in22 = theta2 >= pi/4   && theta2 <= 3*pi/4;
         in23 = theta2 >= 3*pi/4 && theta2 <= 5*pi/4;
         in24 = theta2 >= 5*pi/4 && theta2 <= 7*pi/4;
         
         reloop = in11&&in21 || in12&&in22 || in13&&in23 || in14&&in24;
      end
   end
   theta1rec(iter+1) = theta1;
   theta2rec(iter+1) = theta2;
end
disp(['Iterative method terminated ', num2str(toc) ]);

Solu1 = Solu1B; Solu2 = Solu2B; Solu3 = Solu3B; Solu4 = Solu4B;

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

%% Segments visu
%toplot = fsolu2(1:2:end-1);
%%toplot = usolu1(2:2:end);
%figure; hold on;
%minn1 = min(toplot); maxn1 = max(toplot);
%for i=1:nboun2
%   no1 = boundary2(i,2); no2 = boundary2(i,3);
%   x1 = nodesu(no1,1); y1 = nodesu(no1,2);
%   x2 = nodesu(no2,1); y2 = nodesu(no2,2);
%   
%   x = (toplot(i)-minn1)/(maxn1-minn1);
%   rgb = rgbmap(x);
%   plot( [x1,x2], [y1,y2], 'Color', rgb, 'LineWidth', 3 );
%end
%%axis equal;
%colormap('default')
%h = colorbar();
%ytick = get (h, 'ytick');
%set (h, 'yticklabel', sprintf ( '%g|', (maxn1-minn1)*ytick+minn1 ));

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

try
figure;
hold on;
plot(theta1pi,'Color','red');
plot(theta2pi,'Color','blue');
plot(theta1ref/pi*theta1pi./theta1pi,'Color','black'); % I love disgusting hacks
plot(theta2ref/pi*theta1pi./theta1pi,'Color','black');
legend('theta1/pi', 'theta2/pi','theta1ref/pi', 'theta2ref/pi');
end

%% Recover the reference
step = (xy2r-xy1r)/30;
n = [-step(2);step(1)]; n = n/norm(n); % Normal
nodes3r = [ xy1r(1):step(1):xy2r(1) ; xy1r(2):step(2):xy2r(2) ];
nodes3r = nodes3r'; nnodes3r = size(nodes3r,1);
% Check if need to reverse the nodes (for visu purposes)
if norm(nodes3r(2,:)-nodes3b(1,:)) < norm(nodes3r(1,:)-nodes3b(1,:))
   nodes3r = nodes3r(end:-1:1,:);
end
nodes3s = nodes3r + 1e-3*ones(size(nodes3r,1),1)*n';
nodes3i = nodes3r - 1e-3*ones(size(nodes3r,1),1)*n';
urs = passMesh2D( nodes, elements, nodes3s, [], [u1,u2,u3,u4] );
uri = passMesh2D( nodes, elements, nodes3i, [], [u1,u2,u3,u4] );
urg = uri-urs;  % Vectorial gap
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
curv = sqrt( (nodes3b(:,1)-nodes3b(1,1)).^2 + (nodes3b(:,2)-nodes3b(1,2)).^2 );
%toplot1  = sqrt( ucrsol1(1:2:end-1).^2 + ucrsol1(2:2:end).^2 );
%toplot2  = sqrt( ucrsol2(1:2:end-1).^2 + ucrsol2(2:2:end).^2 );
%toplot3  = sqrt( ucrsol3(1:2:end-1).^2 + ucrsol3(2:2:end).^2 );
%toplot4  = sqrt( ucrsol4(1:2:end-1).^2 + ucrsol4(2:2:end).^2 );
%toplotr1 = sqrt( urg(1:2:end-1,1).^2 + urg(2:2:end,1).^2 );
%toplotr2 = sqrt( urg(1:2:end-1,2).^2 + urg(2:2:end,2).^2 );
%toplotr3 = sqrt( urg(1:2:end-1,3).^2 + urg(2:2:end,3).^2 );
%toplotr4 = sqrt( urg(1:2:end-1,4).^2 + urg(2:2:end,4).^2 );

toplot1  = ucrsol1(1:2:end-1);
toplot2  = ucrsol2(1:2:end-1);
toplot3  = ucrsol3(1:2:end-1);
toplot4  = ucrsol4(1:2:end-1);
toplotr1 = urg(1:2:end-1,1);
toplotr2 = urg(1:2:end-1,2);
toplotr3 = urg(1:2:end-1,3);
toplotr4 = urg(1:2:end-1,4);

try
figure;
hold on;
%plot( curv, toplot1,'Color','green' );
%plot( curv, toplot2,'Color','black' );
%plot( curv, toplot3,'Color','blue' );
plot( curv, toplot4,'Color','red' );
%plot( curvr, toplotr1,'Color','green' );
%plot( curvr, toplotr2,'Color','black' );
%plot( curvr, toplotr3,'Color','blue' );
plot( curvr, toplotr4,'Color','red' );
legend('[[u]] identified (4)', '[[u]] reference (4)');
%legend('[[u]] identified (1)', '[[u]] identified (2)', '[[u]] identified (3)', '[[u]] identified (4)',...
%       '[[u]] reference (1)', '[[u]] reference (2)', '[[u]] reference (3)', '[[u]] reference (4)');
end

%erroru = norm(usolu1(b2nodesnoD)-ur1(b2nodesnoD))   / norm(ur1(b2nodesnoD));
%errorf = norm(fsolu1no(b2nodesnoN)-fr1(b2nodesnoN)) / norm(fr1(b2nodesnoN));
