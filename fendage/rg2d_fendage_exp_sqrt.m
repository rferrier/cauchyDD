% 27/06/2018
% Essai de fendage, 1 seule donnée, forme en racine

tic
close all;
clear all;

addpath(genpath('./tools'));

% Parameters
E           = 17000; % MPa : Young modulus
nu          = 0.2;    % Poisson ratio
fscalar     = 1;    % N.mm-1 : Loading on the plate
mat         = [0, E, nu];
mur         = 1e8;%2e3;    % Regularization parameter
regular     = 1;      % Use the derivative regularization matrix (0 : Id, 1 : derivative)
derivative  = 1;    % Use H1 or H1/2 regularization
froreg      = 1;      % Frobenius preconditioner
x0          = 10;     % Crack Length (to optimize)
anglestep   = 1;%pi/1000;  % Step in angle for Finite Differences anglestep = 0 means auto-adaptation
kauto       = 10;     % Coefficient for the auto-adaptation
nbstep      = 100;     % Nb of Newton Iterations
Npg         = 2;      % Nb Gauss points
%ordertest   = 20;     % Order of test fonctions (by default 20)
ndofcrack   = 20;      % Nb of elements on the crack
ratio       = 0.062;   % mm/pix ratio
order       = 1;       % Order of the mesh
th          = 72.5;   % Thickness (unused as it is only for the loadings)
t0          = 0;%263;     % Time for the computation (263 is good) 0 means all
frf         = 0; % Give dof to the forces (modified RG)
fru         = 0; % Give dof to the displacement (modified RG)
gamma       = 1e1; % Weight of the data (modified RG)
regularized = 1; % Use the regularized data
keep_disgus = 1; % Keep the disgusting values near of the crack (for not regularized only)

% Manually disable some nodes for Dirichlet
%d_nodes = [37,59,60,68,10,74,11,75,12,13,76,14,77,15,78,79];
d_nodes = [];

% Material properties
lambda = nu*E/((1+nu)*(1-2*nu));
mu     = E/(2*(1+nu));

% Boundary conditions
dirichlet  = [ 12,1,0 ; 12,2,0 ];
neumann1   = [ 4,1,fscalar ; 4,2,-fscalar ; 9,1,-fscalar ; 9,2,-fscalar ];

dirichlet0 = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ];
neumann0   = [ 1,1,0 ; 1,2,0 ; 3,1,0 ; 3,2,0 ]; % 3 are the discarded Dirichlet
%neumann0   = [ 1,1,0 ; 1,2,0 ];

%% Import the experimental data
if regularized==1
   Exp = load('fendage/dic_e9f_mechreg_200_coarser_relaxGV.mat');
else
   Exp = load('fendage/dic_e9f_cohesive_coarse.mat');
end

Mesh = Exp.Mesh;

nodes    = Mesh.coordinates(:,[1,2]); nnodes = size(nodes,1);
elements = Mesh.connectivity(:,[2:4]);
Ume      = Mesh.U; nt = size(Ume,2);

% Gap
gap = load('fendage/delta_U_distance.mat');
gap = ratio * gap.delta_U; lmax = size(gap,1);
ixe = 1:nt;
ygrec = 0:lmax-1; ygrec = ygrec*62/(lmax-1);

% Change the formate of Umes
Umes = zeros(2*nnodes, nt);
for i=1:nt
   Umee = Ume{i};
   Umes(1:2:2*nnodes-1,i) = Umee(:,1);
   Umes(2:2:2*nnodes,i) = Umee(:,2);
end

M = rigidModes(nodes);
Umes = Umes - M*((M'*M)\M')*Umes;

% Get the right size
Umes = ratio*Umes; nodes = ratio*nodes;

%if regularized==0 % Truncate (for visu purposes)
%   for i=1:nt
%      Umes(find(Umes(:,i)>7e-2), i)  = 7e-2;
%      Umes(find(Umes(:,i)<-7e-2), i) = -7e-2;
%   end
%end

% Vectorial gap from the measurements directly
if regularized==1
   gapv = Umes([2*14-1;2*14],:) - Umes([2*11-1;2*11],:);
   figure; hold on;
   set(gca, 'fontsize', 20);
   surf(ixe,ygrec,gap);
   shading interp;
   colorbar();
   %legend('Gap in function of the time')
else
   gapv = Umes([2*13-1;2*13],:) - Umes([2*12-1;2*12],:);
   indp = [13,128:142,3]; indm = [12,113:127,4];
   gato = zeros(2*size(indp,2),nt);
   gato(1:2:end-1,:) = - Umes(2*indp-1,:) + Umes(2*indm-1,:);
   gato(2:2:end,:)   = - Umes(2*indp,:) + Umes(2*indm,:);

   gato = min(gato,.08); gato = max(gato,-.01);

   lmax = size(indp,2);
   ixe = 1:nt;
   ygrec = 0:lmax-1; ygrec = ygrec*62/(lmax-1);
   figure; hold on;
   set(gca, 'fontsize', 20);
   surf(ixe,ygrec,gato(1:2:end-1,:));
   shading interp;
   colorbar();
   %legend('Gap in function of the time')
end

% Get the crack's size
crlen = load('fendage/crack_tip_2methods.mat');
lenFEMU = ratio*crlen.crack_FEMU;
lenIDIC = ratio*crlen.crack_IDIC;

%try
%figure;
%hold on;
%plot(lenFEMU,'Color','red');
%plot(lenIDIC,'Color','blue');
%legend('Crack length FEMU', 'Crack length IDIC');
%end

if regularized==1
   %% Manually create the boundaries
   n101 = [6,51:-1:35];         b101 = [n101(1:end-1)',n101(2:end)'];
   n102 = [34:-1:19,1];         b102 = [n102(1:end-1)',n102(2:end)'];
   n2   = [1,198:-1:164,18];    b2 = [n2(1:end-1)',n2(2:end)'];
   n3   = [18,163:-1:153,17];   b3 = [n3(1:end-1)',n3(2:end)'];
   n4   = [17,152:-1:146,16];   b4 = [n4(1:end-1)',n4(2:end)'];
   n51  = [16,145:-1:143];      b51 = [n51(1:end-1)',n51(2:end)'];
   n52  = [143,142,15];         b52 = [n52(1:end-1)',n52(2:end)'];
   n6   = [15,141:-1:137,14];   b6 = [n6(1:end-1)',n6(2:end)'];
   b1411 = [14,11]; b43 = [4,3];
   n7   = [11,114:-1:110,10];   b7 = [n7(1:end-1)',n7(2:end)'];
   n82  = [10,109,108];         b82 = [n82(1:end-1)',n82(2:end)'];
   n81  = [108,107,106,9];      b81 = [n81(1:end-1)',n81(2:end)'];
   n9   = [9,105:-1:99,8];      b9 = [n9(1:end-1)',n9(2:end)'];
   n10  = [8,98:-1:87,7];       b10 = [n10(1:end-1)',n10(2:end)'];
   n11  = [7,86:-1:52,6];       b11 = [n11(1:end-1)',n11(2:end)'];
   n121 = [35,5,4];             b121 = [n121(1:end-1)',n121(2:end)'];
   n122 = [3,2,34];             b122 = [n122(1:end-1)',n122(2:end)'];

   %% Add a few elements in the riff
   elemadd = [ 14,115,136 ; 14,115,11 ;...
               136,116,135 ; 136,116,115 ;...
               135,117,134 ; 135,117,116 ;...
               134,118,133 ; 134,118,117 ;...
               133,119,132 ; 133,119,118 ;...
               132,120,131 ; 132,120,119 ;...
               131,121,130 ; 131,121,120 ;...
               130,122,129 ; 130,122,121 ;...
               129,123,128 ; 129,123,122 ;...
               128,124,127 ; 128,124,123 ;...
               127,125,126 ; 127,125,124 ;...
               126,12,13 ; 126,12,125 ;...
               13,199,210 ; 13,199,12 ;...
               210,200,211 ; 210,200,199 ;...
               211,201,212 ; 211,201,200 ;...
               212,202,213 ; 212,202,201 ;...
               213,203,214 ; 213,203,202 ;...
               214,204,215 ; 214,204,203 ;...
               215,205,216 ; 215,205,204 ;...
               216,206,217 ; 216,206,205 ;...
               217,207,218 ; 217,207,206 ;...
               218,208,219 ; 218,208,207 ;...
               219,209,220 ; 219,209,208 ;...
               220,4,3 ; 220,4,209 ];

   elements = [ elements ; elemadd ];

   % Assemble it into 2 boundaries : bo1 : all data, bo2 : no Neumann
   bo1 = [b101;b102;b2;b3;b52;b6;b7;b82;b10;b11];
   bo2 = [b121;b43;b122;b4;b51;b1411;b81;b9];

   boundary = [ [ones(size(bo1,1),1),bo1] ; [2*ones(size(bo2,1),1),bo2] ];
else
   %% Manually create the boundaries
   n101 = [6,38:-1:29];         b101 = [n101(1:end-1)',n101(2:end)'];
   n102 = [28:-1:19,1];         b102 = [n102(1:end-1)',n102(2:end)'];
   n2   = [1,112:-1:91,18];     b2 = [n2(1:end-1)',n2(2:end)'];
   n3   = [18,90:-1:84,17];     b3 = [n3(1:end-1)',n3(2:end)'];
   n4   = [17,83:-1:80,16];     b4 = [n4(1:end-1)',n4(2:end)'];
   n51  = [16,79,78];           b51 = [n51(1:end-1)',n51(2:end)'];
   n52  = [78,15];              b52 = [n52(1:end-1)',n52(2:end)'];
   n6   = [15,77,14,76,13];     b6 = [n6(1:end-1)',n6(2:end)'];
   b1411 = [13,12]; b43 = [4,3];
   n7   = [12,75,11,74,10];     b7 = [n7(1:end-1)',n7(2:end)'];
   n82  = [10,73];              b82 = [n82(1:end-1)',n82(2:end)'];
   n81  = [73,72,9];            b81 = [n81(1:end-1)',n81(2:end)'];
   n9   = [9,71:-1:68,8];       b9 = [n9(1:end-1)',n9(2:end)'];
   n10  = [8,67:-1:61,7];       b10 = [n10(1:end-1)',n10(2:end)'];
   n11  = [7,60:-1:39,6];       b11 = [n11(1:end-1)',n11(2:end)'];
   n121 = [29,5,4];             b121 = [n121(1:end-1)',n121(2:end)'];
   n122 = [3,2,28];             b122 = [n122(1:end-1)',n122(2:end)'];

   %% Add a few elements in the riff
   elemadd = [ 13,12,113 ; 13,113,128 ; ...
               128,113,114 ; 128,114,129 ; ...
               129,114,115 ; 129,115,130 ; ...
               130,115,116 ; 130,116,131 ; ...
               131,116,117 ; 131,117,132 ; ...
               132,117,118 ; 132,118,133 ; ...
               133,118,119 ; 133,119,134 ; ...
               134,119,120 ; 134,120,135 ; ...
               135,120,121 ; 135,121,136 ; ...
               136,121,122 ; 136,122,137 ; ...
               137,122,123 ; 137,123,138 ; ...
               138,123,124 ; 138,124,139 ; ...
               139,124,125 ; 139,125,140 ; ...
               140,125,126 ; 140,126,141 ; ...
               141,126,127 ; 141,127,142 ; ...
               142,127,4 ; 142,4,3 ];

   elements = [ elements ; elemadd ];

   % Assemble it into 2 boundaries : bo1 : all data, bo2 : no Neumann
   if keep_disgus==1
      bo1 = [b101;b102;b2;b3;b52;b6;b7;b82;b10;b11];
      bo2 = [b121;b43;b122;b4;b51;b1411;b81;b9];
   else
      bo1 = [b101;b102;b2;b3;b10;b11];
      bo2 = [b121;b43;b122;b4;b51;b1411;b81;b9];
      bo3 = [b52;b6;b7;b82];
   end

   boundary = [ [ones(size(bo1,1),1),bo1] ; [2*ones(size(bo2,1),1),bo2] ];
end

%mesh2GMSH( nodes, elements, boundary, 'fendage/correlimesh' );

%[K,C,nbloq,node2c,c2node] = Krig2 (nodes,elements,mat,order,boundary,[],1);
%Kinter = K( 1:2*nnodes, 1:2*nnodes );
%Fmes = Kinter*Umes;
Fmes = zeros(size(Umes));
%plotGMSH({Umes,'U_mes';Fmes,'F_mes'}, elements, nodes, 'output/mesure');
%save('fields/F.mat','Fmes');

%% TODO : try to remove the out of plane movements

if t0 == 0
   u1 = Umes; f1 = Fmes;
else
   u1 = Umes(:,t0); f1 = Fmes(:,t0);
   nt = 1; % Only 1 time step
end

%% Find the reference angles of the crack
if regularized==1
   x15 = nodes(14,1); y15 = nodes(14,2); % To see the detail
else
   x15 = nodes(13,1); y15 = nodes(13,2); % To see the detail
end

xmin = min(nodes(:,1)); xmax = max(nodes(:,1)); Lx = xmax-xmin;
ymin = min(nodes(:,2)); ymax = max(nodes(:,2)); Ly = ymax-ymin;
xb = .5*(xmin+xmax); yb = .5*(ymin+ymax);

U = ymax-y15;

xy1r = [(xmax+xmin)/2;y15];
%xy1r = [(xmax+xmin)/2;ymin];
xy2r = [(xmax+xmin)/2;ymax];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[ nodes,elements,ntoelem2,boundary,order] = readmesh( 'meshes/fendage/fendage_direct2.msh' );
%nnodes = size(nodes,1);
%[K2,C2,nbloq2,node2c2,c2node2] = Krig2 (nodes,elements,mat,order,boundary,dirichlet,1);
%Kinter2 = K2( 1:2*nnodes, 1:2*nnodes );


%%% Pass meshes
%UF = [u1,f1];
%UR = passMesh2D(nodes, elements, nodes, elements, UF, 0);

%ur1 = UR(:,1); fr1 = UR(:,2);

%u1n = ur1;

%plotGMSH({ur1,'U1'},elements, nodes, 'output/bound');

ur1 = u1; fr1 = f1;

%% No discrete measure point stuff
nmin = 1:nnodes;

%% Build the maps of the Dirichlet and Neumann bounds
%b2nodesD = cell(4,2); % nb of boundaries * nb of dimensions
%b2nodesN = cell(4,2);
b2nodesD = []; b2nodesN = []; b2nodesTOT = [];
for i=[1:5,8:12]
   [~, b2node] = mapBound(i, boundary, nnodes);
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
nboun2 = size(boundary,1); nelem2 = size(elements,1);
boun2vol2 = zeros( nboun2, 1 ); extnorm2 = zeros( nboun2, 2 );
for i=1:nboun2
   % Volumic element
   no1 = boundary(i,2); no2 = boundary(i,3); % only with 2 nodes even if order > 1
   cand1 = rem( find(elements==no1)-1,nelem2 )+1; % find gives line + column*size
   cand2 = rem( find(elements==no2)-1,nelem2 )+1;
   boun2vol2(i) = intersect(cand1, cand2); % If everything went well, there is only one
   
   % Exterior normal
   elt = boun2vol2(i); no3 = setdiff( elements( elt, 1:3 ), [no1,no2]);
   x1 = nodes(no1,1); y1 = nodes(no1,2);
   x2 = nodes(no2,1); y2 = nodes(no2,2);
   x3 = nodes(no3,1); y3 = nodes(no3,2);
   extnorm2(i,:) = [ y1-y2 , -(x1-x2) ];
   extnorm2(i,:) = extnorm2(i,:)/norm(extnorm2(i,:));
   if (x3-x2)*(y1-y2) - (y3-y2)*(x1-x2) > 0 % Check that the normal is exterior
      extnorm2(i,:) = -extnorm2(i,:);
   end
end

nboun1 = size(boundary,1); nelem1 = size(elements,1);
boun2vol1 = zeros( nboun1, 1 ); extnorm1 = zeros( nboun1, 2 );
sigr1 = zeros( nboun1, 3*order ); sigr2 = zeros( nboun1, 3*order );
urr1  = zeros( nboun1, 2+2*order );
for i=1:nboun1 % TODO : rationalize the stuff (don't do 2 times the same work)
   % Volumic element
   no1 = boundary(i,2); no2 = boundary(i,3); % only with 2 nodes even if order > 1
   cand1 = rem( find(elements==no1)-1,nelem1 )+1; % find gives line + column*size
   cand2 = rem( find(elements==no2)-1,nelem1 )+1;
   boun2vol1(i) = intersect(cand1, cand2); % If everything went well, there is only one
   
   % Exterior normal
   elt = boun2vol1(i); no3 = setdiff( elements( elt, 1:3 ), [no1,no2]);
   x1 = nodes(no1,1); y1 = nodes(no1,2);
   x2 = nodes(no2,1); y2 = nodes(no2,2);
   x3 = nodes(no3,1); y3 = nodes(no3,2);
   extnorm1(i,:) = [ y1-y2 , -(x1-x2) ];
   extnorm1(i,:) = extnorm1(i,:)/norm(extnorm1(i,:));
   if (x3-x2)*(y1-y2) - (y3-y2)*(x1-x2) > 0 % Check that the normal is exterior
      extnorm1(i,:) = -extnorm1(i,:);
   end
   
   % ur
   urr1(i,1:4) = u1( [2*no1-1,2*no1,2*no2-1,2*no2] );
   if order == 2
      no4 = boundary(i,4);
      urr1(i,5:6) = u1( [2*no4-1,2*no4] );
   end
end

% Give a dof to each node...
nodes2dof = zeros(2*nnodes,1); dof2nodes = []; szD = 0;
for i=1:2*nnodes
   if min(size(find(b2nodesD==i)))==0
      nodes2dof(i) = szD+1;
      dof2nodes(szD+1) = i;
      szD = szD+1;
   end
end

% ...and then give a dof to each element
boun2dof = zeros(nboun2,1);
for i=1:nboun2
   boun2dof(i,1) = nodes2dof(2*boundary(i,2)-1);
   boun2dof(i,2) = nodes2dof(2*boundary(i,2));
   boun2dof(i,3) = nodes2dof(2*boundary(i,3)-1);
   boun2dof(i,4) = nodes2dof(2*boundary(i,3));
end

%% Order the dofs : force and displacement on the boundary
nodesbound2 = unique(boundary(:,[2,3])); % Nodes on the boundary
nnbound2    = size(nodesbound2,1);
boun2doftot = zeros(nboun2,2);

for i=1:nboun2 % Stores the nodes associated to each element of bound2 (in boundary numbering)
   boun2doftot(i,1) = find(nodesbound2==boundary(i,2)); % Rem : boundary(i,1) is the no of boundary
   boun2doftot(i,2) = find(nodesbound2==boundary(i,3));
end

% Local - total dof correspondancy
dofloc2tot = unique(boun2doftot(:));
doftot2loc = zeros(nnodes,1);
for i=1:nnodes
   a = find(dofloc2tot==i);
   if min(size(a))>0
      doftot2loc(i) = a;
   end
end

boun2dofloc = doftot2loc(boun2doftot);

knownD = [];
knownN = [];
for i=1:nboun2 % Chooses the known displacement dofs
   index = boundary(i,1);
   isimposed = find( dirichlet0(:,1) == index );
   dofs = dirichlet0(isimposed,2);
   knownD = [ knownD ; 2*boun2dofloc(i,1)-2+dofs ; 2*boun2dofloc(i,2)-2+dofs ];
   
   isimposed = find( neumann0(:,1) == index );
   dofs = neumann0(isimposed,2);
   knownN = [ knownN ; 2*i-2+dofs ];
end

knownD = unique(knownD); knownN = unique(knownN); % Remove redondnacy
knownD = intersect( knownD, [2*nmin-1;2*nmin] ); % Discrete measurement    /!\ There is a global/local confusion /!\

tofindD = setdiff( 1:2*nnbound2 , knownD );
tofindN = setdiff( 1:2*nboun2 , knownN );
%%

%[ nodesu,elementsu,ntoelemu,boundaryu,orderu] = readmesh( 'meshes/fendage/fendage_direct2.msh' );
%nnodesu = size(nodesu,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test-functions shit
% Build the polynomial test functions.
%load('conditions20_2dc_nu2.mat','-ascii');
%M = spconvert(conditions20_2dc_nu2); clear('conditions20_2dc');
nmax = 20;

% Build the condition matrix by hand
Mc = zeros(2*(nmax+1)^2);
for i=0:nmax
   for j=0:nmax
      index = (nmax+1)*i + j+1; nadd = (nmax+1)^2;

      if i>1
         Mc( (nmax+1)*(i-2) + j+1, index ) = (lambda+2*mu)*i*(i-1);       % a
         Mc( nadd + (nmax+1)*(i-2) + j+1, nadd + index ) = mu*i*(i-1);    % b
      end
      if j>1
         Mc( (nmax+1)*i + j-1, index ) = mu*j*(j-1);                         % a
         Mc( nadd + (nmax+1)*i + j-1, nadd + index ) = (lambda+2*mu)*j*(j-1); % b
      end
      if i>0 && j>0
         Mc( (nmax+1)*(i-1) + j, nadd + index ) = (lambda+mu)*i*j;    % b
         Mc( nadd + (nmax+1)*(i-1) + j, index ) = (lambda+mu)*i*j;    % a
      end
   end
end

ncoef =2*(nmax+1)^2; neq = ncoef;
M = Mc;
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
nelemu = size(elements,1); nn = size(elements,2); %nn=3

% Build the list of segment elements
nseg = nn*nelemu; segel = zeros(nseg,2); nbu = size(boundary,1);
for j=1:nelemu
   segel(nn*(j-1)+1,1) = elements(j,2);
   segel(nn*(j-1)+1,2) = elements(j,3);
   segel(nn*(j-1)+2,1) = elements(j,1);
   segel(nn*(j-1)+2,2) = elements(j,3);
   segel(nn*(j-1)+3,1) = elements(j,2);
   segel(nn*(j-1)+3,2) = elements(j,1);
end
j = 1;
while j <= nseg % Remove redundancy (I don't use unique because of inverted values)
   lis1 = find( segel(1:j-1,:) == segel(j,1) );
   lis1( find(lis1>j-1) ) = lis1( find(lis1>j-1) ) - j + 1;
   lis2 = find( segel(1:j-1,:) == segel(j,2) );
   lis2( find(lis2>j-1) ) = lis2( find(lis2>j-1) ) - j + 1;
   
   % also remove boundary elements
   lis3 = find( boundary(:,2:3) == segel(j,1) );
   lis3( find(lis3>nbu) ) = lis3( find(lis3>nbu) ) - nbu;
   lis4 = find( boundary(:,2:3) == segel(j,2) );
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
f_known1 = zeros(2*nboun2,nt); u_known1 = zeros(2*nnbound2,nt);
   
for i=1:nboun2
   bonod = boundary(i,1);
   fer1 = [0;0]; % All 0 for the known part
     
   indicesLoc = [2*boun2dofloc(i,:)-1,2*boun2dofloc(i,:)]; % Local Displacement dofs
   indicesGlo = [2*boundary(i,[2,3])-1,2*boundary(i,[2,3])]; % Global Displacement dofs
      
   u_known1(indicesLoc,:) = ur1(indicesGlo,:);  
   f_known1([2*i-1,2*i],:) = fer1*ones(1,nt);
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
   bonod = boundary(i,:); exno = extnorm2(i,:)';
   
   no1 = bonod(2); no2 = bonod(3); boname = bonod(1);
   x1 = nodes(no1,1); y1 = nodes(no1,2);
   x2 = nodes(no2,1); y2 = nodes(no2,2);
      
   len = sqrt( (x1-x2)^2 + (y1-y2)^2 );

   indDtot = [2*boun2dofloc(i,1)-1,2*boun2dofloc(i,1),...
              2*boun2dofloc(i,2)-1,2*boun2dofloc(i,2)]; % U dofs the element is associated to
   indNtot = [2*i-1,2*i]; % F dofs the element is associated to
                
   for j=1:Ndots
      xg = Xg(j,:); wg = Wg(j);
      xgr  = ( (1-xg)*[x1;y1] + xg*[x2;y2] ); % abscissae

      %X = xgr(1); Y = xgr(2);
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
            sloc11a = 1/Lx*(lambda+2*mu)*XY;
            sloc22a = 1/Lx*lambda*XY;
            sloc12b = 1/Lx*mu*XY;
         end
         if jj>0
            XY = jj.*Xxg.^ii.*Yyg.^(jj-1);
            sloc11b = 1/Lx*lambda*XY;
            sloc22b = 1/Lx*(lambda+2*mu)*XY;
            sloc12a = 1/Lx*mu*XY;
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
f_known1(tofindN,:) = 0;
u_known1(tofindD,:) = 0;
Piu = eye(size(u_known1,1)); Piu(tofindD,tofindD) = 0; % Projectors of measurements
Piu([2*d_nodes-1,2*d_nodes],[2*d_nodes-1,2*d_nodes]) = 0; % Remove the disgusting nodes (global/local confusion as usual)
Pif = eye(size(f_known1,1)); Pif(tofindN,tofindN) = 0;

ndofs = 2*(ndofcrack+1);
if frf==0, ndofs=ndofs+size(tofindN,2); else ndofs=ndofs+size(f_known1,1); end
if fru==0, ndofs=ndofs+size(tofindD,2); else ndofs=ndofs+size(u_known1,1); end
%ndofs = size(f_known1(tofindN,:),1) + size(u_known1(tofindD,:),1) + 2*(ndofcrack+1);

%% Restrict the operators
Rfm = Rf(:,tofindN); Rfr = Rf(:,knownN);
Rum = Ru(:,tofindD); Rur = Ru(:,knownD);

if fru==1, LhsA = -Ru; else, LhsA = -Rum; end % If there is freedom, more dof appear
if frf==1, LhsB =  Rf; else, LhsB =  Rfm; end
Rhs1 = zeros(nftest,nt);
if frf==0, Rhs1=Rhs1-Rfr*f_known1(knownN,:); end
if fru==0, Rhs1=Rhs1+Rur*u_known1(knownD,:); end
%Rhs1 = Rur*u_known1(knownD,:) - Rfr*f_known1(knownN,:);
   
% Build the matrix that passes f on the nodes from the bound
Fntob = zeros( 2*nnbound2, 2*nboun2 );
nodeMass = zeros(2*nnbound2); elemMass = zeros(2*nboun2);
for i=1:nboun2
   coef1 = [ boun2dofloc(i,1), boun2dofloc(i,2) ];
   len = norm(nodes(boundary(i,2),:) - nodes(boundary(i,3),:));
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
   Df = zeros(2*nboun2); % No derivative for f
   for i=1:nboun2
      coefU = [ boun2dofloc(i,1) , boun2dofloc(i,2) ];
      len = norm(nodes(boundary(i,2),:) - nodes(boundary(i,3),:));
      Du( 2*coefU-1, 2*coefU-1 ) = Du( 2*coefU-1, 2*coefU-1 ) + 1/len*[1,-1;-1,1];
      Du( 2*coefU, 2*coefU )     = Du( 2*coefU, 2*coefU ) + 1/len*[1,-1;-1,1];
      Df( 2*i-1, 2*i-1 ) = Df( 2*i-1, 2*i-1 ) + len^2;
      Df( 2*i, 2*i )     = Df( 2*i, 2*i ) + len^2;
   end
   %Du = eye(2*nnbound2);
   Duu  = Du(tofindD,tofindD); Dfu  = Df(tofindN,tofindN);
   Duk  = Du(knownD,knownD);   Dfk  = Df(knownN,knownN);
   Duuk = Du(tofindD,knownD);  Dfuk = Df(tofindN,knownN);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the crack RHS (not from angle)

xy1 = xy1r; xy2 = xy2r;

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
   
curv = sqrt( (nodes3(:,1)-nodes3(1,1)).^2 + (nodes3(:,2)-nodes3(1,2)).^2 );

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
         sloc11a = 1/Lx*(lambda+2*mu)*XY;
         sloc22a = 1/Lx*lambda*XY;
         sloc12b = 1/Lx*mu*XY;
      end
      if jj>0
         XY = jj.*Xxg.^ii.*Yyg.^(jj-1);
         sloc11b = 1/Lx*lambda*XY;
         sloc22b = 1/Lx*(lambda+2*mu)*XY;
         sloc12a = 1/Lx*mu*XY;
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
if froreg == 1 && min(size(LhsB))>0
   kB = sqrt(norm(LhsA,'fro')^2+norm(Ruc,'fro')^2)/norm(LhsB,'fro'); % In order to regularize the stuff
else
   kB = 1;
end

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
      if derivative == 1
   for i=1:nboun3
      coefU = [ i , i+1 ];
      len = norm(nodes3(i,:) - nodes3(i+1,:));
      D3( 2*coefU-1, 2*coefU-1 ) = D3( 2*coefU-1, 2*coefU-1 ) + 1/len*[1,-1;-1,1];
      D3( 2*coefU, 2*coefU )     = D3( 2*coefU, 2*coefU )     + 1/len*[1,-1;-1,1];
   end
else % derivative == 1/2
   for i=1:nboun3
      xg  = .5 * (nodes3(i,:) + nodes3(i+1,:));
      le1 = norm(nodes3(i,:) - nodes3(i+1,:));
      for j=1:nboun3
         if i==j, continue; end

         yg = .5 * (nodes3(j,:) + nodes3(j+1,:));
         le2 = norm(nodes3(j,:) - nodes3(j+1,:));
         len = norm(xg-yg);
         D3( 2*[i,j]-1, 2*[i,j]-1 ) = D3( 2*[i,j]-1, 2*[i,j]-1 ) + le1*le2/len^2*[1,-1;-1,1];
         D3( 2*[i,j], 2*[i,j] )     = D3( 2*[i,j], 2*[i,j] )     + le1*le2/len^2*[1,-1;-1,1];

         D3( 2*[i+1,j+1]-1, 2*[i+1,j+1]-1 ) = ...
           D3( 2*[i+1,j+1]-1, 2*[i+1,j+1]-1 ) + le1*le2/len^2*[1,-1;-1,1];
         D3( 2*[i+1,j+1], 2*[i+1,j+1] ) = ...
               D3( 2*[i+1,j+1], 2*[i+1,j+1] ) + le1*le2/len^2*[1,-1;-1,1];

         if i+1 != j
            D3( 2*[i+1,j]-1, 2*[i+1,j]-1 ) = D3( 2*[i+1,j]-1, 2*[i+1,j]-1 ) + le1*le2/len^2*[1,-1;-1,1];
            D3( 2*[i+1,j], 2*[i+1,j] )     = D3( 2*[i+1,j], 2*[i+1,j] )     + le1*le2/len^2*[1,-1;-1,1];
         end
         if j+1 != i
            D3( 2*[i,j+1]-1, 2*[i,j+1]-1 ) = D3( 2*[i,j+1]-1, 2*[i,j+1]-1 ) + le1*le2/len^2*[1,-1;-1,1];
            D3( 2*[i,j+1], 2*[i,j+1] )     = D3( 2*[i,j+1], 2*[i,j+1] )     + le1*le2/len^2*[1,-1;-1,1];
         end
      end
   end
   D3 = .25*(D3+D3');
end
D3 = norm(Dfu,'fro')/norm(D3,'fro')*D3;
Dsmall = [ 0*Mum, Z12, Z13  ; Z12', 0*Mfm, Z23 ; Z13', Z23', D3 ];

if regular == 1
   if frf==1 && fru==1
      Z12 = zeros( size(Du,1) , size(Df,2) ); Z13 = zeros( size(Du,1), size(D3,2) );
      Z23 = zeros( size(Df,1), size(D3,2) );
      Dtot = [ Du, Z12, Z13 ; Z12', Df, Z23 ; Z13', Z23', D3 ];
      L121 = zeros(ndofs,nt); L21 = zeros(nt);
   elseif frf==1
      Z12 = zeros( size(Duu,1) , size(Df,2) ); Z13 = zeros( size(Duu,1), size(D3,2) );
      Z23 = zeros( size(Df,1), size(D3,2) );
      Dtot = [ Duu, Z12, Z13 ; Z12', Df, Z23 ; Z13', Z23', D3 ];
      L12 = [ Du(tofindD,knownD) ; zeros(size(Df,2),size(Duk,2)) ; zeros(2*nboun3+2,size(Duk,2)) ];
      L121 = L12*u_known1(knownD,:); L21 = u_known1(knownD,:)'*Du(knownD,knownD)*u_known1(knownD,:);
   elseif fru==1
      Z12 = zeros( size(Du,1) , size(Dfu,2) ); Z13 = zeros( size(Du,1), size(D3,2) );
      Z23 = zeros( size(Dfu,1), size(D3,2) );
      Dtot = [ Du, Z12, Z13 ; Z12', Dfu, Z23 ; Z13', Z23', D3 ];
      L121 = zeros(ndofs,nt); L21 = zeros(nt);
   else
      Z12 = zeros( size(Duu,1) , size(Dfu,2) ); Z13 = zeros( size(Duu,1), size(D3,2) );
      Z23 = zeros( size(Dfu,1), size(D3,2) );
      Dtot = [ Duu, Z12, Z13 ; Z12', Dfu, Z23 ; Z13', Z23', D3 ];
      L12 = [ Du(tofindD,knownD) ; zeros(size(Dfu,2),size(Duk,2)) ; zeros(2*nboun3+2,size(Duk,2)) ];
      L121 = L12*u_known1(knownD,:); L21 = u_known1(knownD,:)'*Du(knownD,knownD)*u_known1(knownD,:);
   end
   L = Dtot;
else
%      L = eye(size(A));
   L = Mtot;
   L21 = 0; L22 = 0; L23 = 0; L24 = 0;
   L121 = zeros(ndofs,1);
end
ninfty = 0;
%L = eye(size(A));

sL = real(L^(1/2));   % Square root of L (as usual, real should not be)

% Assemble the projector of measurements
if frf==1 && fru==1
   Pi = [ Piu, Z12, Z13 ; Z12', Pif, Z23 ; Z13', Z23', zeros(size(D3)) ];
   Ix0 = [ u_known1 ; f_known1 ; zeros(size(D3),nt) ];% Measurements
elseif frf==1
   Pi = [ Piu, Z12, Z13 ; Z12', zeros(size(Pif)), Z23 ; Z13', Z23', zeros(size(D3)) ];
   Ix0 = [ u_known1(tofindD,:) ; f_known1 ; zeros(size(D3),nt) ];
elseif fru==1
   Pi = [ zeros(size(Piu)), Z12, Z13 ; Z12', Pif, Z23 ; Z13', Z23', zeros(size(D3)) ];
   Ix0 = [ u_known1 ; f_known1(tofindN,:) ; zeros(size(D3),nt) ];
else
   Pi = zeros(ndofs);
   Ix0 = [ u_known1(tofindD,:) ; f_known1(tofindN,:) ; zeros(size(D3),nt) ];
end
Ix0PiIx0 = Ix0'*Pi*Ix0;

x0rec = zeros(nbstep,nt);%x0*ones(1,nt);
x0b   = x0;
x0s   = zeros(nt,1);

if anglestep == 0 % Initialize the step
   anglestep1 = 100/kauto;
else
   anglestep1 = anglestep;
end

Solu1  = zeros(ndofs,nt);
Solu10 = zeros(ndofs,nt);
phi = zeros(ndofcrack,nt);

for tt = 1:nt
   for iter = 1:ndofcrack%nbstep % Newton loop
      pbmax = 1;
%   if iter == nbstep
%      pbmax == 0;
%   end
      if anglestep == 0 && iter > 1
         anglestep1 = dtheta/kauto;
      end

      for pb = 0:0%pbmax % Construct all the problems

         if pb == 0
            x0c = curv(iter+1);%x0;
         else %if pb == 1
            x0c = x0+anglestep1;
         end

         dx = x0c-curv; indp = find(dx>=0);
         s0 = zeros(nboun3+1,1);
         s0(indp) = dx(indp).^1; % Only the positive part
         s1 = zeros(ndofs,1); s2 = zeros(ndofs,1);
         s1(end-2*nboun3-1:2:end-1) = s0;
         s2(end-2*nboun3:2:end)     = s0;
         B = [ eye(ndofs,ndofs-2*nboun3-2) , s1 , s2 ];

         MAT = B'*Lhs'*Lhs*B + mur*B'*L*B;
         VE1 = B'*Lhs'*Rhs1(:,tt)-mur*B'*L121(:,tt) + gamma*B'*Pi*Ix0(:,tt);

         tokeep = 1:ndofs-2*nboun3;
         [Ll, Uu, Pp] = lu (MAT(tokeep,tokeep)); % No more zerobound

         SoluRG = zeros(size(B,2),1);
         SoluRG(tokeep) = Uu \ ( Ll \ ( Pp * ( VE1(tokeep) ) ) );

         SoluRG = B*SoluRG; % Pass onto the physical basis
         Solu1(:,tt) = SoluRG;

         if pb == 0
            nor1 = Solu1(:,tt)'*Lhs'*Lhs*Solu1(:,tt) - 2*Solu1(:,tt)'*Lhs'*Rhs1(:,tt) + Rhs1(:,tt)'*Rhs1(:,tt);
            nor2 = mur*Solu1(:,tt)'*L*Solu1(:,tt) + 2*mur*Solu1(:,tt)'*L121(:,tt) + mur*L21(tt,tt);
            nor3 = gamma*Solu1(:,tt)'*Pi*Solu1(:,tt) - 2*gamma*Solu1(:,tt)'*Pi*Ix0(:,tt) + gamma*Ix0PiIx0(tt,tt);
            phi( iter, tt )  = nor1+nor2+nor3;
            phi0( iter, tt ) = Rhs1(:,tt)'*Rhs1(:,tt) + mur*L21(tt,tt) + gamma*Ix0PiIx0(tt,tt) ; % Worst residual possible

            regu = nor2/mur; resu = nor1; % Regularization norm and residual
            remu = nor3/gamma;

            res1 = Lhs*Solu1(:,tt) - Rhs1(:,tt); % Residuals
            rel1 = sL*Solu1(:,tt);

            res = res1;
            rel = rel1;
            Ax  = Lhs*Solu1(:,tt);
            xt  = Solu1(:,tt);
            sLx = sL*Solu1(:,tt);
            L1x = Solu1(:,tt)'*L121(:,tt);

            if phi(iter,tt) == min(phi(1:iter,tt))
               Solu10(:,tt) = Solu1(:,tt); % Store the values
               x0s(tt) = x0c;
               nodes3b = nodes3;
            end

         elseif pb == 1
            D1    = (Lhs*Solu1(:,tt) - Ax)/anglestep1;
            Dx1   = (Solu1(:,tt) - xt)/anglestep1;
            DL1   = (sL*Solu1(:,tt) - sLx)/anglestep1;
            DL121 = (Solu1(:,tt)'*L121(:,tt) - L1x)/anglestep1;
         end
      end
%      D = D1; DL = DL1; DL12 = DL121;

%      dtheta = - ( D'*D + mur*DL'*DL ) \ ( D'*res + mur*DL'*rel + mur*DL12'*ones(size(DL12,1),1) );

%      x0 = x0 + dtheta;
%      x0rec(iter,tt) = x0;
   end
   %x0
end
disp(['Iterative method terminated ', num2str(toc) ]);

%figure; plot(phi);

Solu1 = Solu10;

% Post-process : reconstruct u and f from Solu
fsolu1 = zeros(2*nboun2,nt); usolu1 = zeros(2*nnodes,nt);
szB = 0;
szD = size(LhsA,2);
for i=1:szD
   usolu1(dof2nodes(i),:) = Solu1(i);
end
for i=1:nboun2
   boname = boundary(i,1);
   
   if max(size(neumann0))>0
      frdof = neumann0(find(neumann0(:,1)==boname),2) ;  % Dof of the given f
   else
      frdof = [];
   end
   fmdof = setdiff([1;2],frdof);  % Dof of the unknown f

   indicesB = szD + (szB+1:szB+size(fmdof,1));
   fsolu1(2*i-2+fmdof,:) = kB*Solu1(indicesB,:);
   szB = szB + size(fmdof,1);
end

ucrsol1 = Solu1(end-2*nnodes3+1:end,:);

% Compute the values of f at the dofs
fsolu1no = zeros(2*nnodes,nt);
for i=1:size(b2nodesnoN)
   b2nodesnoNi = b2nodesnoN(i);
   noode  = floor((b2nodesnoNi+1)/2); % node
   boound = rem( find(boundary(:,2:3)==noode)-1, nboun2 )+1; % boundaries
   
   % Compute lengths
   n1 = boundary(boound(1),2); n2 = boundary(boound(1),3);
   len1 = sqrt( (nodes(n1,1)-nodes(n2,1))^2 + (nodes(n1,2)-nodes(n2,2))^2 );

   % Case when there is no other bound with this node
   if max(size(boound))<2
      if rem(b2nodesnoNi,2)==0
         fsolu1no(b2nodesnoNi,:) = fsolu1(2*boound(1),:)*len1;
      else
         fsolu1no(b2nodesnoNi,:) = fsolu1(2*boound(1)-1,:)*len1;
      end
   else
      n1 = boundary(boound(2),2); n2 = boundary(boound(2),3);
      len2 = sqrt( (nodes(n1,1)-nodes(n2,1))^2 + (nodes(n1,2)-nodes(n2,2))^2 );
      if rem(b2nodesnoNi,2)==0
         fsolu1no(b2nodesnoNi,:) = .5*( fsolu1(2*boound(1),:)*len1 + fsolu1(2*boound(2),:)*len2 );
      else
         fsolu1no(b2nodesnoNi,:) = .5*( fsolu1(2*boound(1)-1,:)*len1 + fsolu1(2*boound(2)-1,:)*len2 );
      end
   end
end

%% Graph for u
%toplot = usolu1(b2nodesnoD); toplot2 = ur1(b2nodesnoD);
%figure;
%hold on;
%plot(toplot(2:2:end),'Color','red');
%plot(toplot2(2:2:end),'Color','blue');
%legend('Uy identified', 'Uy reference');

%% Recover the reference
step = (xy2r-xy1r)/30;
n = [-step(2);step(1)]; n = n/norm(n); % Normal
if step(1)==0
   nodes3r = [ xy1r(1)*ones(1,31) ; xy1r(2):step(2):xy2r(2) ];
elseif step(2)==0
   nodes3r = [ xy1r(1):step(1):xy2r(1) ; xy1r(2)*ones(1,31) ];
else
   nodes3r = [ xy1r(1):step(1):xy2r(1) ; xy1r(2):step(2):xy2r(2) ];
end
nodes3r = nodes3r'; nnodes3r = size(nodes3r,1);
% Check if need to reverse the nodes (for visu purposes)
if norm(nodes3r(end,:)-nodes3b(1,:)) < norm(nodes3r(1,:)-nodes3b(1,:))
   nodes3r = nodes3r(end:-1:1,:);
end
nodes3s = nodes3r + 3*ones(size(nodes3r,1),1)*n';
nodes3i = nodes3r - 3*ones(size(nodes3r,1),1)*n';
urs = passMesh2D( nodes, elements, nodes3s, [], u1 );
uri = passMesh2D( nodes, elements, nodes3i, [], u1 );
urg = -(uri-urs)*(extnorm3(1,:)*n);  % Vectorial gap
%urg([1,2,end-1,end],:) = 0; % Overwrite the strange stuff that comes from the fact that we get out of the domain
curvr = sqrt( (nodes3r(:,1)-nodes3r(1,1)).^2 + (nodes3r(:,2)-nodes3r(1,2)).^2 );

%if nt==1
%   % Graph for f
%   if min(size(b2nodesnoN))>0
%      toplot = fsolu1no(b2nodesnoN,:); toplot2 = fr1(b2nodesnoN,:);
%      try
%      figure;
%      hold on;
%      plot(toplot(1:2:end-1),'Color','red');
%      plot(toplot2(1:2:end-1),'Color','blue');
%      legend('Fy identified (1)', 'Fy reference (1)');
%      end
%   end
%   % Graph for u
%   if min(size(b2nodesnoD))>0
%      toplot = usolu1(b2nodesnoD); toplot2 = ur1(b2nodesnoD);
%      try
%      figure;
%      hold on;
%      plot(toplot(1:2:end-1),'Color','red');
%      plot(toplot2(1:2:end-1),'Color','blue');
%      legend('Uy identified (1)', 'Uy reference (1)');
%      end
%   end
%else
%   % Graph for f
%   if min(size(b2nodesnoN))>0
%      toplot = fsolu1no(b2nodesnoN,:); toplot2 = fr1(b2nodesnoN,:);
%      try
%      figure;
%      hold on;
%      surf(toplot(1:2:end-1));
%      shading interp;
%      colorbar();
%      end
%   end
%   % Graph for u
%   if min(size(b2nodesnoD))>0
%      toplot = usolu1(b2nodesnoD,:); toplot2 = ur1(b2nodesnoD,:);
%      try
%      figure;
%      hold on;
%      surf(toplot(1:2:end-1,:));
%      shading interp;
%      colorbar();
%      end
%   end
%end

% Graph for [[u]]
curv = sqrt( (nodes3b(:,1)-nodes3b(1,1)).^2 + (nodes3b(:,2)-nodes3b(1,2)).^2 );

n = extnorm3(1,:)';
toplot1  = ucrsol1(1:2:end-1,:)*n(1) + ucrsol1(2:2:end,:)*n(2);
toplotr1 = urg(1:2:end-1,:)*n(1) + urg(2:2:end,:)*n(2);

if nt==1
   try
   figure;
   hold on;
   set(gca, 'fontsize', 20);
   plot( curv, toplot1,'Color','red','LineWidth',3 );
   plot( curvr, toplotr1,'Color','blue','LineWidth',3 );
   legend('[[u]] identified (1)', '[[u]] reference (1)');
   end
else
   try
   figure;
   hold on;
   set(gca, 'fontsize', 20);
   surf( 1:nt, curv, toplot1 );
   shading interp;
   colorbar();
   end
   try
   figure;
   hold on;
   set(gca, 'fontsize', 20);
   plot( curv, toplot1(:,263),'Color','red','LineWidth',3 );
   plot( curv, toplot1(:,215),'Color','magenta','LineWidth',3 );
   plot( curv, toplot1(:,176),'Color','blue','LineWidth',3 );
   plot( curv, toplot1(:,100),'Color','green','LineWidth',3 );
   plot( curv, toplot1(:,45),'Color','black','LineWidth',3 );
%   plot( curv, toplot1(:,215)*toplot1(1,263)/toplot1(1,215),'Color','magenta','LineWidth',3 );
%   plot( curv, toplot1(:,176)*toplot1(1,263)/toplot1(1,176),'Color','blue','LineWidth',3 );
%   plot( curv, toplot1(:,100)*toplot1(1,263)/toplot1(1,100),'Color','green','LineWidth',3 );
%   plot( curv, toplot1(:,45)*toplot1(1,263)/toplot1(1,45),'Color','black','LineWidth',3 );
   legend('[[u]] 263', '[[u]] 215', '[[u]] 176', '[[u]] 100', '[[u]] 45');
   end
   try
   figure;
   hold on;
   set(gca, 'fontsize', 20);
   plot( curvr, toplotr1(:,263),'Color','red','LineWidth',3 );
   plot( curvr, toplotr1(:,215),'Color','magenta','LineWidth',3 );
   plot( curvr, toplotr1(:,176),'Color','blue','LineWidth',3 );
   plot( curvr, toplotr1(:,100),'Color','green','LineWidth',3 );
   plot( curvr, toplotr1(:,45),'Color','black','LineWidth',3 );
%   plot( curvr, toplotr1(:,215)*toplotr1(1,263)/toplotr1(1,215),'Color','magenta','LineWidth',3 );
%   plot( curvr, toplotr1(:,176)*toplotr1(1,263)/toplotr1(1,176),'Color','blue','LineWidth',3 );
%   plot( curvr, toplotr1(:,100)*toplotr1(1,263)/toplotr1(1,100),'Color','green','LineWidth',3 );
%   plot( curvr, toplotr1(:,45)*toplotr1(1,263)/toplotr1(1,45),'Color','black','LineWidth',3 );
   legend('[[u]] 263', '[[u]] 215', '[[u]] 176', '[[u]] 100', '[[u]] 45');
   end
end

%% Compute the derivatives and such of toplot1
i=1:ndofcrack+1; im1 = 1:ndofcrack; ip1 = 2:ndofcrack+1;
M1 = sparse(im1,im1,-1,ndofcrack,ndofcrack+1);
M2 = sparse(im1,ip1,1,ndofcrack,ndofcrack+1);
Md1 = M1+M2; % Derivative matrix

toplod = Md1*toplot1;

i=1:ndofcrack; im1 = 1:ndofcrack-1; ip1 = 2:ndofcrack;
M1 = sparse(im1,im1,-1,ndofcrack-1,ndofcrack);
M2 = sparse(im1,ip1,1,ndofcrack-1,ndofcrack);
Md2 = M1+M2; % Derivative matrix

toplodd = Md2*toplod;

%if nt==1
%   try
%   figure;
%   hold on;
%   set(gca, 'fontsize', 20);
%   plot( curv(1:end-1), toplod,'Color','red','LineWidth',3 );
%   %plot( curvr, toplotr1,'Color','blue','LineWidth',3 );
%   legend('Derivative of [[u]] identified (1)');
%   end
%else
%   try
%   figure;
%   hold on;
%   set(gca, 'fontsize', 20);
%   surf( 1:nt, curv(1:end-1), toplod ); % /!\ curv(1:end-1) is not very clean
%   shading interp;
%   colorbar();
%   end
%end

%if nt==1
%   try
%   figure;
%   hold on;
%   set(gca, 'fontsize', 20);
%   plot( curv(2:end-1), toplodd,'Color','red','LineWidth',3 );
%   legend('Derivative of [[u]] identified (1)');
%   end
%else
%   try
%   figure;
%   hold on;
%   set(gca, 'fontsize', 20);
%   surf( 1:nt, curv(2:end-1), toplodd ); % /!\ curv(2:end-1) is not very clean
%   shading interp;
%   colorbar();
%   end
%end

% Find the last negative second derivative point / ceil data
lcrack = x0s;

%% The same according to the reference (gap function)
%lcrackr = zeros(nt,1);
%for i=1:nt
%%   indexmaxr = max(find(toplodd(:,i)<0));
%   epscrit = max(gap(:))/3;
%   %epscrit = max(gap(:,i))/100;
%   indexmaxr = min( find( gap(:,i) < epscrit ) );
%   if min(size(indexmaxr))==0
%      lcrackr(i) = 0;
%   else
%%   lcrack(i) = curvdd(indexmax);
%      lcrackr(i) = ygrec(indexmaxr);
%   end
%end

try
figure;
hold on;
plot(lenFEMU,'Color','red');
plot(lenIDIC,'Color','blue');
plot(lcrack,'Color','black');
%plot(lcrackr,'Color','green');
%legend('Crack length FEMU', 'Crack length IDIC', 'Crack length RG', 'Crack length gap');
legend('Crack length FEMU', 'Crack length IDIC', 'Crack length RG');
end


%erroru = norm(usolu1(b2nodesnoD)-ur1(b2nodesnoD))   / norm(ur1(b2nodesnoD));
%errorf = norm(fsolu1no(b2nodesnoN)-fr1(b2nodesnoN)) / norm(fr1(b2nodesnoN));
