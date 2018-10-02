% 25/09/2018
% Essai de fendage, 1 seule donnée, code simplifié : L-curve

tic
close all;
clear all;

addpath(genpath('./tools'));

% Parameters
E           = 17000; % MPa : Young modulus
nu          = 0.2;    % Poisson ratio
fscalar     = 1;    % N.mm-1 : Loading on the plate
mat         = [0, E, nu];
regular     = 1;      % Use the derivative regularization matrix (0 : Id, 1 : derivative)
froreg      = 1;      % Frobenius preconditioner
Npg         = 2;      % Nb Gauss points
nuzawa      = 1;     % (nuzawa = 1 means no Uzawa (and nuzawa = 0 means)
kuzawa      = 0;%1e2;     % Parameters of the Uzawa algorithm (kuzawa = 0 means best parameter)
ndofcrack   = 50;      % Nb of elements on the crack
ratio       = 0.062;   % mm/pix ratio
order       = 1;       % Order of the mesh
th          = 72.5;   % Thickness (unused as it is only for the loadings)
t0          = 0;%263;     % Time for the computation (263 is good) 0 means all
gapbound    = 0; % Take the gap on the known extremity into account

murs = 10.^[6,7,8,9,10,11,12,13]; % Regularization parameters
%murs = [2e8,4e8,6e8,8e8,1e9,2e9,4e9,6e9,8e9,1e10,2e10,4e10];
%murs = [9e8,1e9,1.5e9,2e9,3e9,4e9,5e9,6e9,7e9,8e9,9e9,1e10];
%murs = [1e5,1e6,1e7,1e8,2e8,4e8,6e8,8e8,1e9,2e9,4e9,6e9,8e9,1e10,2e10,4e10,1e11,1e12,1e13];
%murs = [1e5,1e6,1e7,1e8,2e8,4e8,6e8,8e8,9e8,1e9,2e9,3e9,4e9,5e9,6e9,7e9,8e9,9e9,1e10,2e10,4e10,1e11,1e12,1e13];

% Material properties
lambda = nu*E/((1+nu)*(1-2*nu));
mu     = E/(2*(1+nu));

% Boundary conditions
dirichlet  = [ 12,1,0 ; 12,2,0 ];
neumann1   = [ 4,1,fscalar ; 4,2,-fscalar ; 9,1,-fscalar ; 9,2,-fscalar ];

dirichlet0 = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ];
neumann0   = [ 1,1,0 ; 1,2,0 ];

%% Import the experimental data
Exp = load('fendage/dic_e9f_mechreg_200_coarser_relaxGV.mat');
Mesh = Exp.Mesh;

nodes    = Mesh.coordinates(:,[1,2]); nnodes = size(nodes,1);
elements = Mesh.connectivity(:,[2:4]);
Ume      = Mesh.U; nt = size(Ume,2);

% Gap
gap = load('fendage/delta_U_distance.mat');
gap = ratio * gap.delta_U; lmax = size(gap,1);
ixe = 1:nt;
ygrec = 0:lmax-1; ygrec = ygrec*62/(lmax-1);
figure;
surf(ixe,ygrec,gap);
shading interp;
colorbar();
legend('Gap in function of the time')

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

% Vectorial gap from the measurements directly
gapv = Umes([2*14-1;2*14],:) - Umes([2*11-1;2*11],:);

% Get the crack's size
crlen = load('fendage/crack_tip_2methods.mat');
lenFEMU = ratio*crlen.crack_FEMU;
lenIDIC = ratio*crlen.crack_IDIC;

try
figure;
hold on;
plot(lenFEMU,'Color','red');
plot(lenIDIC,'Color','blue');
legend('Crack length FEMU', 'Crack length IDIC');
end

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
x15 = nodes(14,1); y15 = nodes(14,2); % To see the detail

xmin = min(nodes(:,1)); xmax = max(nodes(:,1)); Lx = xmax-xmin;
ymin = min(nodes(:,2)); ymax = max(nodes(:,2)); Ly = ymax-ymin;
xb = .5*(xmin+xmax); yb = .5*(ymin+ymax);

U = ymax-y15;

xy1r = [(xmax+xmin)/2;y15];
xy2r = [(xmax+xmin)/2;ymax];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
   cand1 = rem( find(elements==no1),nelem1 ); % find gives line + column*size
   cand2 = rem( find(elements==no2),nelem1 );
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

ndofs = size(f_known1(tofindN,:),1) + size(u_known1(tofindD,:),1) + 2*(ndofcrack+1);

%% Restrict the operators
Rfm = Rf(:,tofindN); Rfr = Rf(:,knownN);
Rum = Ru(:,tofindD); Rur = Ru(:,knownD);

LhsA = -Rum;
LhsB = Rfm;
Rhs1 = Rur*u_known1(knownD,:) - Rfr*f_known1(knownN,:);
   
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
   Df = eye(2*nboun2); % No derivative for f
   for i=1:nboun2
      coefU = [ boun2dofloc(i,1) , boun2dofloc(i,2) ];
      len = norm(nodes(boundary(i,2),:) - nodes(boundary(i,3),:));
      Du( 2*coefU-1, 2*coefU-1 ) = Du( 2*coefU-1, 2*coefU-1 ) + 1/len*[1,-1;-1,1];
      Du( 2*coefU, 2*coefU )     = Du( 2*coefU, 2*coefU ) + 1/len*[1,-1;-1,1];
   end

   Duu  = Du(tofindD,tofindD); Dfu  = Df(tofindN,tofindN);
   Duk  = Du(knownD,knownD);   Dfk  = Df(knownN,knownN);
   Duuk = Du(tofindD,knownD);  Dfuk = Df(tofindN,knownN);
end

nbnotwork = 0; % Records the number of failures
previous = []; % Stores the index of the previous solutions

nbstep = max(size(murs));
regu = zeros(nbstep,1); resu = zeros(nbstep,1);

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
   L = Dtot;% + Mtot * norm(Dtot,'fro')/(10*norm(Mtot,'fro'));

   L12 = [ Du(tofindD,knownD) ; zeros(size(Dfu,2),size(Duk,2)) ; zeros(2*nboun3+2,size(Duk,2)) ];
   L2 = Du(knownD,knownD);
   L121 = L12*u_known1(knownD,:); L21 = u_known1(knownD,:)'*Du(knownD,knownD)*u_known1(knownD,:);
   if gapbound == 1 % Add gapv as reference on the left
      totreg = [ndofs-2*nnodes3+1, ndofs-2*nnodes3+2];
      L121o = L121; L21o = L21; % Store those terms for the computation of the cost-function
      L121 = L121 + [ zeros(size(Duu,2),2) ; zeros(size(Dfu,2),2) ; D3(:,[1,2]) ] * gapv;
      L21 = L21 + gapv'*D3([1,2],[1,2])*gapv;
   end
else

   L = Mtot;
   L21 = 0; L22 = 0; L23 = 0; L24 = 0;
   L121 = zeros(ndofs,1);
end
ninfty = 0;
sL = real(L^(1/2));   % Square root of L (as usual, real should not be)

% Build the contact inequation condition
C = zeros(nnodes3, 2*nnodes3);
C(1,2*1-1) = extnorm3(1,1); C(1,2*1) = extnorm3(1,2);
C(nnodes3,2*nnodes3-1) = extnorm3(end,1); C(nnodes3,2*nnodes3) = extnorm3(end,2);
for i=2:nnodes3-1
   C(i,2*i-1) = .5*(extnorm3(i-1,1) + extnorm3(i,1));
   C(i,2*i)   = .5*(extnorm3(i-1,2) + extnorm3(i,2));
end
C = [ zeros(nnodes3,ndofs-2*nnodes3), C ];

Solu10 = zeros(size(L,1),nt,nbstep);

for iter = 1:nbstep % L-curve loop
   mur = murs(iter);
   for pb = 0:0

      MAT = Lhs'*Lhs + mur*L;
      VE1 = Lhs'*Rhs1 - mur*L121;

      f = zeros(nnodes3,nt); Ctf = C'*f;
      respos = zeros(nuzawa,1); df = zeros(nuzawa);

      if gapbound == 1
         toremove = [ ndofs-2*nnodes3+2, ndofs-2*nnodes3+1, ndofs, ndofs-1 ]; % Remove first and last node
      else
         toremove = [ ndofs, ndofs-1 ]; % Remove the last node
      end

      tokeep = setdiff(1:size(L,1),toremove);
      [Ll, Uu, Pp] = lu (MAT(tokeep,tokeep)); % No more zerobound
      kuzawa1 = .999*min(eig(MAT(tokeep,tokeep)));

      SoluRG = zeros(size(L,1),nt);

      for i=1:nuzawa % Uzawa for the contact problems
         SoluRG(tokeep,:) = Uu \ ( Ll \ ( Pp * ( VE1(tokeep,:) + Ctf(tokeep,:) ) ) );
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

      if gapbound == 1
         SoluRG(totreg,:) = gapv;
      end

      Solu1 = SoluRG;

      if pb == 0
         nor1 = diag( Solu1'*Lhs'*Lhs*Solu1 - 2*Solu1'*Lhs'*Rhs1 + Rhs1'*Rhs1);
         if gapbound == 1 % DOn't count 2 times the additionnal terms
            nor2 = diag( mur*Solu1'*L*Solu1 + 2*mur*Solu1'*L121o + mur*L21o  );
         else
            nor2 = diag( mur*Solu1'*L*Solu1 + 2*mur*Solu1'*L121 + mur*L21  );
         end
         phi( iter )  = sum(nor1+nor2);

         regu(iter) = sum(nor2)/mur;
         resu(iter) = sum(nor1); % Regularization norm and residual

         Solu10(:,:,iter) = Solu1; % Store the values

         if phi(iter) == min(phi)
            nodes3b = nodes3;
         end

      end
   end
end
disp(['Iterative method terminated ', num2str(toc) ]);

%% L-curve
zebest = findCorner2 (resu, regu, 3, 1); bestmu = murs(zebest);
try
figure; hold on;
loglog(resu,regu,'-*');
loglog(resu(zebest),regu(zebest),'o');
end
Solu1 = Solu10(:,:,zebest);

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
   plot( curv, toplot1(:,176),'Color','blue','LineWidth',3 );
   plot( curv, toplot1(:,100),'Color','green','LineWidth',3 );
   plot( curv, toplot1(:,45),'Color','black','LineWidth',3 );
%   plot( curv, toplot1(:,176)*toplot1(1,263)/toplot1(1,176),'Color','blue','LineWidth',3 );
%   plot( curv, toplot1(:,100)*toplot1(1,263)/toplot1(1,100),'Color','green','LineWidth',3 );
%   plot( curv, toplot1(:,45)*toplot1(1,263)/toplot1(1,45),'Color','black','LineWidth',3 );
   legend('[[u]] 263', '[[u]] 176', '[[u]] 100', '[[u]] 45');
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

% Find the last negative second derivative point / ceil data
lcrack = zeros(nt,1);
curvdd = curv(2:end-1); % /!\ curv(1:end-1) is not very clean
for i=1:nt
%   indexmax = max(find(toplodd(:,i)<0));
   epscrit = max(toplot1(:))/3;
   %epscrit = max(toplot1(:,i))/3;
   indexmax = min( find( toplot1(:,i) < epscrit ) );
   if min(size(indexmax))==0
      lcrack(i) = 0;
   else
%   lcrack(i) = curvdd(indexmax);
      lcrack(i) = curv(indexmax);
   end
end

% The same according to the reference (gap function)
lcrackr = zeros(nt,1);
for i=1:nt
%   indexmaxr = max(find(toplodd(:,i)<0));
   epscrit = max(gap(:))/3;
   %epscrit = max(gap(:,i))/100;
   indexmaxr = min( find( gap(:,i) < epscrit ) );
   if min(size(indexmaxr))==0
      lcrackr(i) = 0;
   else
%   lcrack(i) = curvdd(indexmax);
      lcrackr(i) = ygrec(indexmaxr);
   end
end

try
figure;
hold on;
plot(lenFEMU,'Color','red');
plot(lenIDIC,'Color','blue');
plot(lcrack,'Color','black');
plot(lcrackr,'Color','green');
legend('Crack length FEMU', 'Crack length IDIC', 'Crack length RG', 'Crack length gap');
end


%erroru = norm(usolu1(b2nodesnoD)-ur1(b2nodesnoD))   / norm(ur1(b2nodesnoD));
%errorf = norm(fsolu1no(b2nodesnoN)-fr1(b2nodesnoN)) / norm(fr1(b2nodesnoN));
