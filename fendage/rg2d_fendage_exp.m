% 27/06/2018
% Essai de fendage, 1 seule donnée, code simplifié

tic
close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E           = 17000; % MPa : Young modulus
nu          = 0.2;    % Poisson ratio
fscalar     = 1;    % N.mm-1 : Loading on the plate
mat         = [0, E, nu];
mur         = 1e4;%2e3;    % Regularization parameter
regular     = 1;      % Use the derivative regularization matrix (0 : Id, 1 : derivative)
froreg      = 1;      % Frobenius preconditioner
theta1      = pi/2;
theta2      = 3*pi/2;   % angles of the crack
anglestep   = 0;%pi/1000;  % Step in angle for Finite Differences anglestep = 0 means auto-adaptation
kauto       = 10;     % Coefficient for the auto-adaptation
nbstep      = 1;     % Nb of Newton Iterations
Npg         = 2;      % Nb Gauss points
%ordertest   = 20;     % Order of test fonctions (by default 20)
nuzawa      = 100;     % (nuzawa = 1 means no Uzawa (and nuzawa = 0 means)
kuzawa      = 0;%1e2;     % Parameters of the Uzawa algorithm (kuzawa = 0 means best parameter)
ndofcrack   = 20;      % Nb of elements on the crack
ratio       = 0.062;   % mm/pix ratio
order       = 1;       % Order of the mesh

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
gap = gap.delta_U;
figure;
plot(gap(1,:))
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

mesh2GMSH( nodes, elements, boundary, 'fendage/correlimesh' );

[K,C,nbloq,node2c,c2node] = Krig2 (nodes,elements,mat,order,boundary,[],1);
Kinter = K( 1:2*nnodes, 1:2*nnodes );
Fmes = Kinter*Umes;

plotGMSH({Umes,'U_mes';Fmes,'F_mes'}, elements, nodes, 'output/mesure');

% No 263 is good
u1 = Umes(:,263); f1 = Fmes(:,263);

%% Find the reference angles of the crack
x15 = nodes(14,1); y15 = nodes(14,2); % To see the detail

xmin = min(nodes(:,1)); xmax = max(nodes(:,1)); Lx = xmax-xmin;
ymin = min(nodes(:,2)); ymax = max(nodes(:,2)); Ly = ymax-ymin;
xb = .5*(xmin+xmax); yb = .5*(ymin+ymax);

U = ymax-y15;

xy1r = [(xmax+xmin)/2;y15];
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
   cand1 = rem( find(elements==no1),nelem2 ); % find gives line + column*size
   cand2 = rem( find(elements==no2),nelem2 );
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
%% Test-functions shit
% Build the polynomial test functions.
load('conditions20_2dc_nu2.mat','-ascii');
M = spconvert(conditions20_2dc_nu2); clear('conditions20_2dc');
nmax = 20;

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
f_known1 = zeros(2*nboun2,1); u_known1 = zeros(2*nnbound2,1);
   
for i=1:nboun2
   bonod = boundary(i,1);
   fer1 = [0;0]; % All 0 for the known part
     
   indicesLoc = [2*boun2dofloc(i,:)-1,2*boun2dofloc(i,:)]; % Local Displacement dofs
   indicesGlo = [2*boundary(i,[2,3])-1,2*boundary(i,[2,3])]; % Global Displacement dofs
      
   u_known1(indicesLoc) = ur1(indicesGlo);  
   f_known1([2*i-1,2*i]) = fer1;
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
f_known1(tofindN) = 0;
u_known1(tofindD) = 0;

ndofs = size(f_known1(tofindN),1) + size(u_known1(tofindD),1) + 2*(ndofcrack+1);

%% Restrict the operators
Rfm = Rf(:,tofindN); Rfr = Rf(:,knownN);
Rum = Ru(:,tofindD); Rur = Ru(:,knownD);

LhsA = -Rum;
LhsB = Rfm;
Rhs1 = Rur*u_known1(knownD) - Rfr*f_known1(knownN);
   
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

theta1rec = theta1;
theta2rec = theta2;
theta1b   = theta1;
theta2b   = theta2;
nbnotwork = 0; % Records the number of failures
previous = []; % Stores the index of the previous solutions

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

   for pb = 0:0%pbmax % Construct all the problems

      if pb == 0
         theta1c = theta1; theta2c = theta2;
      elseif pb == 1
         theta1c = theta1+anglestep1; theta2c = theta2;
      else
         theta1c = theta1; theta2c = theta2+anglestep2;
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

            %X = xgr(1); Y = xgr(2);
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
%      if iter == 1 && pb == 0 % kB must be defined only once
         if froreg == 1 && min(size(LhsB))>0
            kB = sqrt(norm(LhsA,'fro')^2+norm(Ruc,'fro')^2)/norm(LhsB,'fro'); % In order to regularize the stuff
         %   kB = norm(LhsA,'fro')/norm(LhsB,'fro');
         else
            kB = 1;
         end
%      end
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
         L121 = L12*u_known1(knownD); L21 = u_known1(knownD)'*Du(knownD,knownD)*u_known1(knownD);
      else
%      L = eye(size(A));
         L = Mtot;
         L21 = 0; L22 = 0; L23 = 0; L24 = 0;
         L121 = zeros(ndofs,1);
      end
      ninfty = 0;
   %L = eye(size(A));

      sL = real(L^(1/2));   % Square root of L (as usual, real should not be)

      MAT = Lhs'*Lhs + mur*L;
      VE1 = Lhs'*Rhs1-mur*L121;

      % Build the contact inequation condition
      C = zeros(nnodes3, 2*nnodes3);
      C(1,2*1-1) = extnorm3(1,1); C(1,2*1) = extnorm3(1,2);
      C(nnodes3,2*nnodes3-1) = extnorm3(end,1); C(nnodes3,2*nnodes3) = extnorm3(end,2);
      for i=2:nnodes3-1
         C(i,2*i-1) = .5*(extnorm3(i-1,1) + extnorm3(i,1));
         C(i,2*i)   = .5*(extnorm3(i-1,2) + extnorm3(i,2));
      end
      C = [ zeros(nnodes3,ndofs-2*nnodes3), C ];
      f = zeros(nnodes3,1); Ctf = C'*f;
      respos = zeros(nuzawa,1); df = zeros(nuzawa,1);

%      if zerobound == 1
         %toremove = [ ndofs, ndofs-1, ndofs-2*nnodes3+2, ndofs-2*nnodes3+1 ];
         toremove = [ ndofs, ndofs-1 ]; % Remove the last node
         tokeep = setdiff(1:size(L,1),toremove);
         [Ll, Uu, Pp] = lu (MAT(tokeep,tokeep)); % No more zerobound
         kuzawa1 = .999*min(eig(MAT(tokeep,tokeep)));
%      else
%         [Ll, Uu, Pp] = lu (MAT);  % Pp*MAT = Ll*Uu
%         kuzawa1 = .999*min(eig(MAT));
%      end

      SoluRG = zeros(size(L,1),1);

      for i=1:nuzawa % Uzawa for the contact problems
         SoluRG(tokeep,:) = Uu \ ( Ll \ ( Pp * ( VE1(tokeep) + Ctf(tokeep,:) ) ) );
         %SoluRG = Uu \ ( Ll \ ( Pp * ( VE1 + Ctf ) ) );
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

      Solu1 = SoluRG(:,1);

      if pb == 0
         nor1 = Solu1'*Lhs'*Lhs*Solu1 - 2*Solu1'*Lhs'*Rhs1 + Rhs1'*Rhs1 + mur*Solu1'*L*Solu1 + 2*mur*Solu1'*L121 + mur*L21;
         phi( iter )  = nor1;

%         res1 = Lhs*Solu1 - Rhs1; % Residuals

%         rel1 = sL*Solu1;

%         res = [ res1 ];
%         rel = [ rel1 ];
%         Ax  = [ Lhs*Solu1 ];
%         xt  = [ Solu1 ];
%         sLx = [ sL*Solu1 ];
%         L1x = [ Solu1'*L121 ];
%         Lh0 = Lhs; sL0 = sL; L1210 = L121;

         Solu10 = Solu1; % Store the values

         if phi(iter) == min(phi)
            nodes3b = nodes3;
         end

%      elseif pb == 1
%         D1    = ([ Lhs*Solu1 ] - Ax)/anglestep1;
%         Dx1   = ([ Solu1 ] - xt)/anglestep1;
%         DL1   = ([ sL*Solu1 ] - sLx)/anglestep1;
%         DL121 = ([ Solu1'*L121 ] - L1x)/anglestep1;
%         Lh1   = (Lhs-Lh0)/anglestep1;

%%         D1    = ((Lhs-Lh0) * [Solu10,Solu20,Solu30,Solu40] )/anglestep1; D1 = D1(:);
%%         DL1   = ((sL -sL0) * [Solu10,Solu20,Solu30,Solu40] )/anglestep1; DL1 = DL1(:);
%%         DL121 = [ Solu10'*(L121-L1210) , Solu20'*(L122-L1220) , Solu30'*(L123-L1230) , Solu40'*(L124-L1240) ]/anglestep1; DL121 = DL121(:);
%      elseif pb == 2
%         D2    = ([ Lhs*Solu1 ] - Ax)/anglestep2;
%         Dx2   = ([ Solu1 ] - xt)/anglestep2;
%         DL2   = ([ sL*Solu1 ] - sLx)/anglestep2;
%         DL122 = ([ Solu1'*L121 ] - L1x)/anglestep2;
%         Lh2   = (Lhs-Lh0)/anglestep2;

%%         D2    = ((Lhs-Lh0) * [Solu10,Solu20,Solu30,Solu40] )/anglestep2; D2 = D2(:);
%%         DL2   = ((sL -sL0) * [Solu10,Solu20,Solu30,Solu40] )/anglestep2; DL2 = DL2(:);
%%         DL122 = [ Solu10'*(L121-L1210) , Solu20'*(L122-L1220) , Solu30'*(L123-L1230) , Solu40'*(L124-L1240) ]/anglestep2; DL122 = DL122(:);
      end
   end
%   D = [D1,D2]; DL = [DL1,DL2]; DL12 = [DL121,DL122];% Dx = [Dx1,Dx2];
%%   Dd = [ Lh1*Solu10 , Lh2*Solu10 ; Lh1*Solu20 , Lh2*Solu20 ;...
%%          Lh1*Solu30 , Lh2*Solu30 ; Lh1*Solu40 , Lh2*Solu40 ];

%%   ZL = zeros(size(L));
%%   Aile = [ L, ZL, ZL, ZL ; ...
%%            ZL, L, ZL, ZL ; ...
%%            ZL, ZL, L, ZL ; ...
%%            ZL, ZL, ZL, L ];

%   %dtheta = - ( D'*D + mur*Dx'*Aile*Dx ) \ ( D'*res - mur*Dx'*[L121;L122;L123;L124] );
%   dtheta = - ( D'*D + mur*DL'*DL ) \ ( D'*res + mur*DL'*rel + mur*DL12'*ones(size(DL12,1),1) );
%   %dtheta = - ( D'*D ) \ ( D'*res );
%   %dtheta = - ( Dd'*Dd ) \ ( Dd'*res );

%   theta1 = theta1 + dtheta(1); theta1 = mod(theta1,2*pi);
%   theta2 = theta2 + dtheta(2); theta2 = mod(theta2,2*pi);

%   theta1rec(iter+1) = theta1;
%   theta2rec(iter+1) = theta2;
end
disp(['Iterative method terminated ', num2str(toc) ]);

Solu1 = Solu10;

% Post-process : reconstruct u and f from Solu
fsolu1 = zeros(2*nboun2,1); usolu1 = zeros(2*nnodes,1);
szB = 0;
szD = size(LhsA,2);
for i=1:szD
   usolu1(dof2nodes(i)) = Solu1(i);
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
   fsolu1(2*i-2+fmdof) = kB*Solu1(indicesB);
   szB = szB + size(fmdof,1);
end

ucrsol1 = Solu1(end-2*nnodes3+1:end);

% Compute the values of f at the dofs
fsolu1no = zeros(2*nnodes);
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
         fsolu1no(b2nodesnoNi) = fsolu1(2*boound(1))*len1;
      else
         fsolu1no(b2nodesnoNi) = fsolu1(2*boound(1)-1)*len1;
      end
   else
      n1 = boundary(boound(2),2); n2 = boundary(boound(2),3);
      len2 = sqrt( (nodes(n1,1)-nodes(n2,1))^2 + (nodes(n1,2)-nodes(n2,2))^2 );
      if rem(b2nodesnoNi,2)==0
         fsolu1no(b2nodesnoNi) = .5*( fsolu1(2*boound(1))*len1 + fsolu1(2*boound(2))*len2 );
      else
         fsolu1no(b2nodesnoNi) = .5*( fsolu1(2*boound(1)-1)*len1 + fsolu1(2*boound(2)-1)*len2 );
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

theta1pi = theta1rec/pi;
theta2pi = theta2rec/pi;

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

% Graph for f
if min(size(b2nodesnoN))>0
   toplot = fsolu1no(b2nodesnoN); toplot2 = fr1(b2nodesnoN);
   try
   figure;
   hold on;
   plot(toplot(1:2:end-1),'Color','red');
   plot(toplot2(1:2:end-1),'Color','blue');
   legend('Fy identified (1)', 'Fy reference (1)');
   end
end
% Graph for u
if min(size(b2nodesnoD))>0
   toplot = usolu1(b2nodesnoD); toplot2 = ur1(b2nodesnoD);
   try
   figure;
   hold on;
   plot(toplot(1:2:end-1),'Color','red');
   plot(toplot2(1:2:end-1),'Color','blue');
   legend('Uy identified (1)', 'Uy reference (1)');
   end
end

% Graph for [[u]]
curv = sqrt( (nodes3b(:,1)-nodes3b(1,1)).^2 + (nodes3b(:,2)-nodes3b(1,2)).^2 );

n = extnorm3(1,:)';
toplot1  = ucrsol1(1:2:end-1)*n(1) + ucrsol1(2:2:end)*n(2);
toplotr1 = urg(1:2:end-1,1)*n(1) + urg(2:2:end,1)*n(2);

try
figure;
hold on;
set(gca, 'fontsize', 20);
plot( curv, toplot1,'Color','red','LineWidth',3 );
plot( curvr, toplotr1,'Color','blue','LineWidth',3 );
legend('[[u]] identified (1)', '[[u]] reference (1)');
end

%erroru = norm(usolu1(b2nodesnoD)-ur1(b2nodesnoD))   / norm(ur1(b2nodesnoD));
%errorf = norm(fsolu1no(b2nodesnoN)-fr1(b2nodesnoN)) / norm(fr1(b2nodesnoN));
