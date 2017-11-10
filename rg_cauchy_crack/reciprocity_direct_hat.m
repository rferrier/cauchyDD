% 02/11/2017
% Problèmes directs et de Cauchy par écart à la réciprocité, 
% identification par fonctions chapeau

tic
close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E          = 210000; % MPa : Young modulus
nu         = 0.3;    % Poisson ratio
fscalar    = 250;    % N.mm-1 : Loading on the plate
mat        = [0, E, nu];
br         = .0;      % Noise level
jmax       = 40;     % Eigenvalues truncation number
regular    = 0;      % Use the regularization matrix
upper_term = 0;      % 1 : use i=0:10, j=0:10, 0 : use i>=0,j>=0,i+j<=10
froreg     = 1;      % frobenius preconditioner
recompute  = 1;      % Recompute the operators

% Boundary conditions
dirichlet  = [0,1,0 ; 0,2,0 ; 0,3,0];
neumann1   = [3,2,fscalar ; 1,2,-fscalar];
neumann2   = [2,1,fscalar ; 4,1,-fscalar];
neumann3   = [3,1,fscalar ; 3,2,fscalar ; 2,1,fscalar ; 2,2,fscalar ; ...
              1,1,-fscalar ; 1,2,-fscalar ; 4,1,-fscalar ; 4,2,-fscalar];
neumann4   = [3,1,-fscalar ; 3,2,fscalar ; 2,1,fscalar ; 2,2,-fscalar ; ...
              1,1,fscalar ; 1,2,-fscalar ; 4,1,-fscalar ; 4,2,fscalar];
              
%dirichlet0 = [1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ; 3,1,0 ; 3,2,0 ; 4,1,0 ; 4,2,0];
%neumann0 = []
dirichlet0 = [1,1,0 ; 1,2,0 ; 3,1,0 ; 3,2,0];
neumann0   = [2,1,0 ; 2,2,0 ; 4,1,0 ; 4,2,0];
%dirichlet0 = [1,1,0 ; 1,2,0 ; 3,1,0 ; 3,2,0 ; 4,1,0 ; 4,2,0];
%neumann0   = [4,1,0 ; 4,2,0];

% First, import the mesh
[ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nu.msh' );

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ nodes2,elements2,ntoelem2,boundary2,order] = readmesh( 'meshes/rg_refined/plate_nu.msh' );
nnodes2 = size(nodes2,1);
[K2,C2,nbloq2,node2c2,c2node2] = Krig2 (nodes2,elements2,mat,order,boundary2,dirichlet);
Kinter2 = K2( 1:2*nnodes2, 1:2*nnodes2 );

% Mapbounds
[ node2b12, b2node12 ] = mapBound( 1, boundary2, nnodes2 );
[ node2b22, b2node22 ] = mapBound( 2, boundary2, nnodes2 );
[ node2b32, b2node32 ] = mapBound( 3, boundary2, nnodes2 );
[ node2b42, b2node42 ] = mapBound( 4, boundary2, nnodes2 );

indexbound  = [2*b2node1-1 ; 2*b2node1 ; 2*b2node2-1 ; 2*b2node2 ;...
               2*b2node3-1 ; 2*b2node3 ; 2*b2node4-1 ; 2*b2node4];

 % Sort the nodes (in case order > 1)
[~,order5]  = sort( nodes( b2node5, 1 ) );
[~,order6]  = sort( nodes( b2node6, 1 ) );

icrack5x    = [ 2*b2node5(order5)-1 ];
icrack6x    = [ 2*b2node6(order6)-1 ];
icrack5y    = [ 2*b2node5(order5) ];
icrack6y    = [ 2*b2node6(order6) ];

indexbound2 = [2*b2node12-1 ; 2*b2node12 ; 2*b2node22-1 ; 2*b2node22 ;...
               2*b2node32-1 ; 2*b2node32 ; 2*b2node42-1 ; 2*b2node42];
% Pass f and u on the uncracked mesh
ur1 = zeros( 2*nnodes2, 1 ); fr1 = zeros( 2*nnodes2, 1 );
ur1(indexbound2) = u1(indexbound);
fr1(indexbound2) = f1(indexbound);
% Same for the second one
ur2 = zeros( 2*nnodes2, 1 ); fr2 = zeros( 2*nnodes2, 1 );
ur2(indexbound2) = u2(indexbound);
fr2(indexbound2) = f2(indexbound);
% 3
ur3 = zeros( 2*nnodes2, 1 ); fr3 = zeros( 2*nnodes2, 1 );
ur3(indexbound2) = u3(indexbound);
fr3(indexbound2) = f3(indexbound);
% 4
ur4 = zeros( 2*nnodes2, 1 ); fr4 = zeros( 2*nnodes2, 1 );
ur4(indexbound2) = u4(indexbound);
fr4(indexbound2) = f4(indexbound);

% Add the noise
u1n = u1; u2n = u2; u3n = u3; u4n = u4;
br1 = randn(2*nnodes,1); br2 = randn(2*nnodes,1);
br3 = randn(2*nnodes,1); br4 = randn(2*nnodes,1);
% noise = load('noises/105.mat'); br1 = noise.br1; br2 = noise.br2;
u1 = ( 1 + br*br1 ) .* u1; u2 = ( 1 + br*br2 ) .* u2;
u3 = ( 1 + br*br3 ) .* u3; u4 = ( 1 + br*br4 ) .* u4;

plotGMSH({ur1,'U_bound'}, elements2, nodes2, 'bound');

%% Build the maps of the Dirichlet and Neumann bounds
%b2nodesD = cell(4,2); % nb of boundaries * nb of dimensions
%b2nodesN = cell(4,2);
b2nodesD = []; b2nodesN = []; b2nodesTOT = [];
for i=1:4
   [~, b2node] = mapBound(i, boundary, nnodes);
   dofD = dirichlet0(find(dirichlet0(:,1)==i),2); % dof of the given dirichlet
   dofN = neumann0(find(neumann0(:,1)==i),2); % dof of the given neumann
   for j=1:size(dofD,1)
      b2nodesD = [b2nodesD ; 2*b2node-2+dofD(j)];
   end
   b2nodesTOT = [b2nodesTOT ; 2*b2node-1 ; 2*b2node];
%   if max(size(dofN)) > 0
   for j=1:size(dofN,1)
      b2nodesN = [b2nodesN ; 2*b2node-2+dofN(j)];
   end
end
b2nodesnoD = setdiff(b2nodesTOT,b2nodesD);  % Missing Dirichlet conditions
b2nodesnoN = setdiff(b2nodesTOT,b2nodesN);  % Missing Neumann conditions (unused except for visualization)
% Put u to 0 at those missing BCs
u1(b2nodesnoD) = 0; u2(b2nodesnoD) = 0; u3(b2nodesnoD) = 0; u4(b2nodesnoD) = 0;

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

nboun1 = size(boundary,1); nelem1 = size(elements,1);
boun2vol1 = zeros( nboun1, 1 ); extnorm1 = zeros( nboun1, 2 );
sigr1 = zeros( nboun1, 3*order ); sigr2 = zeros( nboun1, 3*order );
urr1  = zeros( nboun1, 2+2*order ); urr2 = zeros( nboun1, 2+2*order );
urr3  = zeros( nboun1, 2+2*order ); urr4 = zeros( nboun1, 2+2*order );
for i=1:nboun1
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

[ nodesu,elementsu,ntoelemu,boundaryu,orderu] = readmesh( 'meshes/rg_refined/plate_nu.msh' );
nnodesu = size(nodesu,1);

disp([ 'Direct problem solved and data management ', num2str(toc) ]);
if recompute == 1
   tic
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% Solve via reciprocity gap
   % Build the polynomial test functions.
   load('conditions20_2d.mat','-ascii');
   M = spconvert(conditions20_2d); clear('conditions20_2d');
   nmax = 20;
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
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% Compute the RG and the Lhs
   Rhs1  = zeros(nftest,1); Rhs2  = zeros(nftest,1);
   Rhs3  = zeros(nftest,1); Rhs4  = zeros(nftest,1);
   LhsA  = zeros(nftest,szD); LhsB  = zeros(nftest,0);
   
   for k=1:nftest
      Rhs1(k) = 0; Rhs2(k) = 0;
      coefa = coef(1:(nmax+1)^2,k);
      coefb = coef((nmax+1)^2+1:2*(nmax+1)^2,k);
      
      szA = 0; szB = 0;
      for i=1:nboun2
         bonod = boundary2(i,:); exno = extnorm2(i,:)';
      
         no1 = bonod(2); no2 = bonod(3); boname = bonod(1);
         x1 = nodes2(no1,1); y1 = nodes2(no1,2);
         x2 = nodes2(no2,1); y2 = nodes2(no2,2);
         len = sqrt( (x1-x2)^2 + (y1-y2)^2 );
         
         dirichletDof = boun2dof(i,:);
         suppressInDD = find(dirichletDof==0); % remove the zeros
         keepInDD = setdiff([1;2;3;4],suppressInDD);
         
         if max(size(dirichlet0))>0
            urdof = dirichlet0(find(dirichlet0(:,1)==boname),2) ;  % Dof of the given u
         else
            urdof = [];
         end
         if max(size(neumann0))>0
            frdof = neumann0(find(neumann0(:,1)==boname),2) ;  % Dof of the given f
         else
            frdof = [];
         end
         umdof = setdiff([1;2],urdof);  % Dof of the unknown u
         fmdof = setdiff([1;2],frdof);  % Dof of the unknown f
         
         indicesA = dirichletDof( keepInDD );
         indicesB = szB+1:szB+size(fmdof,1);
         %LhsA(k,indicesA) = zeros(1,size(umdof,1));
         LhsB(k,indicesB) = zeros(1,size(fmdof,1));
         
         if order==1
            Ng = 2; % (anyway, it's inexact)
         elseif order==2
            Ng  = 2; no3 = bonod(4);
            x3  = nodes2(no3,1); y3 = nodes2(no3,2);
         end
         [ Xg, Wg ] = gaussPt1d( Ng );
                   
         for j=1:Ng
            xg = Xg(j); wg = Wg(j);
            
            % Interpolations
            if order == 1
               uer1 = transpose( (1-xg)*urr1(i,1:2) + xg*urr1(i,3:4) ); % [ux;uy] on the
               uer2 = transpose( (1-xg)*urr2(i,1:2) + xg*urr2(i,3:4) ); % Gauss point
               uer3 = transpose( (1-xg)*urr3(i,1:2) + xg*urr3(i,3:4) ); % [ux;uy] on the
               uer4 = transpose( (1-xg)*urr4(i,1:2) + xg*urr4(i,3:4) ); % Gauss point
               xgr  = ( (1-xg)*[x1;y1] + xg*[x2;y2] ); % abscissae of the Gauss point in the physical space
            elseif order == 2
               uer1 = transpose( urr1(i,1:2) + ...
                      xg*(4*urr1(i,5:6)-3*urr1(i,1:2)-urr1(i,3:4)) + ...   % [ux;uy] on the
                      xg^2*(2*urr1(i,3:4)+2*urr1(i,1:2)-4*urr1(i,5:6)) ); % Gauss point
               uer2 = transpose( urr2(i,1:2) + ...
                      xg*(4*urr2(i,5:6)-3*urr2(i,1:2)-urr2(i,3:4)) + ...  
                      xg^2*(2*urr2(i,3:4)+2*urr2(i,1:2)-4*urr2(i,5:6)) );
               uer3 = transpose( urr3(i,1:2) + ...
                      xg*(4*urr3(i,5:6)-3*urr3(i,1:2)-urr3(i,3:4)) + ...   % [ux;uy] on the
                      xg^2*(2*urr3(i,3:4)+2*urr3(i,1:2)-4*urr3(i,5:6)) ); % Gauss point
               uer4 = transpose( urr4(i,1:2) + ...
                      xg*(4*urr4(i,5:6)-3*urr4(i,1:2)-urr4(i,3:4)) + ...  
                      xg^2*(2*urr4(i,3:4)+2*urr4(i,1:2)-4*urr4(i,5:6)) );
               xgr  = ( [x1;y1] + xg*(4*[x3;y3]-3*[x1;y1]-[x2;y2]) + ...
                      xg^2*(2*[x2;y2]+2*[x1;y1]-4*[x3;y3]) );
            end
      
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
   
            % Remove missing info
            fer1(fmdof) = 0; fer2(fmdof) = 0; fer3(fmdof) = 0; fer4(fmdof) = 0;
            
            Rhs1(k) = Rhs1(k) + len * wg * ( fer1'*vpa - fpa'*uer1 );
            Rhs2(k) = Rhs2(k) + len * wg * ( fer2'*vpa - fpa'*uer2 );
            Rhs3(k) = Rhs3(k) + len * wg * ( fer3'*vpa - fpa'*uer3 );
            Rhs4(k) = Rhs4(k) + len * wg * ( fer4'*vpa - fpa'*uer4 );
            
            fpaTimesPhitest = [fpa(1)*(1-xg);fpa(2)*(1-xg);fpa(1)*xg;fpa(2)*xg]; 
            fpaTimesPhitest = fpaTimesPhitest( keepInDD );
   
            LhsA(k,indicesA) = LhsA(k,indicesA) + len * wg * fpaTimesPhitest';
            LhsB(k,indicesB) = LhsB(k,indicesB) - len * wg * vpa(fmdof)';
         end
         szA = szA + size(umdof,1);  szB = szB + size(fmdof,1);
      end
   end
   disp([ 'Right hand side generated ', num2str(toc) ]);
else
   Anb  = load('rg_cauchy_crack/reciprocity_direct_hat.mat');
   LhsA  = Anb.LhsA; LhsB  = Anb.LhsB;
   Rhs1 = Anb.Rhs1; Rhs2 = Anb.Rhs2; Rhs3 = Anb.Rhs3; Rhs4 = Anb.Rhs4;
end
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve the linear system and recover the unknowns
if froreg == 1
   kB = norm(LhsA,'fro')/norm(LhsB,'fro'); % In order to regularize the stuff
else
   kB = 1;
end
Lhs = [LhsA,kB*LhsB];

Rhs = Rhs1;
A = Lhs'*Lhs; b = Lhs'*Rhs;

L = eye(size(A));
ninfty = 0;

[Q,Theta] = eig(A,L); Q = Q*real((Q'*L*Q)^(-1/2)); % real should not be, but you know, numerical shit...
thetas = diag(Theta);
[thetas,Ind] = sort( thetas,'descend' );
Q = Q(:,Ind);
Thetas = diag(thetas); 
Theta = Theta(Ind,Ind);

disp([ 'Rectangular system pinversed ', num2str(toc) ]);
% Plot the Picard stuff
imax = min( find(thetas/thetas(ninfty+1)<1e-16) );
if size(imax,1) == 0
    imax = size(thetas,1);
end

tplo = thetas(ninfty+1:imax); bplo1 = Q'*b; bplo1 = bplo1(ninfty+1:imax);
rplo1 = (Q'*Lhs'*Rhs1)./thetas; rplo1 = rplo1(1:imax);
bplo2 = Q'*Lhs'*Rhs2; bplo2 = bplo2(ninfty+1:imax);
rplo2 = (Q'*Lhs'*Rhs2)./thetas; rplo2 = rplo2(1:imax);
figure
hold on;
plot(log10(abs(tplo)),'Color','green');
plot(log10(abs(bplo1)),'Color','red');
plot(log10(abs(rplo1)),'Color','black');
plot(log10(abs(bplo2)),'Color','magenta');
plot(log10(abs(rplo2)),'Color','blue');
legend('Singular values','Rhs1','sol1','Rhs2','sol2');

% Filter eigenvalues
if jmax == 0
   jmax = size(Thetas,1);
end
ThetaT = Thetas( ninfty+1:jmax , ninfty+1:jmax );

bT1 = Q'*Lhs'*Rhs1; bT1 = bT1(ninfty+1:jmax);
bT2 = Q'*Lhs'*Rhs2; bT2 = bT2(ninfty+1:jmax);
bT3 = Q'*Lhs'*Rhs1; bT3 = bT3(ninfty+1:jmax);
bT4 = Q'*Lhs'*Rhs2; bT4 = bT4(ninfty+1:jmax);

Solu1 = Q(:,ninfty+1:jmax) * (ThetaT\bT1);
Solu2 = Q(:,ninfty+1:jmax) * (ThetaT\bT2);
Solu3 = Q(:,ninfty+1:jmax) * (ThetaT\bT3);
Solu4 = Q(:,ninfty+1:jmax) * (ThetaT\bT4);

% Post-process : reconstruct u and f from Solu
fsolu1 = zeros(2*nboun2,1); usolu1 = zeros(2*nnodes,1);
fsolu2 = zeros(2*nboun2,1); usolu2 = zeros(2*nnodes,1);
fsolu3 = zeros(2*nboun2,1); usolu3 = zeros(2*nnodes,1);
fsolu4 = zeros(2*nboun2,1); usolu4 = zeros(2*nnodes,1);
szB = 0;
for i=1:szD
   usolu1(dof2nodes(i)) = Solu1(i);
   usolu2(dof2nodes(i)) = Solu2(i);
   usolu3(dof2nodes(i)) = Solu3(i);
   usolu4(dof2nodes(i)) = Solu4(i);
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
   fsolu2(2*i-2+fmdof) = kB*Solu2(indicesB);
   fsolu3(2*i-2+fmdof) = kB*Solu3(indicesB);
   fsolu4(2*i-2+fmdof) = kB*Solu4(indicesB);
   szB = szB + size(fmdof,1);
end

% Segments visu
toplot = fsolu1(2:2:end);
%toplot = usolu1(2:2:end);
figure; hold on;
minn1 = min(toplot); maxn1 = max(toplot);
for i=1:nboun2
   no1 = boundary(i,2); no2 = boundary(i,3);
   x1 = nodesu(no1,1); y1 = nodesu(no1,2);
   x2 = nodesu(no2,1); y2 = nodesu(no2,2);
   
   x = (toplot(i)-minn1)/(maxn1-minn1);
   rgb = rgbmap(x);
   plot( [x1,x2], [y1,y2], 'Color', rgb, 'LineWidth', 3 );
end
%axis equal;
colormap('default')
h = colorbar();
ytick = get (h, 'ytick');
set (h, 'yticklabel', sprintf ( '%g|', (maxn1-minn1)*ytick+minn1 ));

% Compute the values of f at the dofs
fsolu1no = zeros(2*nnodes);
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
   else
      fsolu1no(b2nodesnoNi) = .5*( fsolu1(2*boound(1)-1)*len1 + fsolu1(2*boound(2)-1)*len2 );
   end
end

% Graph for u
toplot = usolu1(b2nodesnoD); toplot2 = u1ref(b2nodesnoD);
figure;
hold on;
plot(toplot(2:2:end),'Color','red');
plot(toplot2(2:2:end),'Color','blue');

% Graph for f
toplot = fsolu1no(b2nodesnoN); toplot2 = f1ref(b2nodesnoN);
figure;
hold on;
plot(toplot(2:2:end),'Color','red');
plot(toplot2(2:2:end),'Color','blue');

erroru = norm(usolu1(b2nodesnoD)-u1ref(b2nodesnoD))   / norm(u1ref(b2nodesnoD));
errorf = norm(fsolu1no(b2nodesnoN)-f1ref(b2nodesnoN)) / norm(f1ref(b2nodesnoN));
