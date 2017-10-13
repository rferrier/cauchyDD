% 09/05/2017
% Détection de fissure quelconque par écart à la réciprocité

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
jmax       = 35;     % Eigenvalues truncation number
regular    = 0;      % Use the regularization matrix
upper_term = 0;      % 1 : use i=0:10, j=0:10, 0 : use i>=0,j>=0,i+j<=10

thetareg   = pi/3;  % Regularization angle
dreg       = .2;    % Regularization distance
alreg      = 1e-5;     % Regularization parameter

% Boundary conditions
dirichlet  = [0,1,0 ; 0,2,0 ; 0,3,0];
neumann1   = [3,2,fscalar ; 1,2,-fscalar];
neumann2   = [2,1,fscalar ; 4,1,-fscalar];
neumann3   = [3,1,fscalar ; 3,2,fscalar ; 2,1,fscalar ; 2,2,fscalar ; ...
              1,1,-fscalar ; 1,2,-fscalar ; 4,1,-fscalar ; 4,2,-fscalar];
neumann4   = [3,1,-fscalar ; 3,2,fscalar ; 2,1,fscalar ; 2,2,-fscalar ; ...
              1,1,fscalar ; 1,2,-fscalar ; 4,1,-fscalar ; 4,2,fscalar];

% First, import the mesh
[ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc.msh' );

nnodes = size(nodes,1);

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

% Compute stress :
sigma = stress(u1,E,nu,nodes,elements,order,1,ntoelem);

% Output :
plotGMSH({u1,'U1';u2,'U2';u3,'U3';u4,'U4'}, elements, nodes, 'output/reference');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import the uncracked domain
[ nodes2,elements2,ntoelem2,boundary2,order2] = readmesh( 'meshes/rg_refined/plate_nn.msh' );

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

disp([ 'Direct problem solved and data management ', num2str(toc) ]);
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Identify the crack
% Build the polynomial test functions.
load('conditions20_2d.mat','-ascii');
M = spconvert(conditions20_2d); clear('conditions20_2d');
nmax = 20;
ncoef =2*(nmax+1)^2; neq = ncoef;

% Suppress the superior terms
if upper_term == 0 % Rem : there will be loads of 0 in coef
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
[ nodesu,elementsu,ntoelemu,boundaryu,orderu] = readmesh( 'meshes/rg_refined/plate_nu.msh' );
nnodesu = size(nodesu,1);

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
   lis1( find(lis1>j-1) ) = lis1( find(lis1>j-1) ) - j+1;
   lis2 = find( segel(1:j-1,:) == segel(j,2) );
   lis2( find(lis2>j-1) ) = lis2( find(lis2>j-1) ) - j+1;
   
   % also remove boundary elements
   lis3 = find( boundaryu(:,2:3) == segel(j,1) );
   lis3( find(lis3>nbu) ) = lis3( find(lis3>nbu) ) - nbu;
   lis4 = find( boundaryu(:,2:3) == segel(j,2) );
   lis4( find(lis4>nbu) ) = lis4( find(lis4>nbu) ) - nbu;
   
   if size( intersect(lis1,lis2), 1 ) > 0 || size( intersect(lis3,lis4), 1 ) > 0% || (size(lis3,1) > 0 || size(lis4,1) > 0)%
      segel(j,:) = [];
      j = j-1;
      nseg = nseg-1;
   end
   j = j+1;
end

Lhs = zeros(nftest,2*nseg); % Decomposition functions are constant per segment
nor = zeros(nseg,2); % Normal
for j=1:nseg % Compute the 1D integrals
   bonod = segel(j,:);

   no1 = bonod(1); no2 = bonod(2);
   x1  = nodesu(no1,1); y1 = nodesu(no1,2);
   x2  = nodesu(no2,1); y2 = nodesu(no2,2);
   nor(j,:) = [y1-y2 ; x2-x1]'; nor(j,:) = nor(j,:)/norm(nor(j,:));% Normal
   % Orientate the normals
   if nor(j,1) == 0
      nor(j,:) = sign(nor(j,2))*nor(j,:);
   else
      nor(j,:) = sign(nor(j,1))*nor(j,:);
   end
%   if nor(j,2) == 0
%      nor(j,:) = sign(nor(j,1))*nor(j,:);
%   else
%      nor(j,:) = sign(nor(j,2))*nor(j,:);
%   end

   len = sqrt( (x1-x2)^2 + (y1-y2)^2 );
   
   Ng = 1;
   [ Xg, Wg ] = gaussPt1d( Ng );
             
   for k=1:size(Wg,1)
      xg = Xg(k,:); wg = Wg(k);
      
      % Interpolation
      xgr  = ( (1-xg)*[x1;y1] + xg*[x2;y2] ); ... % abscissae

      X = xgr(1)/Lx; Y = xgr(2)/Lx;
      
      for i=1:nftest
         if norm(coef(1:(nmax+1)^2,i)) == 0 % Zero column : no need to continue
            continue;
         end
         
         coefa = coef(1:(nmax+1)^2,i);
         coefb = coef((nmax+1)^2+1:2*(nmax+1)^2,i);
         
         jj = 0;
         for ii=1:nmax
            aij = coefa( (nmax+1)*ii + jj+1 );
            bij = coefb( (nmax+1)*ii + jj+1 );
            % Compute sigma.n
            s11 = E/(1-nu^2)*aij*ii*X^(ii-1)*Y^jj;
            s22 = nu*E/(1-nu^2)*aij*ii*X^(ii-1)*Y^jj;
            s12 = E/(2*(1+nu))*bij*ii*X^(ii-1)*X^jj;
            
            st = [s11,s12 ; s12,s22];
            ft = st*nor(j,:)';
            Lhs(i,2*j-1) = Lhs(i,2*j-1) + len * wg * ft(1); % increment the integral
            Lhs(i,2*j)   = Lhs(i,2*j)   + len * wg * ft(2);
         end
         
         ii = 0;
         for jj=1:nmax
            aij = coefa( (nmax+1)*ii + jj+1 );
            bij = coefb( (nmax+1)*ii + jj+1 );
            % Compute sigma.n
            s11 = nu*E/(1-nu^2)*bij*jj*X^ii*Y^(jj-1);
            s22 = E/(1-nu^2)*bij*jj*X^ii*Y^(jj-1);
            s12 = E/(2*(1+nu))*aij*jj*X^ii*Y^(jj-1);
            
            st = [s11,s12 ; s12,s22];
            ft = st*nor(j,:)';
            Lhs(i,2*j-1) = Lhs(i,2*j-1) + len * wg * ft(1); % increment the integral
            Lhs(i,2*j)   = Lhs(i,2*j)   + len * wg * ft(2);
         end
         
         for ii=1:nmax
            for jj=1:nmax
               aij = coefa( (nmax+1)*ii + jj+1 );
               bij = coefb( (nmax+1)*ii + jj+1 );
               % Compute sigma.n
               s11 = E/(1-nu^2)*aij*ii*X^(ii-1)*Y^jj + nu*E/(1-nu^2)*bij*jj*X^ii*Y^(jj-1);
               s22 = nu*E/(1-nu^2)*aij*ii*X^(ii-1)*Y^jj + E/(1-nu^2)*bij*jj*X^ii*Y^(jj-1);
               s12 = E/(2*(1+nu)) * (aij*jj*X^ii*Y^(jj-1) + bij*ii*X^(ii-1)*X^jj);
 
               st = [s11,s12 ; s12,s22];
               ft = st*nor(j,:)';
               Lhs(i,2*j-1) = Lhs(i,2*j-1) + len * wg * ft(1); % increment the integral
               Lhs(i,2*j)   = Lhs(i,2*j)   + len * wg * ft(2);
            end
         end
      end

   end
end
disp([ 'Left hand side generated ', num2str(toc) ]);
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the RG
Rhs1  = zeros(nftest,1);
Rhs2  = zeros(nftest,1);
Rhs3  = zeros(nftest,1);
Rhs4  = zeros(nftest,1);

for k=1:nftest

   Rhs1(k) = 0; Rhs2(k) = 0;
   coefa = coef(1:(nmax+1)^2,k);
   coefb = coef((nmax+1)^2+1:2*(nmax+1)^2,k);
   
   for i=1:nboun2
      bonod = boundary2(i,:); exno = extnorm2(i,:)';
   
      no1 = bonod(2); no2 = bonod(3);
      x1 = nodes2(no1,1); y1 = nodes2(no1,2);
      x2 = nodes2(no2,1); y2 = nodes2(no2,2);
      len = sqrt( (x1-x2)^2 + (y1-y2)^2 );
      
      if order==1
         Ng = 1; % Max (anyway, it's inexact)
      elseif order==2
         Ng  = 1; no3 = bonod(4);
         x3  = nodes2(no3,1); y3 = nodes2(no3,2);
      end
      [ Xg, Wg ] = gaussPt1d( Ng );
                
      for j=1:Ng
         xg = Xg(j); wg = Wg(j);
         
         % Interpolations
         if order==1
            uer1 = transpose( (1-xg)*urr1(i,1:2) + xg*urr1(i,3:4) ); % [ux;uy] on the
            uer2 = transpose( (1-xg)*urr2(i,1:2) + xg*urr2(i,3:4) ); % Gauss point
            uer3 = transpose( (1-xg)*urr3(i,1:2) + xg*urr3(i,3:4) ); % [ux;uy] on the
            uer4 = transpose( (1-xg)*urr4(i,1:2) + xg*urr4(i,3:4) ); % Gauss point
            xgr  = ( (1-xg)*[x1;y1] + xg*[x2;y2] ); % abscissae of the Gauss point in the physical space
         elseif order==2
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

         Rhs1(k) = Rhs1(k) + len * wg * ( fer1'*vpa - fpa'*uer1 );
         Rhs2(k) = Rhs2(k) + len * wg * ( fer2'*vpa - fpa'*uer2 );
         Rhs3(k) = Rhs3(k) + len * wg * ( fer3'*vpa - fpa'*uer3 );
         Rhs4(k) = Rhs4(k) + len * wg * ( fer4'*vpa - fpa'*uer4 );
      end
   end
   
end
disp([ 'Right hand side generated ', num2str(toc) ]);
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve the linear system and recover the volumic loading
Rhs = Rhs1;
A = Lhs'*Lhs; b = Lhs'*Rhs;

%[Q,Theta,P] = svd(A);
%thetas = diag(Theta);
%% figure; plot(log10(diag(Theta)));
%[thetas,Ind] = sort( thetas,'descend' );
%Q = Q(:,Ind); P = P(:,Ind);
%Thetas = diag(thetas); 
%Theta = Theta(Ind,Ind);
%
%disp([ 'Rectangular system pinversed ', num2str(toc) ]);
%% Plot the Picard stuff
%imax = min( find(thetas/thetas(1)<1e-16) );
%if size(imax,1) == 0
%    imax = size(thetas,1);
%end
%
%tplo = thetas(1:imax); bplo = Q'*b; bplo = bplo(1:imax);
%rplo = (Q'*b)./thetas; rplo = rplo(1:imax);
%figure
%hold on;
%plot(log10(abs(tplo)),'Color','blue');
%plot(log10(abs(bplo)),'Color','red');
%plot(log10(abs( rplo )),'Color','black');
%
%% Filter eigenvalues
%ThetaT = Thetas( 1:jmax , 1:jmax );
%bT     = Q'*b; bT = bT(1:jmax);
%
%Solu = P(:,1:jmax) * (ThetaT\bT);

%[Q,Theta,P] = svd(Lhs);
%thetas = diag(Theta);
%% figure; plot(log10(diag(Theta)));
%[thetas,Ind] = sort( thetas,'descend' );
%Q = Q(:,Ind); P = P(:,Ind);
%Thetas = diag(thetas); 
%Theta = Theta(Ind,Ind);

%% Build the regularization matrix
if regular == 1
   L = zeros(size(A));
   for i = 1:nnodesu
      % Find the segments
      seg = find(segel==i);
      seg(find(seg>nseg)) = seg(find(seg>nseg)) - nseg;
      
      sseg = size(seg,1);
      L(2*seg-1,2*seg-1) = L(2*seg-1,2*seg-1) -1;
      L(2*seg,2*seg)     = L(2*seg,2*seg) -1;
      for j=1:sseg
         L(2*seg(j)-1,2*seg(j)-1) = L(2*seg(j)-1,2*seg(j)-1) + sseg; % Le = sseg-1 on the diagonal
         L(2*seg(j),2*seg(j))     = L(2*seg(j),2*seg(j)) + sseg;
      end
   end
   ninfty = size(L,1)-rank(L); % To remove infinitys
elseif regular == 2
   L = zeros(size(A));
   for i = 1:nnodesu
      seg = find(segel==i);
      seg(find(seg>nseg)) = seg(find(seg>nseg)) - nseg;
      sseg = size(seg,1);
      L(2*seg-1,2*seg-1) = L(2*seg-1,2*seg-1) -1;
      L(2*seg,2*seg)     = L(2*seg,2*seg) -1;
      for j=1:sseg
         L(2*seg(j)-1,2*seg(j)-1) = L(2*seg(j)-1,2*seg(j)-1) + sseg;
         L(2*seg(j),2*seg(j))     = L(2*seg(j),2*seg(j)) + sseg;
      end
   end
   L = pinv(L);
   ninfty = size(L,1)-rank(L); % I love theoreme of rank
elseif regular == 3
   L = zeros(size(A));
   for i = 1:nseg
      for j = i+1:nseg
         n0 = intersect( segel(i,:), segel(j,:) );
         if size( n0 ) > 0 % They have a point in common
            n1 = setdiff( segel(i,:), n0 ); n2 = setdiff( segel(j,:), n0 );
            x1 = nodesu(n1,1)-nodesu(n0,1); y1 = nodesu(n1,2)-nodesu(n0,2);
            x2 = nodesu(n2,1)-nodesu(n0,1); y2 = nodesu(n2,2)-nodesu(n0,2);
            alpha = atan2(y2,x2) - atan2(y1,x1);
            
            L(2*i-1,2*i-1) = L(2*i-1,2*i-1) + 1;
            L(2*i,2*i)     = L(2*i,2*i) + 1;
            L(2*j-1,2*j-1) = L(2*j-1,2*j-1) + 1;
            L(2*j,2*j)     = L(2*j,2*j) + 1;
            
            L(2*i-1,2*j-1) = L(2*i-1,2*j-1) + cos(alpha);
            L(2*i,2*j)     = L(2*i,2*j) + cos(alpha);
            L(2*j-1,2*i-1) = L(2*j-1,2*i-1) + cos(alpha);
            L(2*j,2*i)     = L(2*j,2*i) + cos(alpha);
         end
      end
   end
   ninfty = size(L,1)-rank(L); % To remove infinitys
elseif regular == 4
   L = zeros(size(A));
   for i = 1:nseg
      for j = i+1:nseg
%         n0 = intersect( segel(i,:), segel(j,:) );
%         if size( n0 ) > 0 % They have a point in common
%            n1 = setdiff( segel(i,:), n0 ); n2 = setdiff( segel(j,:), n0 );
%            x1 = nodesu(n1,1)-nodesu(n0,1); y1 = nodesu(n1,2)-nodesu(n0,2);
%            x2 = nodesu(n2,1)-nodesu(n0,1); y2 = nodesu(n2,2)-nodesu(n0,2);
%            alpha = atan2(y2,x2) - atan2(y1,x1);
%            
%            L(2*i-1,2*i-1) = L(2*i-1,2*i-1) + 1;
%            L(2*i,2*i)     = L(2*i,2*i) + 1;
%            L(2*j-1,2*j-1) = L(2*j-1,2*j-1) + 1;
%            L(2*j,2*j)     = L(2*j,2*j) + 1;
%            valuebe = (pi-alpha)^2/thetareg^2;
%            L(2*i-1,2*j-1) = L(2*i-1,2*j-1) + valuebe;
%            L(2*i,2*j)     = L(2*i,2*j) + valuebe;
%            L(2*j-1,2*i-1) = L(2*j-1,2*i-1) + valuebe;
%            L(2*j,2*i)     = L(2*j,2*i) + valuebe;
%         else % Look for quasi-parallel elements

            nmali = nor(i,:)'; % Local normal
            nmalj = nor(j,:)'; % Local normal

            n1 = segel(i,1); n2 = segel(i,2); n3 = segel(j,1); n4 = segel(j,2);
            x1 = nodesu(n1,1); x2 = nodesu(n2,1); x3 = nodesu(n3,1); x4 = nodesu(n4,1);
            y1 = nodesu(n1,2); y2 = nodesu(n2,2); y3 = nodesu(n3,2); y4 = nodesu(n4,2);
            t1 = [x2-x1;y2-y1]; t2 = [x4-x3;y4-y3]; t1 = t1/norm(t1); t2 = t2/norm(t2);
            dist = sqrt( (x1+x2-x3-x4)^2 + (y1+y2-y3-y4)^2 ) / 2;
            l1   = sqrt( (x1-x2)^2 + (y1-y2)^2 ) / 2;
            l2   = sqrt( (x3-x4)^2 + (y3-y4)^2 ) / 2;
            
            % Projection on s1
            a1 = [x1,y1]*t1; a2 = [x2,y2]*t1; a3 = [x3,y3]*t1; a4 = [x4,y4]*t1;
            % Test if the projection of the center of s2 is in s1
            l2s1 = 0;
            if ((a3+a4)/2-a1) * ((a3+a4)/2-a2) < 0
               l2s1 = (a4-a3);
            end
            
            % Projection on s2
            a1 = [x1,y1]*t2; a2 = [x2,y2]*t2; a3 = [x3,y3]*t2; a4 = [x4,y4]*t2;
            % Test if the projection of the center of s1 is in s2
            l1s2 = 0;
            if ((a1+a2)/2-a3) * ((a1+a2)/2-a4) < 0
               l1s2 = (a1-a2);
            end

            valuebe = dreg^2/dist^2*l1s2^2*l2s1^2/l1^2/l2^2*sign(nmali'*nmalj); % Question of sign of the normal
            
            if dreg^2/dist^2 > 1
               L(2*i-1,2*j-1) = L(2*i-1,2*j-1) + valuebe;
               L(2*i,2*j)     = L(2*i,2*j) + valuebe;
               L(2*j-1,2*i-1) = L(2*j-1,2*i-1) + valuebe;
               L(2*j,2*i)     = L(2*j,2*i) + valuebe;
            end
      end
      
%      L(2*i-1,2*i-1) = L(2*i-1,2*i-1) + sign( nmali'*[1;0] ); % Again, the sign
%      L(2*i,2*i)     = L(2*i,2*i) + sign( nmali'*[0;1] );
      
   end
%   L = L + norm(L)*eye(size(L));

   % Every raw of L must hav the same norm (in order not to penalize particular values)
   nL = 0;
   for i=1:size(L,2)
      cand = norm(L(:,i));
      if nL < cand
         nL = cand;
      end
   end
   
   for i=1:size(L,2)
      cand = sqrt( nL^2 - norm(L(:,i))^2 );
      L(i,i) = L(i,i) + cand;
   end

   ninfty = size(L,1)-rank(L); % To remove infinitys
elseif regular == 5
   L = zeros(size(A));
   for i = 1:nseg
      nmali = nor(i,:)'; % Local normal
      L(2*i-1,2*i-1) = L(2*i-1,2*i-1) + sign( nmali'*[1;0] ); % Again, the sign
      L(2*i,2*i)     = L(2*i,2*i) + sign( nmali'*[0;1] );
   end
   ninfty = 0;
else
   L = eye(size(A));
   ninfty = 0;
end

%A = A+alreg*norm(A)/norm(L)*L;
%Solu1 = A\(Lhs'*Rhs1);
%Solu2 = A\(Lhs'*Rhs2);

[Q,Theta] = eig(A,L); Q = Q*(Q'*L*Q)^(-1/2);
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

%tplo = thetas(1:imax); bplo1 = Q'*Rhs1; bplo1 = bplo1(1:imax);
%rplo1 = (Q'*Rhs1)./thetas; rplo1 = rplo1(1:imax);
%bplo2 = Q'*Rhs2; bplo2 = bplo2(1:imax);
%rplo2 = (Q'*Rhs2)./thetas; rplo2 = rplo2(1:imax);
%figure
%hold on;
%plot(log10(abs(tplo)),'Color','green');
%plot(log10(abs(bplo1)),'Color','red');
%plot(log10(abs(rplo1)),'Color','black');
%plot(log10(abs(bplo2)),'Color','magenta');
%plot(log10(abs(rplo2)),'Color','blue');
%legend('Singular values','Rhs1','sol1','Rhs2','sol2');
%
%% Filter eigenvalues
%if jmax == 0
%   jmax = size(Thetas,1);
%end
%ThetaT = Thetas( 1:jmax , 1:jmax );
%bT     = Q'*Rhs; bT = bT(1:jmax);
%
%bT1 = Q'*Rhs1; bT1 = bT1(1:jmax);
%bT2 = Q'*Rhs2; bT2 = bT2(1:jmax);
%
%Solu = P(:,1:jmax) * (ThetaT\bT);
%Solu1 = P(:,1:jmax) * (ThetaT\bT1);
%Solu2 = P(:,1:jmax) * (ThetaT\bT2);

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

Solu1 = Q(:,ninfty+1:jmax) * (ThetaT\bT1);
Solu2 = Q(:,ninfty+1:jmax) * (ThetaT\bT2);

% Re-loop over segments to build normal gap
nogap  = zeros(nseg,1); nogap1 = zeros(nseg,1); nogap2 = zeros(nseg,1);
ntoseg = zeros(nnodesu,1);
Sgap   = zeros(nnodesu,1);
for j=1:nseg
   no1 = segel(j,1); no2 = segel(j,2);
   nogap1(j) = norm( Solu1( [2*j-1,2*j] ) ); 
   nogap2(j) = norm( Solu2( [2*j-1,2*j] ) );
   nogap(j)  = sqrt(nogap1(j)^2*nogap2(j)^2);
%   ntoseg(no1) = ntoseg(no1) + 1; ntoseg(no2) = ntoseg(no2) + 1;
%   Sgap(no1) = Sgap(no1) + nogap(j); Sgap(no2) = Sgap(no2) + nogap(j);
end
Sno = Sgap./ntoseg;
%Sno = .5*(Sno-abs(Sno)); % Negative part
%Sno = Red*Solu;

%% Output
%plotGMSH({Sno,'Normal gap'}, elementsu, nodesu, 'output/identification');

% Segments visu
figure;
hold on;
plot( [nodes(5,1),nodes(6,1)], [nodes(5,2),nodes(6,2)], 'Color', [.6,.6,.6], 'LineWidth', 5 );
nogapp1 = abs(nogap1)-min(abs(nogap1)); maxn1 = max(nogapp1);
for i=1:nseg
   no1 = segel(i,1); no2 = segel(i,2);
   x1 = nodesu(no1,1); y1 = nodesu(no1,2);
   x2 = nodesu(no2,1); y2 = nodesu(no2,2);
   
   x = nogapp1(i)/maxn1;
   rgb = rgbmap(x);
   plot( [x1,x2], [y1,y2], 'Color', rgb, 'LineWidth', 3 );
end
%axis equal;
colormap("default")
h = colorbar();
ytick = get (h, "ytick");
set (h, "yticklabel", sprintf ( "%g|", maxn1*ytick+min(nogap1) ));

figure;
hold on;
plot( [nodes(5,1),nodes(6,1)], [nodes(5,2),nodes(6,2)], 'Color', [.6,.6,.6], 'LineWidth', 5 );
nogapp2 = abs(nogap2)-min(abs(nogap2)); maxn2 = max(nogapp2);
for i=1:nseg
   no1 = segel(i,1); no2 = segel(i,2);
   x1 = nodesu(no1,1); y1 = nodesu(no1,2);
   x2 = nodesu(no2,1); y2 = nodesu(no2,2);
   
   x = nogapp2(i)/maxn2;
   rgb = rgbmap(x);
   plot( [x1,x2], [y1,y2], 'Color', rgb, 'LineWidth', 3 );
end
%axis equal;
colormap("default")
h = colorbar();
ytick = get (h, "ytick");
set (h, "yticklabel", sprintf ( "%g|", maxn2*ytick+min(nogap2) ));

figure;
hold on;
plot( [nodes(5,1),nodes(6,1)], [nodes(5,2),nodes(6,2)], 'Color', [.6,.6,.6], 'LineWidth', 5 );
nogapp = abs(nogap)-min(abs(nogap)); maxn = max(nogapp);
for i=1:nseg
   no1 = segel(i,1); no2 = segel(i,2);
   x1 = nodesu(no1,1); y1 = nodesu(no1,2);
   x2 = nodesu(no2,1); y2 = nodesu(no2,2);
   
   x = nogapp(i)/maxn;
   rgb = rgbmap(x);
   plot( [x1,x2], [y1,y2], 'Color', rgb, 'LineWidth', 3 );
end
%axis equal;
colormap("default")
h = colorbar();
ytick = get (h, "ytick");
set (h, "yticklabel", sprintf ( "%g|", maxn2*ytick+min(nogap2) ));

%% Per node visu
%figure;
%patch('Faces',elementsu(:,1:3),'Vertices',nodesu,'FaceVertexCData',Sno,'FaceColor','interp');
%colorbar; axis equal;

%% Reference visu (force on the displacement basis, it's ugly)
%fxy = reshape(f,2,[])';  fx = fxy(:,1);  fy = fxy(:,2);
%figure;
%patch('Faces',elements(:,1:3),'Vertices',nodes,'FaceVertexCData',fx,'FaceColor','interp');
%colorbar;
%figure;
%patch('Faces',elements(:,1:3),'Vertices',nodes,'FaceVertexCData',fy,'FaceColor','interp');
%colorbar;

% On line visu
%Y = Ymin:Ly/100:Ymax; X = Xmoy*ones(1,size(Y,2)); XY = [X;Y];
%Sline = passMesh2D(nodesu, elementsu, XY', [], Sno);
%Slinexy = reshape(Sline,2,[])';  Slinex = Slinexy(:,1);  Sliney = Slinexy(:,2);
%
%ref = loadV(1,1) + loadV(2,1)*Y + loadV(3,1)*Y.^2 + loadV(4,1)*Y.^3;
%figure;
%hold on;
%plot(Y,Sliney);
%plot(Y,ref,'Color','red')