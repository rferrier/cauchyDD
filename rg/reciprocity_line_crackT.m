% 23/05/2017
% Détection de fissure quelconque par écart à la réciprocité
% Régularisation non-linéaire topologique

tic
close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E          = 210000; % MPa : Young modulus
nu         = 0.3;    % Poisson ratio
fscalar    = 250;    % N.mm-1 : Loading on the plate
mat        = [0, E, nu];
br         = 0;      % Noise level
jmax       = 0;     % Eigenvalues truncation number (if 0, automatic Picard choose)
%niter      = 4;    %26-49 Number of regularization iterations
ncrack     = 10001;    % nb of cracks (odd : 1 crack, even : 2 cracks), 5 : 1% noise, 7 : 10% noise, 9 : corner crack, 11 : U crack
                    % 51  : refined postpro mesh
                    % 101 : direct & integrals refined, basic crack
                    % 103 : idem for the small crack
                    % 1001 : more test-functions
                    % 10001 : analysis on a mesh marking the crack
co         = [1,1,1,1]; % Coefficients for each RG test-case
recompute  = 0; % Recompute A and b

% List of the operations :
% 1 : Computation without any suppression
% 2 : 2 elts per node algorithm
% 3 : Building of chains and suppression of monoms
% 4 : Proximity-based suppression of shortest chains
%operations = [1,2,3,4,4,4,4,4];
%regD       = [0,0,0,.1,.2,.3,.4,.5];   
% Minimal distance between 2 chains (0 means only one) 

%operations = [1,2,3,4,4,4,4];
%regD       = [0,0,0,.2,.4,.6,.8];

operations = [1,2,3,4,4,4,4,4,4,4,4];
regD       = [0,0,0,.1,.2,.3,.4,.5,.6,.7,.8];

%operations = [1,2,3,4];
%regD       = [0,0,0,0];
%operations = [1,2];
%regD       = [0,0];
%operations = 1;
%regD       = 0;

niter = size(operations,2);

%% In order to build the (Petrov-Galerkin) Left Hand Side
if ncrack == 10001
   [ nodesu,elementsu,ntoelemu,boundaryu,orderu] = readmesh( 'meshes/rg_refined/plate_ncu.msh' );
elseif ncrack < 50 || ncrack > 99
   [ nodesu,elementsu,ntoelemu,boundaryu,orderu] = readmesh( 'meshes/rg_refined/plate_nu.msh' );
else
   [ nodesu,elementsu,ntoelemu,boundaryu,orderu] = readmesh( 'meshes/rg_refined/plate_nu_r.msh' );
end
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

% Build the shortcut list of ccordinates of segment elements
coorseg = [ nodesu(segel(:,1),1), nodesu(segel(:,1),2), ...
            nodesu(segel(:,2),1), nodesu(segel(:,2),2), ...
            (nodesu(segel(:,1),1)+nodesu(segel(:,2),1))/2, ...
            (nodesu(segel(:,1),2)+nodesu(segel(:,2),2))/2 ];

if recompute == 1
   % Boundary conditions
   dirichlet  = [0,1,0 ; 0,2,0 ; 0,3,0];
   neumann1   = [3,2,fscalar ; 1,2,-fscalar];
   neumann2   = [2,1,fscalar ; 4,1,-fscalar];
   neumann3   = [3,1,fscalar ; 3,2,fscalar ; 2,1,fscalar ; 2,2,fscalar ; ...
                 1,1,-fscalar ; 1,2,-fscalar ; 4,1,-fscalar ; 4,2,-fscalar];
   neumann4   = [3,1,-fscalar ; 3,2,fscalar ; 2,1,fscalar ; 2,2,-fscalar ; ...
                 1,1,fscalar ; 1,2,-fscalar ; 4,1,-fscalar ; 4,2,fscalar];
   
   % First, import the mesh
   if ncrack == 2
      [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc2.msh' );
   elseif ncrack == 1
      [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc.msh' );
   elseif ncrack == 3
      [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc3.msh' );
   elseif ncrack == 5
      [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc.msh' );
   elseif ncrack == 7
      [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc.msh' );
   elseif ncrack == 9
      [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc9.msh' );
   elseif ncrack == 11
      [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc11.msh' );
   elseif ncrack == 51
      [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc.msh' );
   elseif ncrack == 101
      [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc_r.msh' );
   elseif ncrack == 103
      [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc_r_3.msh' );
   elseif ncrack == 111
      [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc_r_11.msh' );
   elseif ncrack == 109
      [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc_r_9.msh' );
   elseif ncrack == 1001
      [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc.msh' );
   elseif ncrack == 10001
      [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc.msh' );
   else
      warning('This test-case was not meshed')
   end
   
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
   uin = K\f1;
   u1 = uin(1:2*nnodes,1);
   f1 = Kinter*u1;
   
   f2  = loading(nbloq,nodes,boundary,neumann2);
   uin = K\f2;
   u2 = uin(1:2*nnodes,1);
   f2 = Kinter*u2;
   
   f3  = loading(nbloq,nodes,boundary,neumann3);
   uin = K\f3;
   u3 = uin(1:2*nnodes,1);
   f3 = Kinter*u3;
   
   f4  = loading(nbloq,nodes,boundary,neumann4);
   uin = K\f4;
   u4 = uin(1:2*nnodes,1);
   f4 = Kinter*u4;
   
   ui = reshape(u1,2,[])';  ux = ui(:,1);  uy = ui(:,2);
   
   % Compute stress :
   sigma = stress(u1,E,nu,nodes,elements,order,1,ntoelem);
   
   % Output :
   plotGMSH({u1,'U1';u2,'U2';u3,'U3';u4,'U4'}, elements, nodes, 'output/reference');
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Import the uncracked domain
   if ncrack < 100 || ncrack > 999
      [ nodes2,elements2,ntoelem2,boundary2,order2] = readmesh( 'meshes/rg_refined/plate_nn.msh' );
   else
      [ nodes2,elements2,ntoelem2,boundary2,order2] = readmesh( 'meshes/rg_refined/plate_nn_r.msh' );
   end
   
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
   % third
   ur3 = zeros( 2*nnodes2, 1 ); fr3 = zeros( 2*nnodes2, 1 );
   ur3(indexbound2) = u3(indexbound);
   fr3(indexbound2) = f3(indexbound);
   % fourth
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
   
%   plotGMSH({ur1,'U_bound'}, elements2, nodes2, 'bound');
   
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
%   sigr1 = zeros( nboun1, 3*order ); sigr2 = zeros( nboun1, 3*order );
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
   
   disp([ 'Direct problem solved and data management ', num2str(toc) ]);
   tic
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% Identify the crack
   % Build the polynomial test functions.
   if ncrack<1000 || ncrack > 9999
      load('conditions10_2d.mat','-ascii');
      M = spconvert(conditions10_2d); clear('conditions10_2d');
   else
      load('conditions20_2d.mat','-ascii');
      M = spconvert(conditions20_2d); clear('conditions20_2d');
   end
   
   nmax = 10;
   ncoef =2*(nmax+1)^2; neq = ncoef;
   
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
   
      Rhs1(k) = 0; Rhs2(k) = 0; Rhs3(k) = 0; Rhs4(k) = 0;
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

else
   if ncrack == 1
      Anb = load('fields/matrix.mat');
      % The mesh is still needed
      [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc.msh' );
   elseif ncrack == 2
      Anb = load('fields/matrix2.mat');
      [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc2.msh' );
   elseif ncrack == 3
      Anb = load('fields/matrix3.mat');
      [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc3.msh' );
   elseif ncrack == 5
      Anb = load('fields/matrix5.mat');
      [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc.msh' );
   elseif ncrack == 7
      Anb = load('fields/matrix7.mat');
      [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc.msh' );
   elseif ncrack == 9
      Anb = load('fields/matrix9.mat');
      [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc9.msh' );
   elseif ncrack == 11
      Anb = load('fields/matrix11.mat');
      [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc11.msh' );
   elseif ncrack == 51
      Anb = load('fields/matrix51.mat');
      [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc.msh' );
   elseif ncrack == 101
      Anb = load('fields/matrix101.mat');
      [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc_r.msh' );
   elseif ncrack == 103
      Anb = load('fields/matrix103.mat');
      [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc_r_3.msh' );
   elseif ncrack == 109
      Anb = load('fields/matrix109.mat');
      [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc_r_9.msh' );
   elseif ncrack == 111
      Anb = load('fields/matrix111.mat');
      [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc_r_11.msh' );
   elseif ncrack == 1001
      Anb = load('fields/matrix1001.mat');
      [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc.msh' );
   elseif ncrack == 10001
      Anb = load('fields/matrix10001.mat');
      [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nc.msh' );
   end
   Lhs = Anb.Lhs; Rhs1 = Anb.Rhs1; Rhs2 = Anb.Rhs2;
   Rhs3 = Anb.Rhs3; Rhs4 = Anb.Rhs4;
   
end
%tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve the linear system and recover the volumic loading
Rhs = Rhs1;
A = Lhs'*Lhs; b = Lhs'*Rhs;

res1 = zeros(niter,1); res2 = zeros(niter,1);
res3 = zeros(niter,1); res4 = zeros(niter,1);

%% Build the regularization matrix
%L = eye(size(A));
%ninfty = 0;

% Initialize the vectors
nogap1  = ones(nseg,1); nogap2  = ones(nseg,1);
nogap3  = ones(nseg,1); nogap4  = ones(nseg,1);

oldauthorized = 1:nseg;

ind1 = zeros(niter,1); % Stores the stopping index
ind2 = zeros(niter,1);
ind3 = zeros(niter,1);
ind4 = zeros(niter,1);

for i = 1:niter

   % Build the non-local criterion
   criterion1 = co(1)*nogap1 + co(2)*nogap2 + co(3)*nogap3 + co(4)*nogap4;
   criterion  = criterion1;

   if operations(i) == 1
      authorized = find(criterion>=-1); % Take everybody at the first iteration
   elseif operations(i) == 2 % Keep only chains and identify them
      forbitten = [];
      for j=1:size(nodesu,1)
         segs = find( segel == j );
         segs = mod(segs-1,nseg)+1; % Because of the 2 columns
         
         if size(segs,1) > 1
            m1 = max(criterion(segs));
            s1 = segs ( find( criterion(segs) == m1 ) ); s1 = s1(1); % Case there are 2 equal values
            segs2 = setdiff(segs,s1);
            m2 = max(criterion(segs2));
            s2 = segs2 ( find( criterion(segs2) == m2 ) ); s2 = s2(1);
            
            forbitten = [forbitten ; setdiff(segs,[s1;s2]) ];
         end
      end
      authorized = setdiff( (1:nseg)', forbitten);
      
      % Identify the chains
      chains = {};
      toexplore = authorized;  % Next elements that remain to explore
      sanity1 = 0;
      while size(toexplore)>0 && sanity1 < 1000
         active = toexplore(1);
         curchain = active; % Current chain to explore
         
         sanity2 = 0;
         while size(active) > 0 && sanity2 < 1000
            % Find the neighours of active elements
            active2 = [];
            for j=1:size(active,2)
               act = active(j);
               no1 = segel(act,1); no2 = segel(act,2);
               segs1 = find( segel(toexplore,:) == no1 );
               segs1 = mod(segs1-1,size(toexplore,1))+1;
               seg1  = toexplore(segs1);
               seg1 = setdiff(seg1,act); % We're looking only for the other one
               segs2 = find( segel(toexplore,:) == no2 );
               segs2 = mod(segs2-1,size(toexplore,1))+1;
               seg2  = toexplore(segs2);
               seg2 = setdiff(seg2,act); %
               
               if size(seg1) > 0 % There must be an other segment AND it mussn't be in curchain
                  if size(find(curchain==seg1),2) == 0
                     % Last condition : inner product (angle)
                     noa = segel(act,:);
                     nos = segel(seg1,:);
                     nC1 = intersect(noa, nos); nC2 = setdiff(noa, nC1); nC3 = setdiff(nos, nC1);
                     C1 = nodesu(nC1,:); C2 = nodesu(nC2,:); C3 = nodesu(nC3,:);
                     C1C2 = C2-C1; C1C3 = C3-C1;
                     if C1C2*C1C3' < 0
                        active2(end+1) = seg1;
                     end
                  end
               end
               if size(seg2) > 0
                  if size(find(curchain==seg2),2) == 0
                     noa = segel(act,:);
                     nos = segel(seg2,:);
                     nC1 = intersect(noa, nos); nC2 = setdiff(noa, nC1); nC3 = setdiff(nos, nC1);
                     C1 = nodesu(nC1,:); C2 = nodesu(nC2,:); C3 = nodesu(nC3,:);
                     C1C2 = C2-C1; C1C3 = C3-C1;
                     if C1C2*C1C3' < 0
                        active2(end+1) = seg2;
                     end
                  end
               end
            end
            curchain = [curchain,active2];
            active = active2;
            sanity2 = sanity2 + 1;
         end
         chains{end+1} = curchain;
         toexplore = setdiff( toexplore, curchain );
         sanity1 = sanity1+1;

      end
      
   elseif operations(i) == 3 % Kill monoms
      forbitten = [];
      authorized = oldauthorized;
      
      % Suppress chains of size 1
      toremove = [];
%      size(chains,2)
      for j=1:size(chains,2)
         curchain = chains{j};
         if size(curchain,2) == 1
            forbitten(end+1) = curchain;
            toremove(end+1) = j;
         end
      end
      chains(toremove) = [];%
%      size(chains,2)
%      for j1=1:size(authorized,1)
%         j = authorized(j1);
%         no1 = segel(j,1); no2 = segel(j,2);
%         segs1 = find( segel(authorized,:) == no1 );% segs1 = mod(segs1-1,size(authorized,1))+1;
%         segs2 = find( segel(authorized,:) == no2 );% segs2 = mod(segs2-1,size(authorized,1))+1;
%         
%         if size(segs1) == 1 && size(segs2) == 1
%            forbitten = [forbitten ; j ];
%         end
%      end
      authorized = setdiff( authorized, forbitten);
      
   elseif operations(i) == 4 % Look for the largest chain
      authorized = oldauthorized;
      gaps = zeros( size(chains,2), 1);
      
      for j=1:size(chains,2)
         curchain = chains{j};
         for k = 1:size(curchain,2)
            no = segel( curchain(k), : );
            len = sqrt( (nodesu(no(1),1)-nodesu(no(2),1))^2 + ...
                        (nodesu(no(1),2)-nodesu(no(2),2))^2 );
            gaps(j) = gaps(j) + len*criterion( curchain(k) );
         end
      end
      
      if regD(i) == 0
         [m,nimozq] = max(gaps);
         % All the other chains are sealed (#totalitarism)
         authorized = unique( chains{nimozq}' );
      else % Each chain suppress smallest chains in its sourroundings
      
         % Preliminary : find the nodes of a chain
         nodechain = cell( 1, size(chains,2) );
         for j=1:size(chains,2)
            curchain = chains{j};
            curnodes = segel(curchain,:);
            curnodes = unique( reshape(curnodes,[],1) );
            nodechain{j} = curnodes;
         end
      
         forbitten = [];
         for j=1:size(chains,2)
         
            if gaps(j) == 0
               continue;
            end
         
            curchain = nodechain{j};
            for k=1:size(chains,2) % Look for an other chain to kill
               if j==k % Don't autokill
                  continue;
               end
               if gaps(j) < gaps(k) % Small beast doesn't eat the bigger
                  continue;
               end
               
               secchain = nodechain{k};
%               nbin = 0; % Count the number of nodes of the prey inside the
                          % radius of at least one node of the predator
               iskillable = 1;
               for n=1:size(secchain,1)
                  Cs = nodesu(secchain(n),:);
                  iskillable = 0;
                  for m=1:size(curchain,1)
                     Cc = nodesu(curchain(m),:);
                     if norm(Cs-Cc) < regD(i) % One more inside !
%                        nbin = nbin + 1;
                        iskillable = 1;
                        break;
                     end
                  end
                  if iskillable == 0 % One node isn't in the radius
                     break;
                  end
               end
%               nbin
%               size(secchain,1)
               if iskillable == 1 % Kill him
                  forbitten = [ forbitten , chains{k} ];
               end
%               if nbin == size(secchain,1) % The prey is close to the predator
%                  forbitten = [ forbitten , chains{k} ];
%%               elseif nbin > size(secchain,1)
%%                  disp('totoPaKontan');
%               end
               
            end
         end
      end
      authorized = setdiff( authorized, forbitten);
   else
      authorized = oldauthorized;
   end
   
   oldauthorized = authorized;
   authorized = [ 2*authorized-1 ; 2*authorized ];

   [Q,Theta] = eig( A(authorized,authorized) );%,L); Q = Q*(Q'*L*Q)^(-1/2);
   thetas = diag(Theta);
   [thetas,Ind] = sort( thetas,'descend' );
   Q = Q(:,Ind);
   Thetas = diag(thetas); 
   Theta = Theta(Ind,Ind);
   
%   disp([ 'Rectangular system pinversed ', num2str(toc) ]);
   % Plot the Picard stuff
   imax = min( find(thetas/thetas(1)<1e-16) );
   if size(imax,1) == 0
       imax = size(thetas,1);
   end
   
   tplo = thetas(1:imax); b1 = Lhs'*Rhs1; bplo1 = Q'*b1(authorized);
   bplo1 = bplo1(1:imax); rplo1 = (Q'*b1(authorized))./thetas;
   rplo1 = rplo1(1:imax); b2 = Lhs'*Rhs2; bplo2 = Q'*b2(authorized);
   bplo2 = bplo2(1:imax); rplo2 = (Q'*b2(authorized))./thetas;
   rplo2 = rplo2(1:imax); b3 = Lhs'*Rhs3; bplo3 = Q'*b3(authorized);
   bplo3 = bplo3(1:imax); rplo3 = (Q'*b3(authorized))./thetas;
   rplo3 = rplo3(1:imax); b4 = Lhs'*Rhs4; bplo4 = Q'*b4(authorized);
   bplo4 = bplo4(1:imax); rplo4 = (Q'*b4(authorized))./thetas;
   rplo4 = rplo4(1:imax);
   
   % Remove Zeros in rploi
   me1 = mean(abs(rplo1))/1e5; arplo1 = max(me1,abs(rplo1));
   me2 = mean(abs(rplo2))/1e5; arplo2 = max(me2,abs(rplo2));
   me3 = mean(abs(rplo3))/1e5; arplo3 = max(me3,abs(rplo3));
   me4 = mean(abs(rplo4))/1e5; arplo4 = max(me4,abs(rplo4));   

   % Compute Picard stopping indices (induces singular matrix warning...)
   ind1(i) = findPicard2 (log10(arplo1), 7, 1, 2);
   ind2(i) = findPicard2 (log10(arplo2), 7, 1, 2);
   ind3(i) = findPicard2 (log10(arplo3), 7, 1, 2);
   ind4(i) = findPicard2 (log10(arplo4), 7, 1, 2);
   
   % Postpro
%   if i>1
%      ind1(i) = min(ind1(i),ind1(i-1));
%      ind2(i) = min(ind2(i),ind2(i-1));
%      ind3(i) = min(ind3(i),ind3(i-1));
%      ind4(i) = min(ind4(i),ind4(i-1));
%   end
   
   if i==1
      figure
      hold on;
      plot(log10(abs(tplo)),'Color','green');
      plot(log10(abs(bplo1)),'Color','red');
      plot(log10(abs(rplo1)),'Color','black');
      plot(log10(abs(bplo2)),'Color','magenta');
      plot(log10(abs(rplo2)),'Color','blue');
      legend('Singular values','Rhs1','sol1','Rhs2','sol2');
      
%      [indm,p] = findPicard2 (log10(abs(rplo3)), 7, 1, 2);
%      n = size(p,1);
%      t = 1:.05:imax; tt = zeros(n,20*(imax-1)+1);
%      for j=1:n
%         tt(j,:) = t.^(n-j);
%      end
%      px = p'*tt;
%      
%      figure
%      hold on;
%      plot(log10(abs(rplo3)),'Color','black');
%      plot(t,px,'Color','red');
%      legend('sol3','poly3');

      figure
      hold on;
      plot(log10(abs(tplo)),'Color','green');
      plot(log10(abs(bplo3)),'Color','red');
      plot(log10(abs(rplo3)),'Color','black');
      plot(log10(abs(bplo4)),'Color','magenta');
      plot(log10(abs(rplo4)),'Color','blue');
      legend('Singular values','Rhs3','sol3','Rhs4','sol4');

%      figure
%      hold on;
%      plot(log10(abs(rplo1)),'Color','black');
%      plot(log10(abs(rplo2)),'Color','blue');
%      plot(log10(abs(rplo3)),'Color','red');
%      plot(log10(abs(rplo4)),'Color','green');
%      legend('sol1','sol2','sol3','sol4');
   end
   
   % Filter eigenvalues
   if jmax == 0
      jmax = size(Thetas,1);
      jmax1 = ind1; jmax2 = ind2; jmax3 = ind3; jmax4 = ind4;
   else
      jmax = min( size(Thetas,1) , jmax );
      jmax1 = jmax; jmax2 = jmax; jmax3 = jmax; jmax4 = jmax;
   end
   
   ThetaT  = Thetas( 1:jmax , 1:jmax );
   ThetaT1 = Thetas( 1:jmax1 , 1:jmax1 );
   ThetaT2 = Thetas( 1:jmax2 , 1:jmax2 );
   ThetaT3 = Thetas( 1:jmax3 , 1:jmax3 );
   ThetaT4 = Thetas( 1:jmax4 , 1:jmax4 );
   
   bT1 = Q'*b1(authorized); bT1 = bT1(1:jmax1);
   bT2 = Q'*b2(authorized); bT2 = bT2(1:jmax2);
   bT3 = Q'*b3(authorized); bT3 = bT3(1:jmax3);
   bT4 = Q'*b4(authorized); bT4 = bT4(1:jmax4);
   
   Solu1 = zeros(size(b)); Solu2 = zeros(size(b));
   Solu3 = zeros(size(b)); Solu4 = zeros(size(b));
   
   Solu1(authorized) = Q(:,1:jmax1) * (ThetaT1\bT1);
   Solu2(authorized) = Q(:,1:jmax2) * (ThetaT2\bT2);
   Solu3(authorized) = Q(:,1:jmax3) * (ThetaT3\bT3);
   Solu4(authorized) = Q(:,1:jmax4) * (ThetaT4\bT4);
   
   res1(i) = norm( Lhs'*Lhs*Solu1 - Lhs'*Rhs1 ) / norm( Lhs'*Rhs1 );
   res2(i) = norm( Lhs'*Lhs*Solu2 - Lhs'*Rhs2 ) / norm( Lhs'*Rhs2 );
   res3(i) = norm( Lhs'*Lhs*Solu3 - Lhs'*Rhs3 ) / norm( Lhs'*Rhs3 );
   res4(i) = norm( Lhs'*Lhs*Solu4 - Lhs'*Rhs4 ) / norm( Lhs'*Rhs4 );
   
   % Re-loop over segments to build normal gap
   nogap  = zeros(nseg,1); nogap1 = zeros(nseg,1); nogap2 = zeros(nseg,1);
   nogap3 = zeros(nseg,1); nogap4 = zeros(nseg,1);
   ntoseg = zeros(nnodesu,1);
   Sgap   = zeros(nnodesu,1);
   for j=1:nseg
      no1 = segel(j,1); no2 = segel(j,2);
      nogap1(j) = norm( Solu1( [2*j-1,2*j] ) );
      nogap2(j) = norm( Solu2( [2*j-1,2*j] ) );
      nogap3(j) = norm( Solu3( [2*j-1,2*j] ) );
      nogap4(j) = norm( Solu4( [2*j-1,2*j] ) );
      nogap(j)  = sqrt(nogap1(j)^2*nogap2(j)^2);
   end

end

%% Segments visu
%figure;
%hold on;
%plot( [nodes(5,1),nodes(6,1)], [nodes(5,2),nodes(6,2)], 'Color', [.6,.6,.6], 'LineWidth', 5 );
%if ncrack == 2
%   plot( [nodes(7,1),nodes(8,1)], [nodes(7,2),nodes(8,2)], 'Color', [.6,.6,.6], 'LineWidth', 5 );
%end
%nogapp1 = abs(nogap1)-min(abs(nogap1)); maxn1 = max(nogapp1);
%%      for i=1:nseg
%for j=1:size(oldauthorized,1)
%   i = oldauthorized(j);
%   no1 = segel(i,1); no2 = segel(i,2);
%   x1 = nodesu(no1,1); y1 = nodesu(no1,2);
%   x2 = nodesu(no2,1); y2 = nodesu(no2,2);
%   
%   x = nogapp1(i)/maxn1;
%   rgb = rgbmap(x);
%   plot( [x1,x2], [y1,y2], 'Color', rgb, 'LineWidth', 3 );
%end
%%axis equal;
%axis([0 1 0 1]);
%colormap("default")
%h = colorbar();
%ytick = get (h, "ytick");
%set (h, "yticklabel", sprintf ( "%g|", maxn1*ytick+min(nogap1) ));

% Segments visu
figure;
hold on;
plot( [nodes(5,1),nodes(6,1)], [nodes(5,2),nodes(6,2)], 'Color', [.6,.6,.6], 'LineWidth', 5 );
if ncrack == 2
   plot( [nodes(7,1),nodes(8,1)], [nodes(7,2),nodes(8,2)], 'Color', [.6,.6,.6], 'LineWidth', 5 );
elseif ncrack == 9 || ncrack == 109 || ncrack == 11 || ncrack == 111
   plot( [nodes(7,1),nodes(6,1)], [nodes(7,2),nodes(6,2)], 'Color', [.6,.6,.6], 'LineWidth', 5 );
end
if ncrack == 11 || ncrack == 111
   plot( [nodes(7,1),nodes(8,1)], [nodes(7,2),nodes(8,2)], 'Color', [.6,.6,.6], 'LineWidth', 5 );
end
nogapp1 = co(1)*nogap1 + co(2)*nogap2 + co(3)*nogap3 + co(4)*nogap4;
maxn1 = max(nogapp1);
%      for i=1:nseg
for j=1:size(oldauthorized,1)
   i = oldauthorized(j);
   no1 = segel(i,1); no2 = segel(i,2);
   x1 = nodesu(no1,1); y1 = nodesu(no1,2);
   x2 = nodesu(no2,1); y2 = nodesu(no2,2);
   
   x = nogapp1(i)/maxn1;
   rgb = rgbmap(x);
   plot( [x1,x2], [y1,y2], 'Color', rgb, 'LineWidth', 3 );
end
%axis equal;
axis([0 1 0 1]);
colormap("default")
h = colorbar();
ytick = get (h, "ytick");
set (h, "yticklabel", sprintf ( "%g|", maxn1*ytick+min(nogap1) ));

%% Segments visu
%figure;
%hold on;
%plot( [nodes(5,1),nodes(6,1)], [nodes(5,2),nodes(6,2)], 'Color', [.6,.6,.6], 'LineWidth', 5 );
%if ncrack == 2
%   plot( [nodes(7,1),nodes(8,1)], [nodes(7,2),nodes(8,2)], 'Color', [.6,.6,.6], 'LineWidth', 5 );
%end
%%nogapp1 = abs(nogap1+nogap2); maxn1 = max(nogapp1);
%nogapp1 = criterion1; maxn1 = max(criterion1);
%%      for i=1:nseg
%for j=1:size(oldauthorized,1)
%   i = oldauthorized(j);
%   no1 = segel(i,1); no2 = segel(i,2);
%   x1 = nodesu(no1,1); y1 = nodesu(no1,2);
%   x2 = nodesu(no2,1); y2 = nodesu(no2,2);
%   
%   x = nogapp1(i)/maxn1;
%   rgb = rgbmap(x);
%   plot( [x1,x2], [y1,y2], 'Color', rgb, 'LineWidth', 3 );
%end
%%axis equal;
%axis([0 1 0 1]);
%colormap("default")
%h = colorbar();
%ytick = get (h, "ytick");
%set (h, "yticklabel", sprintf ( "%g|", maxn1*ytick+min(nogap1) ));

%figure;
%hold on;
%plot(res1);
%plot(res2,'Color','black');
%
figure;
hold on;
plot(ind1);
plot(ind2,'Color','black');
plot(ind3,'Color','red');
plot(ind4,'Color','green');
legend('Picard index 1','Picard index 2','Picard index 3','Picard index 4');