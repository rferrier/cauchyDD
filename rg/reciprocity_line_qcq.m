% 21/04/2017
% Détection de source par écart à la réciprocité
% Intégrations par PG

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E          = 210000; % MPa : Young modulus
nu         = 0.3;    % Poisson ratio
fscalar    = 250;    % N.mm-1 : Loading on the plate
mat        = [0, E, nu];
br         = .0;      % Noise level

% Boundary conditions
% first index  : index of the boundary
% second index : 1=x, 2=y
% third        : value
% [0,1,value] marks a dirichlet regularization therm on x
dirichlet  = [0,1,0 ; 0,2,0 ; 0,3,0];

% First, import the mesh
[ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_nn.msh' );

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
L = Xmax-Xmin;

f = volumicLoad( 3, nodes, elements, 2, [-fscalar*Xmoy,0;fscalar,0] );
uin = K\f;
u = uin(1:2*nnodes,1);
f = Kinter*u;

ui = reshape(u,2,[])';  ux = ui(:,1);  uy = ui(:,2);

% Compute stress :
sigma = stress(u,E,nu,nodes,elements,order,1,ntoelem);

% Output :
plotGMSH({ux,'U_x';uy,'U_y';u,'U_vect';sigma,'stress'}, elements, nodes, 'output/reference');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create the data vectors
ur = zeros( 2*nnodes, 1 ); fr = zeros( 2*nnodes, 1 );
ur(indexbound) = u(indexbound);
%fr(indexbound) = f(indexbound); % All the boundaries are free Neumann

% Add the noise
%u1n = u1; u2n = u2;
%br1 = randn(2*nnodes,1); br2 = randn(2*nnodes,1);
%% noise = load('noises/105.mat'); br1 = noise.br1; br2 = noise.br2;
%u1 = ( 1 + br*br1 ) .* u1; u2 = ( 1 + br*br2 ) .* u2;

%plotGMSH({ur1,'U_bound'}, elements2, nodes2, 'bound');

%% Preliminary stuff : find the volumic elements corresponding to the boundaries
nboun = size(boundary,1); nelem = size(elements,1);
boun2vol = zeros( nboun, 1 ); extnorm = zeros( nboun, 2 );
urr = zeros( nboun, 2+2*order );
for i=1:nboun
   % Volumic element
   no1 = boundary(i,2); no2 = boundary(i,3); % only with 2 nodes even if order > 1
   cand1 = rem( find(elements==no1),nelem ); % find gives line + column*size
   cand2 = rem( find(elements==no2),nelem );
   boun2vol(i) = intersect(cand1, cand2); % If everything went well, there is only one
   
   % Exterior normal
   elt = boun2vol(i); no3 = setdiff( elements( elt, 1:3 ), [no1,no2]);
   x1 = nodes(no1,1); y1 = nodes(no1,2);
   x2 = nodes(no2,1); y2 = nodes(no2,2);
   x3 = nodes(no3,1); y3 = nodes(no3,2);
   extnorm(i,:) = [ y1-y2 , -(x1-x2) ];
   extnorm(i,:) = extnorm(i,:)/norm(extnorm(i,:));
   if (x3-x2)*(y1-y2) - (y3-y2)*(x1-x2) > 0 % Check that the normal is exterior
      extnorm(i,:) = -extnorm(i,:);
   end
   
   % ur
   urr(i,1:4) = u( [2*no1-1,2*no1,2*no2-1,2*no2] );
   if order == 2
      no4 = boundary(i,4);
      urr(i,5:6) = u( [2*no4-1,2*no4] );
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Identify the volumic effort
% Build the polynomial test functions.
load('conditions30_2d.mat');%,'ascii');
M = spconvert(conditions30_2d); clear('conditions30_2d');
nmax = 30;
ncoef =2*(nmax+1)^2; neq = ncoef;

% Suppress zeros rows
toremove = find( norm(M,'inf','rows') == 0 );
M(toremove,:) = [];

coef   = null(full(M));
nftest = size(coef,2); % Nb of avaliable test functions

%% Build the (Petrov-Galerkin) Left Hand Side
[ nodesu,elementsu,ntoelemu,boundaryu,orderu] = readmesh( 'meshes/rg_refined/plate_nu.msh' );
nnodesu = size(nodesu,1);

nelemu = size(elementsu,1);
Lhs = zeros(nftest,2*nelemu); % Decomposition functions are constant per element

for j=1:nelemu % Compute the integrals
   bonod = elementsu(j,:);

   no1 = bonod(1); no2 = bonod(2); no3 = bonod(3);
   x1 = nodesu(no1,1); y1 = nodesu(no1,2);
   x2 = nodesu(no2,1); y2 = nodesu(no2,2);
   x3 = nodesu(no3,1); y3 = nodesu(no3,2);

   S = .5 * ((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
   
   Ng = 1;
   [ Xg, Wg ] = gaussPt( Ng );
             
   for k=1:size(Wg,1)
      xg = Xg(k,:); wg = Wg(k);
      
      % Interpolation
      xgr  = (1-xg(1)-xg(2))*[x1;y1] + xg(1)*[x2;y2] + xg(2)*[x3;y3] ; ... % abscissae

      X = xgr(1)/L; Y = xgr(2)/L;
      
      for i=1:nftest
         for ii=0:nmax
            for jj=0:nmax
               coefa = coef(1:(nmax+1)^2,i);
               coefb = coef((nmax+1)^2+1:2*(nmax+1)^2,i);
               % Distinguish x and y components 
               aij = coefa( (nmax+1)*ii + jj+1 );
               bij = coefb( (nmax+1)*ii + jj+1 );
               Lhs(i,2*j-1) = Lhs(i,2*j-1) + S * wg * aij*X^ii*Y^jj; % increment the integral
               Lhs(i,2*j)   = Lhs(i,2*j)   + S * wg * bij*X^ii*Y^jj;
            end
         end
      end

   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the RG
Rhs  = zeros(nftest,1);

for k=1:nftest

   Rhs(k) = 0;
   coefa = coef(1:(nmax+1)^2,k);
   coefb = coef((nmax+1)^2+1:2*(nmax+1)^2,k);
   
   for i=1:nboun
      bonod = boundary(i,:); exno = extnorm(i,:)';
   
      no1 = bonod(2); no2 = bonod(3);
      x1 = nodes(no1,1); y1 = nodes(no1,2);
      x2 = nodes(no2,1); y2 = nodes(no2,2);
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
            uer = transpose( (1-xg)*urr(i,1:2) + xg*urr(i,3:4) ); % [ux;uy] on the
            xgr  = ( (1-xg)*[x1;y1] + xg*[x2;y2] ); % abscissae of the Gauss point in the physical space
         elseif order==2
            uer = transpose( urr(i,1:2) + ...
                   xg*(4*urr(i,5:6)-3*urr(i,1:2)-urr(i,3:4)) + ...   % [ux;uy] on the
                   xg^2*(2*urr(i,3:4)+2*urr(i,1:2)-4*urr(i,5:6)) ); % Gauss point
            xgr  = ( [x1;y1] + xg*(4*[x3;y3]-3*[x1;y1]-[x2;y2]) + ...
                   xg^2*(2*[x2;y2]+2*[x1;y1]-4*[x3;y3]) );
         end
   
         % Free Neumann boundaries
         fer = [0;0];
         
%         ixigrec = Q'*[xgr(1);xgr(2)];
%         X = (ixigrec(1)-offset)/L-.5; Y = (ixigrec(2)+Cte2)/L;
         X = xgr(1)/L; Y = xgr(2)/L; % Use the standard basis
         
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

         sloc = 1/L*[sloc11,sloc12;sloc12,sloc22];
         spa = sloc;
         fpa = spa*exno;

         Rhs(k) = Rhs(k) + len * wg * ( fer'*vpa - fpa'*uer );
      end
   end
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve the linear system and recover the volumic loading
Solu = (Lhs'*Lhs)\Lhs'*Rhs;

% Redistribute the identified loading on the nodes (for visu)
Sno = zeros(2*nnodesu,1);
for i=1:size(elementsu,1)
   Xloc = nodes(elementsu(i,:),:);    % Extract and adapt coords
   nno = size(Xloc,1);

   ne = size(elementsu,2);
   Se = []; maps = []; %TODEBUG
   for j=1:ne
      Se = [ Se, Solu(2*i-1)/ntoelem(elementsu(i,j),1), ...
                 Solu(2*i)/ntoelem(elementsu(i,j),1) ]; % Pass to nodes
      maps = [maps,2*elementsu(i,j)-1,2*elementsu(i,j)];
   end
 
   Sno(maps,1) = Sno(maps,1) + Se';
end

% Output
Si = reshape(Sno,2,[])';  Sx = Si(:,1);  Sy = Si(:,2);
plotGMSH({Sx,'F_x';Sy,'F_y';Sno,'F_vect'}, elementsu, nodesu, 'output/identification');