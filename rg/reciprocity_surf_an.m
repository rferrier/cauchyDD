%05/04/2017
% Détection de fissure 3D plane par écart à la réciprocité
% Intégrations par PG, matrice analytique

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 210000; % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 250;    % N.mm-2 : Loading
mat = [0, E, nu];
regmu   = 0;      % Regularization parameter

ordp      = 2; % Order of polynom
loadfield = 2; % If 0 : recompute the reference problem and re-pass mesh
               % If 2 : meshes are conformal, do everything
               % If 3 : meshes are conformal, store the u field
               % If 4 : meshes are conformal, read the u field

usepolys   = 1;
plotref    = 1;
comperror  = 1;

% Load the polynom's matrix
load('../analytique/conditions20.mat','-ascii');
Expression1 = spconvert(conditions20);
M = Expression1;
% load('../analytique/conditions.mat');
% Expression1 = sparse(Expression1);
% M = Expression1;

nmax = 20;
ncoef = 3*(nmax+1)^3;
neq = 3*((nmax+1)^3+(nmax+1)^2);

% Sanity check
if 2*ordp > nmax
    warning('Interpolation order is too high wrt computed polynoms');
end

% cracked_mesh = 'meshes/rg3dpp/plate_c_710t10u.msh';
% uncracked_mesh = 'meshes/rg3dpp/plate710t10u.msh';
% cracked_mesh = 'meshes/rg3dpp/plate_c_710t10.msh';
% uncracked_mesh = 'meshes/rg3dpp/plate710t10.msh';
cracked_mesh = 'meshes/rg3dm/platem_c.msh';
uncracked_mesh = 'meshes/rg3dm/platem.msh';

centCrack = [4;3;1]; % Point on the crack (for reference)

if loadfield ~= 1 && loadfield ~= 4
   tic
   % Boundary conditions
   % first index  : index of the boundary
   % second index : 1=x, 2=y
   % third        : value
   % [0,1,value] marks a dirichlet regularization therm on x
   dirichlet  = [0,1,0 ; 0,2,0 ; 0,3,0 ; 0,4,0 ; 0,5,0 ; 0,6,0];
   neumann1   = [2,3,fscalar ; 1,3,-fscalar];
   neumann2   = [4,1,fscalar ; 6,1,-fscalar];
   
   % First, import the mesh
   [ nodes,elements,ntoelem,boundary,order] = readmesh3D( cracked_mesh );
   nnodes = size(nodes,1);
   
   % mapBounds
   [ node2b1, b2node1 ]   = mapBound3D( 1, boundary, nnodes );
   [ node2b2, b2node2 ]   = mapBound3D( 2, boundary, nnodes );
   [ node2b3, b2node3 ]   = mapBound3D( 3, boundary, nnodes );
   [ node2b4, b2node4 ]   = mapBound3D( 4, boundary, nnodes );
   [ node2b5, b2node5 ]   = mapBound3D( 5, boundary, nnodes );
   [ node2b6, b2node6 ]   = mapBound3D( 6, boundary, nnodes );
   
   indexbound  = [3*b2node1-2 ; 3*b2node1-1 ; 3*b2node1 ;...
                  3*b2node2-2 ; 3*b2node2-1 ; 3*b2node2 ;...
                  3*b2node3-2 ; 3*b2node3-1 ; 3*b2node3 ;...
                  3*b2node4-2 ; 3*b2node4-1 ; 3*b2node4 ;...
                  3*b2node5-2 ; 3*b2node5-1 ; 3*b2node5 ;...
                  3*b2node6-2 ; 3*b2node6-1 ; 3*b2node6 ];
   
   % Then, build the stiffness matrix :
   [K,C,nbloq,node2c,c2node] = Krig3D (nodes,elements,mat,order,boundary,dirichlet);
   Kinter = K( 1:3*nnodes, 1:3*nnodes );
   
   f1  = loading3D(nbloq,nodes,boundary,neumann1);
   uin = K\f1;
   u1 = uin(1:3*nnodes,1);
   f1 = Kinter*u1;
   
   f2  = loading3D(nbloq,nodes,boundary,neumann2);
   uin = K\f2;
   u2 = uin(1:3*nnodes,1);
   f2 = Kinter*u2;
   
   ui = reshape(u1,3,[])'; ux = ui(:,1); uy = ui(:,2); uz = ui(:,3);
   
   % Compute stress :
   sigma = stress3D(u1,mat,nodes,elements,order,1,ntoelem);
   
   % Save the field
   if loadfield == 3
      save('fields/U.mat','u1','u2','f1','f2');
   end
   
   % Output :
   plotGMSH3D({ux,'U_x';uy,'U_y';uz,'U_z';u1,'U_vect';sigma,'stress'}, elements, nodes, 'reference');
   disp([ 'Direct problem solved ', num2str(toc) ]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import the uncracked domain /!\ MUST BE THE SAME !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (except for the crack)
[ nodes2,elements2,ntoelem2,boundary2,order2] = readmesh3D( uncracked_mesh );
nnodes2 = size(nodes2,1);
[K2,C2,nbloq2,node2c2,c2node2] = Krig3D (nodes2,elements2,mat,order2,boundary2,[]);
Kinter2 = K2( 1:3*nnodes2, 1:3*nnodes2 );

% Mapbounds
[ node2b12, b2node12 ] = mapBound3D( 1, boundary2, nnodes2 );
[ node2b22, b2node22 ] = mapBound3D( 2, boundary2, nnodes2 );
[ node2b32, b2node32 ] = mapBound3D( 3, boundary2, nnodes2 );
[ node2b42, b2node42 ] = mapBound3D( 4, boundary2, nnodes2 );
[ node2b52, b2node52 ] = mapBound3D( 5, boundary2, nnodes2 );
[ node2b62, b2node62 ] = mapBound3D( 6, boundary2, nnodes2 );

indexbound2 = [3*b2node12-2 ; 3*b2node12-1 ; 3*b2node12 ;...
               3*b2node22-2 ; 3*b2node22-1 ; 3*b2node22 ;...
               3*b2node32-2 ; 3*b2node32-1 ; 3*b2node32 ;...
               3*b2node42-2 ; 3*b2node42-1 ; 3*b2node42 ;...
               3*b2node52-2 ; 3*b2node52-1 ; 3*b2node52 ;...
               3*b2node62-2 ; 3*b2node62-1 ; 3*b2node62 ];
               
if loadfield == 0
   %% Pass f and u on the uncracked mesh
   UFr = passMesh3D (nodes, elements, nodes2, elements2, [u1,u2,f1,f2], boundary2);
   ur1 = UFr(:,1); ur2 = UFr(:,2); fr1 = UFr(:,3); fr2 = UFr(:,4);
   plotGMSH3D({ur1,'u1';ur2,'u2'}, elements2, nodes2, 'us');
elseif loadfield == 1
   %load the field
   UFR = load('fields/UFr107.mat'); UFr = UFR.UFr;
   ur1 = UFr(:,1); ur2 = UFr(:,2); fr1 = UFr(:,3); fr2 = UFr(:,4);
   plotGMSH3D({ur1,'u1';ur2,'u2'}, elements2, nodes2, 'us');
elseif loadfield == 2 || loadfield == 3 %: conformal mesh /!\
   ur1 = zeros( 3*nnodes2, 1 );   ur2 = zeros( 3*nnodes2, 1 );
   fr1 = zeros( 3*nnodes2, 1 );   fr2 = zeros( 3*nnodes2, 1 );
   ur1(indexbound2) = u1(indexbound);   ur2(indexbound2) = u2(indexbound);
   fr1(indexbound2) = f1(indexbound);   fr2(indexbound2) = f2(indexbound);
elseif loadfield == 4
   % For the reference stuff
   [ nodes,elements,ntoelem,boundary,order] = readmesh3D( cracked_mesh );
   
   UU = load('fields/U.mat'); u1 = UU.u1; u2 = UU.u2; f1 = UU.f1; f2 = UU.f2;
   ui = reshape(u1,3,[])'; ux = ui(:,1); uy = ui(:,2); uz = ui(:,3);
   ur1 = zeros( 3*nnodes2, 1 );   ur2 = zeros( 3*nnodes2, 1 );
   fr1 = zeros( 3*nnodes2, 1 );   fr2 = zeros( 3*nnodes2, 1 );
   ur1(indexbound2) = u1(indexbound2);   ur2(indexbound2) = u2(indexbound2);
   fr1(indexbound2) = f1(indexbound2);   fr2(indexbound2) = f2(indexbound2);
end

%% Preliminary stuff : find the volumic elements corresponding to the boundaries
boundary2( find( boundary2(:,1)==7 ),: ) = []; % Suppress internal fictive boundary 7
nboun2 = size(boundary2,1); nelem2 = size(elements2,1);
boun2vol2 = zeros( nboun2, 1 ); extnorm2 = zeros( nboun2, 3 );
for i=1:nboun2
   % Volumic element
   no1 = boundary2(i,2); no2 = boundary2(i,3);
   no3 = boundary2(i,4); % only with 3 nodes even if order > 1
   cand1 = rem( find(elements2==no1),nelem2 ); % find gives line + column*size
   cand2 = rem( find(elements2==no2),nelem2 );
   cand3 = rem( find(elements2==no3),nelem2 );
%   cand2 = rem( find(elements2(cand1,:)==no2),size(cand1,1) );
%   cand3 = rem( find(elements2(cand2,:)==no3),size(cand2,1) );
   cand4 = intersect(cand1, cand2);
   cand5 = intersect(cand4, cand3);
   if cand5 == 0
      cand5 = nelem2; % Ok, it's not optimized
   end
   boun2vol2(i) = cand5; % If everything went well, there is only one

   % Exterior normal
   elt = boun2vol2(i); no4 = setdiff( elements2( elt, 1:4 ), [no1,no2,no3]);
   x1 = nodes2(no1,1); y1 = nodes2(no1,2); z1 = nodes2(no1,3);
   x2 = nodes2(no2,1); y2 = nodes2(no2,2); z2 = nodes2(no2,3);
   x3 = nodes2(no3,1); y3 = nodes2(no3,2); z3 = nodes2(no3,3);
   x4 = nodes2(no4,1); y4 = nodes2(no4,2); z4 = nodes2(no4,3);
   
   extnorm2(i,:) = [ (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1),...
                     -(x2-x1)*(z3-z1) + (x3-x1)*(z2-z1),...
                     (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)]; % Vectorial product
   extnorm2(i,:) = extnorm2(i,:)/norm(extnorm2(i,:));
   if extnorm2(i,:) * [x4-x1;y4-y1;z4-z1] > 0 % Check that the normal is exterior
      extnorm2(i,:) = -extnorm2(i,:);
   end
end

nboun1 = size(boundary,1); nelem1 = size(elements,1);
boun2vol1 = zeros( nboun1, 1 ); extnorm1 = zeros( nboun1, 2 );
urr1  = zeros( nboun1, 3+3*order ); urr2 = zeros( nboun1, 3+3*order );
for i=1:nboun2
%   % Volumic element
   no1 = boundary2(i,2); no2 = boundary2(i,3);
   no3 = boundary2(i,4); % only with 3 nodes even if order > 1

   % ur
   map = [3*no1-2,3*no1-1,3*no1,3*no2-2,3*no2-1,3*no2,3*no3-2,3*no3-1,3*no3];
   urr1(i,1:9) = u1( map );
   urr2(i,1:9) = u2( map );
   if order == 2
      % Rem : no4 is the last summit of the tetraedron
      no5 = boundary(i,5); no6 = boundary(i,6); no7 = boundary(i,7);
      map = [3*no5-2,3*no5-1,3*no5,3*no6-2,3*no6-1,3*no6,3*no7-2,3*no7-1,3*no7];
      urr1(i,10:18) = u1( map );
      urr2(i,10:18) = u2( map );
   end
end

% Reconstruct the reference efforts
neumann1   = [2,3,fscalar ; 1,3,-fscalar];
neumann2   = [4,1,fscalar ; 6,1,-fscalar];
fr1 = loading3D(nbloq2,nodes2,boundary2,neumann1);
fr2 = loading3D(nbloq2,nodes2,boundary2,neumann2);

tic
%% First determine the crack's normal.

% Clean redondant stuff in indexbound2
i = 1;
while i <= size(indexbound2,1)
   if find( indexbound2(i) == indexbound2(1:i-1) )
      indexbound2(i) = [];
      i = i-1;
   end
   i = i+1;
end

% Loop over boundaries to compute the integrals
R11 = 0; R21 = 0; R31 = 0; R41 = 0; R51 = 0; R61 = 0;
R12 = 0; R22 = 0; R32 = 0; R42 = 0; R52 = 0; R62 = 0;
for i=1:nboun2 % boundary1 and boundary 2 are supposed to be the same
   bonod = boundary2(i,:); exno = extnorm2(i,:)';

   no1 = bonod(2); no2 = bonod(3); no3 = bonod(4);
   x1 = nodes2(no1,1); y1 = nodes2(no1,2); z1 = nodes2(no1,3);
   x2 = nodes2(no2,1); y2 = nodes2(no2,2); z2 = nodes2(no2,3);
   x3 = nodes2(no3,1); y3 = nodes2(no3,2); z3 = nodes2(no3,3);
   vecprod = [ (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1),...
               -(x2-x1)*(z3-z1) + (x3-x1)*(z2-z1),...
               (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)];
   S = .5*norm(vecprod);
   
   if order==1
      Ng = 1;
   elseif order==2
      Ng  = 2; 
      no4 = bonod(5); no5 = bonod(6); no6 = bonod(7);
      x4  = nodes2(no4,1); y4 = nodes2(no4,2); z4 = nodes2(no4,3);
      x5  = nodes2(no5,1); y5 = nodes2(no5,2); z5 = nodes2(no5,3);
      x6  = nodes2(no6,1); y6 = nodes2(no6,2); z6 = nodes2(no6,3);
   end
   [ Xg, Wg ] = gaussPt( Ng );
             
   for j=1:size(Wg,1)
      xg = Xg(j,:); wg = Wg(j);
      
      % Interpolations
      if order==1
         uer1 = transpose( (1-xg(1)-xg(2))*urr1(i,1:3) + ... % [ux;uy] on the
                           xg(1)*urr1(i,4:6) + ...           % Gauss point
                           xg(2)*urr1(i,7:9) );        
         uer2 = transpose( (1-xg(1)-xg(2))*urr2(i,1:3) + ... % [ux;uy] on the
                           xg(1)*urr2(i,4:6) + ...           % Gauss point
                           xg(2)*urr2(i,7:9) );
         xgr  = (1-xg(1)-xg(2))*[x1;y1;z1] + xg(1)*[x2;y2;z2] + xg(2)*[x3;y3;z3] ; ... % abscissae

      elseif order==2
         uer1 = transpose( -(1-xg(1)-xg(2))*(1-2*(1-xg(1)-xg(2)))*urr1(i,1:3) + ...
                -xg(1)*(1-2*xg(1))*urr1(i,4:6) + ...           % [ux;uy] on the
                -xg(2)*(1-2*xg(2))*urr1(i,7:9) + ...          % Gauss point
                4*xg(1)*(1-xg(1)-xg(2))*urr1(i,10:12) + ...
                4*xg(1)*xg(2)*urr1(i,13:15) + ...
                4*xg(2)*(1-xg(1)-xg(2))*urr1(i,16:18) );                      
         uer2 = transpose( -(1-xg(1)-xg(2))*(1-2*(1-xg(1)-xg(2)))*urr2(i,1:3) + ...
                -xg(1)*(1-2*xg(1))*urr2(i,4:6) + ...
                -xg(2)*(1-2*xg(2))*urr2(i,7:9) + ...
                4*xg(1)*(1-xg(1)-xg(2))*urr2(i,10:12) + ...
                4*xg(1)*xg(2)*urr2(i,13:15) + ...
                4*xg(2)*(1-xg(1)-xg(2))*urr2(i,16:18) );  
         xgr  = ( -(1-xg(1)-xg(2))*(1-2*(1-xg(1)-xg(2)))*[x1;y1;z1] + ...
                -xg(1)*(1-2*xg(1))*[x2;y2;z2] + ...
                -xg(2)*(1-2*xg(2))*[x3;y3;z3] + ...
                4*xg(1)*(1-xg(1)-xg(2))*[x4;y4;z4] + ...
                4*xg(1)*xg(2)*[x5;y5;z5] + ...
                4*xg(2)*(1-xg(1)-xg(2))*[x6;y6;z6] ); 
      end

      % Reference force from the BC's
      if exno(1) == 1
         fer1 = [0;0;0]; fer2 = [fscalar;0;0];
      elseif exno(1) == -1
         fer1 = [0;0;0]; fer2 = -[fscalar;0;0];
      elseif exno(2) == 1
         fer1 = [0;0;0]; fer2 = [0;0;0];
      elseif exno(2) == -1
         fer1 = [0;0;0]; fer2 = [0;0;0];
      elseif exno(3) == 1
         fer1 = [0;0;fscalar]; fer2 = [0;0;0];
      elseif exno(3) == -1
         fer1 = -[0;0;fscalar]; fer2 = [0;0;0];
      end

       % Test fields
      s1 = [1,0,0;0,0,0;0,0,0];
      s2 = [0,0,0;0,1,0;0,0,0];
      s3 = [0,0,0;0,0,0;0,0,1];
      s4 = [0,.5,0;.5,0,0;0,0,0];
      s5 = [0,0,.5;0,0,0;.5,0,0];
      s6 = [0,0,0;0,0,.5;0,.5,0];
      
      f1 = s1*exno; f2 = s2*exno; f3 = s3*exno;
      f4 = s4*exno; f5 = s5*exno; f6 = s6*exno;

      v1 = 1/E*[ xgr(1) ; -nu*xgr(2) ; -nu*xgr(3) ];
      v2 = 1/E*[ -nu*xgr(1) ; xgr(2) ; -nu*xgr(3) ];
      v3 = 1/E*[ -nu*xgr(1) ; -nu*xgr(2) ; xgr(3) ];
      v4 = (1+nu)/(2*E)*[ xgr(2) ; xgr(1) ; 0 ];
      v5 = (1+nu)/(2*E)*[ xgr(3) ; 0 ; xgr(1) ];
      v6 = (1+nu)/(2*E)*[ 0 ; xgr(3) ; xgr(2) ];

      R11 = R11 + S * wg * ( fer1'*v1 - f1'*uer1 ); % increment the integral
      R21 = R21 + S * wg * ( fer1'*v2 - f2'*uer1 );
      R31 = R31 + S * wg * ( fer1'*v3 - f3'*uer1 );
      R41 = R41 + S * wg * ( fer1'*v4 - f4'*uer1 );
      R51 = R51 + S * wg * ( fer1'*v5 - f5'*uer1 );
      R61 = R61 + S * wg * ( fer1'*v6 - f6'*uer1 );
      
      R12 = R12 + S * wg * ( fer2'*v1 - f1'*uer2 );
      R22 = R22 + S * wg * ( fer2'*v2 - f2'*uer2 );
      R32 = R32 + S * wg * ( fer2'*v3 - f3'*uer2 );
      R42 = R42 + S * wg * ( fer2'*v4 - f4'*uer2 );
      R52 = R52 + S * wg * ( fer2'*v5 - f5'*uer2 );
      R62 = R62 + S * wg * ( fer2'*v6 - f6'*uer2 );
   end
end
RR1 = [R11, R41, R51 ; R41, R21, R61 ; R51, R61, R31];
RR2 = [R12, R42, R52 ; R42, R22, R62 ; R52, R62, R32];

% Normalize R1
normR1 = sqrt( 2*(R11^2+R21^2+R31^2+2*(R41^2+R51^2+R61^2)) - (R11+R21+R31)^2 );
RR1 = RR1/normR1;
[phi1,Lambda1] = eig( RR1 );

% Normalize R2
normR2 = sqrt( 2*(R12^2+R22^2+R32^2+2*(R42^2+R52^2+R62^2)) - (R12+R22+R32)^2 );
RR2 = RR2/normR2;
[phi2,Lambda2] = eig( RR2 );

% Vectorial product (provided the null eigenvalue is the second one)
normal = [phi1(2,2)*phi2(3,2) - phi2(2,2)*phi1(3,2)
          phi1(3,2)*phi2(1,2) - phi2(3,2)*phi1(1,2)
          phi1(1,2)*phi2(2,2) - phi2(1,2)*phi1(2,2)];

normal = normal/norm(normal);
          
% Build the base-change matrix : [x;y;z] = Q*[X;Y;Z], [X;Y;Z] = Q'*[x;y;z]
t = sign(normal(3))*[1 ; 0 ; -(normal(1))/normal(3)]; t = t/norm(t); % Total random t
v = [normal(2)*t(3) - t(2)*normal(3)
     normal(3)*t(1) - t(3)*normal(1)
     normal(1)*t(2) - t(1)*normal(2)]; % v = n^t
Q = [ t , v , normal ];

%% Then the constant

% First, find the minimal point (in order to have Cte < 0)
norep = Q'*nodes2'; K = min(norep(3,:));

vt = zeros(3*nnodes2, 1);
vv = zeros(3*nnodes2, 1);

Xs = zeros(nnodes2, 1);
Ys = zeros(nnodes2, 1);
Zs = zeros(nnodes2, 1);
for i=1:nnodes2
   x = nodes2(i,1);
   y = nodes2(i,2);
   z = nodes2(i,3);
   % Change base (for coordiantes)
   ixigrec = Q'*[x;y;z]; X = ixigrec(1); Y = ixigrec(2); Z = ixigrec(3)-K;
   Xs(i) = X; Ys(i) = Y; Zs(i) = Z;

   vloct = [ -X^2/(2*E) - nu*Y^2/(2*E) + (2+nu)*Z^2/(2*E) ;...
             nu*X*Y/E ;...
             nu*X*Z/E ];
   vlocv = [ nu*X*Y/E ;...
             -Y^2/(2*E) - nu*X^2/(2*E) + (2+nu)*Z^2/(2*E) ;...
             nu*Y*Z/E ];
   % Change base (for vector), in the other direction
   vxyt = Q*vloct; vt(3*i-2) = vxyt(1); vt(3*i-1) = vxyt(2); vt(3*i) = vxyt(3);
   vxyv = Q*vlocv; vv(3*i-2) = vxyv(1); vv(3*i-1) = vxyv(2); vv(3*i) = vxyv(3);
end

% Norm of [[ut]] (case 1 and 2)
normT  = sqrt( abs( 2*(R11^2+R21^2+R31^2+2*(R41^2+R51^2+R61^2)) ...
                    - 2*(R11+R21+R31)^2 ) );
normT2 = sqrt( abs( 2*(R12^2+R22^2+R32^2+2*(R42^2+R52^2+R62^2)) ...
                    - 2*(R12+R22+R32)^2 ) );

Rt = 0; Rv = 0; Rt2 = 0; Rv2 = 0;
for i=1:nboun2 % boundary1 and boundary 2 are supposed to be the same
   bonod = boundary2(i,:); exno = extnorm2(i,:)';

   no1 = bonod(2); no2 = bonod(3); no3 = bonod(4);
   x1 = nodes2(no1,1); y1 = nodes2(no1,2); z1 = nodes2(no1,3);
   x2 = nodes2(no2,1); y2 = nodes2(no2,2); z2 = nodes2(no2,3);
   x3 = nodes2(no3,1); y3 = nodes2(no3,2); z3 = nodes2(no3,3);
   vecprod = [ (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1),...
               -(x2-x1)*(z3-z1) + (x3-x1)*(z2-z1),...
               (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)];
   S = .5*norm(vecprod);
   
   if order==1
      Ng = 2;
   elseif order==2
      Ng  = 3; 
      no4 = bonod(5); no5 = bonod(6); no6 = bonod(7);
      x4  = nodes2(no4,1); y4 = nodes2(no4,2); z4 = nodes2(no4,3);
      x5  = nodes2(no5,1); y5 = nodes2(no5,2); z5 = nodes2(no5,3);
      x6  = nodes2(no6,1); y6 = nodes2(no6,2); z6 = nodes2(no6,3);
   end
   [ Xg, Wg ] = gaussPt( Ng );
             
   for j=1:size(Wg,1)
      xg = Xg(j,:); wg = Wg(j);
      
      % Interpolations
      if order==1
         uer1 = transpose( (1-xg(1)-xg(2))*urr1(i,1:3) + ... % [ux;uy] on the
                           xg(1)*urr1(i,4:6) + ...           % Gauss point
                           xg(2)*urr1(i,7:9) );        
         uer2 = transpose( (1-xg(1)-xg(2))*urr2(i,1:3) + ...
                           xg(1)*urr2(i,4:6) + ...
                           xg(2)*urr2(i,7:9) );
         xgr  = (1-xg(1)-xg(2))*[x1;y1;z1] + xg(1)*[x2;y2;z2] + xg(2)*[x3;y3;z3] ; ... % abscissae

      elseif order==2
         uer1 = transpose( -(1-xg(1)-xg(2))*(1-2*(1-xg(1)-xg(2)))*urr1(i,1:3) + ...
                -xg(1)*(1-2*xg(1))*urr1(i,4:6) + ...           % [ux;uy] on the
                -xg(2)*(1-2*xg(2))*urr1(i,7:9) + ...          % Gauss point
                4*xg(1)*(1-xg(1)-xg(2))*urr1(i,10:12) + ...
                4*xg(1)*xg(2)*urr1(i,13:15) + ...
                4*xg(2)*(1-xg(1)-xg(2))*urr1(i,16:18) );                      
         uer2 = transpose( -(1-xg(1)-xg(2))*(1-2*(1-xg(1)-xg(2)))*urr2(i,1:3) + ...
                -xg(1)*(1-2*xg(1))*urr2(i,4:6) + ...
                -xg(2)*(1-2*xg(2))*urr2(i,7:9) + ...
                4*xg(1)*(1-xg(1)-xg(2))*urr2(i,10:12) + ...
                4*xg(1)*xg(2)*urr2(i,13:15) + ...
                4*xg(2)*(1-xg(1)-xg(2))*urr2(i,16:18) );  
         xgr  = ( -(1-xg(1)-xg(2))*(1-2*(1-xg(1)-xg(2)))*[x1;y1;z1] + ...
                -xg(1)*(1-2*xg(1))*[x2;y2;z2] + ...
                -xg(2)*(1-2*xg(2))*[x3;y3;z3] + ...
                4*xg(1)*(1-xg(1)-xg(2))*[x4;y4;z4] + ...
                4*xg(1)*xg(2)*[x5;y5;z5] + ...
                4*xg(2)*(1-xg(1)-xg(2))*[x6;y6;z6] ); 
      end

      % Reference force from the BC's
      if exno(1) == 1
         fer1 = [0;0;0]; fer2 = [fscalar;0;0];
      elseif exno(1) == -1
         fer1 = [0;0;0]; fer2 = -[fscalar;0;0];
      elseif exno(2) == 1
         fer1 = [0;0;0]; fer2 = [0;0;0];
      elseif exno(2) == -1
         fer1 = [0;0;0]; fer2 = [0;0;0];
      elseif exno(3) == 1
         fer1 = [0;0;fscalar]; fer2 = [0;0;0];
      elseif exno(3) == -1
         fer1 = -[0;0;fscalar]; fer2 = [0;0;0];
      end

      ixigrec = Q'*[xgr(1);xgr(2);xgr(3)];
      X = ixigrec(1); Y = ixigrec(2); Z = ixigrec(3)-K;
      
       % Test fields
      sloct = [-X,0,Z;0,0,0;Z,0,0];
      slocv = [0,0,0;0,-Y,Z;0,Z,0];
      st = Q*sloct*Q'; sv = Q*slocv*Q';
      
      ft = st*exno; fv = sv*exno;

      vloct = [ -X^2/(2*E) - nu*Y^2/(2*E) + (2+nu)*Z^2/(2*E) ;...
                nu*X*Y/E ;...
                nu*X*Z/E ];
      vlocv = [ nu*X*Y/E ;...
                -Y^2/(2*E) - nu*X^2/(2*E) + (2+nu)*Z^2/(2*E) ;...
                nu*Y*Z/E ];
      vt = Q*vloct; vv = Q*vlocv;

      Rt  = Rt + S * wg * ( fer1'*vt - ft'*uer1 ); % increment the integral
      Rv  = Rv + S * wg * ( fer1'*vv - fv'*uer1 );
      Rt2 = Rt2 + S * wg * ( fer2'*vt - ft'*uer2 );
      Rv2 = Rv2 + S * wg * ( fer2'*vv - fv'*uer2 );
   end
end

Cte  = -sqrt(Rt^2+Rv^2)/normT - K; % K was chosen so that Cte+K is negative
Cte2 = -sqrt(Rt2^2+Rv2^2)/normT2 - K;

% Reference : we know that the point P belongs to the plane.
Pt = centCrack; QPt = Q'*Pt; CteR = -QPt(3);

%% And now, the crack itself
% Compute L ( (not so) provisionnal formula assuming that we have a particular case)
b = max(nodes2(:,2)) - min(nodes2(:,2)) ;
bt = max(nodes2(:,1)) - min(nodes2(:,1));
a = t(3)/t(1)*bt; Lx = sqrt(a^2+bt^2); Ly = b; Lx = Lx/2; % Interval is [-1;1]

disp([ 'Plane found ', num2str(toc) ]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Polynomial interpolation
if usepolys == 1
   tic
   
   % Build the xyz vectors in the crack's base
   xyz = nodes2';
   Xyz = Q'*xyz; Xs = Xyz(1,:)/Lx; Ys = Xyz(2,:)/Lx; Zs = (Xyz(3,:)+Cte)/Lx;
   X0 = (max(Xs)+min(Xs))/2; Y0 = (max(Ys)+min(Ys))/2;
   Xso = Xs-X0; Yso = Ys-Y0;
   
   nodes2b = transpose(Q*nodes2'); %[Xso',Yso',Zs'];
   [K2b,~,~,~,~] = Krig3D (nodes2b,elements2,mat,order2,boundary2,[]);
   Kinter2b = K2b( 1:3*nnodes2, 1:3*nnodes2 );
   
   % First, determine the coefficients
   lam = nu*E/((1+nu)*(1-2*nu)) ;
   mu = E/(2*(1+nu)) ;
   
   Kr = eye(ncoef);
    
%     % Physical rotated stiffness matrix
%     nodes2b = transpose(Q*nodes2');
%     [K2b,~,~,~,~] = Krig3D (nodes2b,elements2,mat,order2,boundary2,[]);
%     Kinter2b = K2b( 1:3*nnodes2, 1:3*nnodes2 );
%     Kr = zeros(ncoef);
%     
%     % Compute the diagonal energy operator
%     for ii=0:nmax
%         for jj=0:nmax
%             for pp=0:nmax
%                 index = (nmax+1)^2*ii + (nmax+1)*jj + pp+1;
%                 indey = index + (nmax+1)^3;
%                 indez = index + 2*(nmax+1)^3;
%                 u = Xso.^ii.*Yso.^jj.*Zs.^pp;
%                 uxx = reshape([u;u-u;u-u],1,[])'; % Because size(K) = 3*nnodes
%                 uyy = reshape([u-u;u;u-u],1,[])';
%                 uzz = reshape([u-u;u-u;u],1,[])';
%                 fx = Kinter2b*uxx; fy = Kinter2b*uyy; fz = Kinter2b*uzz;
%                 uxx(~indexbound2) = 0; fx(~indexbound2) = 0; 
%                 uyy(~indexbound2) = 0; fy(~indexbound2) = 0; 
%                 uzz(~indexbound2) = 0; fz(~indexbound2) = 0; 
%                 Kr(index,index) = uxx'*fx; Kr(indey,indey) = uyy'*fy; 
%                 Kr(indez,indez) = uzz'*fz;
%             end
%         end
%     end

%     Krd = diag(Kr); % debug stuff
    
   kmax = ordp;% Simplicity
    
   % Rhs : find (sigma.ez).ez
   n0 = 3*(nmax+1)^3+2*(nmax+1)^2; % Basic offset
   b = zeros(neq,(kmax+1)^2);
   num = 1;
   for k=0:kmax
       for l=0:kmax
          b(n0 + (nmax+1)*k + l+1,num) = 1;
          num = num+1;
       end
   end
    
   % Pre-process M and b to remove the 0=0 equations.
   toremove = [];
   for k=1:size(M,1)
      if norm(M(k,:)) == 0
         toremove(end+1) = k;
      end
   end
   M(toremove,:) = []; b(toremove,:) = [];
    
   % Find the coefficients thanks to the Kernel.
   %R = null(full(M));
   % First, find the kernel of M
   epsi=1.e-14; % criterion for the kernel
   [Ll, Uu, Pp, Qq, Rr] = lu (M);  % P * (R \ M) * Q = L * U
   if (norm(diag(Ll)-1,'inf')~=0), warning('diag L has 0 values'); end
   if (size(Uu,1)~=size(Uu,2)), warning('U is not squared'); end
   % Find the kernel of U 
   z1 = find(abs(diag(Uu))<=epsi); 
   z1p = setdiff([1:size(Uu,1)]',z1);
   % z1p id invertible, 
   % Schur complement on z1
   UU = Uu(z1p,z1p)\Uu(z1p,z1);
   U2 = Uu(z1,z1)-Uu(z1,z1p)*UU;
   N2 = null(full(U2));
   Rk  = zeros(size(Uu,1),size(N2,2));
   Rk(z1,:) = N2;
   Rk(z1p,:) = -UU*N2;
   Rk = Qq*Rk; % operate the permutations
   
   c0 = M\b;
   RTKR = Rk'*Kr*Rk;
   alpha = -RTKR\Rk'*Kr*c0;
   coef = c0 + Rk*alpha;
    
   % Double Kernel
%    Mm = M(1:3*(nmax+1),:);
%    Mc = M(3*(nmax+1)+1:end,:);
%    bc = b(3*(nmax+1)+1:end,:);
%    N = null(Mm);
%    CN = Mc*N;
%    R = null(CN);
%    al0 = CN\bc;
%    RTNTKNR = R'*N'*K*N*R;
%    alpha = -RTNTKNR\R'*N'*K*N*al0;
%    al = al0 + R*alpha;
%    coef = N*al;

    % Check the residual
    rescoefi = zeros(size(b,2),1) ;
    for k=1:size(b,2)
       rescoefi(k) = norm(b(:,k)-M*coef(:,k))^2 / norm(b(:,k))^2;
    end
    rescoef = norm(rescoefi);
    if rescoef > 1e-10
        disp(['Residual : ',num2str(rescoef)]);
    end

   disp([ 'Coefficients found ', num2str(toc) ]);
   tic
    
   % compute the RG
   Rhs  = zeros((ordp+1)^2,1); % Rhs
   Rhse = zeros(3*(nmax+1)^3,1); % Elementary terms of the rhs wrt basis polynomials
   vpa  = zeros(3*nnodes2, 1);
   
   for i=1:nboun2 % boundary1 and boundary 2 are supposed to be the same
      bonod = boundary2(i,:); exno = extnorm2(i,:)';

      no1 = bonod(2); no2 = bonod(3); no3 = bonod(4);
      x1 = nodes2(no1,1); y1 = nodes2(no1,2); z1 = nodes2(no1,3);
      x2 = nodes2(no2,1); y2 = nodes2(no2,2); z2 = nodes2(no2,3);
      x3 = nodes2(no3,1); y3 = nodes2(no3,2); z3 = nodes2(no3,3);
      vecprod = [ (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1),...
                  -(x2-x1)*(z3-z1) + (x3-x1)*(z2-z1),...
                  (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)];
      S = .5*norm(vecprod);

      if order==1
         Ng = min( 12, max(1, ceil((2*ordp+2)/2)) );
      elseif order==2
         Ng  = min( 12, max(1, ceil((2*ordp+3)/2)) ); 
         no4 = bonod(5); no5 = bonod(6); no6 = bonod(7);
         x4  = nodes2(no4,1); y4 = nodes2(no4,2); z4 = nodes2(no4,3);
         x5  = nodes2(no5,1); y5 = nodes2(no5,2); z5 = nodes2(no5,3);
         x6  = nodes2(no6,1); y6 = nodes2(no6,2); z6 = nodes2(no6,3);
      end
      [ Xg, Wg ] = gaussPt( Ng );

      for j=1:size(Wg,1)
         xg = Xg(j,:); wg = Wg(j);

         % Interpolations
         if order==1
            uer1 = transpose( (1-xg(1)-xg(2))*urr1(i,1:3) + ... % [ux;uy] on the
                              xg(1)*urr1(i,4:6) + ...           % Gauss point
                              xg(2)*urr1(i,7:9) );        
            xgr  = (1-xg(1)-xg(2))*[x1;y1;z1] + xg(1)*[x2;y2;z2] + xg(2)*[x3;y3;z3] ; ... % abscissae

         elseif order==2
            uer1 = transpose( -(1-xg(1)-xg(2))*(1-2*(1-xg(1)-xg(2)))*urr1(i,1:3) + ...
                   -xg(1)*(1-2*xg(1))*urr1(i,4:6) + ...           % [ux;uy] on the
                   -xg(2)*(1-2*xg(2))*urr1(i,7:9) + ...          % Gauss point
                   4*xg(1)*(1-xg(1)-xg(2))*urr1(i,10:12) + ...
                   4*xg(1)*xg(2)*urr1(i,13:15) + ...
                   4*xg(2)*(1-xg(1)-xg(2))*urr1(i,16:18) );                       
           xgr  = ( -(1-xg(1)-xg(2))*(1-2*(1-xg(1)-xg(2)))*[x1;y1;z1] + ...
                   -xg(1)*(1-2*xg(1))*[x2;y2;z2] + ...
                   -xg(2)*(1-2*xg(2))*[x3;y3;z3] + ...
                   4*xg(1)*(1-xg(1)-xg(2))*[x4;y4;z4] + ...
                   4*xg(1)*xg(2)*[x5;y5;z5] + ...
                   4*xg(2)*(1-xg(1)-xg(2))*[x6;y6;z6] ); 
         end

         % Reference force from the BC's
         if exno(1) == 1
            fer1 = [0;0;0];
         elseif exno(1) == -1
            fer1 = [0;0;0];
         elseif exno(2) == 1
            fer1 = [0;0;0];
         elseif exno(2) == -1
            fer1 = [0;0;0];
         elseif exno(3) == 1
            fer1 = [0;0;fscalar];
         elseif exno(3) == -1
            fer1 = -[0;0;fscalar];
         end

         ixigrec = Q'*[xgr(1);xgr(2);xgr(3)];
         X = (ixigrec(1)/Lx-X0); Y = (ixigrec(2)/Lx-Y0); Z = (ixigrec(3)+Cte)/Lx;

         for ii=0:nmax
             for jj=0:nmax
                 for pp=0:nmax
                     s11a = (lam+2*mu)*ii*X^(ii-1)*Y^jj*Z^pp;
                     s22a = lam*ii*X^(ii-1)*Y^jj*Z^pp;
                     s33a = lam*ii*X^(ii-1)*Y^jj*Z^pp;
                     s12a = mu*jj*X^ii*Y^(jj-1)*Z^pp;
                     s13a = mu*pp*X^ii*Y^jj*Z^(pp-1);
                     s23a = 0;
                           
                     s11b = lam*jj*X^ii*Y^(jj-1)*Z^pp;
                     s22b = (lam+2*mu)*jj*X^ii*Y^(jj-1)*Z^pp;
                     s33b = lam*jj*X^ii*Y^(jj-1)*Z^pp;
                     s12b = mu*ii*X^(ii-1)*Y^jj*Z^pp;
                     s13b = 0;
                     s23b = mu*pp*X^ii*Y^jj*Z^(pp-1);
                           
                     s11c = lam*pp*X^ii*Y^jj*Z^(pp-1);
                     s22c = lam*pp*X^ii*Y^jj*Z^(pp-1);
                     s33c = (lam+2*mu)*pp*X^ii*Y^jj*Z^(pp-1);
                     s12c = 0;
                     s13c = mu*ii*X^(ii-1)*Y^jj*Z^pp;
                     s23c = mu*jj*X^ii*Y^(jj-1)*Z^pp;

                     v1a = X^ii*Y^jj*Z^pp; v2a = 0; v3a = 0;
                     v1b = 0; v2b = X^ii*Y^jj*Z^pp; v3b = 0;
                     v1c = 0; v2c = 0; v3c = X^ii*Y^jj*Z^pp;
                     

                     % In this version, the field v is multiplied by Lx
                     slocp = [s11a,s12a,s13a;s12a,s22a,s23a;s13a,s23a,s33a];
                     sp = Q*slocp*Q';
                     fpa = sp*exno;
                     vlocp = Lx*[ v1a ; v2a ; v3a ];
                     vpa = Q*vlocp;

                     slocp = [s11b,s12b,s13b;s12b,s22b,s23b;s13b,s23b,s33b];
                     sp = Q*slocp*Q';
                     fpb = sp*exno;
                     vlocp = Lx*[ v1b ; v2b ; v3b ];
                     vpb = Q*vlocp;

                     slocp = [s11c,s12c,s13c;s12c,s22c,s23c;s13c,s23c,s33c];
                     sp = Q*slocp*Q';
                     fpc = sp*exno;
                     vlocp = Lx*[ v1c ; v2c ; v3c ];
                     vpc = Q*vlocp;

                     % increment the integrals
                     Rhse( (nmax+1)^2*ii + (nmax+1)*jj + pp+1 ) = ...
                         Rhse( (nmax+1)^2*ii + (nmax+1)*jj + pp+1 ) + ...
                           S * wg * ( fer1'*vpa - fpa'*uer1 ); 
                     Rhse( (nmax+1)^3 + (nmax+1)^2*ii + (nmax+1)*jj + pp+1 ) = ...
                         Rhse( (nmax+1)^3 + (nmax+1)^2*ii + (nmax+1)*jj + pp+1 ) + ...
                           S * wg * ( fer1'*vpb - fpb'*uer1 );
                     Rhse( 2*(nmax+1)^3 + (nmax+1)^2*ii + (nmax+1)*jj + pp+1 ) = ...
                         Rhse( 2*(nmax+1)^3 + (nmax+1)^2*ii + (nmax+1)*jj + pp+1 ) + ...
                           S * wg * ( fer1'*vpc - fpc'*uer1 );
                     
                 end
             end
         end

      end
   end

         
   for k = 0:ordp
      for l = 0:ordp
         
         coefabc = coef(:,(kmax+1)*k+l+1);
         coefa = coefabc(1:(nmax+1)^3);
         coefb = coefabc((nmax+1)^3+1:2*(nmax+1)^3);
         coefc = coefabc(2*(nmax+1)^3+1:end);
         
         Rhs(1+l+(ordp+1)*k) = 0;
         
         for ii=0:nmax
             for jj=0:nmax
                 for pp=0:nmax
                     aijp = coefa( (nmax+1)^2*ii + (nmax+1)*jj + pp+1 );
                     bijp = coefb( (nmax+1)^2*ii + (nmax+1)*jj + pp+1 );
                     cijp = coefc( (nmax+1)^2*ii + (nmax+1)*jj + pp+1 );
                     
                     Rhs(1+l+(ordp+1)*k) = Rhs(1+l+(ordp+1)*k) + ...
                         aijp * Rhse( (nmax+1)^2*ii + (nmax+1)*jj + pp+1 ) + ...
                         bijp * Rhse( (nmax+1)^3 + (nmax+1)^2*ii + (nmax+1)*jj + pp+1 ) + ...
                         cijp * Rhse( 2*(nmax+1)^3 + (nmax+1)^2*ii + (nmax+1)*jj + pp+1 ) ;
                         
                 end
             end
         end

      end
   end

%% Old way
%    for k = 0:ordp
%       for l = 0:ordp
%          
%          coefabc = coef(:,(kmax+1)*k+l+1);
%          coefa = coefabc(1:(nmax+1)^3);
%          coefb = coefabc((nmax+1)^3+1:2*(nmax+1)^3);
%          coefc = coefabc(2*(nmax+1)^3+1:end);
%          
%          Rhs(1+l+(ordp+1)*k) = 0;
%          for i=1:nboun2 % boundary1 and boundary 2 are supposed to be the same
%             bonod = boundary2(i,:); exno = extnorm2(i,:)';
%          
%             no1 = bonod(2); no2 = bonod(3); no3 = bonod(4);
%             x1 = nodes2(no1,1); y1 = nodes2(no1,2); z1 = nodes2(no1,3);
%             x2 = nodes2(no2,1); y2 = nodes2(no2,2); z2 = nodes2(no2,3);
%             x3 = nodes2(no3,1); y3 = nodes2(no3,2); z3 = nodes2(no3,3);
%             vecprod = [ (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1),...
%                         -(x2-x1)*(z3-z1) + (x3-x1)*(z2-z1),...
%                         (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)];
%             S = .5*norm(vecprod);
%             
%             if order==1
%                Ng = max(1,ceil((ordp)/2)+1);
%             elseif order==2
%                Ng  = max(1,ceil((ordp+1)/2)+1); 
%                no4 = bonod(5); no5 = bonod(6); no6 = bonod(7);
%                x4  = nodes2(no4,1); y4 = nodes2(no4,2); z4 = nodes2(no4,3);
%                x5  = nodes2(no5,1); y5 = nodes2(no5,2); z5 = nodes2(no5,3);
%                x6  = nodes2(no6,1); y6 = nodes2(no6,2); z6 = nodes2(no6,3);
%             end
%             [ Xg, Wg ] = gaussPt( Ng );
%                       
%             for j=1:size(Wg,1)
%                xg = Xg(j,:); wg = Wg(j);
%                
%                % Interpolations
%                if order==1
%                   uer1 = transpose( (1-xg(1)-xg(2))*urr1(i,1:3) + ... % [ux;uy] on the
%                                     xg(1)*urr1(i,4:6) + ...           % Gauss point
%                                     xg(2)*urr1(i,7:9) );        
%                   xgr  = (1-xg(1)-xg(2))*[x1;y1;z1] + xg(1)*[x2;y2;z2] + xg(2)*[x3;y3;z3] ; ... % abscissae
%          
%                elseif order==2
%                   uer1 = transpose( -(1-xg(1)-xg(2))*(1-2*(1-xg(1)-xg(2)))*urr1(i,1:3) + ...
%                          -xg(1)*(1-2*xg(1))*urr1(i,4:6) + ...           % [ux;uy] on the
%                          -xg(2)*(1-2*xg(2))*urr1(i,7:9) + ...          % Gauss point
%                          4*xg(1)*(1-xg(1)-xg(2))*urr1(i,10:12) + ...
%                          4*xg(1)*xg(2)*urr1(i,13:15) + ...
%                          4*xg(2)*(1-xg(1)-xg(2))*urr1(i,16:18) );                       
%                   xgr  = ( -(1-xg(1)-xg(2))*(1-2*(1-xg(1)-xg(2)))*[x1;y1;z1] + ...
%                          -xg(1)*(1-2*xg(1))*[x2;y2;z2] + ...
%                          -xg(2)*(1-2*xg(2))*[x3;y3;z3] + ...
%                          4*xg(1)*(1-xg(1)-xg(2))*[x4;y4;z4] + ...
%                          4*xg(1)*xg(2)*[x5;y5;z5] + ...
%                          4*xg(2)*(1-xg(1)-xg(2))*[x6;y6;z6] ); 
%                end
%          
%                % Reference force from the BC's
%                if exno(1) == 1
%                   fer1 = [0;0;0];
%                elseif exno(1) == -1
%                   fer1 = [0;0;0];
%                elseif exno(2) == 1
%                   fer1 = [0;0;0];
%                elseif exno(2) == -1
%                   fer1 = [0;0;0];
%                elseif exno(3) == 1
%                   fer1 = [0;0;fscalar];
%                elseif exno(3) == -1
%                   fer1 = -[0;0;fscalar];
%                end
%          
%                ixigrec = Q'*[xgr(1);xgr(2);xgr(3)];
%                X = ixigrec(1)/Lx-X0; Y = ixigrec(2)/Lx-Y0; Z = (ixigrec(3)+Cte)/Lx;
%                
%                 % Test fields
%                s11 = 0; s22 = 0; s33 = 0; s12 = 0; s13 = 0; s23 = 0;
%                v1 = 0; v2 = 0; v3 = 0;
% 
%                for ii=0:nmax
%                    for jj=0:nmax
%                        for pp=0:nmax
%                            aijp = coefa( (nmax+1)^2*ii + (nmax+1)*jj + pp+1 );
%                            bijp = coefb( (nmax+1)^2*ii + (nmax+1)*jj + pp+1 );
%                            cijp = coefc( (nmax+1)^2*ii + (nmax+1)*jj + pp+1 );
%                            s11 = s11 + (lam+2*mu)*aijp*ii*X^(ii-1)*Y^jj*Z^pp ...
%                                       + lam*bijp*jj*X^ii*Y^(jj-1)*Z^pp ...
%                                       + lam*cijp*pp*X^ii*Y^jj*Z^(pp-1);
%                            s22 = s22 + lam*aijp*ii*X^(ii-1)*Y^jj*Z^pp ...
%                                      + (lam+2*mu)*bijp*jj*X^ii*Y^(jj-1)*Z^pp ...
%                                      + lam*cijp*pp*X^ii*Y^jj*Z^(pp-1);
%                            s33 = s33 + lam*aijp*ii*X^(ii-1)*Y^jj*Z^pp ...
%                                      + lam*bijp*jj*X^ii*Y^(jj-1)*Z^pp ...
%                                      + (lam+2*mu)*cijp*pp*X^ii*Y^jj*Z^(pp-1);
%                            s12 = s12 + mu*aijp*jj*X^ii*Y^(jj-1)*Z^pp ...
%                                      + mu*bijp*ii*X^(ii-1)*Y^jj*Z^pp;
%                            s13 = s13 + mu*aijp*pp*X^ii*Y^jj*Z^(pp-1) ...
%                                      + mu*cijp*ii*X^(ii-1)*Y^jj*Z^pp;
%                            s23 = s23 + mu*bijp*pp*X^ii*Y^jj*Z^(pp-1) ...
%                                      + mu*cijp*jj*X^ii*Y^(jj-1)*Z^pp;
%                                  
%                            v1 = v1 + aijp*X^ii*Y^jj*Z^pp;
%                            v2 = v2 + bijp*X^ii*Y^jj*Z^pp;
%                            v3 = v3 + cijp*X^ii*Y^jj*Z^pp;
%                        end
%                    end
%                end
%                
%                % In this version, the field v is multiplied by Lx
%                slocp = [s11,s12,s13;s12,s22,s23;s13,s23,s33];
%                sp = Q*slocp*Q';
%                fp = sp*exno;
% 
%                vlocp = Lx*[ v1 ; v2 ; v3 ];
%                vp = Q*vlocp;
%                
%                Rhs(1+l+(ordp+1)*k) = Rhs(1+l+(ordp+1)*k) + ...
%                      S * wg * ( fer1'*vp - fp'*uer1 ); % increment the integral
%             end
%          end
% 
%       end
%    end

   %% Build the Lhs : /!\ you must have the same odd numerotation as previously
   L1x = min(Xs)-X0; L2x = max(Xs)-X0;
   L1y = min(Ys)-Y0; L2y = max(Ys)-Y0;
   Lhs = zeros( (ordp+1)^2 );
   for i=0:ordp
      for j=0:ordp
         for k=0:ordp
            for l=0:ordp
               ordx = i+k+1;
               ordy = j+l+1;
               Lhs(j+1+(ordp+1)*i,l+1+(ordp+1)*k) = ...
                    Lx^2*(L2x^ordx - L1x^ordx)/ordx * (L2y^ordy - L1y^ordy)/ordy;
                    % Lx* beacuse there is a variable change x' = Lx*x and y'=Lx*y
            end
         end
      end
   end
   
   McCoef = -Lhs\Rhs;  % - in order to have positive gap on the crack (sign is arbitrary)
   
   Xs = Xs*Lx; Ys = Ys*Lx; % use the standard basis (reverse homotetical dilatation)
   % plot the identified normal gap
   nxs = (max(Xs)-min(Xs))/100; nys = (max(Ys)-min(Ys))/100;
   X = min(Xs):nxs:max(Xs); Y = min(Ys):nys:max(Ys); 
   solup = zeros(101,101);
   for k=0:ordp
      for l=0:ordp
         solup = solup + McCoef(1+l+(ordp+1)*k) .* (X/Lx-X0)'.^k * (Y/Lx-Y0).^l;
      end
   end
   solup = solup'; % prepare for plot
   % Center of the Circle
%    Cc = centCrack; Cc = Q'*Cc; Rad = 2; zed = max(max(solup));
   
   figure;
   hold on;
   surf(X,Y,solup);
   shading interp;
   colorbar();
%    teta = 0:.1:2*pi+.1; ixe = Rad*cos(teta)+Cc(1); igrec = Rad*sin(teta)+Cc(2);
%    plot3( ixe, igrec, zed*ones(1,size(ixe,2)) , ...
%                                   'Color', 'black',  'LineWidth', 3 );
%   drawCircle ( Cc(1), Cc(2), 2, 'Color', 'black', 'LineWidth', 3 );
   axis('equal');
   
   csvwrite('fields/rg3d_poly.csv',solup);
   csvwrite('fields/rg3d_X.csv',X);
   csvwrite('fields/rg3d_Y.csv',Y);
   
   % "Zoom"
   figure;
   hold on;
   surf(X(7:end-6),Y(7:end-6),solup(7:end-6,7:end-6)); % Again
   shading interp;
   colorbar();
%    teta = 0:.1:2*pi+.1; ixe = Rad*cos(teta)+Cc(1); igrec = Rad*sin(teta)+Cc(2);
%    plot3( ixe, igrec, zed*ones(1,size(ixe,2)) , ...
%                                   'Color', 'black',  'LineWidth', 3 );
   axis('equal');
   
   disp([ 'Polynomial method ', num2str(toc) ]);
   
   if comperror == 1
      Xp = X; Yp = Y;
      X = X'*ones(size(Yp)); Y = Y'*ones(size(Xp));
      X = transpose(reshape(X',[],1)); Y = transpose(reshape(Y,[],1));
      
      Z1 = (-CteR-1e-8)*ones(size(X)); Z2 = (-CteR+1e-8)*ones(size(X));
      XYZ1 = Q*[X;Y;Z1]; % Use the physical base for abscissa
      XYZ2 = Q*[X;Y;Z2];
      Uxyz = transpose( Q'*[ux,uy,uz]' ); % Use the normal base for U
      Uxyz = reshape(Uxyz',[],1);% Re-stick the components together

      uplo = passMesh3D(nodes, elements, [XYZ1';XYZ2'], [], Uxyz);
      uplo = uplo(3*size(XYZ1,2)+1:end)-uplo(1:3*size(XYZ1,2)); % Compute the gap
      
      solup2 = reshape(solup,[],1);
      poly_error = norm( solup2-uplo(3:3:end) ) / norm(uplo(3:3:end)); % the mesh is regular
      
      % DEBUG plots
     uplo1 = reshape(uplo(3:3:end),[],size(Yp,2));
     figure;
     hold on;
     surf(Xp,Yp,uplo1);
     shading interp;
     colorbar();
%      zed = max(max(uplo));
%      plot3( ixe, igrec, zed*ones(1,size(ixe,2)) , ...
%                               'Color', 'black',  'LineWidth', 3 );
     axis('equal');
%      
%      udiff = solup2-uplo(3:3:end); udiff = reshape(udiff,[],size(Yp,2));
%      figure;
%      hold on;
%      surf(Xp,Yp,udiff);
%      shading interp;
%      colorbar();
%      zed = max(max(uplo));
%      plot3( ixe, igrec, zed*ones(1,size(ixe,2)) , ...
%                               'Color', 'black',  'LineWidth', 3 );
%      axis('equal');
      
   end
   
end

% Plot on the line X = 4
figure;
hold on;
nys = (max(Ys)-min(Ys))/100;
Y = min(Ys):nys:max(Ys); X = 4;

if usepolys == 1
   solup = 0;
   for k=0:ordp
      for l=0:ordp
         solup = solup + McCoef(1+l+(ordp+1)*k) .* (X/Lx-X0)'.^k * (Y/Lx-Y0).^l;
      end
   end
   
   plot( Y, solup, 'Color', 'black' );
end

% And the reference
if plotref == 1
   X = X*ones(size(Y));
   Z1 = (-CteR-1e-8)*ones(size(X)); Z2 = (-CteR+1e-8)*ones(size(X));
   XYZ1 = Q*[X;Y;Z1]; % Use the physical base for abscissa
   XYZ2 = Q*[X;Y;Z2];

   Uxyz = transpose( Q'*[ux,uy,uz]' ); % Use the normal base for U
   Uxyz = reshape(Uxyz',[],1); % Re-stick the components together
   uplo = passMesh3D(nodes, elements, [XYZ1';XYZ2'], [], Uxyz);

   uplo = uplo(304:end)-uplo(1:303);%  % Compute the gap
   plot( Y, uplo(3:3:end,1), 'Color', 'red' );
   csvwrite('fields/rg3d_poly2d.csv',[Y',uplo(3:3:end,1),solup']);
   %legend
end