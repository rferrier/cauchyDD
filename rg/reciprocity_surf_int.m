%01/03/2017
% Détection de fissure 3D plane par écart à la réciprocité
% Intégrations par PG

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 210000; % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 250;    % N.mm-1 : Loading on the plate
mat = [0, E, nu];
regmu   = 0;      % Regularization parameter

nbase     = 2; % Number of Fourier basis functions
ordp      = 6; % Order of polynom
loadfield = 2; % If 0 : recompute the reference problem and re-pass mesh
               % If 2 : meshes are conformal

usefourier = 0;
usepolys   = 1;
plotref    = 0;

centCrack = [4;3;1]; % Center of the crack (for reference)

if loadfield ~= 1
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
   [ nodes,elements,ntoelem,boundary,order] = readmesh3D( 'meshes/rg3dpp/plate_c_710t10.msh' );
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
   
   % Output :
   plotGMSH3D({ux,'U_x';uy,'U_y';uz,'U_z';u1,'U_vect';sigma,'stress'}, elements, nodes, 'reference');
   disp([ 'Direct problem solved ', num2str(toc) ]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import the uncracked domain /!\ MUST BE THE SAME !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (except for the crack)
[ nodes2,elements2,ntoelem2,boundary2,order2] = readmesh3D( 'meshes/rg3dpp/plate710t10.msh' );
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
else  % loadfield == 2 : conformal mesh /!\
   ur1 = zeros( 3*nnodes2, 1 );   ur2 = zeros( 3*nnodes2, 1 );
   fr1 = zeros( 3*nnodes2, 1 );   fr2 = zeros( 3*nnodes2, 1 );
   ur1(indexbound2) = u1(indexbound);   ur2(indexbound2) = u2(indexbound);
   fr1(indexbound2) = f1(indexbound);   fr2(indexbound2) = f2(indexbound);
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

%% Compute and apply v fields (for plane constraint)
%v1 = zeros(3*nnodes2, 1);
%v2 = zeros(3*nnodes2, 1);
%v3 = zeros(3*nnodes2, 1);
%v4 = zeros(3*nnodes2, 1);
%v5 = zeros(3*nnodes2, 1);
%v6 = zeros(3*nnodes2, 1);
%
%for i=1:nnodes2
%   x = nodes2(i,1);
%   y = nodes2(i,2);
%   z = nodes2(i,3);
%   
%   v1(3*i-2) = 1/E*x;
%   v2(3*i-2) = -nu/E*x;
%   v3(3*i-2) = -nu/E*x;
%   v1(3*i-1) = -nu/E*y;
%   v2(3*i-1) = 1/E*y;
%   v3(3*i-1) = -nu/E*y;
%   v1(3*i)   = -nu/E*z;
%   v2(3*i)   = -nu/E*z;
%   v3(3*i)   = 1/E*z;
%   
%   v4(3*i-2) = (1+nu)/(2*E)*y;
%   v4(3*i-1) = (1+nu)/(2*E)*x;
%   v5(3*i-2) = (1+nu)/(2*E)*z;
%   v5(3*i)   = (1+nu)/(2*E)*x;
%   v6(3*i-1) = (1+nu)/(2*E)*z;
%   v6(3*i)   = (1+nu)/(2*E)*y;
%end
%%
%f1 = Kinter2*v1;
%f2 = Kinter2*v2;
%f3 = Kinter2*v3;
%f4 = Kinter2*v4;
%f5 = Kinter2*v5;
%f6 = Kinter2*v6;

% Clean redondant stuff in indexbound2
i = 1;
while i <= size(indexbound2,1)
   if find( indexbound2(i) == indexbound2(1:i-1) )
      indexbound2(i) = [];
      i = i-1;
   end
   i = i+1;
end

%R11 = (fr1(indexbound2)'*v1(indexbound2) - f1(indexbound2)'*ur1(indexbound2));
%R21 = (fr1(indexbound2)'*v2(indexbound2) - f2(indexbound2)'*ur1(indexbound2));
%R31 = (fr1(indexbound2)'*v3(indexbound2) - f3(indexbound2)'*ur1(indexbound2));
%R41 = (fr1(indexbound2)'*v4(indexbound2) - f4(indexbound2)'*ur1(indexbound2));
%R51 = (fr1(indexbound2)'*v5(indexbound2) - f5(indexbound2)'*ur1(indexbound2));
%R61 = (fr1(indexbound2)'*v6(indexbound2) - f6(indexbound2)'*ur1(indexbound2));
%RR1e = [R11, R41, R51 ; R41, R21, R61 ; R51, R61, R31];
%
%R12 = (fr2(indexbound2)'*v1(indexbound2) - f1(indexbound2)'*ur2(indexbound2));
%R22 = (fr2(indexbound2)'*v2(indexbound2) - f2(indexbound2)'*ur2(indexbound2));
%R32 = (fr2(indexbound2)'*v3(indexbound2) - f3(indexbound2)'*ur2(indexbound2));
%R42 = (fr2(indexbound2)'*v4(indexbound2) - f4(indexbound2)'*ur2(indexbound2));
%R52 = (fr2(indexbound2)'*v5(indexbound2) - f5(indexbound2)'*ur2(indexbound2));
%R62 = (fr2(indexbound2)'*v6(indexbound2) - f6(indexbound2)'*ur2(indexbound2));
%RR2e = [R12, R42, R52 ; R42, R22, R62 ; R52, R62, R32];

% Second method : loop over boundaries to compute the integrals
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
t = [1 ; 0 ; -(normal(1))/normal(3)]; t = t/norm(t); % Total random t
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

%ft = Kinter2*vt;
%fv = Kinter2*vv;
%
%Rte  = (fr1(indexbound2)'*vt(indexbound2) - ft(indexbound2)'*ur1(indexbound2));
%Rve  = (fr1(indexbound2)'*vv(indexbound2) - fv(indexbound2)'*ur1(indexbound2));
%Rt2e = (fr2(indexbound2)'*vt(indexbound2) - ft(indexbound2)'*ur2(indexbound2));
%Rv2e = (fr2(indexbound2)'*vv(indexbound2) - fv(indexbound2)'*ur2(indexbound2));

% Alternative way
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
a = t(3)/t(1)*bt; Lx = sqrt(a^2+bt^2); Ly = b;

disp([ 'Plane found ', num2str(toc) ]);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourier sum
if usefourier == 1
   tic
   Rp       = zeros(nbase+1,1);
   Rm       = zeros(nbase+1,1);
   lambdax  = zeros(nbase+1,1);
   lambday  = zeros(nbase+1,1);
   fournpp  = zeros(nbase+1,nbase+1);
   fournpm  = zeros(nbase+1,nbase+1);
   fournmp  = zeros(nbase+1,nbase+1);
   fournmm  = zeros(nbase+1,nbase+1);
   akan     = zeros(nbase,nbase);
   bkan     = zeros(nbase,nbase);
   ckan     = zeros(nbase,nbase);
   dkan     = zeros(nbase,nbase);
   
   for kpx=2:nbase+1
      for kpy=2:nbase+1
         kx = kpx-1; ky = kpy-1;
         vp = zeros(3*nnodes2,1);
         lambdax(kpx) = 2*kx*pi/Lx; lambday(kpy) = 2*ky*pi/Ly;
         lambda = sqrt(lambdax(kpx)^2+lambday(kpy)^2);
         
         for sx = [1,-1]
            lambdaxe = sx*lambdax(kpx);
            for sy = [1,-1]
               lambdaye = sy*lambday(kpy);
%               for i=1:nnodes2
%                  x = nodes2(i,1);
%                  y = nodes2(i,2);
%                  z = nodes2(i,3);
%                  % Change base (for coordiantes)
%                  ixigrec = Q'*[x;y;z]; X = ixigrec(1); Y = ixigrec(2); Z = ixigrec(3)+Cte; %%%%%%%%%%!!!!!!!! + or - ?
%                  Xs(i) = X; Ys(i) = Y; Zs(i) = Z;
%                  
%                  v1 = -I*lambdaxe*exp(-I*lambdaxe*X-I*lambdaye*Y)* ...
%                                     ( exp(lambda*Z)+exp(-lambda*Z) );
%                  v2 = -I*lambdaye*exp(-I*lambdaxe*X-I*lambdaye*Y)* ...
%                                     ( exp(lambda*Z)+exp(-lambda*Z) );
%                  v3 = lambda*exp(-I*lambdaxe*X-I*lambdaye*Y)* ...
%                                     ( exp(lambda*Z)-exp(-lambda*Z) );
%                  vloc = [ v1 ; v2 ; v3 ];
%                  % Change base (for vector), in the other direction
%                  vxy = Q*vloc; vp(3*i-2) = vxy(1); vp(3*i-1) = vxy(2); vp(3*i) = vxy(3);
%               end
%               fp = Kinter2*vp;
%               Rpe = (fr1(indexbound2)'*vp(indexbound2) - fp(indexbound2)'*ur1(indexbound2));
           
               Rp = 0;
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
                     Ng = 8;
                  elseif order==2
                     Ng  = 8; 
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
                     X = ixigrec(1); Y = ixigrec(2); Z = ixigrec(3)+Cte;
                     
                      % Test fields
                     s11 = -E/(1+nu)*lambdaxe^2*exp(-I*lambdaxe*X-I*lambdaye*Y)* ...
                                        ( exp(lambda*Z)+exp(-lambda*Z) );
                     s12 = -E/(1+nu)*lambdaxe*lambdaye*exp(-I*lambdaxe*X-I*lambdaye*Y)* ...
                                        ( exp(lambda*Z)+exp(-lambda*Z) );
                     s13 = -I*E/(1+nu)*lambdaxe*lambda*exp(-I*lambdaxe*X-I*lambdaye*Y)* ...
                                        ( exp(lambda*Z)-exp(-lambda*Z) );
                     s22 = -E/(1+nu)*lambdaye^2*exp(-I*lambdaxe*X-I*lambdaye*Y)* ...
                                        ( exp(lambda*Z)+exp(-lambda*Z) );
                     s23 = -I*E/(1+nu)*lambdaye*lambda*exp(-I*lambdaxe*X-I*lambdaye*Y)* ...
                                        ( exp(lambda*Z)-exp(-lambda*Z) );
                     s33 = -E/(1+nu)*lambda^2*exp(-I*lambdaxe*X-I*lambdaye*Y)* ...
                                        ( exp(lambda*Z)+exp(-lambda*Z) );
                     slocp = [s11,s12,s13;s12,s22,s23;s13,s23,s33];
                     sp = Q*slocp*Q';
                     fp = sp*exno;
               
                     v1 = -I*lambdaxe*exp(-I*lambdaxe*X-I*lambdaye*Y)* ...
                                        ( exp(lambda*Z)+exp(-lambda*Z) );
                     v2 = -I*lambdaye*exp(-I*lambdaxe*X-I*lambdaye*Y)* ...
                                        ( exp(lambda*Z)+exp(-lambda*Z) );
                     v3 = lambda*exp(-I*lambdaxe*X-I*lambdaye*Y)* ...
                                        ( exp(lambda*Z)-exp(-lambda*Z) );
                     vlocp = [ v1 ; v2 ; v3 ];
                     vp = Q*vlocp;
                     
                     Rp = Rp + S * wg * ( fer1'*vp - fp'*uer1 ); % increment the integral
                  end
               end
                  
               % Fourier coefficients
               if sx==1 && sy==1
                  fournpp(kpx,kpy) = (1+nu)/(2*E*Lx*Ly*lambda^2)*Rp;
               elseif sx==1 && sy==-1
                  fournpm(kpx,kpy) = (1+nu)/(2*E*Lx*Ly*lambda^2)*Rp;
               elseif sx==-1 && sy==1
                  fournmp(kpx,kpy) = (1+nu)/(2*E*Lx*Ly*lambda^2)*Rp;
               elseif sx==-1 && sy==-1
                  fournmm(kpx,kpy) = (1+nu)/(2*E*Lx*Ly*lambda^2)*Rp;
               end

            end
         end
         % u = acoscos+bcossin+csincos+dsinsin
         akan(kx,ky)   = real(fournpp(kpx,kpy) + fournpm(kpx,kpy) +...
                              fournmp(kpx,kpy) + fournmm(kpx,kpy));
         bkan(kx,ky)   = imag(fournpp(kpx,kpy) - fournpm(kpx,kpy) +...
                              fournmp(kpx,kpy) - fournmm(kpx,kpy));
         ckan(kx,ky)   = imag(fournpp(kpx,kpy) + fournpm(kpx,kpy) -...
                              fournmp(kpx,kpy) - fournmm(kpx,kpy));
         dkan(kx,ky)   = -real(fournpp(kpx,kpy) - fournpm(kpx,kpy) -...
                               fournmp(kpx,kpy) + fournmm(kpx,kpy));
      end
   end
   
   %% DEBUG :
   %ui = reshape(imag(vp),2,[])';  ux = ui(:,1);  uy = ui(:,2);
%   plotGMSH({ Xs,'x' ; Ys,'y' ; Zs,'z' ; real(vp),'U_vect'}, elements2, nodes2, 'test field');
   
   % The constant term
   fournpp(1,1) = -(R11+R21+R31)/(Lx*Ly);
   
   % plot the identified normal gap
   nxs = (max(Xs)-min(Xs))/100; nys = (max(Ys)-min(Ys))/100;
   X = min(Xs):nxs:max(Xs); Y = min(Ys):nys:max(Ys); 
   solu = fournpp(1,1);
   for kpx=2:nbase+1
      for kpy=2:nbase+1
         solu = solu + akan(kpx-1,kpy-1)*cos(lambdax(kpx)*X')*cos(lambday(kpy)*Y)...
                     + bkan(kpx-1,kpy-1)*cos(lambdax(kpx)*X')*sin(lambday(kpy)*Y)...
                     + ckan(kpx-1,kpy-1)*sin(lambdax(kpx)*X')*cos(lambday(kpy)*Y)...
                     + dkan(kpx-1,kpy-1)*sin(lambdax(kpx)*X')*sin(lambday(kpy)*Y);
      end
   end
   solu = solu';  % Very imporant and no hack : comes from the function surf
   % Center of the Circle
   Cc = [4;3;1]; Cc = Q'*Cc; Rad = 2; zed = max(max(solu));
   
   figure;
   hold on;
   surf(X,Y,solu);
   shading interp;
   colorbar();
%   drawCircle ( Cc(1), Cc(2), 2, 'Color', 'black', 'LineWidth', 3 );
   teta = 0:.1:2*pi+.1; ixe = Rad*cos(teta)+Cc(1); igrec = Rad*sin(teta)+Cc(2);
   plot3( ixe, igrec, zed*ones(1,size(ixe,2)) , ...
                                  'Color', 'black',  'LineWidth', 3 );
   disp([ 'Fourier method ', num2str(toc) ]);
   axis('equal');
end

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
   
   coef = cell(ordp+1); % coefficients of the polynoms
   
   for k = 0:ordp
      for l = 0:ordp
         sze = floor((k+3)/2)*floor((l+2)/2) + floor((k+2)/2)*floor((l+3)/2)...
               + floor((k+2)/2)*floor((l+2)/2); % Number of unknowns
         Rhsco = zeros(sze,1);
         Lhsca = zeros( sze, floor((k+3)/2)*floor((l+2)/2) );
         Lhscb = zeros( sze, floor((k+2)/2)*floor((l+3)/2) );
         Lhscc = zeros( sze, floor((k+2)/2)*floor((l+2)/2) );
         
         neq = 1;  % index of the no of the equation
         % div sigma.x
         for i=1:floor((k-1)/2)
            for j=1:floor(l/2)
               inda = (j+1) + (floor(l/2)+1)*i; % index of aij in coef
               indb = (j+1) + (floor((l+1)/2)+1)*i; % index of bij in coef
               indc = (j+1) + (floor(l/2)+1)*i; % index of cij in coef
               
               Lhsca( neq, inda ) = (lam+2*mu)*(k+1-2*i)*(k-2*i);
               Lhscb( neq, indb ) = (lam+mu)*(l+1-2*j)*(k-2*i);
               Lhscc( neq, indc ) = (lam+mu)*(2*i+2*j+1)*(k-2*i);
               Lhsca( neq, inda-1+floor(l/2)+1 ) = ... % a_{i+1,j-1}
                                    mu*(l-2*j+2)*(l-2*j+1);
               Lhsca( neq, inda+floor(l/2)+1 ) = ... % a_{i+1,j}
                                    mu*(2*i+2*j+2)*(2*i+2*j+1);

               neq = neq+1;
            end
         end
         for j=0:floor( (l-2)/2 )
            Lhsca( neq, j+1 ) = mu*(l-2*j)*(l-2*j-1);
            Lhsca( neq, j+2 ) = mu*(2*j+2)*(2*j+1);
            neq =neq+1;
         end
         for i=0:floor((k-1)/2)
            Lhsca( neq, 1 + (floor(l/2)+1)*i ) = (lam+2*mu)*(k+1-2*i)*(k-2*i);
            Lhscb( neq, 1 + (floor((l+1)/2)+1)*i ) = (lam+mu)*(l+1)*(k-2*i);
            Lhscc( neq, 1 + (floor(l/2)+1)*i ) = (lam+mu)*(k-2*i)*(2*i+1);
            Lhsca( neq, 1 + (floor(l/2)+1)*(i+1) ) = mu*(2*i+2)*(2*i+1);
            neq = neq+1;
         end
         % div sigma.y
         for i=1:floor(k/2)
            for j=1:floor((l-1)/2)
               inda = (j+1) + (floor(l/2)+1)*i; % index of aij in coef
               indb = (j+1) + (floor((l+1)/2)+1)*i; % index of bij in coef
               indc = (j+1) + (floor(l/2)+1)*i; % index of cij in coef
               
               Lhsca( neq, inda ) = (lam+mu)*(k+1-2*i)*(l-2*j);
               Lhscb( neq, indb ) = (lam+2*mu)*(l+1-2*j)*(l-2*j);
               Lhscc( neq, indc ) = (lam+mu)*(2*i+2*j+1)*(l-2*j);
               Lhscb( neq, indb+1-floor((l+1)/2)-1 ) = ... % b_{i-1,j+1}
                                    mu*(k-2*i+2)*(k-2*i+1);
               Lhscb( neq, indb+1 ) = ... % b_{i,j+1}
                                    mu*(2*i+2*j+2)*(2*i+2*j+1);
                                    
               neq = neq+1;
            end
         end
         for i=0:floor( (k-2)/2 )
            Lhscb( neq, (floor((l+1)/2)+1)*i+1 ) = mu*(k-2*i)*(k-2*i-1); % and not *(i+1)
            Lhscb( neq, (floor((l+1)/2)+1)*(i+1)+1 ) = mu*(2*i+2)*(2*i+1);
            neq =neq+1;
         end
         for j=0:floor((l-1)/2)
            Lhsca( neq, j+1 ) = (lam+mu)*(l-2*j)*(k+1);
            Lhscb( neq, j+1 ) = (lam+2*mu)*(l+1-2*j)*(l-2*j);
            Lhscc( neq, j+1 ) = (lam+mu)*(l-2*j)*(2*j+1);
            Lhscb( neq, j+2 ) = mu*(2*j+2)*(2*j+1);
            neq = neq+1;
         end
         % div sigma.z
         for i=1:floor(k/2)
            for j=1:floor(l/2)
               inda = (j+1) + (floor(l/2)+1)*i; % index of aij in coef
               indb = (j+1) + (floor((l+1)/2)+1)*i; % index of bij in coef
               indc = (j+1) + (floor(l/2)+1)*i; % index of cij in coef
               
               Lhsca( neq, inda ) = (lam+mu)*(k+1-2*i)*(2*i+2*j);
               Lhscb( neq, indb ) = (lam+mu)*(l+1-2*j)*(2*i+2*j);
               Lhscc( neq, indc ) = (lam+2*mu)*(2*i+2*j+1)*(2*i+2*j);
               Lhscc( neq, indc-floor(l/2)-1 ) = ... % c_{i-1,j}
                                    mu*(k-2*i+2)*(k-2*i+1);
               Lhscc( neq, indc-1 ) = ... % c_{i,j-1}
                                    mu*(l-2*j+2)*(l-2*j+1);
                                    
               neq = neq+1;
            end
         end
         for j=1:floor(l/2)
            Lhsca( neq, j+1 ) = (lam+mu)*(2*j)*(k+1);
            Lhscb( neq, j+1 ) = (lam+mu)*(l+1-2*j)*(2*j);
            Lhscc( neq, j+1 ) = (lam+2*mu)*(2*j)*(2*j+1);
            Lhscc( neq, j ) = mu*(l-2*j+2)*(l-2*j+1);
            neq = neq+1;
         end
         for i=1:floor(k/2)
            Lhsca( neq, 1 + (floor(l/2)+1)*i ) = (lam+mu)*(2*i)*(k+1-2*i);
            Lhscb( neq, 1 + (floor((l+1)/2)+1)*i ) = (lam+mu)*(l+1)*(2*i);
            Lhscc( neq, 1 + (floor(l/2)+1)*i ) = (lam+2*mu)*(2*i)*(2*i+1);%c0i
            Lhscc( neq, 1 + (floor(l/2)+1)*(i-1) ) = mu*(k-2*i+2)*(k-2*i+1);%c0i-1
            neq = neq+1;
         end
         % arbitrary equations
%         for i=0:floor(k/2) % ci0 = 0
%            Lhscc( neq, 1 + (floor(l/2)+1)*i ) = E;
%            neq = neq+1;
%         end
%         for j=0:floor(l/2) % c0j = 0
%            Lhscc( neq, j+1 ) = E;   % I use E instead of 1 for the condition
%            neq = neq+1;
%         end
%         Lhscc(neq, end) = E; neq = neq+1; % ckl = 0 (greatest term in z)
%         Lhsca( neq, 1) = E; Lhscb( neq, 1) = -E; % a00 = b00
%         neq = neq+1;
%         Lhsca( neq, end) = E; Lhscb( neq, end) = -E; % akl = bkl
%         neq = neq+1;

%         for j=1:floor(l/2) % c0j = 0
%            Lhscc( neq, j+1 ) = E;
%            neq = neq+1;
%         end
%         Lhsca( neq, 1) = E; Lhscb( neq, 1) = -E; % a00 = b00
%         neq = neq+1;
         % Purpose equation
         Lhsca( neq, 1) = lam*(k+1);
         Lhscb( neq, 1) = lam*(l+1);
         Lhscc( neq, 1) = (lam+2*mu);
         Rhsco( neq ) = Lx; % Because of the homotecy
         
         Lhsco = [ Lhsca , Lhscb , Lhscc ];
         
         % It's cleaner to remove the 0=0 equations at the end
         Lhsco(neq+1:end,:) = []; Rhsco(neq+1:end,:) = [];
         
         % Solve the linear problem to find the coefficients
%         coef{ k+1, l+1 } = Lhsco\Rhsco;
         
         % Ker stuff
%         U = E*null(Lhsco);  % E* in order to have a better condition number
%         Lhsto = [ Lhsco ; U' ]; Rhsto = [ Rhsco ; zeros(size(U,2), 1) ];
%         Solpro = Lhsto\Rhsto;
%         coef{ k+1, l+1 } = Solpro(1:sze);
         
         %% Minimization problem to find the best polynom
         ka = zeros(floor((k+1)/2)*(floor(l/2)+1),1);
         kb = zeros(floor(k/2)*(floor((l+1)/2)+1),1);
         kc = zeros(floor(k/2)*(floor(l/2)+1),1);
         
         for ii = 0:floor((k+1)/2)
            for jj = 0:floor(l/2)
               u = Xso.^(k+1-2*ii).*Yso.^(l-2*jj).*Zs.^(2*ii+2*jj);
               u = reshape([u;u-u;;u-u],1,[])'; % Because size(K) = 3*nnodes
               f = Kinter2b*u;
               u(~indexbound2) = 0; f(~indexbound2) = 0; 
               ka( 1+jj+ii*(floor(l/2)+1), 1 ) =...
%                        max( abs( u ) );
%                        sum( ( u ).^2 );
                        u'*f;
            end
         end
         for ii = 0:floor(k/2)
            for jj = 0:floor((l+1)/2)
               u = Xso.^(k-2*ii).*Yso.^(l+1-2*jj).*Zs.^(2*ii+2*jj);
               u = reshape([u;u-u;;u-u],1,[])'; % Because size(K) = 3*nnodes
               f = Kinter2b*u;
               u(~indexbound2) = 0; f(~indexbound2) = 0; 
               kb( 1+jj+ii*(floor((l+1)/2)+1), 1 ) =...
%                       max( abs( u ) );
%                       sum( ( u ).^2 );
                       u'*f;

            end
         end
         for ii = 0:floor(k/2)
            for jj = 0:floor(l/2)
               u = Xso.^(k-2*ii).*Yso.^(l-2*jj).*Zs.^(2*ii+2*jj+1);
               u = reshape([u;u-u;;u-u],1,[])'; % Because size(K) = 3*nnodes
               f = Kinter2b*u;
               u(~indexbound2) = 0; f(~indexbound2) = 0; 
               kc( 1+jj+ii*(floor(l/2)+1), 1 ) =...
%                       max( abs( u ) );
%                       sum( ( u ).^2 );
                       u'*f;
            end
         end
      
         zab = zeros( size(ka,1), size(kb,1) );
         zac = zeros( size(ka,1), size(kc,1) );
         zbc = zeros( size(kb,1), size(kc,1) );
      
         % Mm consists in minimizing the worst amplification for each a, b or c (E^2 just for condition number)
         Mm = [ diag(ka),zab,zac ; zab',diag(kb),zbc ; zac',zbc',diag(kc) ]*E; %*E^2
%         Mm = [ diag(ka).^2,zab,zac ; zab',diag(kb).^2,zbc ; zac',zbc',diag(kc).^2 ];
         
         % Minimize under the constraints
         Lhsto2 = [ Mm , Lhsco' ; Lhsco , zeros(size(Lhsco,1)) ];
         Rhsto2 = [ zeros(size(Mm,1), 1) ; Rhsco ];
         Solpro2 = Lhsto2\Rhsto2;
         coef{ k+1, l+1 } = Solpro2(1:sze);
      
      end
   end

   % compute the RG
   Rhs  = zeros((ordp+1)^2,1);
   Rhse = zeros((ordp+1)^2,1);
   vpa  = zeros(3*nnodes2, 1);
   
   for k = 0:ordp
      for l = 0:ordp
         sze = floor((k+3)/2)*floor((l+2)/2) + floor((k+2)/2)*floor((l+3)/2)...
               + floor((k+2)/2)*floor((l+2)/2); % Number of unknowns
               
         coefabc = coef{ k+1, l+1 }; % extract data
         coefa   = coefabc( 1:floor((k+3)/2)*floor((l+2)/2) );
         coefb   = coefabc( floor((k+3)/2)*floor((l+2)/2) + 1 : ...
                floor((k+3)/2)*floor((l+2)/2) + floor((k+2)/2)*floor((l+3)/2) );
         coefc   = coefabc( floor((k+3)/2)*floor((l+2)/2) +...
                           floor((k+2)/2)*floor((l+3)/2) + 1 : end );
               
%         for no = 1:nnodes2
%         
%            x = nodes2(no,1);
%            y = nodes2(no,2);
%            z = nodes2(no,3);
%            % Change base (for coordiantes)
%            ixigrec = Q'*[x;y;z]; X = (ixigrec(1))/Lx; Y = ixigrec(2)/Lx;
%            Z = (ixigrec(3)+Cte)/Lx;
%            Xs(no) = X; Ys(no) = Y; Zs(no) = Z;

%            % Build the xyz vector
%            GROXa = zeros( floor((k+3)/2)*floor((l+2)/2), 1 );
%            GROXb = zeros( floor((k+2)/2)*floor((l+3)/2), 1 );
%            GROXc = zeros( floor((k+2)/2)*floor((l+2)/2), 1 );
%            
%            % That's not optimized, but it's clearer
%            for i = 0:floor((k+1)/2)
%               for j = 0:floor(l/2)
%                  GROXa( 1+j+i*(floor(l/2)+1) ) =...
%                                        X^(k+1-2*i)*Y^(l-2*j)*Z^(2*i+2*j);
%               end
%            end
%            for i = 0:floor(k/2)
%               for j = 0:floor((l+1)/2)
%                  GROXb( 1+j+i*(floor((l+1)/2)+1) ) =...
%                                        X^(k-2*i)*Y^(l+1-2*j)*Z^(2*i+2*j);
%               end
%            end
%            for i = 0:floor(k/2)
%               for j = 0:floor(l/2)
%                  GROXc( 1+j+i*(floor(l/2)+1) ) =...
%                                        X^(k-2*i)*Y^(l-2*j)*Z^(2*i+2*j+1);
%               end
%            end
%            
%            % Build the test field
%            v1 = coefa'*GROXa;
%            v2 = coefb'*GROXb;
%            v3 = coefc'*GROXc;
%            vloc = [ v1 ; v2 ; v3 ];
%            % Change base (for vector), in the other direction
%            vxy = Q*vloc; vpa(3*no-2) = vxy(1); vpa(3*no-1) = vxy(2); vpa(3*no) = vxy(3);
%         end
%         fp = Kinter2*vpa;
%         Rhse(1+l+(ordp+1)*k) = (fr1(indexbound2)'*vpa(indexbound2) - fp(indexbound2)'*ur1(indexbound2));
         
         % Compute the center of the crack plane
%         X0 = (max(Xs)+min(Xs))/2; Y0 = (max(Ys)+min(Ys))/2;
         
         Rhs(1+l+(ordp+1)*k) = 0;
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
               Ng = max(1,ceil((ordp)/2)+1);
            elseif order==2
               Ng  = max(1,ceil((ordp+1)/2)+1); 
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
               X = ixigrec(1)/Lx-X0; Y = ixigrec(2)/Lx-Y0; Z = (ixigrec(3)+Cte)/Lx;
               
                % Test fields
               s11 = 0; s22 = 0; s33 = 0; s12 = 0; s13 = 0; s23 = 0;
               for ii=0:floor(k/2)
                  for jj=0:floor(l/2)
                     s11 = s11 + ( (lam+2*mu)*(k+1-2*ii)*coefa( 1+jj+ii*(floor(l/2)+1) ) + ...
                           lam*(l+1-2*jj)*coefb( 1+jj+ii*(floor((l+1)/2)+1) ) + ...
                           lam*(2*ii+2*jj+1)*coefc( 1+jj+ii*(floor(l/2)+1) ) ) ...
                           * X^(k-2*ii)*Y^(l-2*jj)*Z^(2*ii+2*jj);
                     s22 = s22 + ( lam*(k+1-2*ii)*coefa( 1+jj+ii*(floor(l/2)+1) ) + ...
                           (lam+2*mu)*(l+1-2*jj)*coefb( 1+jj+ii*(floor((l+1)/2)+1) ) + ...
                           lam*(2*ii+2*jj+1)*coefc( 1+jj+ii*(floor(l/2)+1) ) ) ...
                           * X^(k-2*ii)*Y^(l-2*jj)*Z^(2*ii+2*jj);
                     s33 = s33 + ( lam*(k+1-2*ii)*coefa( 1+jj+ii*(floor(l/2)+1) ) + ...
                           lam*(l+1-2*jj)*coefb( 1+jj+ii*(floor((l+1)/2)+1) ) + ...
                           (lam+2*mu)*(2*ii+2*jj+1)*coefc( 1+jj+ii*(floor(l/2)+1) ) ) ...
                           * X^(k-2*ii)*Y^(l-2*jj)*Z^(2*ii+2*jj);
                     s13 = s13 + mu*(k-2*ii)*coefc( 1+jj+ii*(floor(l/2)+1) ) * X^(k-2*ii-1)*Y^(l-2*jj)*Z^(2*ii+2*jj+1);
                     s23 = s23 + mu*(l-2*jj)*coefc( 1+jj+ii*(floor(l/2)+1) ) * X^(k-2*ii)*Y^(l-1-2*jj)*Z^(2*ii+2*jj+1);
                  end
               end
               
               for ii=0:floor((k+1)/2)
                  for jj=0:floor(l/2)
                     s12 = s12 + mu*(l-2*jj)*coefa( 1+jj+ii*(floor(l/2)+1) ) * X^(k+1-2*ii)*Y^(l-2*jj-1)*Z^(2*ii+2*jj);
                  end
               end
               for ii=0:floor(k/2)
                  for jj=0:floor((l+1)/2)
                     s12 = s12 + mu*(k-2*ii)*coefb( 1+jj+ii*(floor((l+1)/2)+1) ) * X^(k-1-2*ii)*Y^(l-2*jj+1)*Z^(2*ii+2*jj);
                  end
               end
               
               jj = 0; %ii = 0; 
               %if l >= 2
                  for ii=1:floor((k+1)/2)
                     s13 = s13 + mu*(2*ii+2*jj)*coefa( 1+jj+ii*(floor(l/2)+1) ) * X^(k+1-2*ii)*Y^(l-2*jj)*Z^(2*ii+2*jj-1);
                  end
               %end
               %if l >= 1
                  for ii=1:floor(k/2)
                     s23 = s23 + mu*(2*ii+2*jj)*coefb( 1+jj+ii*(floor((l+1)/2)+1) ) * X^(k-2*ii)*Y^(l+1-2*jj)*Z^(2*ii+2*jj-1);
                  end
               %end
               ii = 0;% jj = 0;
               %if k >= 1
                  for jj=1:floor(l/2)
                     s13 = s13 + mu*(2*ii+2*jj)*coefa( 1+jj+ii*(floor(l/2)+1) ) * X^(k+1-2*ii)*Y^(l-2*jj)*Z^(2*ii+2*jj-1);
                  end
               %end
               %if k >= 2
                  for jj=1:floor((l+1)/2)
                     s23 = s23 + mu*(2*ii+2*jj)*coefb( 1+jj+ii*(floor((l+1)/2)+1) ) * X^(k-2*ii)*Y^(l+1-2*jj)*Z^(2*ii+2*jj-1);
                  end
               %end
               for ii=1:floor((k+1)/2)
                  for jj=1:floor(l/2)
                     s13 = s13 + mu*(2*ii+2*jj)*coefa( 1+jj+ii*(floor(l/2)+1) ) * X^(k+1-2*ii)*Y^(l-2*jj)*Z^(2*ii+2*jj-1);
                  end
               end
               for ii=1:floor(k/2)
                  for jj=1:floor((l+1)/2)
                     s23 = s23 + mu*(2*ii+2*jj)*coefb( 1+jj+ii*(floor((l+1)/2)+1) ) * X^(k-2*ii)*Y^(l+1-2*jj)*Z^(2*ii+2*jj-1);
                  end
               end
%               for ii=0:floor(k/2)   % Fusionned with stuff up there
%                  for jj=0:floor(l/2)
%                     s13 = s13 + mu*(k-2*ii)*coefc( 1+jj+ii*(floor(l/2)+1) ) * X^(k-2*ii-1)*Y^(l-2*jj)*Z^(2*ii+2*jj+1);
%                     s23 = s23 + mu*(l-2*jj)*coefc( 1+jj+ii*(floor(l/2)+1) ) * X^(k-2*ii)*Y^(l-1-2*jj)*Z^(2*ii+2*jj+1);
%                  end
%               end
               
               slocp = 1/Lx*[s11,s12,s13;s12,s22,s23;s13,s23,s33]; % 1/Lx because of the homotecy (fog)' = g'*f'og
               sp = Q*slocp*Q';
               fp = sp*exno;
         
               % Build the xyz vector
               GROXa = zeros( floor((k+3)/2)*floor((l+2)/2), 1 );
               GROXb = zeros( floor((k+2)/2)*floor((l+3)/2), 1 );
               GROXc = zeros( floor((k+2)/2)*floor((l+2)/2), 1 );
               
               % That's not optimized, but it's clearer
               for ii = 0:floor((k+1)/2)
                  for jj = 0:floor(l/2)
                     GROXa( 1+jj+ii*(floor(l/2)+1) ) =...
                                           X^(k+1-2*ii)*Y^(l-2*jj)*Z^(2*ii+2*jj);
                  end
               end
               for ii = 0:floor(k/2)
                  for jj = 0:floor((l+1)/2)
                     GROXb( 1+jj+ii*(floor((l+1)/2)+1) ) =...
                                           X^(k-2*ii)*Y^(l+1-2*jj)*Z^(2*ii+2*jj);
                  end
               end
               for ii = 0:floor(k/2)
                  for jj = 0:floor(l/2)
                     GROXc( 1+jj+ii*(floor(l/2)+1) ) =...
                                           X^(k-2*ii)*Y^(l-2*jj)*Z^(2*ii+2*jj+1);
                  end
               end
               
               % Build the test field
               v1 = coefa'*GROXa;
               v2 = coefb'*GROXb;
               v3 = coefc'*GROXc;
               vlocp = [ v1 ; v2 ; v3 ];
               vp = Q*vlocp;
               
               Rhs(1+l+(ordp+1)*k) = Rhs(1+l+(ordp+1)*k) + ...
                     S * wg * ( fer1'*vp - fp'*uer1 ); % increment the integral
            end
         end

      end
   end

   %% Build the Lhs : /!\ you must have the same odd numerotation as previously
   L1x = min(Xs)-X0; L2x = max(Xs)-X0;
   L1y = min(Ys)-Y0; L2y = max(Ys)-Y0;
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
   
   %% Regularisation (separated in order to do L-curves)
   for i=1:ordp
      for j=1:ordp
         for k=1:ordp
            for l=1:ordp
               ordx = i+k+1;
               ordy = j+l+1;
               Lhsr(j+1+(ordp+1)*i,l+1+(ordp+1)*k) = ...
                    i*k* (L2x^(ordx-1) - L1x^(ordx-1))/(ordx-1) *...
                                        (L2y^ordy - L1y^ordy)/ordy + ...
                    j*l* (L2x^ordx - L1x^ordx)/ordx *...
                                        (L2y^(ordy-1) - L1y^(ordy-1))/(ordy-1);
            end
         end
      end
   end
   Lhsr = Lhs-Lhs;
   for i=1:ordp
      for k=1:ordp
         j = 0; l = 0;
         ordx = i+k+1;
         ordy = j+l+1;
         Lhsr(j+1+(ordp+1)*i,l+1+(ordp+1)*k) = ...
              i*k* (L2x^(ordx-1) - L1x^(ordx-1))/(ordx-1) *...
                                  (L2y^ordy - L1y^ordy)/ordy ;
      end
   end
   
   for j=1:ordp
      for l=1:ordp
         i = 0; k = 0;
         ordx = i+k+1;
         ordy = j+l+1;
         Lhsr(j+1+(ordp+1)*i,l+1+(ordp+1)*k) = ...
              j*l* (L2x^ordx - L1x^ordx)/ordx *...
                                  (L2y^(ordy-1) - L1y^(ordy-1))/(ordy-1);
      end
   end
   %% End of regularization terms
   
   McCoef = -(Lhs+regmu*Lhsr)\Rhs;  % - in order to have positive gap on the crack (sign is arbitrary)
   
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
   Cc = [4;3;1]; Cc = Q'*Cc; Rad = 2; zed = max(max(solup));
   
   figure;
   hold on;
%   surf(X(4:end-3),Y(4:end-3),solu(4:end-3,4:end-3));
   surf(X,Y,solup);   % The /Lx is arnaking very much
   shading interp;
   colorbar();
   teta = 0:.1:2*pi+.1; ixe = Rad*cos(teta)+Cc(1); igrec = Rad*sin(teta)+Cc(2);
   plot3( ixe, igrec, zed*ones(1,size(ixe,2)) , ...
                                  'Color', 'black',  'LineWidth', 3 );
%   drawCircle ( Cc(1), Cc(2), 2, 'Color', 'black', 'LineWidth', 3 );
   axis('equal');
   
   % "Zoom"
   figure;
   hold on;
   surf(X(7:end-6),Y(7:end-6),solup(7:end-6,7:end-6)); % Again
   shading interp;
   colorbar();
   teta = 0:.1:2*pi+.1; ixe = Rad*cos(teta)+Cc(1); igrec = Rad*sin(teta)+Cc(2);
   plot3( ixe, igrec, zed*ones(1,size(ixe,2)) , ...
                                  'Color', 'black',  'LineWidth', 3 );
   axis('equal');
   
   disp([ 'Polynomial method ', num2str(toc) ]);
end

% Plot on the line X = 4
figure;
hold on;
nys = (max(Ys)-min(Ys))/100;
Y = min(Ys):nys:max(Ys); X = 4;

if usefourier == 1
   solu = fournpp(1,1);
   for kpx=2:nbase+1
      for kpy=2:nbase+1
         solu = solu + akan(kpx-1,kpy-1)*cos(lambdax(kpx)*X)*cos(lambday(kpy)*Y)...
                     + bkan(kpx-1,kpy-1)*cos(lambdax(kpx)*X)*sin(lambday(kpy)*Y)...
                     + ckan(kpx-1,kpy-1)*sin(lambdax(kpx)*X)*cos(lambday(kpy)*Y)...
                     + dkan(kpx-1,kpy-1)*sin(lambdax(kpx)*X)*sin(lambday(kpy)*Y);
      end
   end
   
   plot( Y, solu, 'Color', 'blue' );
end

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
   Uxyz = Uxyz'(:); % Re-stick the components together
   
   uplo = passMesh3D(nodes, elements, [XYZ1';XYZ2'], [], Uxyz);
   uplo = uplo(304:end)-uplo(1:303);%  % Compute the gap
   plot( Y, uplo(3:3:end,1), 'Color', 'red' );
   %legend
end