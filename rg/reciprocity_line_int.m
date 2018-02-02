%13/10/16
% Détection de fissure 1D plane par écart à la réciprocité
% Intégrations par PG

close all;
clear all;

I = i;  % Store the complex number (i will be erased)

addpath(genpath('./tools'))

% Parameters
E          = 210000; % MPa : Young modulus
nu         = 0.3;    % Poisson ratio
fscalar    = 250;    % N.mm-1 : Loading on the plate
mat        = [0, E, nu];
mu         = 0;%1e-5;%5e-3;     % Regularization coef
%mu         = 3;     % Regularization coef
dolcurve   = 0;      % Do a L-curve or not
br         = .1;      % Noise level

usefourier = 0;
usepolys   = 1;

nbase = 2; % Number of Fourier basis functions
ordp = 4;  % Number of Polynomial basis functions

useorder = 2; % Order of the FE computation

% Boundary conditions
% first index  : index of the boundary
% second index : 1=x, 2=y
% third        : value
% [0,1,value] marks a dirichlet regularization therm on x
dirichlet  = [0,1,0 ; 0,2,0 ; 0,3,0];
neumann1   = [3,2,fscalar ; 1,2,-fscalar];
neumann2   = [2,1,fscalar ; 4,1,-fscalar];
%neumann2   = [1,1,fscalar ; 2,2,-fscalar ; 3,1,-fscalar ; 4,2,fscalar];

% First, import the mesh
if useorder == 1
   [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_c.msh' );
elseif useorder == 2
   [ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/rg_refined/plate_ct6.msh' );
end

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

% Re-arrange Lagrange in order to constrain only the boundary
%C = zeros(2*nnodes, 3);
%% First, find the ~=barycenter of the solid
%x0 = sum( nodes(:,1) ) / size(nodes,1);
%y0 = sum( nodes(:,2) ) / size(nodes,1);
%
%for i=[b2node1;b2node2;b2node3;b2node4]
%   C(2*i-1,1) = 1;% Regularisation/x
%   C(2*i,2) = 1;  % Regularisation/y
%   x = nodes(i,1);
%   y = nodes(i,2);% Regularisation/theta_z
%   C(2*i,3) = (x-x0);
%   C(2*i-1,3) = (y-y0);
%end
%K = [Kinter,C;C',zeros(3)];
%%

f1  = loading(nbloq,nodes,boundary,neumann1);
uin = K\f1;
u1 = uin(1:2*nnodes,1);
f1 = Kinter*u1;

f2  = loading(nbloq,nodes,boundary,neumann2);
uin = K\f2;
u2 = uin(1:2*nnodes,1);
f2 = Kinter*u2;

ui = reshape(u1,2,[])';  ux = ui(:,1);  uy = ui(:,2);

% Compute stress :
sigma = stress(u1,E,nu,nodes,elements,order,1,ntoelem);

% Output :
plotGMSH({ux,'U_x';uy,'U_y';u1,'U_vect';sigma,'stress'}, elements, nodes, 'output/reference');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import the uncracked domain /!\ MUST BE THE SAME !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (except for the crack)
if useorder == 1
   [ nodes2,elements2,ntoelem2,boundary2,order2] = readmesh( 'meshes/rg_refined/plate_n.msh' );
elseif useorder == 2
   [ nodes2,elements2,ntoelem2,boundary2,order2] = readmesh( 'meshes/rg_refined/plate_nt6.msh' );
end

if order ~= useorder || order2 ~= useorder
   warning('Looks like one of the meshes did not have the right order');
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

% Add the noise
u1n = u1; u2n = u2;
%br1 = randn(2*nnodes,1); br2 = randn(2*nnodes,1);
noise = load('noises/rg2d2.mat'); br1 = noise.br1; br2 = noise.br2;
u1 = ( 1 + br*br1 ) .* u1; u2 = ( 1 + br*br2 ) .* u2;

%plotGMSH({ur1,'U_bound'}, elements2, nodes2, 'bound');

% Debug : compute the displacement gap
Mc = bMass_mat (nodes, boundary, 5);
onfi = ones(2*nnodes,1);
udepx = u1; udepx(icrack5x) = u1(icrack6x);
intx = onfi'*Mc*(u1-udepx);
udepy = u1; udepy(icrack5y) = u1(icrack6y);
inty = onfi'*Mc*(u1-udepy);
intt = sqrt(intx^2+inty^2);
%
% Debug : compute n
x6 = nodes(6,1); y6 = nodes(6,2); x5 = nodes(5,1); y5 = nodes(5,2);
n1 = 1; n2 = -n1*(x6-x5)/(y6-y5);
non = sqrt(n1^2+n2^2);
n1 = n1/non; n2 = n2/non;
R1d = n1*intx; R2d = n2*inty; R3d = .5*(n1*inty+n2*intx);
normRd = sqrt( 2*(R1d^2+R2d^2+2*R3d^2) - (R1d+R2d)^2 );
R1dd = [R1d,R3d;R3d,R2d]/normRd;
[phi1d,Lambda1d] = eig( R1dd/normRd );

%% Debug : compute the displacement gap (second one)
Mc = bMass_mat (nodes, boundary, 5);
onfi = ones(2*nnodes,1);
udepx = u2; udepx(icrack5x) = u2(icrack6x);
intx2 = onfi'*Mc*(u2-udepx);
udepy = u2; udepy(icrack5y) = u2(icrack6y);
inty2 = onfi'*Mc*(u2-udepy);
intt2 = sqrt(intx2^2+inty2^2);

% Mean of the tangential displacement
intta1 = intx*n2-inty*n1;
intta2 = intx2*n2-inty2*n1;
x6 = nodes(6,1); y6 = nodes(6,2); x5 = nodes(5,1); y5 = nodes(5,2);
CteR = x6*n1+y6*n2;
Rt1R = CteR*abs(intta1); Rt2R = CteR*abs(intta2);

% Debug : compute R2
R1d2 = n1*intx2; R2d2 = n2*inty2; R3d2 = .5*(n1*inty2+n2*intx2);
normRd2 = sqrt( 2*(R1d2^2+R2d2^2+2*R3d2^2) - (R1d2+R2d2)^2 );
R2dd = [R1d2,R3d2;R3d2,R2d2];
[phi2d,Lambda2d] = eig( R2dd/normRd2 );

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
   if order == 2
      no4 = boundary(i,4);
      urr1(i,5:6) = u1( [2*no4-1,2*no4] );
      urr2(i,5:6) = u2( [2*no4-1,2*no4] );
   end
end

%% First determine the crack's line.
% Compute and apply v fields (for plane constraint)
%v1 = zeros(2*nnodes2, 1);
%v2 = zeros(2*nnodes2, 1);
%v3 = zeros(2*nnodes2, 1);
%for i=1:nnodes2
%   x = nodes2(i,1);
%   y = nodes2(i,2);
%   v1(2*i-1) = 1/E*x;
%   v2(2*i-1) = -nu/E*x;
%   v3(2*i-1) = (1+nu)/(2*E)*y;
%   v1(2*i)   = -nu/E*y;
%   v2(2*i)   = 1/E*y;
%   v3(2*i)   = (1+nu)/(2*E)*x;
%end
%%
%f1 = Kinter2*v1;
%f2 = Kinter2*v2;
%f3 = Kinter2*v3;

% Clean redondant stuff in indexbound2
i = 1;
while i <= size(indexbound2,1)
   if find( indexbound2(i) == indexbound2(1:i-1) )
      indexbound2(i) = [];
      i = i-1;
   end
   i = i+1;
end

%R11e = (fr1(indexbound2)'*v1(indexbound2) - f1(indexbound2)'*ur1(indexbound2));
%R21e = (fr1(indexbound2)'*v2(indexbound2) - f2(indexbound2)'*ur1(indexbound2));
%R31e = (fr1(indexbound2)'*v3(indexbound2) - f3(indexbound2)'*ur1(indexbound2));
%%[R11s, R31s ; R31s, R21s]
%
%R12e = (fr2(indexbound2)'*v1(indexbound2) - f1(indexbound2)'*ur2(indexbound2));
%R22e = (fr2(indexbound2)'*v2(indexbound2) - f2(indexbound2)'*ur2(indexbound2));
%R32e = (fr2(indexbound2)'*v3(indexbound2) - f3(indexbound2)'*ur2(indexbound2));

% Second method : loop over boundaries to compute the integrals
R11 = 0; R21 = 0; R31 = 0; R12 = 0; R22 = 0; R32 = 0;
for i=1:nboun2 % boundary1 and boundary 2 are supposed to be the same
   bonod = boundary2(i,:); exno = extnorm2(i,:)';

   no1 = bonod(2); no2 = bonod(3);
   x1 = nodes2(no1,1); y1 = nodes2(no1,2);
   x2 = nodes2(no2,1); y2 = nodes2(no2,2);
   len = sqrt( (x1-x2)^2 + (y1-y2)^2 );
   
%   % Compute sigma.n
   nomat = [ exno(1), 0, exno(2) ; ...
             0, exno(2), exno(1) ];% Such that sigma.n = nomat*sier
   
   if order==1
      Ng = 1;
   elseif order==2
      Ng  = 2; no3 = bonod(4);
      x3  = nodes2(no3,1); y3 = nodes2(no3,2);
   end
   [ Xg, Wg ] = gaussPt1d( Ng );
             
   for j=1:Ng
      xg = Xg(j); wg = Wg(j);
      
      % Interpolations
      if order==1
         uer1 = transpose( (1-xg)*urr1(i,1:2) + xg*urr1(i,3:4) ); % [ux;uy] on the
         uer2 = transpose( (1-xg)*urr2(i,1:2) + xg*urr2(i,3:4) ); % Gauss point
         xgr  = ( (1-xg)*[x1;y1] + xg*[x2;y2] ); % abscissae of the Gauss point in the physical space
      elseif order==2
         uer1 = transpose( urr1(i,1:2) + ...
                xg*(4*urr1(i,5:6)-3*urr1(i,1:2)-urr1(i,3:4)) + ...   % [ux;uy] on the
                xg^2*(2*urr1(i,3:4)+2*urr1(i,1:2)-4*urr1(i,5:6)) ); % Gauss point
         uer2 = transpose( urr2(i,1:2) + ...
                xg*(4*urr2(i,5:6)-3*urr2(i,1:2)-urr2(i,3:4)) + ...  
                xg^2*(2*urr2(i,3:4)+2*urr2(i,1:2)-4*urr2(i,5:6)) );
         xgr  = ( [x1;y1] + xg*(4*[x3;y3]-3*[x1;y1]-[x2;y2]) + ...
                xg^2*(2*[x2;y2]+2*[x1;y1]-4*[x3;y3]) );
      end

      % Reference force from the BC's
      if exno(1) == 1
         fer1 = [0;0]; fer2 = [fscalar;0];
      elseif exno(1) == -1
         fer1 = [0;0]; fer2 = -[fscalar;0];
      elseif exno(2) == 1
         fer1 = [0;fscalar]; fer2 = [0;0];
      elseif exno(2) == -1
         fer1 = -[0;fscalar]; fer2 = [0;0];
      end

       % Test fields
      s1 = [1;0;0]; s2 = [0;1;0]; s3 = [0;0;.5];
      f1 = nomat*s1; f2 = nomat*s2; f3 = nomat*s3;
      v1 = 1/E*[ xgr(1) ; -nu*xgr(2) ];
      v2 = 1/E*[ -nu*xgr(1) ; xgr(2) ];
      v3 = (1+nu)/(2*E)*[ xgr(2) ; xgr(1) ];

      R11 = R11 + len * wg * ( fer1'*v1 - f1'*uer1 ); % increment the integral
      R21 = R21 + len * wg * ( fer1'*v2 - f2'*uer1 );
      R31 = R31 + len * wg * ( fer1'*v3 - f3'*uer1 );
      R12 = R12 + len * wg * ( fer2'*v1 - f1'*uer2 );
      R22 = R22 + len * wg * ( fer2'*v2 - f2'*uer2 );
      R32 = R32 + len * wg * ( fer2'*v3 - f3'*uer2 );
   end
end
%RRe = [R11e,R31e;R31e,R21e];
%RR  = [R11,R31;R31,R21];

% Normalize R1
%R11 = R1d; R21 = R2d; R31 = R3d;%DEBUG
normR1 = sqrt( 2*(R11^2+R21^2+2*R31^2) - (R11+R21)^2 );
R1b1 = R11/normR1; R2b1 = R21/normR1; R3b1 = R31/normR1;
%[R1b1, R3b1 ; R3b1, R2b1]

lam11 = (1+R1b1+R2b1)/2; lam21 = -(1-R1b1-R2b1)/2; R1 = [R11,R31;R31,R21];
[phi1,Lambda1] = eig( [R1b1,R3b1;R3b1,R2b1] );
dir11 = [ sqrt( abs(Lambda1(1,1)) ) ; sqrt( abs(Lambda1(2,2)) ) ];
dir21 = [ sign(Lambda1(1,1)) * sqrt( abs(Lambda1(1,1)) ) ;...
          sign(Lambda1(2,2)) * sqrt( abs(Lambda1(2,2)) ) ];
%dir21 = [ sign(Lambda1(1,1)) * sqrt( abs(Lambda1(1,1)) ) ; sqrt( abs(Lambda1(2,2)) ) ];
dir31 = [ sign(Lambda1(1,1)) * sqrt( abs(Lambda1(1,1)) ) ;...
          sqrt( abs(Lambda1(2,2)) ) ];
dir41 = [ sqrt( abs(Lambda1(1,1)) ) ;...
          sign(Lambda1(2,2)) * sqrt( abs(Lambda1(2,2)) ) ];
          
% DEBUG : check Lambda
%Rdeb    = .5*([n1;n2]*[intx,inty] + [intx;inty]*[n1,n2]); % OK
%Lamdeb1 = .5/normR1 * phi1'*( [n1;n2]*[intx,inty] + [intx;inty]*[n1,n2] )*phi1; % OK
          
%dir1 = [ sqrt( abs(lam1) ) ; sqrt( abs(lam2) ) ];
%dir2 = [ sqrt( abs(lam1) ) ; -sqrt( abs(lam2) ) ];
dir11 = phi1*dir11; dir11 = dir11/norm(dir11);
dir21 = phi1*dir21; dir21 = dir21/norm(dir21); % Normal candidates (1)
dir31 = phi1*dir31; dir31 = dir31/norm(dir31);
dir41 = phi1*dir41; dir41 = dir41/norm(dir41);

% Normalize R2
%R12 = R1d2; R22 = R2d2; R32 = R3d2;%DEBUG
normR2 = sqrt( 2*(R12^2+R22^2+2*R32^2) - (R12+R22)^2 );
R1b2 = R12/normR2; R2b2 = R22/normR2; R3b2 = R32/normR2;
%[R12, R32 ; R32, R22]

R2 = [R12,R32;R32,R22];
[phi2,Lambda2] = eig( [R1b2,R3b2;R3b2,R2b2] );
dir12 = [ sqrt( abs(Lambda2(1,1)) ) ; sqrt( abs(Lambda2(2,2)) ) ];
dir22 = [ sign(Lambda2(1,1)) * sqrt( abs(Lambda2(1,1)) ) ;...
          sign(Lambda2(2,2)) * sqrt( abs(Lambda2(2,2)) ) ];
%dir22 = [ sign(Lambda2(1,1)) * sqrt( abs(Lambda2(1,1)) ) ; sqrt( abs(Lambda2(2,2)) ) ];
dir12 = phi2*dir12; dir12 = dir12/norm(dir12);
dir22 = phi2*dir22; dir22 = dir22/norm(dir22); % Normal candidates (2)
dir32 = [ sign(Lambda2(1,1)) * sqrt( abs(Lambda2(1,1)) ) ;...
          sqrt( abs(Lambda2(2,2)) ) ];
dir42 = [ sqrt( abs(Lambda2(1,1)) ) ;...
          sign(Lambda2(2,2)) * sqrt( abs(Lambda2(2,2)) ) ];

% Find the real normal (closest candidates)
dist11 = norm(dir11-dir12); dist12 = norm(dir11-dir22);
dist21 = norm(dir21-dir12); dist22 = norm(dir21-dir22);
dist11m = norm(dir11+dir12); dist12m = norm(dir11+dir22);
dist21m = norm(dir21+dir12); dist22m = norm(dir21+dir22);
mindist = min([dist11,dist12,dist21,dist22,dist11m,dist12m,dist21m,dist22m]);

if dist11 == mindist
   normal = .5*(dir11+dir12);
elseif dist22 == mindist
   normal = .5*(dir21+dir22);
elseif dist12 == mindist
   normal = .5*(dir11+dir22);
elseif dist21 == mindist
   normal = .5*(dir21+dir12);
elseif dist11m == mindist
   normal = .5*(dir11-dir12);
elseif dist22m == mindist
   normal = .5*(dir21-dir22);
elseif dist12m == mindist
   normal = .5*(dir11-dir22);
elseif dist21m == mindist
   normal = .5*(dir21-dir12);
end

% Build the base-change matrix : [x;y] = Q*[X;Y], [X;Y] = Q'*[x;y]
Q = [ normal(2), normal(1) ; - normal(1), normal(2) ];
%Qref = [ n2, n1 ;-n1, n2 ]; % Reference matrix
Qref = Q;

% Plot the mesh
%figure
%ret = patch('Faces',elements(:,1:3),'Vertices',nodes,'FaceAlpha',0);
%axis('equal');

%% Plot the mesh with the normal
%figure
%hold on;
%ret = patch('Faces',elements2(:,1:3),'Vertices',nodes2,'FaceAlpha',0);
%x6 = nodes(6,1); y6 = nodes(6,2); x5 = nodes(5,1); y5 = nodes(5,2); 
%xc = .5*(max(nodes(:,1)) + min(nodes(:,1)));
%yc = .5*(max(nodes(:,2)) + min(nodes(:,2)));
%
%x1 = xc + dir11(1); y1 = yc + dir11(2);
%x2 = xc + dir21(1); y2 = yc + dir21(2);
%xn = xc + n1; yn = yc + n2;
%plot( [xc,x1], [yc,y1] ,'Color', 'red', 'LineWidth',3);
%plot( [xc,x2], [yc,y2] ,'Color', 'green', 'LineWidth',3);
%plot( [xc,xn], [yc,yn] ,'Color', 'cyan', 'LineWidth',2);
%plot( [x5,x6], [y5,y6] ,'Color', 'magenta', 'LineWidth',3);
%axis('equal');
%
%% And the other one
%figure
%hold on;
%ret = patch('Faces',elements2(:,1:3),'Vertices',nodes2,'FaceAlpha',0);
%x6 = nodes(6,1); y6 = nodes(6,2); x5 = nodes(5,1); y5 = nodes(5,2); 
%xc = .5*(max(nodes(:,1)) + min(nodes(:,1)));
%yc = .5*(max(nodes(:,2)) + min(nodes(:,2)));
%
%x1 = xc + dir12(1); y1 = yc + dir12(2);
%x2 = xc + dir22(1); y2 = yc + dir22(2);
%xn = xc + n1; yn = yc + n2;
%plot( [xc,x1], [yc,y1] ,'Color', 'red', 'LineWidth',3);
%plot( [xc,x2], [yc,y2] ,'Color', 'green', 'LineWidth',3);
%plot( [xc,xn], [yc,yn] ,'Color', 'cyan', 'LineWidth',2);
%plot( [x5,x6], [y5,y6] ,'Color', 'magenta', 'LineWidth',3);
%axis('equal');

% Plot the normal
figure
hold on;
%ret = patch('Faces',elements2(:,1:3),'Vertices',nodes2,'FaceAlpha',0);
x6 = nodes(6,1); y6 = nodes(6,2); x5 = nodes(5,1); y5 = nodes(5,2); 
Xmin = min(nodes(:,1)); Xmax = max(nodes(:,1));
Ymin = min(nodes(:,2)); Ymax = max(nodes(:,2));

xc = .5*(max(nodes(:,1)) + min(nodes(:,1)));
yc = .5*(max(nodes(:,2)) + min(nodes(:,2)));

% Rectangle
plot( [Xmin,Xmax], [Ymin,Ymin], 'Color', 'black' );
plot( [Xmax,Xmax], [Ymin,Ymax], 'Color', 'black' );
plot( [Xmax,Xmin], [Ymax,Ymax], 'Color', 'black' );
plot( [Xmin,Xmin], [Ymax,Ymin], 'Color', 'black' );

x2 = xc + 3*normal(1); y2 = yc + 3*normal(2);
xn = xc + 3*n1; yn = yc + 3*n2;
plot( [xc,x2], [yc,y2] ,'Color', 'red', 'LineWidth',3);
plot( [xc,xn], [yc,yn] ,'Color', 'cyan', 'LineWidth',2);
plot( [x5,x6], [y5,y6] ,'Color', 'magenta', 'LineWidth',3);
axis('equal');

%% Then the constant

% First, find the minimal point (in order to have Cte < 0)
norep = Q'*nodes'; K = min(norep(2,:));
%
%vt = zeros(2*nnodes2, 1);
%%vn = zeros(2*nnodes2, 1);
%Xs = zeros(nnodes2, 1);
%Ys = zeros(nnodes2, 1);
%for i=1:nnodes2
%   x = nodes2(i,1);
%   y = nodes2(i,2);
%   % Change base (for coordiantes)
%   ixigrec = Q'*[x;y]; X = ixigrec(1); Y = ixigrec(2)-K;
%   Xs(i) = X; Ys(i) = Y;
%
%   vloc = [ -X^2/(2*E) + (2+nu)*Y^2/(2*E) ; nu*X*Y/E ];
%   % Change base (for vector), in the other direction
%   vxy = Q*vloc; vt(2*i-1) = vxy(1); vt(2*i) = vxy(2);
%end

% Norm of [[ut]] (case 1)
normT  = sqrt( abs( 2*(R11^2+R21^2+2*R31^2) - 2*(R11+R21)^2 ) );
normT2 = sqrt( abs( 2*(R12^2+R22^2+2*R32^2) - 2*(R12+R22)^2 ) );

%ft = Kinter2*vt;
%
%Rte = (fr1(indexbound2)'*vt(indexbound2) - ft(indexbound2)'*ur1(indexbound2));
%Rt2e = (fr2(indexbound2)'*vt(indexbound2) - ft(indexbound2)'*ur2(indexbound2));

Rt = 0; Rt2 = 0;
for i=1:nboun2 % boundary1 and boundary 2 are supposed to be the same
   bonod = boundary2(i,:); exno = extnorm2(i,:)';

   no1 = bonod(2); no2 = bonod(3);
   x1 = nodes2(no1,1); y1 = nodes2(no1,2);
   x2 = nodes2(no2,1); y2 = nodes2(no2,2);
   len = sqrt( (x1-x2)^2 + (y1-y2)^2 );
   
   if order==1
      Ng = 2;
   elseif order==2
      Ng  = 2; no3 = bonod(4);  % Yes, 2 and not 3 (because order = 2Ng-1)
      x3  = nodes2(no3,1); y3 = nodes2(no3,2);
   end
   [ Xg, Wg ] = gaussPt1d( Ng );
             
   for j=1:Ng
      xg = Xg(j); wg = Wg(j);
      
      % Interpolations
      if order==1
         uer1 = transpose( (1-xg)*urr1(i,1:2) + xg*urr1(i,3:4) ); % [ux;uy] on the
         uer2 = transpose( (1-xg)*urr2(i,1:2) + xg*urr2(i,3:4) ); % Gauss point
         xgr  = ( (1-xg)*[x1;y1] + xg*[x2;y2] ); % abscissae of the Gauss point in the physical space
      elseif order==2
         uer1 = transpose( urr1(i,1:2) + ...
                xg*(4*urr1(i,5:6)-3*urr1(i,1:2)-urr1(i,3:4)) + ...   % [ux;uy] on the
                xg^2*(2*urr1(i,3:4)+2*urr1(i,1:2)-4*urr1(i,5:6)) ); % Gauss point
         uer2 = transpose( urr2(i,1:2) + ...
                xg*(4*urr2(i,5:6)-3*urr2(i,1:2)-urr2(i,3:4)) + ...  
                xg^2*(2*urr2(i,3:4)+2*urr2(i,1:2)-4*urr2(i,5:6)) );
         xgr  = ( [x1;y1] + xg*(4*[x3;y3]-3*[x1;y1]-[x2;y2]) + ...
                xg^2*(2*[x2;y2]+2*[x1;y1]-4*[x3;y3]) );
      end

      % Reference force from the BC's
      if exno(1) == 1
         fer1 = [0;0]; fer2 = [fscalar;0];
      elseif exno(1) == -1
         fer1 = [0;0]; fer2 = -[fscalar;0];
      elseif exno(2) == 1
         fer1 = [0;fscalar]; fer2 = [0;0];
      elseif exno(2) == -1
         fer1 = -[0;fscalar]; fer2 = [0;0];
      end

      % Test fields
      ixigrec = Q'*[xgr(1);xgr(2)]; X = ixigrec(1); Y = ixigrec(2)-K;
       
      sloc = [-X,Y;Y,0];
      st = Q*sloc*Q';
      ft = st*exno;
      
      vloc = [ -X^2/(2*E) + (2+nu)*Y^2/(2*E) ; nu*X*Y/E ];
      vt = Q*vloc;

      Rt  = Rt + len * wg * ( fer1'*vt - ft'*uer1 ); % increment the integral
      Rt2 = Rt2 + len * wg * ( fer2'*vt - ft'*uer2 );
   end
end
%Rt = Rte; Rt2 = Rt2e;
Cte  = min( Rt/normT, -Rt/normT) - K;       % Select the negative one
Cte2 = min( Rt2/normT2, -Rt2/normT2 ) - K;  %  /!\ The sign depends on the test case
% Plot the crack, and its estimated lines (there are Y +/- Cte)
figure
hold on;
%ret = patch('Faces',elements2(:,1:3),'Vertices',nodes2,'FaceAlpha',0);
Xmin = min(nodes(:,1)); Xmax = max(nodes(:,1));
Ymin = min(nodes(:,2)); Ymax = max(nodes(:,2));
x6 = nodes(6,1); y6 = nodes(6,2); x5 = nodes(5,1); y5 = nodes(5,2); 

Vp1 = [20;-Cte]; Vp2 = [-20;-Cte];
Vm1 = [20;-Cte2]; Vm2 = [-20;-Cte2];
vp1 = Q*Vp1; vm1 = Q*Vm1; vp2 = Q*Vp2; vm2 = Q*Vm2;

% Rectangle
plot( [Xmin,Xmax], [Ymin,Ymin], 'Color', 'black' );
plot( [Xmax,Xmax], [Ymin,Ymax], 'Color', 'black' );
plot( [Xmax,Xmin], [Ymax,Ymax], 'Color', 'black' );
plot( [Xmin,Xmin], [Ymax,Ymin], 'Color', 'black' );

plot( [vp1(1), vp2(1)], [vp1(2), vp2(2)] ,'Color', 'red', 'LineWidth',3);
plot( [vm1(1), vm2(1)], [vm1(2), vm2(2)] ,'Color', 'green', 'LineWidth',3);
plot( [x5,x6], [y5,y6] ,'Color', 'magenta', 'LineWidth',4);
axis('equal');

%% And now, the crack itself

% Compute L (provisionnal formula assuming that something<Cte<0
tangent = [normal(2) ; -normal(1)]; b = max(nodes(:,1)) - min(nodes(:,1));
a = tangent(2)/tangent(1)*b; L = sqrt(a^2+b^2);

% Convenient values for plot
udepx  = u1n(icrack5x)-u1n(icrack6x);
udepy  = u1n(icrack5y)-u1n(icrack6y);
Xx     = nodes(b2node5,1); Yy = nodes(b2node5,2);
ubase  = Qref'*[udepx,udepy]'; ubase = ubase';
XY     = Qref'*[Xx,Yy]'; XY = XY'; %XY(:,2) = -XY(:,2);
[newX, orderX] = sort(XY(:,1));          % In case order > 1, need to sort
offset = Qref(1,2)/Qref(1,1)*CteR;   % /!\ Only in case rectangle

left  = offset : (-offset+newX(1))/(size(newX,1)) : newX(1) ;
right = newX(end) : (offset + L - newX(end))/(size(newX,1)) : offset + L ;
newXo = newX;
newX  = [left(1:end-1)';
         newX;
         right(2:end)']; % Add a few points

if usefourier == 1
   Rp     = zeros(nbase+1,1);
   Rm     = zeros(nbase+1,1);
   Rpe    = zeros(nbase+1,1);
   Rme    = zeros(nbase+1,1);
   
   lambda = zeros(nbase+1,1);
   fourn  = zeros(nbase+1,1);
   fournm = zeros(nbase+1,1);
   akan   = zeros(nbase,1);
   bkan   = zeros(nbase,1);
   
   for kp=2:nbase+1
      k = kp-1;
      vp = zeros(2*nnodes2,1);
      vm = zeros(2*nnodes2,1);
      lambda(kp) = 2*k*pi/L;
      for sx = [1,-1]  % sx=-1 is not really used, but Debug stuff
         lambdae = sx*lambda(kp);
%         for i=1:nnodes2
%            x = nodes2(i,1);
%            y = nodes2(i,2);
%            % Change base (for coordinates)
%            ixigrec = Q'*[x;y]; X = ixigrec(1); Y = ixigrec(2)+Cte2;
%            Xs(i) = X; Ys(i) = Y;
%            
%            v1 = -I*lambda(kp)*exp(-I*lambdae*X)* ...
%                               ( exp(lambda(kp)*Y)+exp(-lambda(kp)*Y) );
%            v2 = lambda(kp)*exp(-I*lambdae*X)* ...
%                               ( exp(lambda(kp)*Y)-exp(-lambda(kp)*Y) );
%            vloc = [ v1 ; v2 ];
%            % Change base (for vector), in the other direction
%            vxy = Q*vloc; vp(2*i-1) = vxy(1); vp(2*i) = vxy(2);
%         end
%         fp = Kinter2*vp;
%         
%         % Fourier coefficient
%         if sx == 1
%            Rpe(kp) = (fr1(indexbound2)'*vp(indexbound2) - fp(indexbound2)'*ur1(indexbound2));
%            %Rp(kp) = (fr2(indexbound2)'*vp(indexbound2) - fp(indexbound2)'*ur2(indexbound2));
%            if dolcurve == 0 % No L-curve stuff
%               fourn(kp) = - 1/(1+mu*k^2) *(1+nu)/(2*E*L*lambda(kp)^2)*Rpe(kp);
%               akan(k) = 2*real(fourn(kp));
%               bkan(k) = 2*imag(fourn(kp));
%            end
%         else
%            Rpme = (fr1(indexbound2)'*vp(indexbound2) - fp(indexbound2)'*ur1(indexbound2));
%            fournm(kp) = -(1+nu)/(2*E*L*lambda(kp)^2)*Rpme;
%         end
         
         %% Alternative way
         if sx == 1 Rp(kp) = 0; else Rpm = 0; end
         for i=1:nboun2 % boundary1 and boundary 2 are supposed to be the same
            bonod = boundary2(i,:); exno = extnorm2(i,:)';
         
            no1 = bonod(2); no2 = bonod(3);
            x1 = nodes2(no1,1); y1 = nodes2(no1,2);
            x2 = nodes2(no2,1); y2 = nodes2(no2,2);
            len = sqrt( (x1-x2)^2 + (y1-y2)^2 );
            
            if order==1
               Ng = 8; % Give everything you can
            elseif order==2
               Ng  = 8; no3 = bonod(4);
               x3  = nodes2(no3,1); y3 = nodes2(no3,2);
            end
            [ Xg, Wg ] = gaussPt1d( Ng );
                      
            for j=1:Ng
               xg = Xg(j); wg = Wg(j);
               
               % Interpolations
               if order==1
                  uer1 = transpose( (1-xg)*urr1(i,1:2) + xg*urr1(i,3:4) ); % [ux;uy] on the
                  uer2 = transpose( (1-xg)*urr2(i,1:2) + xg*urr2(i,3:4) ); % Gauss point
                  xgr  = ( (1-xg)*[x1;y1] + xg*[x2;y2] ); % abscissae of the Gauss point in the physical space
               elseif order==2
                  uer1 = transpose( urr1(i,1:2) + ...
                         xg*(4*urr1(i,5:6)-3*urr1(i,1:2)-urr1(i,3:4)) + ...   % [ux;uy] on the
                         xg^2*(2*urr1(i,3:4)+2*urr1(i,1:2)-4*urr1(i,5:6)) ); % Gauss point
                  uer2 = transpose( urr2(i,1:2) + ...
                         xg*(4*urr2(i,5:6)-3*urr2(i,1:2)-urr2(i,3:4)) + ...  
                         xg^2*(2*urr2(i,3:4)+2*urr2(i,1:2)-4*urr2(i,5:6)) );
                  xgr  = ( [x1;y1] + xg*(4*[x3;y3]-3*[x1;y1]-[x2;y2]) + ...
                         xg^2*(2*[x2;y2]+2*[x1;y1]-4*[x3;y3]) );
               end
         
               % Reference force from the BC's
               if exno(1) == 1
                  fer1 = [0;0]; fer2 = [fscalar;0];
               elseif exno(1) == -1
                  fer1 = [0;0]; fer2 = -[fscalar;0];
               elseif exno(2) == 1
                  fer1 = [0;fscalar]; fer2 = [0;0];
               elseif exno(2) == -1
                  fer1 = -[0;fscalar]; fer2 = [0;0];
               end
         
               % Test fields
               ixigrec = Q'*[xgr(1);xgr(2)]; X = ixigrec(1); Y = ixigrec(2)+Cte2;
               
               sloc11 = E/(1+nu)*lambdae^2*exp(-I*lambdae*X)* ...
                               ( exp(lambda(kp)*Y)+exp(-lambda(kp)*Y) );
               sloc12 = E/(1+nu)*I*lambdae^2*exp(-I*lambdae*X)* ...
                               ( exp(lambda(kp)*Y)-exp(-lambda(kp)*Y) );

               sloc = [-sloc11,-sloc12;-sloc12,sloc11];
               sp = Q*sloc*Q';
               fp = sp*exno;
               
               v1 = -I*lambda(kp)*exp(-I*lambdae*X)* ...
                               ( exp(lambda(kp)*Y)+exp(-lambda(kp)*Y) );
               v2 = lambda(kp)*exp(-I*lambdae*X)* ...
                               ( exp(lambda(kp)*Y)-exp(-lambda(kp)*Y) );
               vloc = [ v1 ; v2 ];
               vp = Q*vloc;

               if sx == 1
                  Rp(kp) = Rp(kp) + len * wg * ( fer1'*vp - fp'*uer1 );
                  if dolcurve == 0 % No L-curve stuff
                     fourn(kp) = - 1/(1+mu*k^2) *(1+nu)/(2*E*L*lambda(kp)^2)*Rp(kp);
                     akan(k) = 2*real(fourn(kp));
                     bkan(k) = 2*imag(fourn(kp));
                  end
               else
                  Rpm = Rpm + len * wg * ( fer1'*vp - fp'*uer1 );
                  fournm(kp) = -(1+nu)/(2*E*L*lambda(kp)^2)*Rpm;
               end  
            end
         end
            
      end
   end
   
   % The constant term
   fourn(1) = -(R11+R21)/L;
   
   % Invert the operator
   if dolcurve == 1
      listmu = -1:.1:3;
      resid  = zeros( size(listmu) );
      regno  = zeros( size(listmu) );
      rhsk   = zeros(nbase+1,1);
      i = 1;
      for lnmu1 = listmu
         mu1 = 10^lnmu1;
         for kp=2:nbase+1
            %lhsk(kp)  = - (2*E*L*lambda(kp)^2)/(1+nu);
            rhsk(kp)  = - (1+nu)/(2*E*L*lambda(kp)^2)*Rp(kp);
            fourn(kp) = 1/(1+mu1*(kp-1)^2) * rhsk(kp);
            akan(k)   = 2*real(fourn(kp));
            bkan(k)   = 2*imag(fourn(kp));
            resid(i)  = resid(i) + abs(fourn(kp) - rhsk(kp))^2;
            regno(i)  = regno(i) + kp^2 * abs(fourn(kp))^2;
         end
         i = i+1;
      end
      figure;
      loglog(resid,regno);
      xlabel('residual (log)')
      ylabel('norm (log)')
   end
   
   %% DEBUG :
   %ui = reshape(imag(vp),2,[])';  ux = ui(:,1);  uy = ui(:,2);
   %plotGMSH({ux,'U_x';uy,'U_y';imag(vp),'U_vect'}, elements2, nodes2, 'test field');
   
   % Plot the reference normal displacement (first test case)
   soluf = fourn(1) + sum ( [0*newX' ;...  % Hack in case there is 1 element only
              akan.*cos(lambda(2:end)*newX') + bkan.*sin(lambda(2:end)*newX') ] );
   soluf = soluf';
   
   figure
   hold on;
   set(gca, 'fontsize', 20);
   plot(newX, [0*newXo;ubase(:,2);0*newXo],'linewidth',3)
   plot(newX, soluf, 'Color', 'red','linewidth',3)
   legend('Reference gap','reconstructed gap')
   xlabel('X')
   ylabel('[[u]]')

   %% Error computation
   ref = [0*newXo;ubase(:,2);0*newXo]; errorU = ref - soluf;
   nelem = size(newX,1)-1;
   nerrnonrom = 0; nnormaliz = 0; % Well those are ugly names...
   for i=1:nelem
      dlen = newX(i+1)-newX(i);
      Me = dlen/3*[1,.5;.5,1];
      ue = errorU([i,i+1]); uer = ref([i,i+1]);
      nerrnonrom = nerrnonrom + ue'*Me*ue;
      nnormaliz  = nnormaliz  + uer'*Me*uer;
   end
   four_error = nerrnonrom / nnormaliz;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mean square polynoms
if usepolys == 1
   % First, build the polynomial test functions.
   coef = zeros(ordp+2, ordp+1);
   coef( 1:2, 1 ) = L*[-(1-nu^2)/(nu*E) ; 0]; % L* because of the homotecy
   
   for k = 1:ordp
      Rhsco = zeros(k+1,1);
      Lhsco = zeros(k+1);
      
      Axx = zeros(k-1,2);
      Axy = zeros(k-1,2);
      Ayy = zeros(k-1,2);
      Bxx = zeros(k-1,2);
      Bxy = zeros(k-1,2);
      Byy = zeros(k-1,2);
      
      azero = -(1-nu^2)/(nu*E*(k+1));
   
      for i=0:floor( (k-1)/2 )
         Lhsco(2*i+1,2*i+1) = (k+1-2*i)*(k-2*i)*E/(1-nu^2);  % a_i
         Lhsco(2*i+1,2*i+2) = (2*i+1)*(k-2*i) * ( nu*E/(1-nu^2) + E/(2*(1+nu)) ); % b_i
         Lhsco(2*i+1,2*i+3) = (2*i+1)*(2*i+2)*E/(2*(1+nu));  % a_{i+1}
      end
      
      for i=1:floor( k/2 )
         Lhsco(2*i,2*i)   = (k-2*i+2)*(k-2*i+1)*E/(2*(1+nu)) ; % b_{i-1}
         Lhsco(2*i,2*i+1) = 2*i*(k+1-2*i)*( E/(2*(1+nu)) + nu*E/(1-nu^2) ) ; % a_i
         Lhsco(2*i,2*i+2) = 2*i*(2*i+1)*E/(1-nu^2) ; %b_i
      end
   
      C = [ eye(2) , zeros(2,k)]; %impose a0 = a0 and b0 = 0 ...
      Lhsco = [ Lhsco ; C ]; % ... in order to find an unique solution
      Rhsco = L*[ Rhsco ; azero ; 0 ];% Alternatively, you can remove ..
        % the L* on Rhsco and Lhs (it's because of the variable change )
      
      Lhsco(size(Lhsco,1)-2,:) = [];  % 'cause square matrix is life
      Rhsco(size(Rhsco,1)-2,:) = [];
      
      coef( 1:k+2, k+1 ) = Lhsco\Rhsco;
   end

   % Place zeros in coef
   coefa = coef;
   coefb = coef;
   for i=1:size(coefa,1)
      if mod(i,2) == 0
         coefa(i,:) = 0;
      end
      if mod(i,2) == 1
         coefb(i,:) = 0;
      end
   end
   
   %% Compute the RG
   Rhse = zeros(ordp+1,1);
   Rhs  = zeros(ordp+1,1);
   vpa  = zeros(2*nnodes2, 1);

   for k=0:ordp
%      for i=1:nnodes2
%         x = nodes2(i,1);
%         y = nodes2(i,2);
%         ixigrec = Q'*[x;y]; X = (ixigrec(1)-offset)/L; Y = (ixigrec(2)+Cte2)/L;
%         
%         % Build X^k*Y^j
%         GROX = zeros(ordp+2,1);
%         for j = 0:k+1
%            GROX(j+1) = X^(k+1-j)*Y^j;
%         end
%         
%         vloc = [ coefa(:,k+1)'*GROX ; coefb(:,k+1)'*GROX ];
%         vxy = Q*vloc; vpa(2*i-1) = vxy(1); vpa(2*i) = vxy(2);
%      end
%      fpa = Kinter2*vpa;
%      Rhse(k+1) = (fr1'*vpa - fpa'*ur1);

      %% Alternative way
      Rhs(k+1) = 0;
      for i=1:nboun2 % boundary1 and boundary 2 are supposed to be the same
         bonod = boundary2(i,:); exno = extnorm2(i,:)';
      
         no1 = bonod(2); no2 = bonod(3);
         x1 = nodes2(no1,1); y1 = nodes2(no1,2);
         x2 = nodes2(no2,1); y2 = nodes2(no2,2);
         len = sqrt( (x1-x2)^2 + (y1-y2)^2 );
         
         if order==1
            Ng = max(1, ceil((k-3)/2)); % Exact integration ceil((ordp+3)/2);
         elseif order==2
            Ng  = max(1, ceil((k-2)/2)); no3 = bonod(4);
            x3  = nodes2(no3,1); y3 = nodes2(no3,2);
         end
         [ Xg, Wg ] = gaussPt1d( Ng );
                   
         for j=1:Ng
            xg = Xg(j); wg = Wg(j);
            
            % Interpolations
            if order==1
               uer1 = transpose( (1-xg)*urr1(i,1:2) + xg*urr1(i,3:4) ); % [ux;uy] on the
               uer2 = transpose( (1-xg)*urr2(i,1:2) + xg*urr2(i,3:4) ); % Gauss point
               xgr  = ( (1-xg)*[x1;y1] + xg*[x2;y2] ); % abscissae of the Gauss point in the physical space
            elseif order==2
               uer1 = transpose( urr1(i,1:2) + ...
                      xg*(4*urr1(i,5:6)-3*urr1(i,1:2)-urr1(i,3:4)) + ...   % [ux;uy] on the
                      xg^2*(2*urr1(i,3:4)+2*urr1(i,1:2)-4*urr1(i,5:6)) ); % Gauss point
               uer2 = transpose( urr2(i,1:2) + ...
                      xg*(4*urr2(i,5:6)-3*urr2(i,1:2)-urr2(i,3:4)) + ...  
                      xg^2*(2*urr2(i,3:4)+2*urr2(i,1:2)-4*urr2(i,5:6)) );
               xgr  = ( [x1;y1] + xg*(4*[x3;y3]-3*[x1;y1]-[x2;y2]) + ...
                      xg^2*(2*[x2;y2]+2*[x1;y1]-4*[x3;y3]) );
            end
      
            % Reference force from the BC's
            if exno(1) == 1
               fer1 = [0;0]; fer2 = [fscalar;0];
            elseif exno(1) == -1
               fer1 = [0;0]; fer2 = -[fscalar;0];
            elseif exno(2) == 1
               fer1 = [0;fscalar]; fer2 = [0;0];
            elseif exno(2) == -1
               fer1 = -[0;fscalar]; fer2 = [0;0];
            end
            
            ixigrec = Q'*[xgr(1);xgr(2)];
            X = (ixigrec(1)-offset)/L-.5; Y = (ixigrec(2)+Cte2)/L;
            
            % Build X^k+1-j*Y^j
            GROX = zeros(ordp+2,1);
            for j = 0:k+1
               GROX(j+1) = X^(k+1-j)*Y^j;
            end
            
            vloc = [ coefa(:,k+1)'*GROX ; coefb(:,k+1)'*GROX ];
            vpa = Q*vloc;
            
            sloc11 = 0; sloc12 = 0; sloc22 = 0;
            for j=0:floor(k/2)
               sloc11 = sloc11 + ( E/(1-nu^2)*(k+1-2*j)*coefa(2*j+1,k+1) + ...
                        nu*E/(1-nu^2)*(2*j+1)*coefb(2*j+2,k+1) )*X^(k-2*j)*Y^(2*j);
               sloc22 = sloc22 + ( nu*E/(1-nu^2)*(k+1-2*j)*coefa(2*j+1,k+1) + ...
                        E/(1-nu^2)*(2*j+1)*coefb(2*j+2,k+1) )*X^(k-2*j)*Y^(2*j);
            end
            for j=1:floor((k+1)/2)
               sloc12 = sloc12 + ( E/(2*(1+nu))*(2*j)*coefa(2*j+1,k+1)* ...
                                   X^(k+1-2*j)*Y^(2*j-1) );
            end
            for j=0:floor((k-1)/2)
               sloc12 = sloc12 + ( E/(2*(1+nu))*(k-2*j)*coefb(2*j+2,k+1)* ...
                                   X^(k-1-2*j)*Y^(2*j+1) );
            end

            sloc = 1/L*[sloc11,sloc12;sloc12,sloc22];
            spa = Q*sloc*Q';
            fpa = spa*exno;

            Rhs(k+1) = Rhs(k+1) + len * wg * ( fer1'*vpa - fpa'*uer1 );
         end
      end
      
   end

%   L1 = offset; L2 = offset+L;
%   L1 = 0; L2 = 1;
   L1 = -.5; L2 = .5;
   for i=0:ordp
      for j=0:ordp
         ord = i+j+1;
         Lhs(i+1,j+1) = L*(L2^ord - L1^ord)/ord;
         if i>1 && j>1
            Lhs(i+1,j+1) = Lhs(i+1,j+1) + mu*i*j/(i+j-1)*...
                                          (L2^(i+j-1) - L1^(i+j-1));
         end
      end
   end
   
      % Invert the operator
   if dolcurve == 1
      listmu = -6:.1:-2;
      resid  = zeros( size(listmu) );
      regno  = zeros( size(listmu) );
      rhsk   = zeros(nbase+1,1);
      index  = 1;
      for lnmu1 = listmu
         Regop  = Lhs-Lhs;
         mu1 = 10^lnmu1;
         
         % Regularize
         for i=0:ordp
            for j=0:ordp
               if i>1 && j>1
                  Regop(i+1,j+1) = Regop(i+1,j+1) + mu1*i*j/(i+j-1)*...
                                                   (L2^(i+j-1) - L1^(i+j-1));
               end
            end
         end
         
         % Invert
         McCoef = (Lhs+Regop)\Rhs;
         resid(index)  = norm(Lhs*McCoef - Rhs)^2; % Actually, it's resid^2
         regno(index)  = McCoef'*(Regop/mu1)*McCoef;
         index = index+1;
      end
      figure;
      plot(resid,regno);
      xlabel('residual')
      ylabel('norm')
      %res10 = 10.^resid'; reg10 = 10.^regno';
      mu = 10^listmu( findCorner (resid', regno',2,0) );
   end
   
   % Regularize
   for i=0:ordp
      for j=0:ordp
         if i>1 && j>1
            Lhs(i+1,j+1) = Lhs(i+1,j+1) + mu*i*j/(i+j-1)*...
                                          (L2^(i+j-1) - L1^(i+j-1));
         end
      end
   end
   McCoef = Lhs\Rhs;
   nbase = size(McCoef,1);
   
   % Plot the result
   solref = [0*newXo ; ubase(:,2) ; 0*newXo];
   McCoefref = polyfit( newX, solref, nbase-1 );
   McCoefref = McCoefref(end:-1:1);
   
   solu = zeros(size(newX,1),1); solpref = zeros(size(newX,1),1);
   for i=1:nbase
      solu = solu + McCoef(i)*( (newX-offset)./L-.5).^(i-1);
   end
   for i=1:nbase
      solpref = solpref + McCoefref(i)*newX.^(i-1);
   end
   
   figure
   hold on;
   set(gca, 'fontsize', 20);
   plot(newX, solref,'linewidth',3)
   plot(newX, solu, 'Color', 'red','linewidth',3)
   %plot(newX, solpref, 'Color', 'green')
   legend('Reference gap','reconstructed gap')
   xlabel('X')
   ylabel('[[u]]')
   
   % Check the MS :
%   Sref = sum((solpref-solref).^2) / sum(solref.^2);
%   Smc  = sum((solu-solref).^2) / sum(solref.^2);
   %% Error computation
   
   ref = [0*newXo;ubase(:,2);0*newXo]; errorU = ref - solu;
   nelem = size(newX,1)-1;
   nerrnonrom = 0; nnormaliz = 0; % Well those are ugly names...
   for i=1:nelem
      dlen = newX(i+1)-newX(i);
      Me = dlen/3*[1,.5;.5,1];
      ue = errorU([i,i+1]); uer = ref([i,i+1]);
      nerrnonrom = nerrnonrom + ue'*Me*ue;
      nnormaliz  = nnormaliz  + uer'*Me*uer;
   end
   poly_error = nerrnonrom / nnormaliz;

end
