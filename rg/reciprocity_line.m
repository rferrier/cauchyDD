%13/10/16
%Détection de fissure 1D plane par écart à la réciprocité

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E          = 210000; % MPa : Young modulus
nu         = 0.3;    % Poisson ratio
fscalar    = 250;    % N.mm-1 : Loading on the plate
mat        = [0, E, nu];
mu         = 0;%1e-5;%5e-3;     % Regularization coef
%mu         = 3;     % Regularization coef
dolcurve   = 0;      % Do a L-curve or not

usefourier = 0;
usepolys   = 1;

nbase = 2; % Number of Fourier basis functions
ordp = 5;  % Number of Polynomial basis functions

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
plotGMSH({ux,'U_x';uy,'U_y';u1,'U_vect';sigma,'stress'}, elements, nodes, 'reference');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import the uncracked domain /!\ MUST BE THE SAME !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (except for the crack)
if useorder == 1
   [ nodes2,elements2,ntoelem2,boundary2,order2] = readmesh( 'meshes/rg_refined/plate_n.msh' );
elseif useorder == 2
   [ nodes2,elements2,ntoelem2,boundary2,order2] = readmesh( 'meshes/rg_refined/plate_nt6.msh' );
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

%% First determine the crack's line.

% Compute and apply v fields (for plane constraint)
v1 = zeros(2*nnodes2, 1);
v2 = zeros(2*nnodes2, 1);
v3 = zeros(2*nnodes2, 1);
for i=1:nnodes2
   x = nodes2(i,1);
   y = nodes2(i,2);
   v1(2*i-1) = 1/E*x;
   v2(2*i-1) = -nu/E*x;
   v3(2*i-1) = (1+nu)/(2*E)*y;
   v1(2*i)   = -nu/E*y;
   v2(2*i)   = 1/E*y;
   v3(2*i)   = (1+nu)/(2*E)*x;
end
% Debug : check sigma OK
%sigma1 = stress(v1,E,nu,nodes2,elements2,order,1,ntoelem2);
%sigma2 = stress(v2,E,nu,nodes2,elements2,order,1,ntoelem2);
%sigma3 = stress(v3,E,nu,nodes2,elements2,order,1,ntoelem2);
%plotGMSH({sigma1,'sigma1';sigma2,'sigma2';sigma3,'sigma3'}, elements2, nodes2, 'sigmas');
%
f1 = Kinter2*v1;
f2 = Kinter2*v2;
f3 = Kinter2*v3;

% Clean redondant stuff in indexbound2
i = 1;
while i <= size(indexbound2,1)
   if find( indexbound2(i) == indexbound2(1:i-1) )
      indexbound2(i) = [];
      i = i-1;
   end
   i = i+1;
end

R11 = (fr1(indexbound2)'*v1(indexbound2) - f1(indexbound2)'*ur1(indexbound2));
R21 = (fr1(indexbound2)'*v2(indexbound2) - f2(indexbound2)'*ur1(indexbound2));
R31 = (fr1(indexbound2)'*v3(indexbound2) - f3(indexbound2)'*ur1(indexbound2));
%[R11s, R31s ; R31s, R21s]

R12 = (fr2(indexbound2)'*v1(indexbound2) - f1(indexbound2)'*ur2(indexbound2));
R22 = (fr2(indexbound2)'*v2(indexbound2) - f2(indexbound2)'*ur2(indexbound2));
R32 = (fr2(indexbound2)'*v3(indexbound2) - f3(indexbound2)'*ur2(indexbound2));

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

% Plot the mesh with the normal
figure
hold on;
ret = patch('Faces',elements2(:,1:3),'Vertices',nodes2,'FaceAlpha',0);
x6 = nodes(6,1); y6 = nodes(6,2); x5 = nodes(5,1); y5 = nodes(5,2); 
xc = .5*(max(nodes(:,1)) + min(nodes(:,1)));
yc = .5*(max(nodes(:,2)) + min(nodes(:,2)));

x1 = xc + dir11(1); y1 = yc + dir11(2);
x2 = xc + dir21(1); y2 = yc + dir21(2);
xn = xc + n1; yn = yc + n2;
plot( [xc,x1], [yc,y1] ,'Color', 'red', 'LineWidth',3);
plot( [xc,x2], [yc,y2] ,'Color', 'green', 'LineWidth',3);
plot( [xc,xn], [yc,yn] ,'Color', 'cyan', 'LineWidth',2);
plot( [x5,x6], [y5,y6] ,'Color', 'magenta', 'LineWidth',3);
axis('equal');

% And the other one
figure
hold on;
ret = patch('Faces',elements2(:,1:3),'Vertices',nodes2,'FaceAlpha',0);
x6 = nodes(6,1); y6 = nodes(6,2); x5 = nodes(5,1); y5 = nodes(5,2); 
xc = .5*(max(nodes(:,1)) + min(nodes(:,1)));
yc = .5*(max(nodes(:,2)) + min(nodes(:,2)));

x1 = xc + dir12(1); y1 = yc + dir12(2);
x2 = xc + dir22(1); y2 = yc + dir22(2);
xn = xc + n1; yn = yc + n2;
plot( [xc,x1], [yc,y1] ,'Color', 'red', 'LineWidth',3);
plot( [xc,x2], [yc,y2] ,'Color', 'green', 'LineWidth',3);
plot( [xc,xn], [yc,yn] ,'Color', 'cyan', 'LineWidth',2);
plot( [x5,x6], [y5,y6] ,'Color', 'magenta', 'LineWidth',3);
axis('equal');

%% Then the constant

% First, find the minimal point (in order to have Cte < 0)
norep = Q'*nodes'; K = min(norep(2,:));

vt = zeros(2*nnodes2, 1);
%vn = zeros(2*nnodes2, 1);
Xs = zeros(nnodes2, 1);
Ys = zeros(nnodes2, 1);
for i=1:nnodes2
   x = nodes2(i,1);
   y = nodes2(i,2);
   % Change base (for coordiantes)
   ixigrec = Q'*[x;y]; X = ixigrec(1); Y = ixigrec(2)-K;
   Xs(i) = X; Ys(i) = Y;

   vloc = [ -X^2/(2*E) + (2+nu)*Y^2/(2*E) ; nu*X*Y/E ];
   % Change base (for vector), in the other direction
   vxy = Q*vloc; vt(2*i-1) = vxy(1); vt(2*i) = vxy(2);
end

% Norm of [[ut]] (case 1)
normT  = sqrt( abs( 2*(R11^2+R21^2+2*R31^2) - 2*(R11+R21)^2 ) );
normT2 = sqrt( abs( 2*(R12^2+R22^2+2*R32^2) - 2*(R12+R22)^2 ) );

%% Debug : check sigma /!\ it's not straingtforward because diagonalization /!\
%sigmat = stress(vt,E,nu,nodes2,elements2,order,1,ntoelem2);
%scalN = sparse(2*nnodes2, 3*nnodes2);
%scalN( 1:2:2*nnodes2-1 , 1:3:3*nnodes2-2 ) = normal(1);
%scalN( 1:2:2*nnodes2-1 , 3:3:3*nnodes2 ) = normal(2);
%scalN( 2:2:2*nnodes2 , 2:3:3*nnodes2-1 ) = normal(2);
%scalN( 2:2:2*nnodes2 , 3:3:3*nnodes2 ) = normal(1);
%sigNt = scalN*sigmat;
%plotGMSH({sigmat(1:3:3*nnodes2-2),'sigmat';...
%                               Xs,'X';Ys,'Y'}, elements2, nodes2, 'sigmas');

%fn = Kinter2*vn;
ft = Kinter2*vt;

%Rn = fr1(indexbound2)'*vn(indexbound2) - fn(indexbound2)'*ur1(indexbound2);
Rt = (fr1(indexbound2)'*vt(indexbound2) - ft(indexbound2)'*ur1(indexbound2));
%Rn2 = fr2(indexbound2)'*vn(indexbound2) - fn(indexbound2)'*ur2(indexbound2);
Rt2 = (fr2(indexbound2)'*vt(indexbound2) - ft(indexbound2)'*ur2(indexbound2));

%Cte  = 1/normR1*sqrt(Rn^2 + Rt^2);
%Cte2 = 1/normR2*sqrt(Rn2^2 + Rt2^2); % DEBUG
Cte  = min( Rt/normT, -Rt/normT) - K;       % Select the negative one
Cte2 = min( Rt2/normT2, -Rt2/normT2 ) - K;  %  /!\ The sign depends on the test case
% Plot the crack, and its estimated lines (there are Y +/- Cte)
figure
hold on;
ret = patch('Faces',elements2(:,1:3),'Vertices',nodes2,'FaceAlpha',0);
x6 = nodes(6,1); y6 = nodes(6,2); x5 = nodes(5,1); y5 = nodes(5,2); 

Vp1 = [10;-Cte]; Vp2 = [-10;-Cte];
Vm1 = [10;-Cte2]; Vm2 = [-10;-Cte2];
vp1 = Q*Vp1; vm1 = Q*Vm1; vp2 = Q*Vp2; vm2 = Q*Vm2;

plot( [vp1(1), vp2(1)], [vp1(2), vp2(2)] ,'Color', 'red', 'LineWidth',1);
plot( [vm1(1), vm2(1)], [vm1(2), vm2(2)] ,'Color', 'green', 'LineWidth',1);
plot( [x5,x6], [y5,y6] ,'Color', 'magenta', 'LineWidth',2);
axis('equal');

%% And now, the crack itself

% Compute L (provisionnal formula assuming that something<Cte<0
tangent = [normal(2) ; -normal(1)]; b = max(nodes(:,1)) - min(nodes(:,1));
a = tangent(2)/tangent(1)*b; L = sqrt(a^2+b^2);

% Convenient values for plot
udepx  = u1(icrack5x)-u1(icrack6x);
udepy  = u1(icrack5y)-u1(icrack6y);
Xx     = nodes(b2node5,1); Yy = nodes(b2node5,2);
ubase  = Q'*[udepx,udepy]'; ubase = ubase';
XY     = Q'*[Xx,Yy]'; XY = XY';
[newX, orderX] = sort(XY(:,1));          % In case order > 1, need to sort
offset = Q(1,2)/Q(1,1)*CteR;   % /!\ Only in case rectangle

left  = offset : (-offset+newX(1))/(size(newX,1)) : newX(1) ;
right = newX(end) : (offset + L - newX(end))/(size(newX,1)) : offset + L ;
newXo = newX;
newX  = [left(1:end-1)';
         newX;
         right(2:end)']; % Add a few points

if usefourier == 1
   Rp     = zeros(nbase+1,1);
   Rm     = zeros(nbase+1,1);
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
         for i=1:nnodes2
            x = nodes2(i,1);
            y = nodes2(i,2);
            % Change base (for coordinates)
            ixigrec = Q'*[x;y]; X = ixigrec(1); Y = ixigrec(2)+Cte;
            Xs(i) = X; Ys(i) = Y;
            
            v1 = -I*lambda(kp)*exp(-I*lambdae*X)* ...
                               ( exp(lambda(kp)*Y)+exp(-lambda(kp)*Y) );
            v2 = lambda(kp)*exp(-I*lambdae*X)* ...
                               ( exp(lambda(kp)*Y)-exp(-lambda(kp)*Y) );
            vloc = [ v1 ; v2 ];
            % Change base (for vector), in the other direction
            vxy = Q*vloc; vp(2*i-1) = vxy(1); vp(2*i) = vxy(2);
         end
         fp = Kinter2*vp;
         
         % Fourier coefficient
         if sx == 1
            Rp(kp) = (fr1(indexbound2)'*vp(indexbound2) - fp(indexbound2)'*ur1(indexbound2));
            %Rp(kp) = (fr2(indexbound2)'*vp(indexbound2) - fp(indexbound2)'*ur2(indexbound2));
            if dolcurve == 0 % No L-curve stuff
               fourn(kp) = - 1/(1+mu*k^2) *(1+nu)/(2*E*L*lambda(kp)^2)*Rp(kp);
               akan(k) = 2*real(fourn(kp));
               bkan(k) = 2*imag(fourn(kp));
            end
         else
            Rpm = (fr1(indexbound2)'*vp(indexbound2) - fp(indexbound2)'*ur1(indexbound2));
            fournm(kp) = -(1+nu)/(2*E*L*lambda(kp)^2)*Rpm;
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
   solu = fourn(1) + sum ( [0*newX' ;...  % Hack in case there is 1 element only
              akan.*cos(lambda(2:end)*newX') + bkan.*sin(lambda(2:end)*newX') ] );
   solu = solu';
   
   figure
   hold on;
   plot(newX, [0*newXo;ubase(:,2);0*newXo])
   plot(newX, solu, 'Color', 'red')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mean square polynoms
if usepolys == 1
   % First, build the polynomial test functions.
   coef = zeros(ordp+2, ordp+1);
   coef( 1:2, 1 ) = [-(1-nu^2)/(nu*E) ; 0];
   
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
   
      C = [ eye(2) , zeros(2,k)];
      Lhsco = [ Lhsco ; C ]; % find an unique solution
      Rhsco = [ Rhsco ; azero ; 0 ];
      
      Lhsco(size(Lhsco,1)-2,:) = [];  % 'cause square matrix is life
      Rhsco(size(Rhsco,1)-2,:) = [];
      
      coef( 1:k+2, k+1 ) = Lhsco\Rhsco;
   %   Rhsco = [ Rhsco ; azero ; 0 ]  %impose a0 = a0 and b0 = 0 with Lagrange
   %   C = [ eye(2) , zeros(2,k-1)];
   %   Lhsco = [ Lhsco , C' ; C , zeros(2,2) ]
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
   Rhs = zeros(ordp+1,1);
   vpa = zeros(2*nnodes2, 1);
   for k=0:ordp
      for i=1:nnodes2
         x = nodes2(i,1);
         y = nodes2(i,2);
         ixigrec = Q'*[x;y]; X = (ixigrec(1)-offset)/L; Y = (ixigrec(2)+Cte)/L;
         
         % Build X^k*Y^j
         GROX = zeros(ordp+2,1);
         for j = 0:k+1
            GROX(j+1) = X^(k+1-j)*Y^j;
         end
         
         vloc = [ coefa(:,k+1)'*GROX ; coefb(:,k+1)'*GROX ];
         vxy = Q*vloc; vpa(2*i-1) = vxy(1); vpa(2*i) = vxy(2);
      end
      fpa = Kinter2*vpa;
      %Rhs(k+1) = (fr1'*vpa - fpa'*ur1);
      Rhs(k+1) = (fr1(indexbound2)'*vpa(indexbound2) - fpa(indexbound2)'*ur1(indexbound2));
   end
   
   L1 = offset; L2 = offset+L;
   L1 = 0; L2 = 1;
   for i=0:ordp
      for j=0:ordp
         ord = i+j+1;
         Lhs(i+1,j+1) = (L2^ord - L1^ord)/ord;
         if i>1 && j>1
            Lhs(i+1,j+1) = Lhs(i+1,j+1) + mu*i*j/(i+j-1)*...
                                          (L2^(i+j-1) - L1^(i+j-1));
         end
      end
   end
   
   %% Manual stuff for DEBUG
   
   %vp3 = zeros(2*nnodes2, 1);
   %vp2 = zeros(2*nnodes2, 1);
   %vp1 = zeros(2*nnodes2, 1);
   %vp0 = zeros(2*nnodes2, 1);
   %
   %for i=1:nnodes2
   %   x = nodes2(i,1);
   %   y = nodes2(i,2);
   %   % Change base (for coordiantes)
   %   ixigrec = Q'*[x;y]; X = ixigrec(1); Y = ixigrec(2)+Cte;
   %
   %   %vloc = [nu/E*X ; -1/E*Y];
   %   vloc = -[ (1-nu^2)/(nu*E)*X ; 0 ];
   %   vxy = Q*vloc; vp0(2*i-1) = vxy(1); vp0(2*i) = vxy(2);
   %   
   %   %vloc = [ - (1+nu^2)/(2*nu*E)*X^2 + (1+nu)/(nu*E)*Y^2 ; 0 ];
   %   vloc = -[ (1-nu^2)/(2*nu*E)*X^2 - (1+nu)/(nu*E)*Y^2 ; 0 ];
   %   vxy = Q*vloc; vp1(2*i-1) = vxy(1); vp1(2*i) = vxy(2);
   %   
   %   a1 = (1-nu^2)/(3*nu*E); b1 = -2*(1+nu)/(nu*E);
   %   c1 = -(1-nu^2)/(6*E)*( 2*nu*E/(1-nu^2) + E/(1+nu) )*b1;
   %   vloc = -[ a1*X^3 + b1*X*Y^2 ; c1*Y^3 ];
   %   vxy = Q*vloc; vp2(2*i-1) = vxy(1); vp2(2*i) = vxy(2);
   %   
   %   a1 = (1-nu^2)/(4*nu*E); b1 = -3*(1+nu)/(nu*E);
   %   d1 = -(1-nu^2)/(6*E)*( 4*nu*E/(1-nu^2) + 2*E/(1+nu) )*b1;
   %   c1 = -(1+nu)/(6*E) * ( 2*E/(1-nu^2)*b1 + 3*nu*E/(1-nu^2)*d1 + 3*E/(2*(1+nu))*d1 );
   %   vloc = -[ a1*X^4 + b1*X^2*Y^2 + c1*Y^4 ;...
   %            d1*X*Y^3 ];
   %   vxy = Q*vloc; vp3(2*i-1) = vxy(1); vp3(2*i) = vxy(2);
   %end
   %fp3 = Kinter2*vp3;
   %fp2 = Kinter2*vp2;
   %fp1 = Kinter2*vp1;
   %fp0 = Kinter2*vp0;
   %
   %Rp3 = (fr1(indexbound2)'*vp3(indexbound2) - fp3(indexbound2)'*ur1(indexbound2));
   %Rp2 = (fr1(indexbound2)'*vp2(indexbound2) - fp2(indexbound2)'*ur1(indexbound2));
   %Rp1 = (fr1(indexbound2)'*vp1(indexbound2) - fp1(indexbound2)'*ur1(indexbound2));
   %Rp0 = (fr1(indexbound2)'*vp0(indexbound2) - fp0(indexbound2)'*ur1(indexbound2));
   %
   %%% /!\ savage zone /!\
   %%Xs(100)
   %%Ys(100)
   %%sigm = stress(vp1,E,nu,nodes2,elements2,order,1,ntoelem2);
   %%S = [sigm(3*100-2),sigm(300);sigm(300),sigm(3*100-1)]
   %%Q'*S*Q
   %
   %Rhs =  [Rp0 ;Rp1 ; Rp2 ; Rp3];
   %L1 = offset; L2 = offset+L;
   %for i=1:4
   %   for j=1:4
   %      ord = i+j-1;
   %      Lhs(i,j) = (L2^ord - L1^ord)/ord;
   %   end
   %end
   %%Lhs = [L2-L1, (L2^2 - L1^2)/2, (L2^3 - L1^3)/3 ;...
   %%       (L2^2 - L1^2)/2, (L2^3 - L1^3)/3, (L2^4 - L1^4)/4 ;...
   %%       (L2^3 - L1^3)/3, (L2^4 - L1^4)/4, (L2^5 - L1^5)/5 ];
   %%Rhs = -[Rp0]; 
   %%Lhs = [L];
   %%Rhs = Rhs(1:3); Lhs = Lhs(1:3,1:3);
   
      % Invert the operator
   if dolcurve == 1
      listmu = -6:.1:-2;
      resid  = zeros( size(listmu) );
      regno  = zeros( size(listmu) );
      rhsk   = zeros(nbase+1,1);
      index = 1;
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
      solu = solu + McCoef(i)*( (newX-offset)./L).^(i-1);
   end
   for i=1:nbase
      solpref = solpref + McCoefref(i)*newX.^(i-1);
   end
   
   figure
   hold on;
   plot(newX, solref)
   plot(newX, solu, 'Color', 'red')
   %plot(newX, solpref, 'Color', 'green')
   
   % Check the MS :
   Sref = sum((solpref-solref).^2) / sum(solref.^2);
   Smc  = sum((solu-solref).^2) / sum(solref.^2);
end