%02/11/2016
%Détection de fissure 3D plane par écart à la réciprocité

close all;
clear all;

addpath(genpath('./tools'))
addpath(genpath('/usr/share/octave/packages/geometry-2.1.0')) % Geomptry pkg

% Parameters
E       = 210000; % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 250;    % N.mm-1 : Loading on the plate
mat = [0, E, nu];

nbase     = 2; % Number of Fourier basis functions
loadfield = 1; % If 0 : recompute the reference problem and re-pass mesh

if loadfield == 0
   % Boundary conditions
   % first index  : index of the boundary
   % second index : 1=x, 2=y
   % third        : value
   % [0,1,value] marks a dirichlet regularization therm on x
   dirichlet  = [0,1,0 ; 0,2,0 ; 0,3,0 ; 0,4,0 ; 0,5,0 ; 0,6,0];
   neumann1   = [2,3,fscalar ; 1,3,-fscalar];
   neumann2   = [4,1,fscalar ; 6,1,-fscalar];
   
   % First, import the mesh
   [ nodes,elements,ntoelem,boundary,order] = readmesh3D( 'meshes/plate3d_crack0.msh' );
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
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import the uncracked domain /!\ MUST BE THE SAME !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! (except for the crack)
[ nodes2,elements2,ntoelem2,boundary2,order2] = readmesh3D( 'meshes/plate3d.msh' );
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
   UFr = passMesh3D (nodes, elements, nodes2, elements2, [u1,u2,f1,f2]);
   ur1 = UFr(:,1); ur2 = UFr(:,2); fr1 = UFr(:,3); fr2 = UFr(:,4);
   plotGMSH3D({ur1,'u1';ur2,'u2'}, elements2, nodes2, 'us');
else
   %load the field
   UFR = load('fields/UFr.mat'); UFr = UFR.UFr;
   ur1 = UFr(:,1); ur2 = UFr(:,2); fr1 = UFr(:,3); fr2 = UFr(:,4);
   plotGMSH3D({ur1,'u1';ur2,'u2'}, elements2, nodes2, 'us');
end

%fr1p = Kinter2*ur1; fr2p = Kinter2*ur2;  % This is wrong inside the domain,
% but should be right at its boundary provided the crack doesn't touch the boundary
% Well no, actually, ti is wrong everywhere.

neumann1   = [2,3,fscalar ; 1,3,-fscalar];
neumann2   = [4,1,fscalar ; 6,1,-fscalar];
fr1 = loading3D(nbloq2,nodes2,boundary2,neumann1);
fr2 = loading3D(nbloq2,nodes2,boundary2,neumann2);

%% First determine the crack's normal.

% Compute and apply v fields (for plane constraint)
v1 = zeros(3*nnodes2, 1);
v2 = zeros(3*nnodes2, 1);
v3 = zeros(3*nnodes2, 1);
v4 = zeros(3*nnodes2, 1);
v5 = zeros(3*nnodes2, 1);
v6 = zeros(3*nnodes2, 1);

for i=1:nnodes2
   x = nodes2(i,1);
   y = nodes2(i,2);
   z = nodes2(i,3);
   
   v1(3*i-2) = 1/E*x;
   v2(3*i-2) = -nu/E*x;
   v3(3*i-2) = -nu/E*x;
   v1(3*i-1) = -nu/E*y;
   v2(3*i-1) = 1/E*y;
   v3(3*i-1) = -nu/E*y;
   v1(3*i)   = -nu/E*z;
   v2(3*i)   = -nu/E*z;
   v3(3*i)   = 1/E*z;
   
   v4(3*i-2) = (1+nu)/(2*E)*y;
   v4(3*i-1) = (1+nu)/(2*E)*x;
   v5(3*i-2) = (1+nu)/(2*E)*z;
   v5(3*i)   = (1+nu)/(2*E)*x;
   v6(3*i-1) = (1+nu)/(2*E)*z;
   v6(3*i)   = (1+nu)/(2*E)*y;
end
% Debug : check sigma
%sigma1 = stress3D(v1,mat,nodes2,elements2,order,1,ntoelem2);
%sigma2 = stress3D(v2,mat,nodes2,elements2,order,1,ntoelem2);
%sigma3 = stress3D(v3,mat,nodes2,elements2,order,1,ntoelem2);
%plotGMSH3D({sigma1,'sigma1';sigma2,'sigma2';sigma3,'sigma3'}, elements2, nodes2, 'sigmas');
%
f1 = Kinter2*v1;
f2 = Kinter2*v2;
f3 = Kinter2*v3;
f4 = Kinter2*v4;
f5 = Kinter2*v5;
f6 = Kinter2*v6;

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
R41 = (fr1(indexbound2)'*v4(indexbound2) - f4(indexbound2)'*ur1(indexbound2));
R51 = (fr1(indexbound2)'*v5(indexbound2) - f5(indexbound2)'*ur1(indexbound2));
R61 = (fr1(indexbound2)'*v6(indexbound2) - f6(indexbound2)'*ur1(indexbound2));
RR1 = [R11, R41, R51 ; R41, R21, R61 ; R51, R61, R31];

R12 = (fr2(indexbound2)'*v1(indexbound2) - f1(indexbound2)'*ur2(indexbound2));
R22 = (fr2(indexbound2)'*v2(indexbound2) - f2(indexbound2)'*ur2(indexbound2));
R32 = (fr2(indexbound2)'*v3(indexbound2) - f3(indexbound2)'*ur2(indexbound2));
R42 = (fr2(indexbound2)'*v4(indexbound2) - f4(indexbound2)'*ur2(indexbound2));
R52 = (fr2(indexbound2)'*v5(indexbound2) - f5(indexbound2)'*ur2(indexbound2));
R62 = (fr2(indexbound2)'*v6(indexbound2) - f6(indexbound2)'*ur2(indexbound2));
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

% Plot the mesh with the normal
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

% Norm of [[ut]] (case 1)
normT  = sqrt( abs( 2*(R11^2+R21^2+R31^2+2*(R41^2+R51^2+R61^2)) ...
                    - 2*(R11+R21+R31)^2 ) );
normT2 = sqrt( abs( 2*(R12^2+R22^2+R32^2+2*(R42^2+R52^2+R62^2)) ...
                    - 2*(R12+R22+R32)^2 ) );

%% Debug : check sigma /!\ it's not straingtforward because diagonalization /!\
%sigmat = stress3D(vt,mat,nodes2,elements2,order2,1,ntoelem2);
%sigmav = stress3D(vv,mat,nodes2,elements2,order2,1,ntoelem2);
%scalN = sparse(2*nnodes2, 3*nnodes2);
%scalN( 1:2:2*nnodes2-1 , 1:3:3*nnodes2-2 ) = normal(1);
%scalN( 1:2:2*nnodes2-1 , 3:3:3*nnodes2 ) = normal(2);
%scalN( 2:2:2*nnodes2 , 2:3:3*nnodes2-1 ) = normal(2);
%scalN( 2:2:2*nnodes2 , 3:3:3*nnodes2 ) = normal(1);
%sigNt = scalN*sigmat;
%plotGMSH({sigmat(1:3:3*nnodes2-2),'sigmat';...
%                               Xs,'X';Ys,'Y'}, elements2, nodes2, 'sigmas');

ft = Kinter2*vt;
fv = Kinter2*vv;

Rt  = (fr1(indexbound2)'*vt(indexbound2) - ft(indexbound2)'*ur1(indexbound2));
Rv  = (fr1(indexbound2)'*vv(indexbound2) - fv(indexbound2)'*ur1(indexbound2));
Rt2 = (fr2(indexbound2)'*vt(indexbound2) - ft(indexbound2)'*ur2(indexbound2));
Rv2 = (fr2(indexbound2)'*vv(indexbound2) - fv(indexbound2)'*ur2(indexbound2));

Cte  = -sqrt(Rt^2+Rv^2)/normT - K; % K was chosen so that Cte is negative
Cte2 = -sqrt(Rt2^2+Rv2^2)/normT2 - K;

% Reference : we know that the point P belongs to the plane.
Pt = [4;3;1]; QPt = Q'*Pt; CteR = QPt(3);

%% Plot the crack, and its estimated plane
%figure
%hold on;
%ret = patch('Faces',elements2(:,1:3),'Vertices',nodes2,'FaceAlpha',0);
%x6 = nodes(6,1); y6 = nodes(6,2); x5 = nodes(5,1); y5 = nodes(5,2); 
%
%Vp1 = [10;-Cte]; Vp2 = [-10;-Cte];
%Vm1 = [10;-Cte2]; Vm2 = [-10;-Cte2];
%vp1 = Q*Vp1; vm1 = Q*Vm1; vp2 = Q*Vp2; vm2 = Q*Vm2;
%
%plot( [vp1(1), vp2(1)], [vp1(2), vp2(2)] ,'Color', 'red', 'LineWidth',1);
%plot( [vm1(1), vm2(1)], [vm1(2), vm2(2)] ,'Color', 'green', 'LineWidth',1);
%plot( [x5,x6], [y5,y6] ,'Color', 'magenta', 'LineWidth',2);
%axis('equal');

%% And now, the crack itself

% Compute L ( (not so) provisionnal formula assuming that we have a particular case)
b = 7;
a = t(3)/t(1)*b; Lx = sqrt(a^2+b^2); Ly = b;
%
Rp       = zeros(nbase+1,1);
Rm       = zeros(nbase+1,1);
lambdax  = zeros(nbase+1,1);
lambday  = zeros(nbase+1,1);
%lambda  = zeros(nbase+1,1);
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
      %vm = zeros(3*nnodes2,1);
      lambdax(kpx) = 2*kx*pi/Lx; lambday(kpy) = 2*ky*pi/Ly;
      lambda = sqrt(lambdax(kpx)^2+lambday(kpy)^2);
      
      for sx = [1,-1]
         lambdaxe = sx*lambdax(kpx);
         for sy = [1,-1]
            lambdaye = sy*lambday(kpy);
            for i=1:nnodes2
               x = nodes2(i,1);
               y = nodes2(i,2);
               z = nodes2(i,3);
               % Change base (for coordiantes)
               ixigrec = Q'*[x;y;z]; X = ixigrec(1); Y = ixigrec(2); Z = ixigrec(3)+Cte;
               Xs(i) = X; Ys(i) = Y; Zs(i) = Z;
               
               v1 = -I*lambdaxe*exp(-I*lambdaxe*X-I*lambdaye*Y)* ...
                                  ( exp(lambda*Z)+exp(-lambda*Z) );
               v2 = -I*lambdaye*exp(-I*lambdaxe*X-I*lambdaye*Y)* ...
                                  ( exp(lambda*Z)+exp(-lambda*Z) );
               v3 = lambda*exp(-I*lambdaxe*X-I*lambdaye*Y)* ...
                                  ( exp(lambda*Z)-exp(-lambda*Z) );
               vloc = [ v1 ; v2 ; v3 ];
               % Change base (for vector), in the other direction
               vxy = Q*vloc; vp(3*i-2) = vxy(1); vp(3*i-1) = vxy(2); vp(3*i) = vxy(3);
            end
            fp = Kinter2*vp;
      
            Rp = (fr1(indexbound2)'*vp(indexbound2) - fp(indexbound2)'*ur1(indexbound2));
            
            % Fourier coefficient
            if sx==1 && sy==1
               fournpp(kpx,kpy) = -(1+nu)/(2*E*Lx*Ly*lambda^2)*Rp;
            elseif sx==1 && sy==-1
               fournpm(kpx,kpy) = -(1+nu)/(2*E*Lx*Ly*lambda^2)*Rp;
            elseif sx==-1 && sy==1
               fournmp(kpx,kpy) = -(1+nu)/(2*E*Lx*Ly*lambda^2)*Rp;
            elseif sx==-1 && sy==-1
               fournmm(kpx,kpy) = -(1+nu)/(2*E*Lx*Ly*lambda^2)*Rp;
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
      dkan(kx,ky)   = real(-fournpp(kpx,kpy) + fournpm(kpx,kpy) +...
                            fournmp(kpx,kpy) - fournmm(kpx,kpy));
   end
end

%% DEBUG :
%ui = reshape(imag(vp),2,[])';  ux = ui(:,1);  uy = ui(:,2);
plotGMSH({ Xs,'x' ; Ys,'y' ; Zs,'z' ; real(vp),'U_vect'}, elements2, nodes2, 'test field');

% The constant term
fournpp(1,1) = -(R11+R21+R31)/(Lx*Ly);

% plot the identified normal gap
X = 0:Lx/100:Lx; Y = 0:Ly/100:Ly; 
solu = fournpp(1,1);
for kpx=2:nbase+1
   for kpy=2:nbase+1
      solu = solu + akan(kpx-1,kpy-1)*cos(lambdax(kpx)*X')*cos(lambday(kpy)*Y)...
                  + bkan(kpx-1,kpy-1)*cos(lambdax(kpx)*X')*sin(lambday(kpy)*Y)...
                  + ckan(kpx-1,kpy-1)*sin(lambdax(kpx)*X')*cos(lambday(kpy)*Y)...
                  + dkan(kpx-1,kpy-1)*sin(lambdax(kpx)*X')*sin(lambday(kpy)*Y);
   end
end

% Center of the Circle
Cc = [4;3;1]; Cc = Q'*Cc;

figure;
hold on;
surf(X,Y,solu);
colorbar();
drawCircle ( Cc(1), Cc(1), 2, 'Color', 'black', 'LineWidth', 3 );