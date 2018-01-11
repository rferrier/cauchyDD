% 21/04/2017
% Détection de source par écart à la réciprocité
% Intégrations par PG

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
jmax       = 15;     % Eigenvalues truncation number
energy     = 0;      % Use the energy matrix in the GSVD

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
Lx = Xmax-Xmin; Ly = Ymax-Ymin;

% Rigid body modes
r1 = zeros(2*nnodes,1); r2 = r1; r3 = r1;
ind = 2:2:2*nnodes;
r1(ind-1,1) = 1; r2(ind,1) = 1;
moyx = (Xmax-Xmin)/2; moyy = (Ymax-Ymin)/2;
r3(ind-1,1) = -nodes(ind/2,2)+moyy;
r3(ind,1) = nodes(ind/2,1)-moyx;
R = [r1,r2,r3];

% Polynomial, zero mean, loading
% loadV = [ -fscalar/Ly*Ymoy,0,0,0; ...
%           fscalar/Ly,0,0,0;
%           0,0,0,0; ...
%           0,0,0,0];
% loadV = [fscalar/Ly*Ymoy - fscalar/Ly^2*(Ymax^3-Ymin^3) / (3*(Ymax-Ymin)) ,0,0,0; ...
%          -fscalar/Ly,0,0,0; ...
%          fscalar/(Ly^2),0,0,0
%          0,0,0,0];
% loadV = [-fscalar/Ly^2*(Ymax^3-Ymin^3) / (3*(Ymax-Ymin)) ,0,0,0; ...
%          0,0,0,0; ...
%          fscalar/(Ly^2),0,0,0
%          0,0,0,0];
loadV = [-fscalar/Ly^3*(Ymax^4-Ymin^4) / (4*(Ymax-Ymin)) ,0,0,0; ...
         0,0,0,0; ...
         0,0,0,0; ...
         fscalar/Ly^3,0,0,0];
f = volumicLoad( 3, nodes, elements, 2, loadV ); 

% % N ponctual sources
% f = zeros(2*nnodes+nbloq,1);
nnind = zeros(1,1);
% for i=1:1
%     n = randi(2*nnodes);
%     f( n ) = 1;
%     nnind(i) = n;
% end
% % Equilibrate f if needed
% f(1:2*nnodes) = f(1:2*nnodes) - R*( (R'*R) \ (R'*f(1:2*nnodes)) );

MustbeZero = norm(R'*f(1:2*nnodes))

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

disp([ 'Direct problem solved and data management ', num2str(toc) ]);
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Identify the volumic forces
% Build the polynomial test functions.
load('conditions10_2d.mat','-ascii');
M = spconvert(conditions10_2d); clear('conditions10_2d');
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

      X = xgr(1)/Lx; Y = xgr(2)/Lx;
      
      for i=1:nftest
         coefa = coef(1:(nmax+1)^2,i);
         coefb = coef((nmax+1)^2+1:2*(nmax+1)^2,i);
         for ii=0:nmax
            for jj=0:nmax
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
disp([ 'Left hand side generated ', num2str(toc) ]);
tic
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
         x3  = nodes(no3,1); y3 = nodes(no3,2);
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

         Rhs(k) = Rhs(k) + len * wg * ( fer'*vpa - fpa'*uer );
      end
   end
   
end
disp([ 'Right hand side generated ', num2str(toc) ]);
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve the linear system and recover the volumic loading
A = Lhs'*Lhs; b = -Lhs'*Rhs; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PB with (-)

%% Build the force energy matrix

% Matrix such that Red*fe = fn with fe, per elem force and fn per node
% force.
Red = zeros(2*nnodesu,2*nelemu);
for i=1:size(elementsu,1)
   Xloc = nodesu(elementsu(i,:),:);
   nno = size(Xloc,1);

   ne = size(elementsu,2);
   for j=1:ne
      nod = elementsu(i,j);
      Red( 2*nod-1, 2*i-1 ) = Red( 2*nod-1, 2*i-1 ) + 1/ntoelemu(nod,1);
      Red( 2*nod, 2*i ) = Red( 2*nod, 2*i ) + 1/ntoelemu(nod,1);
   end
end

if energy == 0
   L = eye(size(A));
else
%     % Rigid body modes
%     r1u = zeros(2*nnodesu,1); r2u = r1u; r3u = r1u;
%     ind = 2:2:2*nnodesu;
%     r1u(ind-1,1) = 1; r2u(ind,1) = 1;
%     moyx = (Xmax-Xmin)/2;moyy = (Ymax-Ymin)/2;
%     r3u(ind-1,1) = -nodesu(ind/2,2)+moyy;
%     r3u(ind,1) = nodesu(ind/2,1)-moyx;
%     Ru = [r1u,r2u,r3u];
%     Ru = Red'*Ru; % Pass to elements
    
    % Energy inner product : Ep = 1/2 fe'*L*fe
    [K1,~,~,~,~] = Krig2 (nodesu,elementsu,mat,orderu,boundaryu,dirichlet);
    KtRed = K1\[ Red ; zeros(3,size(Red,2)) ]; %size(KtRed)
    KtRed = KtRed( 1:2*nnodesu, : );
    L = Red'*KtRed;
    
%     L1 = [L,Ru;Ru',zeros(3)];
end

% % Rigid body modes
% r1u = zeros(2*nnodesu,1); r2u = r1u; r3u = r1u;
% ind = 2:2:2*nnodesu;
% r1u(ind-1,1) = 1; r2u(ind,1) = 1;
% moyx = (Xmax-Xmin)/2;moyy = (Ymax-Ymin)/2;
% r3u(ind-1,1) = -nodesu(ind/2,2)+moyy;
% r3u(ind,1) = nodesu(ind/2,1)-moyx;
% Ru = [r1u,r2u,r3u];
    
% A = Red*A*Red'; L = Red*L*Red'; b = Red*b;% Pass everything on the nodes
% L1 = [L,1e-7*Ru;1e-7*Ru',zeros(3)];
% A1 = [A,zeros(size(A,1),3) ; zeros(3,size(A,1)),eye(3)];
% b1 = [b;zeros(3,1)];

[Q,Theta,P] = svd(A);
thetas = diag(Theta);
% figure; plot(log10(diag(Theta)));
[thetas,Ind] = sort( thetas,'descend' );
Q = Q(:,Ind); P = P(:,Ind);
Thetas = diag(thetas); 
Theta = Theta(Ind,Ind);

%% A = Q*Theta*V'; L = P*aiS*V'
%[Q,P,V,Theta,aiS] = gsvd( A, L );
%%[Vpp,vpp] = eig(A,L); %Vpp = Vpp*(Vpp'*L*Vpp)^(-1/2);
%thetas = diag(Theta) ./ diag(aiS);
%% figure; plot(log10(diag(Theta)));
%[thetas,Ind] = sort( thetas,'descend' );
%Q = Q(:,Ind); P = P(:,Ind);
%Thetas = diag(thetas); 
%Theta = Theta(Ind,Ind); aiS = aiS(Ind,Ind); V = V(:,Ind);

% [theta,Ind] = sort( diag(Theta),'descend' );
% Q = Q(:,Ind); P = P(:,Ind);
% Theta = diag(theta); aiS = aiS(Ind,Ind); V = V(:,Ind);
disp([ 'Rectangular system pinversed ', num2str(toc) ]);
% Plot the Picard stuff
imax = min( find(thetas/thetas(1)<1e-16) );
if size(imax,1) == 0
    imax = size(thetas,1);
end

% tplo = theta(1:imax); bplo = Q'*b; bplo = bplo(1:imax);
% rplo = (Q'*b)./theta; rplo = rplo(1:imax);
% splo = diag(aiS); splo = splo(1:imax);
tplo = thetas(1:imax); bplo = Q'*b; bplo = bplo(1:imax);
rplo = (Q'*b)./thetas; rplo = rplo(1:imax);
figure
hold on;
plot(log10(abs(tplo)),'Color','blue');
plot(log10(abs(bplo)),'Color','red');
plot(log10(abs( rplo )),'Color','black');

% Filter eigenvalues
% ThetaT = Theta( 1:jmax , 1:jmax );
% bT     = Q'*b; bT = bT(1:jmax);
ThetaT = Thetas( 1:jmax , 1:jmax );
bT     = Q'*b; bT = bT(1:jmax);

% Solu = Q*(Theta\(Q'*b));
% Solu = (V*V')\V(:,1:jmax) * (ThetaT\bT); % (P*P') = eye( whatsoever )
Solu = P(:,1:jmax) * (ThetaT\bT); % (P*P') = eye( whatsoever )   (V*V')^(-1/2)*
% Solu = Red'*Solu;
% % Redistribute the identified loading on the nodes (for visu on GMSH)
% Sno = zeros(2*nnodesu,1);
% for i=1:size(elementsu,1)
%    Xloc = nodes(elementsu(i,:),:);    % Extract coords
%    nno = size(Xloc,1);
% 
%    ne = size(elementsu,2);
%    Se = []; maps = [];
%    for j=1:ne
%       Se = [ Se, Solu(2*i-1)/ntoelemu(elementsu(i,j),1), ...
%                  Solu(2*i)/ntoelemu(elementsu(i,j),1) ]; % Pass to nodes
%       maps = [maps,2*elementsu(i,j)-1,2*elementsu(i,j)];
%    end
%  
%    Sno(maps,1) = Sno(maps,1) + Se';
% end
Sno = Red*Solu;
% Output
Si = reshape(Sno,2,[])';  Sx = Si(:,1);  Sy = Si(:,2);
plotGMSH({Sx,'F_x';Sy,'F_y';Sno,'F_vect'}, elementsu, nodesu, 'output/identification');

% local visu
Seli = reshape(Solu,2,[])';  Selx = Seli(:,1);  Sely = Seli(:,2);
figure;
hold on
patch('Faces',elementsu(:,1:3),'Vertices',nodesu,'FaceVertexCData',Selx,'FaceColor','flat');
colorbar; axis equal;
if nnind > 0 && mod(nnind,2) == 1
    x = nodes(ceil(nnind/2),1); y = nodes(ceil(nnind/2),2); plot(x,y,'+');
end
figure;
hold on
patch('Faces',elementsu(:,1:3),'Vertices',nodesu,'FaceVertexCData',Sely,'FaceColor','flat');
colorbar; axis equal;
if nnind > 0 && mod(nnind,2) == 0
    x = nodes(ceil(nnind/2),1); y = nodes(ceil(nnind/2),2); plot(x,y,'+');
end

% Per node visu
% figure;
% patch('Faces',elementsu(:,1:3),'Vertices',nodesu,'FaceVertexCData',Sy,'FaceColor','interp');
% colorbar;

% % Reference visu (force on the displacement basis, it's ugly)
% fxy = reshape(f,2,[])';  fx = fxy(:,1);  fy = fxy(:,2);
% figure;
% patch('Faces',elements(:,1:3),'Vertices',nodes,'FaceVertexCData',fx,'FaceColor','interp');
% colorbar;
% figure;
% patch('Faces',elements(:,1:3),'Vertices',nodes,'FaceVertexCData',fy,'FaceColor','interp');
% colorbar;

% On line visu
Y = Ymin:Ly/100:Ymax; X = Xmoy*ones(1,size(Y,2)); XY = [X;Y];
Sline = passMesh2D(nodesu, elementsu, XY', [], Sno);
Slinexy = reshape(Sline,2,[])';  Slinex = Slinexy(:,1);  Sliney = Slinexy(:,2);

ref = loadV(1,1) + loadV(2,1)*Y + loadV(3,1)*Y.^2 + loadV(4,1)*Y.^3;
figure;
hold on;
set(gca, 'fontsize', 20);
plot(Y,Sliney,'LineWidth',3);
plot(Y,ref,'Color','red','LineWidth',3);
legend('identification','reference');