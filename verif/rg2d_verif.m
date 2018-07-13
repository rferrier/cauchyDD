% 25/06/2018
% Estimateur d'erreur basé sur RG

tic
close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E           = 210000; % MPa : Young modulus
nu          = 0.3;    % Poisson ratio
fscalar     = 250;    % N.mm-1 : Loading on the plate
mat         = [0, E, nu];
Npg         = 5;      % Nb Gauss points
ordertest   = 10;     % Order of test fonctions
dp          = 0;       % Plane strain

% Boundary conditions
dirichlet  = [0,1,0 ; 0,2,0 ; 0,3,0];
neumann1   = [3,2,fscalar ; 1,2,-fscalar];
neumann1p  = [3,1,-5*fscalar,10*fscalar,0];
neumann2   = [2,1,fscalar ; 4,1,-fscalar];
neumann3   = [3,1,fscalar ; 3,2,fscalar ; 2,1,fscalar ; 2,2,fscalar ; ...
              1,1,-fscalar ; 1,2,-fscalar ; 4,1,-fscalar ; 4,2,-fscalar];
neumann4   = [3,1,-fscalar ; 3,2,fscalar ; 2,1,fscalar ; 2,2,-fscalar ; ...
              1,1,fscalar ; 1,2,-fscalar ; 4,1,-fscalar ; 4,2,fscalar];
%neumann4   = [3,1,-fscalar ; 2,2,-fscalar ; ...
%              1,1,fscalar ; 4,2,fscalar];

dirichlet0 = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ; 3,1,0 ; 3,2,0 ; ... 
               4,1,0 ; 4,2,0 ; 5,1,0 ; 5,2,0 ; 6,1,0 ; 6,2,0 ]; % TODO : remove this useless stuff
neumann0   = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ; 3,1,0 ; 3,2,0 ; ... 
               4,1,0 ; 4,2,0 ; 5,1,0 ; 5,2,0 ; 6,1,0 ; 6,2,0 ];

% First, import the mesh
[ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/verif/plate_c_squared.msh' );
nnodes = size(nodes,1);

% mapBounds
[ node2b1, b2node1 ]   = mapBound( 1, boundary, nnodes );
[ node2b2, b2node2 ]   = mapBound( 2, boundary, nnodes );
[ node2b3, b2node3 ]   = mapBound( 3, boundary, nnodes );
[ node2b4, b2node4 ]   = mapBound( 4, boundary, nnodes );
[ node2b5, b2node5 ]   = mapBound( 5, boundary, nnodes );
[ node2b6, b2node6 ]   = mapBound( 6, boundary, nnodes );
               
% Then, build the stiffness matrix :
[K,C,nbloq,node2c,c2node] = Krig2 (nodes,elements,mat,order,boundary,dirichlet,dp);
Kinter = K( 1:2*nnodes, 1:2*nnodes );

Xmax = max(nodes(:,1)); Xmin = min(nodes(:,1)); Xmoy = (Xmax+Xmin)/2;
Ymax = max(nodes(:,2)); Ymin = min(nodes(:,2)); Ymoy = (Ymax+Ymin)/2;
Lx = Xmax-Xmin; Ly = Ymax-Ymin;

f1  = loading(nbloq,nodes,boundary,neumann1);
f1p = loading(nbloq,nodes,boundary,neumann1p);
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
sigma = stress(u4,E,nu,nodes,elements,order,1,ntoelem,dp);

% Output :
plotGMSH({u1,'U1';u2,'U2';u3,'U3';u4,'U4'}, elements, nodes, 'output/reference');

%% Overkill resolution
[ nodeso,elementso,ntoelemo,boundaryo,ordero] = readmesh( 'meshes/verif/plate_c_squared3.msh' );%/plate_overkill.msh' );
nnodeso = size(nodeso,1);
[Ko,Co,nbloqo,node2co,c2nodeo] = Krig2 (nodeso,elementso,mat,ordero,boundaryo,dirichlet,dp);
f1o  = loading(nbloqo,nodeso,boundaryo,neumann1);
f2o  = loading(nbloqo,nodeso,boundaryo,neumann2);
f3o  = loading(nbloqo,nodeso,boundaryo,neumann3);
f4o  = loading(nbloqo,nodeso,boundaryo,neumann4);
uin = Ko\[f1o,f2o,f3o,f4o];
u1o = uin(1:2*nnodeso,1); u2o = uin(1:2*nnodeso,2);
u3o = uin(1:2*nnodeso,3); u4o = uin(1:2*nnodeso,4);

% Compute the "real" error by the energy
e1 = energy( u1,nodes,elements,mat,order );
e2 = energy( u2,nodes,elements,mat,order );
e3 = energy( u3,nodes,elements,mat,order );
e4 = energy( u4,nodes,elements,mat,order );

e1o = energy( u1o,nodeso,elementso,mat,ordero );
e2o = energy( u2o,nodeso,elementso,mat,ordero );
e3o = energy( u3o,nodeso,elementso,mat,ordero );
e4o = energy( u4o,nodeso,elementso,mat,ordero );

error1 = sqrt(sum(e1o) - sum(e1));
error2 = sqrt(sum(e2o) - sum(e2));
error3 = sqrt(sum(e3o) - sum(e3));
error4 = sqrt(sum(e4o) - sum(e4));

%%% Local error with a passmesh (doesn't work for the order 2 mesh)
%up = passMesh2D(nodes, elements, nodeso, elementso, [u1,u2,u3,u4], 0);
%up1 = up(:,1); up2 = up(:,2); up3 = up(:,3); up4 = up(:,4);

%e1d = energy( u1o-up1,nodeso,elementso,mat,ordero );
%e2d = energy( u2o-up2,nodeso,elementso,mat,ordero );
%e3d = energy( u3o-up3,nodeso,elementso,mat,ordero );
%e4d = energy( u4o-up4,nodeso,elementso,mat,ordero );

%% And display the error
%try
%figure;
%patch('Faces',elementso,'Vertices',nodeso,'FaceVertexCData',e1d,'FaceColor','flat');
%colorbar(); set(colorbar, 'fontsize', 20); axis([0,1,0,1,0,1],'square');
%legend('Overkill error');
%end

%error1d = sqrt(sum(e1d));
%error2d = sqrt(sum(e2d));
%error3d = sqrt(sum(e3d));
%error4d = sqrt(sum(e4d));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ nodes2,elements2,ntoelem2,boundary2,order] = readmesh( 'meshes/verif/plate_c_squared.msh' );
nnodes2 = size(nodes2,1);
%[K2,C2,nbloq2,node2c2,c2node2] = Krig2 (nodes2,elements2,mat,order,boundary2,dirichlet,dp);
%Kinter2 = K2( 1:2*nnodes2, 1:2*nnodes2 );

%% Pass meshes
UF = [u1,u2,u3,u4,f1,f2,f3,f4];
UR = passMesh2D(nodes, elements, nodes2, elements2, UF, 0);

ur1 = UR(:,1); fr1 = UR(:,5);
ur2 = UR(:,2); fr2 = UR(:,6);
ur3 = UR(:,3); fr3 = UR(:,7);
ur4 = UR(:,4); fr4 = UR(:,8);

plotGMSH({ur1,'U1';ur2,'U2';ur3,'U3';ur4,'U4'},elements2, nodes2, 'output/bound');

nmin = 1:nnodes2;

%% Build the maps of the Dirichlet and Neumann bounds
%b2nodesD = cell(4,2); % nb of boundaries * nb of dimensions
%b2nodesN = cell(4,2);
b2nodesD = []; b2nodesN = []; b2nodesTOT = [];
for i=1:6
   [~, b2node] = mapBound(i, boundary2, nnodes2);
   dofD = dirichlet0(find(dirichlet0(:,1)==i),2); % dof of the given dirichlet
   dofN = neumann0(find(neumann0(:,1)==i),2); % dof of the given neumann
   b2nodesDo = []; % Stores the added nodes
   for j=1:size(dofD,1)
      b2nodesD = [b2nodesD ; 2*b2node-2+dofD(j)];
   end
   b2nodesTOT = [b2nodesTOT ; 2*b2node-1 ; 2*b2node];
%   if max(size(dofN)) > 0
   for j=1:size(dofN,1)
      b2nodesN = [b2nodesN ; 2*b2node-2+dofN(j)];
   end
end

% Discrete boundary stuff
b2nodesD = intersect(b2nodesD,[2*nmin-1;2*nmin]);

b2nodesnoD = setdiff(b2nodesTOT,b2nodesD);  % Missing Dirichlet conditions
b2nodesnoN = setdiff(b2nodesTOT,b2nodesN);  % Missing Neumann conditions (unused except for visualization)
% Put u to 0 at those missing BCs
% u1(b2nodesnoD) = 0; u2(b2nodesnoD) = 0; u3(b2nodesnoD) = 0; u4(b2nodesnoD) = 0;

%% Preliminary stuff : find the volumic elements corresponding to the boundaries
nboun2 = size(boundary2,1); nelem2 = size(elements2,1);
boun2vol2 = zeros( nboun2, 1 ); extnorm2 = zeros( nboun2, 2 );
for i=1:nboun2
   % Volumic element
   no1 = boundary2(i,2); no2 = boundary2(i,3); % only with 2 nodes even if order > 1
   cand1 = rem( find(elements2==no1)-1,nelem2 )+1; % find gives line + column*size
   cand2 = rem( find(elements2==no2)-1,nelem2 )+1;
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

nboun1 = size(boundary2,1); nelem1 = size(elements2,1);
boun2vol1 = zeros( nboun1, 1 ); extnorm1 = zeros( nboun1, 2 );
sigr1 = zeros( nboun1, 3*order ); sigr2 = zeros( nboun1, 3*order );
urr1  = zeros( nboun1, 2+2*order ); urr2 = zeros( nboun1, 2+2*order );
urr3  = zeros( nboun1, 2+2*order ); urr4 = zeros( nboun1, 2+2*order );
for i=1:nboun1 % TODO : rationalize the stuff (don't do 2 times the same work)
   % Volumic element
   no1 = boundary2(i,2); no2 = boundary2(i,3); % only with 2 nodes even if order > 1
   cand1 = rem( find(elements2==no1)-1,nelem1 )+1; % find gives line + column*size
   cand2 = rem( find(elements2==no2)-1,nelem1 )+1;
   boun2vol1(i) = intersect(cand1, cand2); % If everything went well, there is only one
   
   % Exterior normal
   elt = boun2vol1(i); no3 = setdiff( elements2( elt, 1:3 ), [no1,no2]);
   x1 = nodes2(no1,1); y1 = nodes2(no1,2);
   x2 = nodes2(no2,1); y2 = nodes2(no2,2);
   x3 = nodes2(no3,1); y3 = nodes2(no3,2);
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
nodes2dof = zeros(2*nnodes2,1); dof2nodes = []; szD = 0;
for i=1:2*nnodes2
   if min(size(find(b2nodesD==i)))==0
      nodes2dof(i) = szD+1;
      dof2nodes(szD+1) = i;
      szD = szD+1;
   end
end

% ...and then give a dof to each element
boun2dof = zeros(nboun2,1);
for i=1:nboun2
   boun2dof(i,1) = nodes2dof(2*boundary2(i,2)-1);
   boun2dof(i,2) = nodes2dof(2*boundary2(i,2));
   boun2dof(i,3) = nodes2dof(2*boundary2(i,3)-1);
   boun2dof(i,4) = nodes2dof(2*boundary2(i,3));
end

%% Order the dofs : force and displacement on the boundary
nodesbound2 = unique(boundary2(:)); % Nodes on the boundary
nnbound2    = size(nodesbound2,1);
boun2doftot = zeros(nboun2,2);

for i=1:nboun2 % Stores the nodes associated to each element of bound2 (in boundary numbering)
   boun2doftot(i,1) = find(nodesbound2==boundary2(i,2)); % Rem : boundary2(i,1) is the no of boundary
   boun2doftot(i,2) = find(nodesbound2==boundary2(i,3));
end

knownD = [];
knownN = [];
for i=1:nboun2 % Chooses the known displacement dofs
   index = boundary2(i,1);
   isimposed = find( dirichlet0(:,1) == index );
   dofs = dirichlet0(isimposed,2);
   knownD = [ knownD ; 2*boun2doftot(i,1)-2+dofs ];
   knownD = [ knownD ; 2*boun2doftot(i,2)-2+dofs ];
   
   isimposed = find( neumann0(:,1) == index );
   dofs = neumann0(isimposed,2);
   knownN = [ knownN ; 2*i-2+dofs ];
end

knownD = unique(knownD); knownN = unique(knownN); % Remove redondnacy
knownD = intersect( knownD, [2*nmin-1;2*nmin] ); % Discrete measurement    /!\ There is a global/local confusion /!\

tofindD = setdiff( 1:2*nnbound2 , knownD );
tofindN = setdiff( 1:2*nnbound2 , knownN );
%%

[ nodesu,elementsu,ntoelemu,boundaryu,orderu] = readmesh( 'meshes/verif/plate_c_squared.msh' );
nnodesu = size(nodesu,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test-functions shit
% Build the polynomial test functions.
if ordertest == 10
   load('conditions10_2d.mat','-ascii');
   M = spconvert(conditions10_2d); clear('conditions10_2d');
   nmax = 10;
elseif ordertest == 20
   if dp == 1 % Plane strain
      load('conditions20_2dc.mat','-ascii');
      M = spconvert(conditions20_2dc); clear('conditions20_2dc');
   else
      load('conditions20_2d.mat','-ascii');
      M = spconvert(conditions20_2d); clear('conditions20_2d');
   end

   nmax = 20;
end
ncoef =2*(nmax+1)^2; neq = ncoef;

% Suppress the superior terms
for i=0:nmax
   for j=0:nmax
      if i+j>nmax
         M((nmax+1)*i+j+1,:) = 0;
         M((nmax+1)*i+j+1,(nmax+1)*i+j+1) = 1;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build the matrix of the energetic inner product

P0 = zeros(2*(nmax+1)^2);
for ii=0:nmax
   nd1 = nmax-ii;
   for jj=0:nd1
      for kk=0:nmax
         nd2 = nmax-kk;
         for ll=0:nd2
            indexij = (nmax+1)*ii + jj + 1;
            indexkl = (nmax+1)*kk + ll + 1;

            e1111 = 0; e1112 = 0; e1211 = 0; e1212 = 0;
            e2212 = 0; e2222 = 0; e2111 = 0; e2112 = 0;

            if ii+kk>1
               e1111 = E/(1-nu^2)*(ii*kk)/((ii+kk-1)*(jj+ll+1));
               e2212 = E/(4*(1+nu))*(ii*kk)/((ii+kk-1)*(jj+ll+1));
            end
            if jj+ll>1
               e1112 = E/(4*(1+nu))*(jj*ll)/((ii+kk+1)*(jj+ll-1));
               e2222 = E/(1-nu^2)*(jj*ll)/((ii+kk+1)*(jj+ll-1));
            end
            if ii+kk>0 && jj+ll>0
               e1211 = nu*E/(1-nu^2)*(ii*ll)/((ii+kk)*(jj+ll));
               e1212 = E/(4*(1+nu))*(jj*kk)/((ii+kk)*(jj+ll));
               e2111 = nu*E/(1-nu^2)*(jj*kk)/((ii+kk)*(jj+ll));
               e2112 = E/(4*(1+nu))*(ii*ll)/((ii+kk)*(jj+ll));
            end

            P0( indexij, indexkl ) = e1111 + 2*e1112;
            P0( (nmax+1)^2 + indexij, indexkl ) = e2111 + 2*e2112;
            P0( indexij, (nmax+1)^2 + indexkl ) = e1211 + 2*e1212;
            P0( (nmax+1)^2 + indexij, (nmax+1)^2 + indexkl ) = e2222 + 2*e2212;
         end
      end
   end
end

P = coef'*P0*coef; % Pass in the test-functions basis
P = .5*(P+P');

[Q,Theta] = eig(P); %Q = real(Q); Theta = real(Theta);% Diagonalize P
[theta,ind] = sort(diag(Theta),'descend'); Q = Q(:,ind);
inc = min(find(theta<1e-12*theta(1))); % Truncate the inner product
sT = diag(sqrt(theta)); Q1 = Q*sT; Q2 = Q*inv(sT); % Normalize
Q1 = Q1(:,1:inc-1); Q2 = Q2(:,1:inc-1); % And truncate
% Q1 passes from the old basis to the new one that is orthogonalized (but smaller because of condition number)

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

%%
% Transform the given data in the proper format
f_known1 = zeros(2*nboun2,1); u_known1 = zeros(2*nnbound2,1);
f_known2 = zeros(2*nboun2,1); u_known2 = zeros(2*nnbound2,1);
f_known3 = zeros(2*nboun2,1); u_known3 = zeros(2*nnbound2,1);
f_known4 = zeros(2*nboun2,1); u_known4 = zeros(2*nnbound2,1);
   
for i=1:nboun2
   bonod = boundary2(i,:); exno = extnorm2(i,:)';
   % Reference force from the BC's
   if exno(1) == 1 % Bound 2
      fer1 = [0;0]; fer2 = [fscalar;0];
      fer3 = [fscalar;fscalar]; fer4 = [fscalar;-fscalar];%fer4 = [0;-fscalar];
   elseif exno(1) == -1 % Bound 4
      fer1 = [0;0]; fer2 = -[fscalar;0];
      fer3 = [-fscalar;-fscalar]; fer4 = [-fscalar;fscalar];%fer4 = [0;fscalar];%
   elseif exno(2) == 1 % Bound 3
      fer1 = [0;fscalar]; fer2 = [0;0];
      fer3 = [fscalar;fscalar]; fer4 = [-fscalar;fscalar];%fer4 = [-fscalar;0];%
   elseif exno(2) == -1 % Bound 1
      fer1 = -[0;fscalar]; fer2 = [0;0];
      fer3 = [-fscalar;-fscalar]; fer4 = [fscalar;-fscalar];%fer4 = [fscalar;0];%
   else % The cracks (no loading)
      fer1 = [0;0]; fer2 = [0;0]; fer3 = [0;0]; fer4 = [0;0];
   end
     
   indicesLoc = [2*boun2doftot(i,:)-1,2*boun2doftot(i,:)]; % Local Displacement dofs
   indicesGlo = [2*boundary2(i,[2,3])-1,2*boundary2(i,[2,3])]; % Global Displacement dofs
      
   u_known1(indicesLoc) = ur1(indicesGlo);
   u_known2(indicesLoc) = ur2(indicesGlo);
   u_known3(indicesLoc) = ur3(indicesGlo);
   u_known4(indicesLoc) = ur4(indicesGlo);
      
   f_known1([2*i-1,2*i]) = fer1;
   f_known2([2*i-1,2*i]) = fer2;
   f_known3([2*i-1,2*i]) = fer3;
   f_known4([2*i-1,2*i]) = fer4;
end

disp([ 'Direct problem solved and data management ', num2str(toc) ]);

%% The operator that passes from the polynomial stuff to the spacial one
Xx  = nodes2(:,1); Yy = nodes2(:,2);
XYp = zeros(2*nnodes2,2*(nmax+1)^2);
for ii=0:nmax
   nd1 = nmax-ii;
   for jj=0:nd1
      indexij = (nmax+1)*ii + jj + 1;
      XYp(1:2:2*nnodes2-1,indexij)              = Xx.^ii.*Yy.^jj;
      XYp(2:2:2*nnodes2  ,(nmax+1)^2 + indexij) = Xx.^ii.*Yy.^jj;
   end
end
Apas = XYp*coef;

%% Debug : get u in the polynomial basis
AAA = Apas'*Apas+1e-10*eye(nftest);
upol1 = AAA\(Apas'*u1);
upol2 = AAA\(Apas'*u2);
upol3 = AAA\(Apas'*u3);
upol4 = AAA\(Apas'*u4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the operators
tic
Ng = Npg;
[ Xg, Wg ] = gaussPt1d( Ng );
Ndots = size(Wg,1); % (= Ng by the way)

%% Build the list of Gauss points, and of construction-functions
Xxg = zeros( nboun2*Ndots,1 ); Yyg = zeros( nboun2*Ndots,1 );
Wwg = zeros( nboun2*Ndots,1 );
Phi = sparse( 2*nboun2*Ndots, 2*nnbound2 );
Psi = sparse( 2*nboun2*Ndots, 2*nboun2 );
exnor = zeros( nboun2*Ndots,2); % Normal

for i=1:nboun2
   bonod = boundary2(i,:); exno = extnorm2(i,:)';
   
   no1 = bonod(2); no2 = bonod(3); boname = bonod(1);
   x1 = nodes2(no1,1); y1 = nodes2(no1,2);
   x2 = nodes2(no2,1); y2 = nodes2(no2,2);
      
   len = sqrt( (x1-x2)^2 + (y1-y2)^2 );

   indDtot = [2*boun2doftot(i,1)-1,2*boun2doftot(i,1),...
              2*boun2doftot(i,2)-1,2*boun2doftot(i,2)]; % U dofs the element is associated to
   indNtot = [2*i-1,2*i]; % F dofs the element is associated to
                
   for j=1:Ndots
      xg = Xg(j,:); wg = Wg(j);
      xgr  = ( (1-xg)*[x1;y1] + xg*[x2;y2] ); % abscissae

      X = xgr(1)/Lx; Y = xgr(2)/Lx; % Normalize
      Swg   = len * wg;

      indexg = Ndots*(i-1) + j;
      Xxg( indexg ) = X; Yyg( indexg ) = Y; Wwg( indexg ) = Swg;
      exnor( indexg, :) = exno';

      Phi( 2*indexg-1, indDtot)  = [ 1-xg, 0, xg, 0 ];
      Phi( 2*indexg  , indDtot)  = [ 0, 1-xg, 0, xg ];
      Psi( 2*indexg-1, indNtot)  = [1,0];
      Psi( 2*indexg  , indNtot)  = [0,1];
   end
end

%% Build the matrix of test-functions
Vv = zeros( size(coef,1), 2*nboun2*Ndots );
Sv = zeros( size(coef,1), 2*nboun2*Ndots );
exnoX = exnor(:,1); exnoY = exnor(:,2);

for ii=0:nmax
   nd1 = nmax-ii;
   for jj=0:nd1
                  
         sloc11a = zeros(nboun2*Ndots,1); sloc12a = zeros(nboun2*Ndots,1);
         sloc22a = zeros(nboun2*Ndots,1); sloc11b = zeros(nboun2*Ndots,1);
         sloc12b = zeros(nboun2*Ndots,1); sloc22b = zeros(nboun2*Ndots,1);
         if ii>0
            XY = ii.*Xxg.^(ii-1).*Yyg.^jj;
            sloc11a = E/(1-nu^2)*XY;
            sloc22a = nu*E/(1-nu^2)*XY;
            sloc12b = E/(2*(1+nu))*XY;
         end
         if jj>0
            XY = jj.*Xxg.^ii.*Yyg.^(jj-1);
            sloc11b = nu*E/(1-nu^2)*XY;
            sloc22b = E/(1-nu^2)*XY;
            sloc12a = E/(2*(1+nu))*XY;
         end

         fpaax = sloc11a.*exnoX + sloc12a.*exnoY;
         fpaay = sloc12a.*exnoX + sloc22a.*exnoY;

         fpabx = sloc11b.*exnoX + sloc12b.*exnoY;
         fpaby = sloc12b.*exnoX + sloc22b.*exnoY;

         XY = Xxg.^ii.*Yyg.^jj;
         vpaax = XY; vpaby = XY;

         index = (nmax+1)*ii + jj + 1;

         Sv( index, 1:2:2*nboun2*Ndots-1 )              = Wwg' .* fpaax';
         Sv( (nmax+1)^2 + index, 1:2:2*nboun2*Ndots-1 ) = Wwg' .* fpabx';

         Sv( index, 2:2:2*nboun2*Ndots )                = Wwg' .* fpaay';
         Sv( (nmax+1)^2 + index, 2:2:2*nboun2*Ndots )   = Wwg' .* fpaby';

         Vv( index, 1:2:2*nboun2*Ndots-1 )              = Wwg' .* vpaax';
         Vv( (nmax+1)^2 + index, 2:2:2*nboun2*Ndots )   = Wwg' .* vpaby';
   end
end

Rupij = Sv*Phi; clear Sv; clear Phi;
Rfpij = Vv*Psi; clear Vv; clear Psi;
Ru = coef'*Rupij; Rf = coef'*Rfpij;

%% Restrict the data on the given part (not necessary, but cleaner)
f_known1(tofindN) = 0; f_known2(tofindN) = 0;
f_known3(tofindN) = 0; f_known4(tofindN) = 0;
u_known1(tofindD) = 0; u_known2(tofindD) = 0;
u_known3(tofindD) = 0; u_known4(tofindD) = 0;

%% Restrict the operators
Rfr = Rf(:,knownN); Rur = Ru(:,knownD);

Rhs1 = Rur*u_known1(knownD) - Rfr*f_known1(knownN);
Rhs2 = Rur*u_known2(knownD) - Rfr*f_known2(knownN);
Rhs3 = Rur*u_known3(knownD) - Rfr*f_known3(knownN);
Rhs4 = Rur*u_known4(knownD) - Rfr*f_known4(knownN);

%% Pass into the orthogonalized basis
%R1 = (Q1'*Q1)\(Q1'*Rhs1); R2 = (Q1'*Q1)\(Q1'*Rhs2);
%R3 = (Q1'*Q1)\(Q1'*Rhs3); R4 = (Q1'*Q1)\(Q1'*Rhs4);
R1 = Q2'*Rhs1; R2 = Q2'*Rhs2;
R3 = Q2'*Rhs3; R4 = Q2'*Rhs4;

try
figure;
hold on;
plot(abs(R1),'Color','red');
plot(abs(R2),'Color','blue');
plot(abs(R3),'Color','black');
plot(abs(R4),'Color','green');
legend('R1','R2','R3','R4');
end

eest1 = norm(R1); % Error estimators
eest2 = norm(R2);
eest3 = norm(R3);
eest4 = norm(R4);

%% Compute the local energetic error
% First, the list of Gauss points and surfaces
Xxg = zeros( nelemu,1 ); Yyg = zeros( nelemu,1 );
Ssg = zeros( nelemu,1); % Surface

% The test-field the energy of which is being computed
v01 = coef*Q2*R1; v02 = coef*Q2*R2;
v03 = coef*Q2*R3; v04 = coef*Q2*R4;
%v01 = coef*Rhs1; v02 = coef*Rhs2;
%v03 = coef*Rhs3; v04 = coef*Rhs4;

Evloc1 = zeros(nelemu,1);
Evloc2 = zeros(nelemu,1);
Evloc3 = zeros(nelemu,1);
Evloc4 = zeros(nelemu,1);

%for i=1:nelemu
%   no1 = elementsu(i,1); no2 = elementsu(i,2); no3 = elementsu(i,3);
%   x1 = nodesu(no1,1); x2 = nodesu(no2,1); x3 = nodesu(no3,1);
%   y1 = nodesu(no1,2); y2 = nodesu(no2,2); y3 = nodesu(no3,2);

%   Sg = .5*abs( (y2-y1)*(x3-x1) - (y3-y1)*(x2-x1) );
%   Xg = (x1+x2+x3)/3;
%   Yg = (y1+y2+y3)/3;

%   %Ssg(i,1) = Sg; Xxg(i,1) = Xg; Yyg(i,1) = Yg;

%   Eloc = zeros(2*(nmax+1)^2);
%   for ii=0:nmax
%      nd1 = nmax-ii;
%      for jj=0:nd1
%         for kk=0:nmax
%            nd2 = nmax-kk;
%            for ll=0:nd2
%               indexij = (nmax+1)*ii + jj + 1;
%               indexkl = (nmax+1)*kk + ll + 1;

%               e1111 = 0; e1112 = 0; e1211 = 0; e1212 = 0;
%               e2212 = 0; e2222 = 0; e2111 = 0; e2112 = 0;

%               if ii+kk>1
%                  e1111 = Sg*E/(1-nu^2)*(ii*kk)*Xg^(ii+kk-2)*Yg^(jj+ll);
%                  e2212 = E/(4*(1+nu))*(ii*kk)*Xg^(ii+kk-2)*Yg^(jj+ll);
%               end
%               if jj+ll>1
%                  e1112 = E/(4*(1+nu))*(jj*ll)*Xg^(ii+kk)*Yg^(jj+ll-2);
%                  e2222 = E/(1-nu^2)*(jj*ll)*Xg^(ii+kk)*Yg^(jj+ll-2);
%               end
%               if ii+kk>0 && jj+ll>0
%                  e1211 = nu*E/(1-nu^2)*(ii*ll)*Xg^(ii+kk-1)*Yg^(jj+ll-1);
%                  e1212 = E/(4*(1+nu))*(jj*kk)*Xg^(ii+kk-1)*Yg^(jj+ll-1);
%                  e2111 = nu*E/(1-nu^2)*(jj*kk)*Xg^(ii+kk-1)*Yg^(jj+ll-1);
%                  e2112 = E/(4*(1+nu))*(ii*ll)*Xg^(ii+kk-1)*Yg^(jj+ll-1);
%               end

%               Eloc( indexij, indexkl ) = e1111 + 2*e1112;
%               Eloc( (nmax+1)^2 + indexij, indexkl ) = e2111 + 2*e2112;
%               Eloc( indexij, (nmax+1)^2 + indexkl ) = e1211 + 2*e1212;
%               Eloc( (nmax+1)^2 + indexij, (nmax+1)^2 + indexkl ) = e2222 + 2*e2212;
%            end
%         end
%      end
%   end

%   Evloc1(i) = v01'*Eloc*v01;
%   Evloc2(i) = v02'*Eloc*v02;
%   Evloc3(i) = v03'*Eloc*v03;
%   Evloc4(i) = v04'*Eloc*v04;
%end

%% Display locally the error
%try
%figure;
%patch('Faces',elementsu,'Vertices',nodesu,'FaceVertexCData',Evloc1 ,'FaceColor','flat');
%colorbar(); set(colorbar, 'fontsize', 20); axis([0,1,0,1,0,1],'square');
%end

%% Display the argmax of the residual
varg01 = XYp*v01; varg02 = XYp*v02; varg03 = XYp*v03; varg04 = XYp*v04;

e01 = energy( varg01,nodes2,elements2,mat,order );
e02 = energy( varg02,nodes2,elements2,mat,order );
e03 = energy( varg03,nodes2,elements2,mat,order );
e04 = energy( varg04,nodes2,elements2,mat,order );

try
figure;
patch('Faces',elementsu,'Vertices',nodesu,'FaceVertexCData',varg01(2:2:2*nnodes2) ,'FaceColor','interp');
colorbar(); set(colorbar, 'fontsize', 20); axis([0,1,0,1,0,1],'square');
end

try
figure;
patch('Faces',elementsu,'Vertices',nodesu,'FaceVertexCData',e01 ,'FaceColor','flat');
colorbar(); set(colorbar, 'fontsize', 20); axis([0,1,0,1,0,1],'square');
end

disp([ 'Error estimation ', num2str(toc) ]);
