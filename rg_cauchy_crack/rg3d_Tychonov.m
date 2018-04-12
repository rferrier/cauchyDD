% 05/03/2018
% Problèmes directs et de Cauchy par écart à la réciprocité, problème 3D, régularisation de Tychonov

tic
close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E          = 210000; % MPa : Young modulus
nu         = 0.3;    % Poisson ratio
fscalar    = 1;      % N.mm-2 : Loading on the plate
mat        = [0, E, nu];
br         = .0;      % Noise level
mur        = 5e4;      % Regularization parameter
%mfr        = 1e4;      % Other regularization parameter
%mur        = 1e9;      % Regularization parameter
regular    = 0;      % Use the derivative regularization matrix (0 : Id, 1 : derivative)
upper_term = 0;      % 1 : use i=0:10, j=0:10, 0 : use i>=0,j>=0,i+j<=10
froreg     = 1;      % frobenius preconditioner
Npg        = 2;      % Nb Gauss points
ordertest  = 20;     % Order of test fonctions
testLcurve = 1;

% Boundary conditions
dirichlet = [3,1,0; 3,2,0 ; 3,3,0
             4,1,0; 4,2,0 ; 4,3,0
             5,1,0; 5,2,0 ; 5,3,0
             6,1,0; 6,2,0 ; 6,3,0 ];
neumann   = [7,3,-fscalar, ];

% BCs for the Cauchy stuff
dirichlet0 = [ 1,1,0 ; 1,2,0 ; 1,3,0 ; 
               3,1,0 ; 3,2,0 ; 3,3,0 ; 
               4,1,0 ; 4,2,0 ; 4,3,0 ; 
               5,1,0 ; 5,2,0 ; 5,3,0 ; 
               6,1,0 ; 6,2,0 ; 6,3,0 ];
neumann0   = [ 1,1,0 ; 1,2,0 ; 1,3,0 ];

%dirichlet0 = [ 2,1,0 ; 2,2,0 ; 2,3,0 ; 
%               3,1,0 ; 3,2,0 ; 3,3,0 ; 
%               4,1,0 ; 4,2,0 ; 4,3,0 ; 
%               5,1,0 ; 5,2,0 ; 5,3,0 ; 
%               6,1,0 ; 6,2,0 ; 6,3,0 ];
%neumann0   = [ 1,1,0 ; 1,2,0 ; 1,3,0 ];

% Import the mesh
%[ nodes,elements,ntoelem,boundary,order ] = readmesh3D( 'meshes/plate3d_charge2.msh' );
[ nodes,elements,ntoelem,boundary,order ] = readmesh3D( 'meshes/rg3d_qcq/plate3d_charge2nt.msh' );
nnodes = size(nodes,1);

[K,C,nbloq] = Krig3D (nodes,elements,mat,order,boundary,dirichlet);
Kinter = K(1:3*nnodes, 1:3*nnodes);
f = loading3D(nbloq,nodes,boundary,neumann);
uin = K\f; uref = uin(1:3*nnodes); fref = Kinter*uref;

sigma = stress3D(uref,mat,nodes,elements,order,1,ntoelem);
plotGMSH3D({uref,'Vect_U';sigma,'stress'}, elements, nodes, 'output/reference');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[ nodes2,elements2,ntoelem2,boundary2,order] = readmesh3D( 'meshes/plate3d_charge2.msh' );
[ nodes2,elements2,ntoelem2,boundary2,order] = readmesh3D( 'meshes/rg3d_qcq/plate3d_charge2_unt.msh' );
%[ nodes2,elements2,ntoelem2,boundary2,order] = readmesh3D( 'meshes/rg3d_qcq/plate3d_charge2_u.msh' );

boundary2( find(boundary2(:,1)==7) , : ) = []; % Hack : remove the double elements in boundary 2, (ie boundary 7)

nnodes2 = size(nodes2,1);
[K2,C2,nbloq2,node2c2,c2node2] = Krig3D(nodes2,elements2,mat,order,boundary2,dirichlet);
Kinter2 = K2( 1:2*nnodes2, 1:2*nnodes2 );

% Slight contraction of the second mesh for the passmesh
mean2x = mean(nodes2(:,1)); mean2y = mean(nodes2(:,2)); mean2z = mean(nodes2(:,3));
nodes2c = nodes2;
nodes2c(:,1) = (1-1e-6)*(nodes2c(:,1)-mean2x) + mean2x;
nodes2c(:,2) = (1-1e-6)*(nodes2c(:,2)-mean2y) + mean2y;
nodes2c(:,3) = (1-1e-6)*(nodes2c(:,3)-mean2z) + mean2z;

UR = passMesh3D( nodes, elements, nodes2c, elements2, [uref,fref] );
ur = UR(:,1); fr = UR(:,2);

plotGMSH3D({ur,'U_bound';fr,'F_bound'}, elements2, nodes2, 'output/bound');

% Add the noise
un = ur;
am = sqrt(mean(ur.^2));
br1 = randn(3*nnodes2,1);
%noise = load('noises/cauchyRG.mat');
%br1 = noise.br1; br2 = noise.br2; br3 = noise.br3; br4 = noise.br4;
%u1 = ( 1 + br*br1 ) .* u1; u2 = ( 1 + br*br2 ) .* u2;
%u3 = ( 1 + br*br3 ) .* u3; u4 = ( 1 + br*br4 ) .* u4;
ur = ur + am*br*br1;

%% Preliminary stuff : find the volumic elements corresponding to the boundaries
nboun2 = size(boundary2,1); nelem2 = size(elements2,1);
boun2vol = zeros( nboun2, 1 ); extnorm = zeros( nboun2, 3 );
frr = zeros( nboun2, 3 );
urr = zeros( nboun2, 9 );
for i=1:nboun2
   % Volumic element
   no1 = boundary2(i,2); no2 = boundary2(i,3); no3 = boundary2(i,4);
   cand1 = rem( find(elements2==no1)-1,nelem2 )+1; % find gives line + column*size
   cand2 = rem( find(elements2==no2)-1,nelem2 )+1;
   cand3 = rem( find(elements2==no3)-1,nelem2 )+1;
   boun2vol(i) = intersect( intersect(cand1, cand2), cand3 ); % If everything went well, there is only one
   
   % Exterior normal
   elt = boun2vol(i); no4 = setdiff( elements2( elt, 1:4 ), [no1,no2,no3]);
   x1 = nodes2(no1,1); y1 = nodes2(no1,2); z1 = nodes2(no1,3);
   x2 = nodes2(no2,1); y2 = nodes2(no2,2); z2 = nodes2(no2,3);
   x3 = nodes2(no3,1); y3 = nodes2(no3,2); z3 = nodes2(no3,3);
   x4 = nodes2(no4,1); y4 = nodes2(no4,2); z4 = nodes2(no4,3);
   
   nx = (y1-y2)*(z1-z3)-(y1-y3)*(z1-z2); %
   ny = (x1-x3)*(z1-z2)-(x1-x2)*(z1-z3); % Vectorial product
   nz = (x1-x2)*(y1-y3)-(x1-x3)*(y1-y2); %
   surfa = .5*sqrt(nx^2+ny^2+nz^2);
   extnorm(i,:) = [ nx ,ny ,nz ] / (2*surfa);

   if nx*(x4-x1)+ny*(y4-y1)+nz*(z4-z1) > 0 % Check that the normal is exterior
      extnorm(i,:) = -extnorm(i,:);
   end
   
   urr(i,:) = ur( [ 3*no1-2,3*no1-1,3*no1,...
                    3*no2-2,3*no2-1,3*no2,...
                    3*no2-2,3*no2-1,3*no2 ] );
   frr(i,:) = 1/(3*surfa)*( fr([3*no1-2,3*no1-1,3*no1]) +...
                           fr([3*no2-2,3*no2-1,3*no2]) +...
                           fr([3*no3-2,3*no3-1,3*no3]) );
   if order == 2
      warning("order 2 unimplemented")
   end
end

knownD = [];
knownN = [];
for i=1:nboun2 % Chooses the known displacement dofs
   index = boundary2(i,1);
   isimposed = find( dirichlet0(:,1) == index );
   dofs = dirichlet0(isimposed,2);
   knownD = [ knownD ; 3*i-3+dofs ];
   
   isimposed = find( neumann0(:,1) == index );
   dofs = neumann0(isimposed,2);
   knownN = [ knownN ; 3*i-3+dofs ];
end

knownD = unique(knownD); knownN = unique(knownN); % Remove redondnacy

% tofindD = setdiff( 1:3*nboun2 , knownD );
tofindN = setdiff( 1:3*nboun2 , knownN );

% Associate numerotation of the nodes on the boundary
buff = boundary2(:,2:4); buff = buff(:);
nodeboun2glob = unique(buff); nnodeboun = size(nodeboun2glob,1);
nodeglob2boun = zeros(nnodes2,1);
for i=1:nnodeboun
   nodeglob2boun(nodeboun2glob(i)) = i;
end
boun2loc = nodeglob2boun(boundary2(:,2:4));

nodekD = [];
for i=1:size(knownD,1)
   dof = mod(knownD(i)+2,3)+1; % No of the dof (1=X, 2=Y, 3=Z)
   B = ceil(knownD(i)/3); % No of the boundary element
   N = boun2loc(B,:);
   nodekD = [ nodekD, 3*N-3+dof ];
end
nodekD  = unique(nodekD); knownD = nodekD; % Erase knownD for name reasons
tofindD = setdiff( 1:3*nnodeboun, knownD );

ndofN = size(tofindN,2); ndofD = size(tofindD,2); % Nbs of dofs to find

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve via reciprocity gap
% Build the polynomial test functions.
if ordertest == 5
   load('conditions5_3D.mat','-ascii');
   M = spconvert(conditions5_3D); clear('conditions5_3D');
   nmax = 5;
   ncoef = 3*(nmax+1)^3; neq = ncoef;
elseif ordertest == 10
   load('conditions10_3D.mat','-ascii');
   M = spconvert(conditions10_3D); clear('conditions10_3D');
   nmax = 10;
   ncoef = 3*(nmax+1)^3; neq = ncoef;
elseif ordertest == 20
   load('conditions20_3D.mat','-ascii');
   M = spconvert(conditions20_3D); clear('conditions20_3D');
   nmax = 20;
   ncoef = 3*(nmax+1)^3; neq = ncoef;
end

% Suppress the superior terms
if upper_term == 0 % Rem : there will be lots of 0 in coef
   for i=0:nmax
      for j=0:nmax
         for k=0:nmax
            if i+j+k>nmax
               index = (nmax+1)^2*i+(nmax+1)*j+k+1;
               M(index,:) = 0;
               M(index,index) = 1;
            end
         end
      end
   end
end

%% Suppress zeros rows
%toremove = [];
%for k=1:size(M,1)
%  if norm(M(k,:)) == 0
%     toremove(end+1) = k;
%  end
%end
%M(toremove,:) = [];

coef  = fatNull(M);    % /!\ Both don't give exactly the same results /!\
%coef = null(full(M));
nftest = size(coef,2); % Nb of avaliable test functions

%%
% Transform the given data in the proper format
f_known = zeros(3*nboun2,1); u_known = zeros(3*nnodeboun,1);
for i=1:nboun2
   ind0       = [ 3*i-2, 3*i-1, 3*i ];
   f_known(ind0)       = frr(i,:);
end
for i=1:nnodeboun
   ind0 = [ 3*i-2, 3*i-1, 3*i ];
   ind1 = [ 3*nodeboun2glob(i)-2, 3*nodeboun2glob(i)-1, 3*nodeboun2glob(i) ];
   u_known(ind0) = ur(ind1);
end

%%
% Restrict the data on the given part (not necessary, but cleaner)
f_refer = f_known; u_refer = u_known;
f_known(tofindN) = 0; u_known(tofindD) = 0;

disp([ 'Direct problem solved and data management ', num2str(toc) ]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%% Compute the operator
Lx = 1;
%   Lx = max(nodes2(:,1)) - min(nodes2(:,1));
   
lambda = nu*E/((1+nu)*(1-2*nu));
mu     = E/(2*(1+nu));
   
if order ~= 1
   warning("order > 1 cases are not implemented");
end
Ng = Npg; nm3   = (nmax+1)^3;
[ Xg, Wg ] = gaussPt( Ng ); Ndots = size(Wg,1);

%% Build the list of Gauss points, and of construction-functions
Xxg = zeros( nboun2*Ndots,1 ); Yyg = zeros( nboun2*Ndots,1 );
Zzg = zeros( nboun2*Ndots,1 ); Wwg = zeros( nboun2*Ndots,1 );
Phi = sparse( 3*nboun2*Ndots, 3*nnodeboun );
Psi = sparse( 3*nboun2*Ndots, 3*nboun2 );
exnor = zeros( nboun2*Ndots,3); % Normal

for i=1:nboun2
   bonod = boundary2(i,:); exno = extnorm(i,:)';
   
   no1 = bonod(2); no2 = bonod(3); no3 = bonod(4); boname = bonod(1);
   x1 = nodes2(no1,1); y1 = nodes2(no1,2); z1 = nodes2(no1,3);
   x2 = nodes2(no2,1); y2 = nodes2(no2,2); z2 = nodes2(no2,3);
   x3 = nodes2(no3,1); y3 = nodes2(no3,2); z3 = nodes2(no3,3);
      
   nsurf = [ (y1-y2)*(z1-z3)-(y1-y3)*(z1-z2);... %
             (x1-x3)*(z1-z2)-(x1-x2)*(z1-z3);... % Vectorial product
             (x1-x2)*(y1-y3)-(x1-x3)*(y1-y2) ];  %
   surfa = .5*norm(nsurf);
     
   % Dofs in the numerotation of the boundary nodes
   N = boun2loc(i,:);
   indDtot = [ 3*N(1)-2, 3*N(1)-1, 3*N(1),...
               3*N(2)-2, 3*N(2)-1, 3*N(2),...
               3*N(3)-2, 3*N(3)-1, 3*N(3) ];
                  
   % Dofs in the numerotation of the boundary elements
   indNtot = [ 3*i-2, 3*i-1, 3*i ];
                
   for j=1:Ndots
      xg = Xg(j,:); wg = Wg(j);
      xgr  = (1-xg(1)-xg(2))*[x1;y1;z1] + xg(1)*[x2;y2;z2] + xg(2)*[x3;y3;z3] ; % abscissae

      X = xgr(1)/Lx; Y = xgr(2)/Lx; Z = xgr(3)/Lx; % Normalize
      Swg   = surfa * wg;

      indexg = Ndots*(i-1) + j;
      Xxg( indexg ) = X; Yyg( indexg ) = Y; Zzg( indexg ) = Z; Wwg( indexg ) = Swg;
      exnor( indexg, :) = exno';

      umxgx = 1-xg(1)-xg(2);
      Phi( 3*indexg-2, indDtot)  = [ umxgx, 0, 0, xg(1), 0, 0, xg(2), 0, 0 ];
      Phi( 3*indexg-1, indDtot)  = [ 0, umxgx, 0, 0, xg(1), 0, 0, xg(2), 0 ];
      Phi( 3*indexg  , indDtot)  = [ 0, 0, umxgx, 0, 0, xg(1), 0, 0, xg(2) ];

      Psi( 3*indexg-2, indNtot) = [1,0,0];
      Psi( 3*indexg-1, indNtot) = [0,1,0];
      Psi( 3*indexg  , indNtot) = [0,0,1];
   end
end

%% Build the matrix of test-functions
Vv = zeros( size(coef,1), 3*nboun2*Ndots );
Sv = zeros( size(coef,1), 3*nboun2*Ndots );
exnoX = exnor(:,1); exnoY = exnor(:,2); exnoZ = exnor(:,3);

for ii=0:nmax
   nd1 = nmax-ii;
   for jj=0:nd1
      nd2 = nmax-ii-jj;
      for kk=0:nd2
                  
         % Build the test field's basis
         if ii>0
            XYZi = ii*Xxg.^(ii-1).*Yyg.^jj.*Zzg.^kk;
            sloc11a = (lambda+2*mu)*XYZi;
            sloc22a = lambda*XYZi;
            sloc33a = lambda*XYZi;
            sloc12b = mu*XYZi;
            sloc13c = mu*XYZi;
         else
            sloc11a = zeros(nboun2*Ndots,1);
            sloc22a = zeros(nboun2*Ndots,1);
            sloc33a = zeros(nboun2*Ndots,1);
            sloc12b = zeros(nboun2*Ndots,1);
            sloc13c = zeros(nboun2*Ndots,1);
         end

         if jj>0
            XYZj = jj*Xxg.^ii.*Yyg.^(jj-1).*Zzg.^kk;
            sloc12a = mu*XYZj; 
            sloc11b = lambda*XYZj; 
            sloc22b = (lambda+2*mu)*XYZj;
            sloc33b = lambda*XYZj;
            sloc23c = mu*XYZj;
         else
            sloc12a = zeros(nboun2*Ndots,1);
            sloc11b = zeros(nboun2*Ndots,1);
            sloc22b = zeros(nboun2*Ndots,1);
            sloc33b = zeros(nboun2*Ndots,1);
            sloc23c = zeros(nboun2*Ndots,1);
         end

         if kk>0
            XYZk = kk*Xxg.^ii.*Yyg.^jj.*Zzg.^(kk-1);
            sloc13a = mu*XYZk; 
            sloc23b = mu*XYZk;
            sloc11c = lambda*XYZk; 
            sloc22c = lambda*XYZk;
            sloc33c = (lambda+2*mu)*XYZk;
         else
            sloc13a = zeros(nboun2*Ndots,1);
            sloc23b = zeros(nboun2*Ndots,1);
            sloc11c = zeros(nboun2*Ndots,1);
            sloc22c = zeros(nboun2*Ndots,1);
            sloc33c = zeros(nboun2*Ndots,1);
         end

         fpaax = sloc11a.*exnoX + sloc12a.*exnoY + sloc13a.*exnoZ;
         fpaay = sloc12a.*exnoX + sloc22a.*exnoY;% + sloc23a.*exnoZ;
         fpaaz = sloc13a.*exnoX + sloc33a.*exnoZ;% + sloc23a.*exnoY 

         fpabx = sloc11b.*exnoX + sloc12b.*exnoY;% + sloc13b.*exnoZ;
         fpaby = sloc12b.*exnoX + sloc22b.*exnoY + sloc23b.*exnoZ;
         fpabz = sloc23b.*exnoY + sloc33b.*exnoZ;% + sloc13b.*exnoX

         fpacx = sloc11c.*exnoX + sloc13c.*exnoZ;% + sloc12c.*exnoY
         fpacy = sloc22c.*exnoY + sloc23c.*exnoZ;% + sloc12c.*exnoX
         fpacz = sloc13c.*exnoX + sloc23c.*exnoY + sloc33c.*exnoZ;

         XYZ = Xxg.^ii.*Yyg.^jj.*Zzg.^kk;
         vpaax = XYZ; vpaby = XYZ; vpacz = XYZ;

         index = (nmax+1)^2*ii + (nmax+1)*jj + kk + 1;
                  
         Sv( index, 1:3:3*nboun2*Ndots-2 )         = Wwg' .* fpaax';
         Sv( nm3 + index, 1:3:3*nboun2*Ndots-2 )   = Wwg' .* fpabx';
         Sv( 2*nm3 + index, 1:3:3*nboun2*Ndots-2 ) = Wwg' .* fpacx';

         Sv( index, 2:3:3*nboun2*Ndots-1 )         = Wwg' .* fpaay';
         Sv( nm3 + index, 2:3:3*nboun2*Ndots-1 )   = Wwg' .* fpaby';
         Sv( 2*nm3 + index, 2:3:3*nboun2*Ndots-1 ) = Wwg' .* fpacy';

         Sv( index, 3:3:3*nboun2*Ndots )           = Wwg' .* fpaaz';
         Sv( nm3 + index, 3:3:3*nboun2*Ndots )     = Wwg' .* fpabz';
         Sv( 2*nm3 + index, 3:3:3*nboun2*Ndots )   = Wwg' .* fpacz';

         Vv( index, 1:3:3*nboun2*Ndots-2 )         = Wwg' .* vpaax';
         Vv( nm3 + index, 2:3:3*nboun2*Ndots-1 )   = Wwg' .* vpaby';
         Vv( 2*nm3 + index, 3:3:3*nboun2*Ndots )   = Wwg' .* vpacz';
      end
   end
end

Ruij = Sv*Phi; clear Sv; clear Phi;
Rfij = Vv*Psi; clear Vv; clear Psi;

Ru = coef'*Ruij; Rf = coef'*Rfij;
clear coef; clear Ruij; clear Rfij; % We can't afford not to clear

% Clear zeros lines
toremove = [];
for i=1:nftest
   if norm(Ru(i,:))==0 && norm(Rf(i,:))==0
      toremove = [toremove,i];
   end
end
Ru(toremove,:) = [];
Rf(toremove,:) = [];

disp([ 'Right hand side generated ', num2str(toc) ]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%% Restrict the operators
Rfm = Rf(:,tofindN); Rfr = Rf(:,knownN);
Rum = Ru(:,tofindD); Rur = Ru(:,knownD);

LhsA = -Rum;
LhsB = Rfm;
Rhs1 = Rur*u_known(knownD);% - Rfr*f_known(knownN); In practice, f_known = 0 on the known bounds

% Build the matrix that passes f on dofs on the nodes from the bound
Fbton = zeros( 3*nnodeboun, 3*nboun2 ); Fntob = zeros( 3*nboun2, 3*nnodeboun ); 
% And the boundary mass matrix
nodeMass = zeros(3*nnodeboun);  % For u
elemMass = zeros(3*nboun2);     % For f
% msur = 10000;  % ~= Infinity
for i=1:nboun2
   coef = [ boun2loc(i,1), boun2loc(i,2), boun2loc(i,3) ];
   
   no1 = boundary2(i,2); no2 = boundary2(i,3); no3 = boundary2(i,4);
   x1 = nodes2(no1,1); y1 = nodes2(no1,2); z1 = nodes2(no1,3);
   x2 = nodes2(no2,1); y2 = nodes2(no2,2); z2 = nodes2(no2,3);
   x3 = nodes2(no3,1); y3 = nodes2(no3,2); z3 = nodes2(no3,3);
   nsurf = [ (y1-y2)*(z1-z3)-(y1-y3)*(z1-z2);... %
             (x1-x3)*(z1-z2)-(x1-x2)*(z1-z3);... % Vectorial product
             (x1-x2)*(y1-y3)-(x1-x3)*(y1-y2) ];  %
   surfa = .5*norm(nsurf);
   
%    if surfa<msur msur = surfa; end % Just a debug asset
   
   Fntob( 3*i-2, 3*coef-2 ) = 1/(3*surfa) * [1,1,1];
   Fntob( 3*i-1, 3*coef-1 ) = 1/(3*surfa) * [1,1,1];
   Fntob( 3*i  , 3*coef   ) = 1/(3*surfa) * [1,1,1];
   Fbton( 3*coef-2, 3*i-2 ) = surfa/3 * [1,1,1];
   Fbton( 3*coef-1, 3*i-1 ) = surfa/3 * [1,1,1];
   Fbton( 3*coef  , 3*i   ) = surfa/3 * [1,1,1];
   
   ico = [ 3*i-2, 3*i-1, 3*i ];
   cco = [ 3*coef-2, 3*coef-1, 3*coef ];
   elemMass(ico,ico) = surfa*eye(3);
   nodeMass(cco,cco) = nodeMass(cco,cco) + ...
                       surfa/3 * [ eye(3), zeros(3), zeros(3) ; ...
                                   zeros(3), eye(3), zeros(3) ; ...
                                   zeros(3), zeros(3), eye(3) ];
end

% Extract interesting parts
Mfm = elemMass(tofindN,tofindN); Mfr = elemMass(knownN,knownN);
Mum = nodeMass(tofindD,tofindD); Mur = nodeMass(knownD,knownD);

% Derivative operators
Du = zeros(9*nboun2,3*nnodeboun);
Df = eye(3*nboun2); % No derivative for f
for i=1:nboun2
   coefU = [ boun2loc(i,1) , boun2loc(i,2) , boun2loc(i,3) ];
      
   no1 = boundary2(i,2); no2 = boundary2(i,3); no3 = boundary2(i,4);
   x1 = nodes2(no1,1); y1 = nodes2(no1,2); z1 = nodes2(no1,3);
   x2 = nodes2(no2,1); y2 = nodes2(no2,2); z2 = nodes2(no2,3);
   x3 = nodes2(no3,1); y3 = nodes2(no3,2); z3 = nodes2(no3,3);

   Jxy = [x2-x1,x3-x1;y2-y1,y3-y1];   %
   Jxz = [x2-x1,x3-x1;z2-z1,z3-z1];   % Jacobians
   Jyz = [y2-y1,y3-y1;z2-z1,z3-z1];   %
      
   D = [-1,1,0;-1,0,1];
%       D = [-1,1,0;-1,0,1;0,0,0];
%       Dexy = D; Dexz = D; Deyz = D; % No respect -- Vae victis
   Dexy = D'*pinv(Jxy); Dexz = D'*pinv(Jxz); Deyz = D'*pinv(Jyz);

%   Be = [ Dexy(1,1)+Dexz(1,1), Dexy(2,1)+Dexz(2,1), Dexy(3,1)+Dexz(3,1) ;
%          Dexy(1,2)+Deyz(1,1), Dexy(2,2)+Deyz(2,1), Dexy(3,2)+Deyz(3,1) ;
%          Dexz(1,2)+Deyz(1,2), Dexz(2,2)+Deyz(2,2), Dexz(3,2)+Deyz(3,2) ]; % (a .5 factor has been ommited)

   Be = [ Dexy(1,1), Dexy(2,1), Dexy(3,1) ;
          Dexy(1,2), Dexy(2,2), Dexy(3,2);
          0, 0, 0 ];

   indelem = [ 3*i-2, 3*i-1, 3*i ];
%   if norm(Du( 3*indelem-2, 3*coefU-2 ),'fro') ~= 0 || norm(Du( 3*indelem-1, 3*coefU-1 ),'fro') ~= 0 || norm(Du( 3*indelem, 3*coefU ),'fro') ~= 0
%      warning('nonul !');
%   end
   Du( 3*indelem-2, 3*coefU-2 ) = Be;
   Du( 3*indelem-1, 3*coefU-1 ) = Be;
   Du( 3*indelem  , 3*coefU   ) = Be;
end
   
Du0 = Du; Df0 = Df;
Du  = Du0(:,tofindD); Df  = Df0(tofindN,tofindN);
Duk = Du0(:,knownD);   Dfk = Df0(knownN,knownN);
% End derivative operator

%% Pure RG : Solve the linear system and recover the unknowns
if froreg == 1
   kB = norm(LhsA,'fro')/norm(LhsB,'fro'); % In order to regularize the stuff
else
   kB = 1;
end
Lhs = [LhsA,kB*LhsB];
   
Rhs = Rhs1;
A = Lhs'*Lhs; b = Lhs'*Rhs; sA = size(A,1);
   
if regular == 1
   Lu = Du'*Du; Lu = Lu/norm(Lu,'fro');
   L = [ Lu, zeros( size(Lu,1) , size(Df,2) ) ;...
         zeros( size(Df,1) , size(Lu,2) ) , Df ];
   L0 = L;
else
   L = [ Mum, zeros(size(Mum,1),size(Mfm,2)) ; zeros(size(Mfm,1),size(Mum,2)) , Mfm ]; % Weighted mass
end

Llu = [ Mum, zeros(size(Mum,1),size(Mfm,2)) ; zeros(size(Mfm,1),size(Mum,2)) , zeros(size(Mfm)) ];
Llf = [ zeros(size(Mum)), zeros(size(Mum,1),size(Mfm,2)) ; zeros(size(Mfm,1),size(Mum,2)) , Mfm ];

if testLcurve == 1
   mumoy = 1000*norm(Lhs'*Lhs, 'fro') / norm(L, 'fro');
   muloc = zeros(20,1); resl = zeros(20,1);
   norl = zeros(20,1); noru = zeros(20,1); norf = zeros(20,1);
   for i=1:30
      muloc(i) = mumoy/3^i;
      Soloc    = (Lhs'*Lhs + muloc(i)*L)\(Lhs'*Rhs1);
      resl(i)  = Soloc'*Lhs'*Lhs*Soloc - 2*Soloc'*Lhs'*Rhs1 + Rhs1'*Rhs1;
      norl(i)  = Soloc'*L*Soloc;
      noru(i)  = Soloc'*Llu*Soloc;
      norf(i)  = Soloc'*Llf*Soloc;
   end

   Solu1RG = (Lhs'*Lhs + mur*L)\(Lhs'*Rhs1);
   res0 = Solu1RG'*Lhs'*Lhs*Solu1RG - 2*Solu1RG'*Lhs'*Rhs1 + Rhs1'*Rhs1;
   norl0 = Solu1RG'*L*Solu1RG;
   norlu = Solu1RG'*Llu*Solu1RG;
   norlf = Solu1RG'*Llf*Solu1RG;

   figure;
   hold on;
   loglog(resl,norl,'Color','red','-*');
   loglog(res0,norl0,'Color','red','o','markersize',15);
   loglog(resl,noru,'Color','blue','-*');
   loglog(res0,norlu,'Color','blue','o','markersize',15);
   loglog(resl,norf,'Color','green','-*');
   loglog(res0,norlf,'Color','green','o','markersize',15);
   legend('L_{tot}','','L_u','','L_f','');

%   figure;
%   hold on;
%   plot(resl,norl,'Color','red','-*');
%   plot(res0,norl0,'Color','red','o','markersize',15);
%   plot(resl,noru,'Color','blue','-*');
%   plot(res0,norlu,'Color','blue','o','markersize',15);
%   plot(resl,norf,'Color','green','-*');
%   plot(res0,norlf,'Color','green','o','markersize',15);
%   legend('L_{tot}','','L_u','','L_f','');
else
   Solu1RG = (Lhs'*Lhs + mur*L)\(Lhs'*Rhs1);
%   Solu1RG = (Lhs'*Lhs + mur*Llu + mfr*Llf)\(Lhs'*Rhs1);
end

Solu1 = Solu1RG;

% Post-process : reconstruct u and f from Solu
fsolu1 = zeros(3*nboun2,1); usolu1 = zeros(3*nnodeboun,1);
usolu1(tofindD) = Solu1( 1:ndofD );
fsolu1(tofindN) = kB*Solu1( ndofD+1:ndofD+ndofN );
fsoluN = Fbton*fsolu1; f_refeN = Fbton*f_refer;

Usol = zeros(3*nnodes2,1); Fsol = zeros(3*nnodes2,1);
Uref = zeros(3*nnodes2,1); Fref = zeros(3*nnodes2,1);
for i=1:nnodeboun
   ind1 = [ 3*nodeboun2glob(i)-2 , 3*nodeboun2glob(i)-1 , 3*nodeboun2glob(i) ];
   ind2 = [ 3*i-2, 3*i-1, 3*i ];
   Usol(ind1) = usolu1(ind2); Fsol(ind1) = fsoluN(ind2);
   Uref(ind1) = u_refer(ind2); Fref(ind1) = f_refeN(ind2);
end

% Add known
FsolT = Fsol; %/!\ this is true because we do know that F=0 on the other points
UsolT = Usol; UsolT(knownD) = Uref(knownD);
Ferro = (FsolT-Fref)/max(abs(Fref)); Ferro = sqrt( Ferro(1:3:end-2).^2 + Ferro(2:3:end-1).^2 + Ferro(3:3:end).^2 );
Uerro = (UsolT-Uref)/max(abs(Uref)); Uerro = sqrt( Uerro(1:3:end-2).^2 + Uerro(2:3:end-1).^2 + Uerro(3:3:end).^2 );
plotGMSH3D({UsolT,'U_sol';FsolT,'F_sol';Uref,'U_ref';Fref,'F_ref';...
            Ferro,'error_f';Uerro,'error_u'}, elements2, nodes2, 'output/solution');
%
erroru = norm(UsolT-Uref) / norm(Uref);
errorf = norm(FsolT-Fref) / norm(Fref);

% Plot on a line
nstep = 100;
xmax = max(nodes2(:,1)); xmin = min(nodes2(:,1)); 
ymax = max(nodes2(:,2)); ymin = min(nodes2(:,2)); 
step = (ymax-ymin)/nstep;
X = .5; Y = ymin:step:ymax; Z = .29999;%.0001;%.29999;
nodesplo = [X*ones(nstep+1,1),Y',Z*ones(nstep+1,1)];
U1 = passMesh3D( nodes2, elements2, nodesplo, [], [UsolT,Uref] );
up = U1(:,1); upref = U1(:,2);

try
figure;
hold on;
set(gca, 'fontsize', 20);
plot(Y,up(3:3:end),'Color','red','LineWidth',3);
plot(Y,upref(3:3:end),'Color','blue','LineWidth',3);
legend('identified','reference');
end

%boun21i = find( boundary2(:,1) == 1 );
boun21i = find( boundary2(:,1) == 2 );
boun22 = boundary2( boun21i, 2:end );

figure;
patch('Faces',boun22,'Vertices',nodes2(:,1:2),'FaceVertexCData',UsolT(3:3:end),'FaceColor','interp');
colorbar(); set(colorbar, 'fontsize', 20);

figure;
patch('Faces',boun22,'Vertices',nodes2(:,1:2),'FaceVertexCData',Uref(3:3:end),'FaceColor','interp');
colorbar(); set(colorbar, 'fontsize', 20);
