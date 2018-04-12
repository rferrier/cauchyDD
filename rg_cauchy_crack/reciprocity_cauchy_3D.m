% 12/01/2018
% Problèmes directs et de Cauchy par écart à la réciprocité, problème 3D

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
jmaxRG     = 40;     % Eigenvalues truncation number
nmaxRG     = 200;    % nb of computed eigenvalues
jmaxRGSP   = 250;
nmaxRGSP   = 400;    % nb of computed eigenvalues
regular    = 0;      % Use the derivative regularization matrix (0 : Id, 1 : derivative, 2 : lumped)
upper_term = 0;      % 1 : use i=0:10, j=0:10, 0 : use i>=0,j>=0,i+j<=10
froreg     = 1;      % frobenius preconditioner
recompute  = 0;      % Recompute the operators
matrixfile = 'fields/rg_cauchy_crack/reciprocity3D_NG8.mat';  % File for the integration matrix
RGorSP     = 1;      % Use RG(1), SP(2) or mix(3)
Npg        = 1;      % Nb Gauss points
ordertest  = 20;     % Order of test fonctions

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
[ nodes2,elements2,ntoelem2,boundary2,order] = readmesh3D( 'meshes/rg3d_qcq/plate3d_charge2nt.msh' );
%[ nodes2,elements2,ntoelem2,boundary2,order] = readmesh3D( 'meshes/rg3d_qcq/plate3d_charge2_u.msh' );

boundary2( find(boundary2(:,1)==7) , : ) = []; % Hack : remove the double elements in boundary 2, (ie boundary 7)

nnodes2 = size(nodes2,1);
[K2,C2,nbloq2,node2c2,c2node2] = Krig3D(nodes2,elements2,mat,order,boundary2,dirichlet);
Kinter2 = K2( 1:2*nnodes2, 1:2*nnodes2 );

UR = passMesh3D( nodes, elements, nodes2, elements2, [uref,fref] );
ur = UR(:,1); fr = UR(:,2);

plotGMSH3D({ur,'U_bound';fr,'F_bound'}, elements2, nodes2, 'output/bound');

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
if recompute == 1
   tic
   
   Lx = 1;
%   Lx = max(nodes2(:,1)) - min(nodes2(:,1));
   
   lambda = nu*E/((1+nu)*(1-2*nu));
   mu     = E/(2*(1+nu));
   
   if order ~= 1
      warning("order > 1 cases are not implemented");
   end
   Ng = Npg;
   [ Xg, Wg ] = gaussPt( Ng );
   Ndots = size(Wg,1);
   nm3   = (nmax+1)^3;
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% Compute the RG
   Ruij = zeros(ncoef,3*nnodeboun); Rfij = zeros(ncoef,3*nboun2);
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
         umxgx = 1-xg(1)-xg(2); % Third component (optimization)
         Swg   = surfa * wg;
                  
         for ii=0:nmax
            nd1 = nmax-ii;
            for jj=0:nd1
               nd2 = nmax-ii-jj;
               for kk=0:nd2
                  %if ii+jj+kk > nmax continue; end % No need to add an other 0
                  
                  % Build the test field's basis
                  if ii>0
                     XYZi = ii*X^(ii-1)*Y^jj*Z^kk;
                     sloc11a = (lambda+2*mu)*XYZi;
                     sloc22a = lambda*XYZi;
                     sloc33a = lambda*XYZi;
                     sloc12b = mu*XYZi;
                     sloc13c = mu*XYZi;
                  else
                     sloc11a = 0; sloc22a = 0; sloc33a = 0; sloc12b = 0; sloc13c = 0;
                  end

                  if jj>0
                     XYZj = jj*X^ii*Y^(jj-1)*Z^kk;
                     sloc12a = mu*XYZj; 
                     sloc11b = lambda*XYZj; 
                     sloc22b = (lambda+2*mu)*XYZj;
                     sloc33b = lambda*XYZj;
                     sloc23c = mu*XYZj;
                  else
                     sloc12a = 0; sloc11b = 0; sloc22b = 0; sloc33b = 0; sloc23c = 0;
                  end

                  if kk>0
                     XYZk = kk*X^ii*Y^jj*Z^(kk-1);
                     sloc13a = mu*XYZk; 
                     sloc23b = mu*XYZk;
                     sloc11c = lambda*XYZk; 
                     sloc22c = lambda*XYZk;
                     sloc33c = (lambda+2*mu)*XYZk;
                  else
                     sloc13a = 0; sloc23b = 0; sloc11c = 0; sloc22c = 0; sloc33c = 0;
                  end

                  sloca = 1/Lx * [ sloc11a, sloc12a, sloc13a ;... % Lx or not Lx ?
                                   sloc12a, sloc22a, 0 ;...
                                   sloc13a, 0, sloc33a ];
                  slocb = 1/Lx * [ sloc11b, sloc12b, 0 ;...
                                   sloc12b, sloc22b, sloc23b ;...
                                   0, sloc23b, sloc33b ];
                  slocc = 1/Lx * [ sloc11c, 0, sloc13c ;...
                                   0, sloc22c, sloc23c ;...
                                   sloc13c, sloc23c, sloc33c ];

                  fpaa = sloca*exno; fpab = slocb*exno; fpac = slocc*exno;

                  XYZ = X^ii*Y^jj*Z^kk;
                  vpaa = [ XYZ ; 0 ; 0 ];
                  vpab = [ 0 ; XYZ ; 0 ];
                  vpac = [ 0 ; 0 ; XYZ ];

                  fpaTimesPhitest0a = [ fpaa(1)*umxgx ; fpaa(2)*umxgx ; fpaa(3)*umxgx ;...
                                        fpaa(1)*xg(1) ; fpaa(2)*xg(1) ; fpaa(3)*xg(1) ;...
                                        fpaa(1)*xg(2) ; fpaa(2)*xg(2) ; fpaa(3)*xg(2) ];
                  fpaTimesPhitest0b = [ fpab(1)*umxgx ; fpab(2)*umxgx ; fpab(3)*umxgx ;...
                                        fpab(1)*xg(1) ; fpab(2)*xg(1) ; fpab(3)*xg(1) ;...
                                        fpab(1)*xg(2) ; fpab(2)*xg(2) ; fpab(3)*xg(2) ];
                  fpaTimesPhitest0c = [ fpac(1)*umxgx ; fpac(2)*umxgx ; fpac(3)*umxgx ;...
                                        fpac(1)*xg(1) ; fpac(2)*xg(1) ; fpac(3)*xg(1) ;...
                                        fpac(1)*xg(2) ; fpac(2)*xg(2) ; fpac(3)*xg(2) ];

                  index = (nmax+1)^2*ii + (nmax+1)*jj + kk + 1;
                  
                  Ruij( index, indDtot ) = Ruij( index, indDtot ) + Swg * fpaTimesPhitest0a';
                  Ruij( nm3 + index, indDtot ) = Ruij( nm3 + index, indDtot ) ...
                                                 + Swg * fpaTimesPhitest0b';
                  Ruij( 2*nm3 + index, indDtot ) = Ruij( 2*nm3 + index, indDtot ) ...
                                                 + Swg * fpaTimesPhitest0c';

                  Rfij(index,indNtot)       = Rfij(index,indNtot) + Swg * vpaa';
                  Rfij(nm3+index,indNtot)   = Rfij(nm3+index,indNtot) + Swg * vpab';
                  Rfij(2*nm3+index,indNtot) = Rfij(2*nm3+index,indNtot) + Swg * vpac';
               end
            end
         end
   
      end
   end
   Ru = coef'*Ruij; Rf = coef'*Rfij;
   disp([ 'Right hand side generated ', num2str(toc) ]);
else
   Anb = load(matrixfile);
   Rf  = Anb.Rf; Ru  = Anb.Ru;
end
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% Derivative operators /!\ Buggy : see rg3d_Tychonov
Du = zeros(3*nnodeboun);
Df = eye(3*nboun2); % No derivative for f
for i=1:nboun2
   coefU = [ boun2loc(i,1) , boun2loc(i,2) , boun2loc(i,3) ];
%       coefV = [ 3*boun2loc(i,1)-2 , 3*boun2loc(i,2)-1 , 3*boun2loc(i,3) ,...
%                 3*boun2loc(i,1)-2 , 3*boun2loc(i,2)-1 , 3*boun2loc(i,3) ,...
%                 3*boun2loc(i,1)-2 , 3*boun2loc(i,2)-1 , 3*boun2loc(i,3) ];
      
   no1 = boundary2(i,2); no2 = boundary2(i,3); no3 = boundary2(i,4);
   x1 = nodes2(no1,1); y1 = nodes2(no1,2); z1 = nodes2(no1,3);
   x2 = nodes2(no2,1); y2 = nodes2(no2,2); z2 = nodes2(no2,3);
   x3 = nodes2(no3,1); y3 = nodes2(no3,2); z3 = nodes2(no3,3);

   Jxy = [x2-x1,x3-x1;y2-y1,y3-y1];   %
   Jxz = [x2-x1,x3-x1;z2-z1,z3-z1];   % Jacobians
   Jyz = [y2-y1,y3-y1;z2-z1,z3-z1];   %
      
   D    = [-1,1,0;-1,0,1];
%       D    = [-1,1,0;-1,0,1;0,0,0];
%       Dexy = D; Dexz = D; Deyz = D; % No respect -- Vae victis
   Dexy = D'*pinv(Jxy); Dexz = D'*pinv(Jxz); Deyz = D'*pinv(Jyz);

%       Be = [ Dexy(1,1)+Dexz(1,1), Dexy(2,1)+Dexz(2,1), Dexy(3,1)+Dexz(3,1) ;
%              Dexy(1,2)+Deyz(1,1), Dexy(2,2)+Deyz(2,1), Dexy(3,2)+Deyz(3,1) ;
%              Dexz(1,2)+Deyz(1,2), Dexz(2,2)+Deyz(2,2), Dexz(3,2)+Deyz(3,2) ]; % (a .5 factor has been ommited)

   Be = [ Dexy(1,1), Dexy(2,1), Dexy(3,1) ;
          Dexy(1,2), Dexy(2,2), Dexy(3,2);
          0, 0, 0 ];

%       Be = D';

   Du( 3*coefU-2, 3*coefU-2 ) = Du( 3*coefU-2, 3*coefU-2 ) + Be;
   Du( 3*coefU-1, 3*coefU-1 ) = Du( 3*coefU-1, 3*coefU-1 ) + Be;
   Du( 3*coefU  , 3*coefU   ) = Du( 3*coefU  , 3*coefU   ) + Be;
end
   
Du0 = Du; Df0 = Df;
Du  = Du0(tofindD,tofindD); Df  = Df0(tofindN,tofindN);
Duk = Du0(knownD,knownD);   Dfk = Df0(knownN,knownN);
% End derivative operator

if RGorSP == 1
   %% Pure RG : Solve the linear system and recover the unknowns
   jmax = jmaxRG;
   if froreg == 1
      kB = norm(LhsA,'fro')/norm(LhsB,'fro'); % In order to regularize the stuff
   else
      kB = 1;
   end
   Lhs = [LhsA,kB*LhsB];
   
   Rhs = Rhs1;
   A = Lhs'*Lhs; b = Lhs'*Rhs; sA = size(A,1);
   
   if regular == 1
      Dtot = [ Du/norm(Du,'fro'), zeros( size(Du,1) , size(Df,2) ) ;... %
               zeros( size(Df,1) , size(Du,2) ) , Df ];
      L = Dtot*Dtot'; L0 = L;
      %L = L + eye(size(A)); % Remove the small eignevalues of L
   else
%       L = eye(size(A));
      L = [ Mum, zeros(size(Mum,1),size(Mfm,2)) ; zeros(size(Mfm,1),size(Mum,2)) , Mfm ]; % Weighted mass
   end
   ninfty = 0;
   
   %L = eye(size(A));
   
   % [~,Theta,Q] = svd(Lhs,L); Q = Q*real((Q'*L*Q)^(-1/2)); Theta = Theta.^2;
   %[Q,Theta] = eig(Lhs'*L*Lhs); Q = Q*real((Q'*L*Q)^(-1/2)); Theta = Theta.^2;
   [Q,Theta] = gradEig(A, L, nmaxRG, 1, Lhs'*Rhs1, 0); Q = real(Q);%Q = Q*real((Q'*L*Q)^(-1/2));
   % [Q,Theta] = eigs(A, L, nmaxRG); Q = real(Q);%Q = Q*real((Q'*L*Q)^(-1/2));
   % [Q,Theta] = eig(A,L); Q = Q*real((Q'*L*Q)^(-1/2)); % real should not be, but you know, numerical shit...
   thetas = diag(Theta);
   [thetas,Ind] = sort( thetas,'descend' );
   Q = Q(:,Ind);
   Thetas = diag(thetas); 
   Theta = Theta(Ind,Ind);
   
   disp([ num2str(size(A,1)), ' - dofs rectangular system pinversed ', num2str(toc) ]);
   % Plot the Picard stuff
   imax = min( find(thetas/thetas(ninfty+1)<1e-16) );
   if size(imax,1) == 0
       imax = size(thetas,1);
   end
   
   tplo = thetas(ninfty+1:imax); bplo1 = Q'*b; bplo1 = bplo1(ninfty+1:imax);
   rplo1 = (Q'*Lhs'*Rhs1)./thetas; rplo1 = rplo1(1:imax);
   figure
   hold on;
   plot(log10(abs(tplo)),'Color','green');
   plot(log10(abs(bplo1)),'Color','red');
   plot(log10(abs(rplo1)),'Color','black');
   legend('Singular values','Rhs1','sol1');
   
   % Filter eigenvalues
   if jmax == 0
      jmax = size(Thetas,1);
   end
   ThetaT = Thetas( ninfty+1:jmax , ninfty+1:jmax );
   bT1 = Q'*Lhs'*Rhs1; bT1 = bT1(ninfty+1:jmax);
   Solu1RG = Q(:,ninfty+1:jmax) * (ThetaT\bT1);

   nori  = zeros(nmaxRG,1); resi  = zeros(nmaxRG,1);
   Dtoti = [ Du/norm(Du,'fro'), zeros( size(Du,1) , size(Df,2) ) ;... %
               zeros( size(Df,1) , size(Du,2) ) , Df ];
   Li = Dtoti*Dtoti';
   for i=1:nmaxRG
      Thetai = Thetas( ninfty+1:i , ninfty+1:i );
      bTi = Q'*Lhs'*Rhs1; bTi = bTi(ninfty+1:i);
      Solu1i = Q(:,ninfty+1:i) * (Thetai\bTi);
      nori(i) = Solu1i'*Li*Solu1i;
      resi(i) = norm(Lhs*Solu1i - Rhs1);
   end

   try
   figure
   hold on;
   loglog(resi,nori,'Color','red','-*');
   loglog(resi(40),nori(40),'Color','red','o','markersize',15);
   end

end
if RGorSP == 3
   % SP+RG : Solve the linear system and recover the unknowns
   tic
   jmax = jmaxRGSP;
   Rhs1 = Rur*u_known(knownD) - Rfr*f_known(knownN);
   if froreg == 1
     kB = norm(Rum,'fro')/norm(Rfm,'fro'); % In order to regularize the stuff
   else
     kB = 1;
   end
   Lhs  = [ -Rum,kB*Rfm,zeros(size(Rur)),zeros(size(Rfr));...
           -Rum,kB*Rfm,zeros(size(Rur)),kB*Rfr ;...
           -Rum,kB*Rfm,-Rur,zeros(size(Rfr)) ];
   %Rhs1 = [ Rur*u_known(knownD) - Rfr*f_known(knownN) ;...
   %         Rur*u_known(knownD) ; -Rfr*f_known(knownN) ];
   Rhs1 = [ Rur*u_known(knownD) ; Rur*u_known(knownD) ; zeros(size(Rfr,1),1) ];
   
   A = Lhs'*Lhs;
   
   if regular == 1 || regular == 2
     Zuf   = zeros( size(Du,1) , size(Df,2) );
     Zuuk  = zeros( size(Du,1) , size(Duk,2) );
     Zufk  = zeros( size(Du,1) , size(Dfk,2) );
     Zfuk  = zeros( size(Df,1) , size(Duk,2) );
     Zffk  = zeros( size(Df,1) , size(Dfk,2) );
     Zukfk = zeros( size(Duk,1) , size(Dfk,2) );
     Dtot = [ Du ,Zuf , Zuuk, Zufk ;...
              Zuf', Df , Zfuk, Zffk ;...
              Zuuk', Zfuk', Duk, Zukfk ;...
              Zufk', Zffk', Zukfk', Dfk ];
     L = Dtot*Dtot';
   else
     L = eye(size(A));
   end
   
   ninfty = 0;
   
   [Q,Theta] = gradEig(A, L, nmaxRGSP, 1, Lhs'*Rhs1, 0 ); %Q = Q*real((Q'*L*Q)^(-1/2));
   %[Q,Theta] = eig(A,L); Q = Q*real((Q'*L*Q)^(-1/2)); % real should not be, but you know, numerical shit...
   thetas = diag(Theta);
   [thetas,Ind] = sort( thetas,'descend' );
   Q = Q(:,Ind);
   Thetas = diag(thetas); 
   Theta = Theta(Ind,Ind);
   
   disp([ num2str(size(A,1)), ' - dofs rectangular system pinversed ', num2str(toc) ]);
   % Plot the Picard stuff
   imax = min( find(thetas/thetas(ninfty+1)<1e-16) );
   if size(imax,1) == 0
      imax = size(thetas,1);
   end
   
   tplo = thetas(ninfty+1:imax); bplo1 = Q'*Lhs'*Rhs1; bplo1 = bplo1(ninfty+1:imax);
   rplo1 = (Q'*Lhs'*Rhs1)./thetas; rplo1 = rplo1(1:imax);
   figure
   hold on;
   plot(log10(abs(tplo)),'Color','green');
   plot(log10(abs(bplo1)),'Color','red');
   plot(log10(abs(rplo1)),'Color','black');
   legend('Singular values','Rhs1','sol1');
   
   % Filter eigenvalues
   if jmax == 0
     jmax = size(Thetas,1);
   end
   ThetaT = Thetas( ninfty+1:jmax , ninfty+1:jmax );
   
   bT1 = Q'*Lhs'*Rhs1; bT1 = bT1(ninfty+1:jmax);
   Solu1RGSP = Q(:,ninfty+1:jmax) * (ThetaT\bT1);
end

%% Choose the solution
if RGorSP == 1
   Solu1 = Solu1RG;
elseif RGorSP == 2
   Solu1 = Solu1SP;
else %RGorSP == 3
   Solu1 = Solu1RGSP;
end
%%

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
X = .5; Y = ymin:step:ymax; Z = .29999;
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

%% Plot on surfaces
%nstep = 100;
%xmax = max(nodes2(:,1)); xmin = min(nodes2(:,1)); 
%ymax = max(nodes2(:,2)); ymin = min(nodes2(:,2)); 
%stepx = (xmax-xmin)/nstep; stepy = (ymax-ymin)/nstep;
%X = xmin:stepx:xmax; Y = ymin:stepx:ymax; Z = .29999;
%Xt = X'*ones(1,nstep+1); Yt = ones(1,nstep+1)'*Y;
%nodesplo = [ Xt(:), Yt(:), Z*ones((nstep+1)^2,1) ];
%U1 = passMesh3D( nodes2, elements2, nodesplo, [], [UsolT,Uref] ); % CPU costly
%up = U1(:,1); upref = U1(:,2);

%upx = up(3:3:end); uprefx = upref(3:3:end);
%upx = reshape(upx,   [nstep+1,nstep+1]);
%% upx = .5*( upx - abs(upx) ); % Negative part
%uprefx = reshape(uprefx,[nstep+1,nstep+1]);

%try
%figure;
%hold on;
%set(gca, 'fontsize', 20);
%surf(X,Y,upx');
%shading interp;
%colorbar(); set(colorbar, 'fontsize', 20);
%end

%try
%figure;
%hold on;
%set(gca, 'fontsize', 20);
%surf(X,Y,uprefx');
%shading interp;
%colorbar(); set(colorbar, 'fontsize', 20);
%end

boun21i = find( boundary2(:,1) == 2 );
boun22 = boundary2( boun21i, 2:end );

figure;
patch('Faces',boun22,'Vertices',nodes2(:,1:2),'FaceVertexCData',UsolT(3:3:end),'FaceColor','interp');
colorbar(); set(colorbar, 'fontsize', 20);

figure;
patch('Faces',boun22,'Vertices',nodes2(:,1:2),'FaceVertexCData',Uref(3:3:end),'FaceColor','interp');
colorbar(); set(colorbar, 'fontsize', 20);
