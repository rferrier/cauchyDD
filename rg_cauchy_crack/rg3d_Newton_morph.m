% 29/10/2018
% ID fissure 3D par Newton -- Morphing de maillage

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
mur        = 1e2;      % Regularization parameter
regular    = 1;      % Use the derivative regularization matrix (0 : Id, 1 : derivative)
froreg     = 1;      % frobenius preconditioner
Npg        = 2;      % Nb Gauss points
ordertest  = 20;     % Order of test fonctions
niter      = 1;      % Nb of iterations in the Newton algorithm
%init       = [0;0;0;0];
init       = [0;0;-1;.5];%  initialization for the plane parameters. If its norm is 0 : use the reference plane
zerobound  = 1;      % Put the boundaries of the crack to 0
step       = 0;%1e-2;   % Step for the finite differences 0 means auto-adptation
nuzawa     = 100;    % Nb of Uzawa iterations
lc1        = .1;     % Size of the mesh of the plane
nelemax    = 4500;   % Max nb of elements in the matrix (for memory usage) (30 Go ram => 5500 max)
schur      = 1;      % Use the Schur complement ?
testkase   = 2;

nbDirichlet = [];
%nbDirichlet = [ 1,20 ; 2,20 ; 3,20 ; 4,20 ; 5,20 ; 6,20 ];

% Boundary conditions
%dirichlet  = [ 0,1,0 ; 0,2,0 ; 0,3,0 ; 0,4,0 ; 0,5,0 ; 0,6,0 ];
dirichlet  = [ 2,1,0 ; 2,2,0 ; 2,3,0 ];

neumann    = {};

neumann{1} = [ 4,1,fscalar ; 6,1,-fscalar ];
neumann{2} = [ 5,2,fscalar ; 3,2,-fscalar ];
neumann{3} = [ 1,3,-fscalar ];

neumann{4} = [ 4,1,fscalar ; 4,2,fscalar ; 5,1,fscalar ; 5,2,fscalar ; ...
               3,1,-fscalar ; 3,2,-fscalar ; 6,1,-fscalar ; 6,2,-fscalar ];
neumann{5} = [ 4,1,fscalar ; 4,2,-fscalar ; 3,1,fscalar ; 3,2,-fscalar ; ...
               5,1,-fscalar ; 5,2,fscalar ; 6,1,-fscalar ; 6,2,fscalar ];
neumann{6} = [ 1,1,-fscalar ; 1,3,-fscalar ; 6,1,-fscalar ; 6,3,-fscalar ];
neumann{7} = [ 4,1,fscalar ; 4,3,-fscalar ; 1,1,fscalar ; 1,3,-fscalar ];
neumann{8} = [ 3,3,-fscalar ; 3,2,-fscalar ; 1,3,-fscalar ; 1,2,-fscalar ];
neumann{9} = [ 1,3,-fscalar ; 1,2,fscalar ; 5,3,-fscalar ; 5,2,fscalar ];

neumann{10} = [ 1,1,-fscalar ; 1,2,-fscalar ; 1,3,-fscalar ;...
                3,1,-fscalar ; 3,2,-fscalar ; 3,3,-fscalar ;...
                6,1,-fscalar ; 6,2,-fscalar ; 6,3,-fscalar ];
neumann{11} = [ 1,1,fscalar ; 1,2,fscalar ; 1,3,-fscalar ;...
                4,1,fscalar ; 4,2,fscalar ; 4,3,-fscalar ;...
                5,1,fscalar ; 5,2,fscalar ; 5,3,-fscalar ];
neumann{12} = [ 1,1,-fscalar ; 1,2,fscalar ; 1,3,-fscalar ;...
                5,1,-fscalar ; 5,2,fscalar ; 5,3,-fscalar ;...
                6,1,-fscalar ; 6,2,fscalar ; 6,3,-fscalar ];
neumann{13} = [ 1,1,fscalar ; 1,2,-fscalar ; 1,3,-fscalar ;...
                3,1,fscalar ; 3,2,-fscalar ; 3,3,-fscalar ;...
                4,1,fscalar ; 4,2,-fscalar ; 4,3,-fscalar ];

%neumann{1} = [ 4,1,fscalar ; 6,1,-fscalar ];
%neumann{2} = [ 5,2,fscalar ; 3,2,-fscalar ];
%neumann{3} = [ 2,3,fscalar ; 1,3,-fscalar ];
%
%neumann{4} = [ 4,1,fscalar ; 4,2,fscalar ; 5,1,fscalar ; 5,2,fscalar ; ...
%               3,1,-fscalar ; 3,2,-fscalar ; 6,1,-fscalar ; 6,2,-fscalar ];
%neumann{5} = [ 4,1,fscalar ; 4,2,-fscalar ; 3,1,fscalar ; 3,2,-fscalar ; ...
%               5,1,-fscalar ; 5,2,fscalar ; 6,1,-fscalar ; 6,2,fscalar ];
%neumann{6} = [ 4,1,fscalar ; 4,3,fscalar ; 2,1,fscalar ; 2,3,fscalar ; ...
%               1,1,-fscalar ; 1,3,-fscalar ; 6,1,-fscalar ; 6,3,-fscalar ];
%neumann{7} = [ 4,1,fscalar ; 4,3,-fscalar ; 1,1,fscalar ; 1,3,-fscalar ; ...
%               2,1,-fscalar ; 2,3,fscalar ; 6,1,-fscalar ; 6,3,fscalar ];
%neumann{8} = [ 2,3,fscalar ; 2,2,fscalar ; 5,3,fscalar ; 5,2,fscalar ; ...
%               3,3,-fscalar ; 3,2,-fscalar ; 1,3,-fscalar ; 1,2,-fscalar ];
%neumann{9} = [ 1,3,-fscalar ; 1,2,fscalar ; 5,3,-fscalar ; 5,2,fscalar ; ...
%               2,3,fscalar ; 2,2,-fscalar ; 3,3,fscalar ; 3,2,-fscalar ];
%
%neumann{10} = [ 2,1,fscalar ; 2,2,fscalar ; 2,3,fscalar ;...
%                4,1,fscalar ; 4,2,fscalar ; 4,3,fscalar ;...
%                5,1,fscalar ; 5,2,fscalar ; 5,3,fscalar ;...
%                1,1,-fscalar ; 1,2,-fscalar ; 1,3,-fscalar ;...
%                3,1,-fscalar ; 3,2,-fscalar ; 3,3,-fscalar ;...
%                6,1,-fscalar ; 6,2,-fscalar ; 6,3,-fscalar ];
%neumann{11} = [ 1,1,fscalar ; 1,2,fscalar ; 1,3,-fscalar ;...
%                4,1,fscalar ; 4,2,fscalar ; 4,3,-fscalar ;...
%                5,1,fscalar ; 5,2,fscalar ; 5,3,-fscalar ;...
%                2,1,-fscalar ; 2,2,-fscalar ; 2,3,fscalar ;...
%                3,1,-fscalar ; 3,2,-fscalar ; 3,3,fscalar ;...
%                6,1,-fscalar ; 6,2,-fscalar ; 6,3,fscalar ];
%neumann{12} = [ 2,1,fscalar ; 2,2,-fscalar ; 2,3,fscalar ;...
%                4,1,fscalar ; 4,2,-fscalar ; 4,3,fscalar ;...
%                3,1,fscalar ; 3,2,-fscalar ; 3,3,fscalar ;...
%                1,1,-fscalar ; 1,2,fscalar ; 1,3,-fscalar ;...
%                5,1,-fscalar ; 5,2,fscalar ; 5,3,-fscalar ;...
%                6,1,-fscalar ; 6,2,fscalar ; 6,3,-fscalar ];
%neumann{13} = [ 2,1,-fscalar ; 2,2,fscalar ; 2,3,fscalar ;...
%                6,1,-fscalar ; 6,2,fscalar ; 6,3,fscalar ;...
%                5,1,-fscalar ; 5,2,fscalar ; 5,3,fscalar ;...
%                1,1,fscalar ; 1,2,-fscalar ; 1,3,-fscalar ;...
%                3,1,fscalar ; 3,2,-fscalar ; 3,3,-fscalar ;...
%                4,1,fscalar ; 4,2,-fscalar ; 4,3,-fscalar ];

% BCs for the stuff
%dirichlet0 = [ 1,1,0 ; 1,2,0 ; 1,3,0 ; 
%               3,1,0 ; 3,2,0 ; 3,3,0 ; 
%               4,1,0 ; 4,2,0 ; 4,3,0 ; 
%               5,1,0 ; 5,2,0 ; 5,3,0 ; 
%               6,1,0 ; 6,2,0 ; 6,3,0 ];
%neumann0   = [ 1,1,0 ; 1,2,0 ; 1,3,0 ; 
%               3,1,0 ; 3,2,0 ; 3,3,0 ; 
%               4,1,0 ; 4,2,0 ; 4,3,0 ; 
%               5,1,0 ; 5,2,0 ; 5,3,0 ; 
%               6,1,0 ; 6,2,0 ; 6,3,0 ];

dirichlet0 = [ 1,1,0 ; 1,2,0 ; 1,3,0 ; 
               2,1,0 ; 2,2,0 ; 2,3,0 ; 
               3,1,0 ; 3,2,0 ; 3,3,0 ; 
               4,1,0 ; 4,2,0 ; 4,3,0 ; 
               5,1,0 ; 5,2,0 ; 5,3,0 ; 
               6,1,0 ; 6,2,0 ; 6,3,0 ];
neumann0   = [ 1,1,0 ; 1,2,0 ; 1,3,0 ; 
               3,1,0 ; 3,2,0 ; 3,3,0 ; 
               4,1,0 ; 4,2,0 ; 4,3,0 ; 
               5,1,0 ; 5,2,0 ; 5,3,0 ; 
               6,1,0 ; 6,2,0 ; 6,3,0 ];

%dirichlet0 = [ 1,1,0 ; 1,2,0 ; 1,3,0 ; 
%               3,1,0 ; 3,2,0 ; 3,3,0 ; 
%               4,1,0 ; 4,2,0 ; 4,3,0 ; 
%               5,1,0 ; 5,2,0 ; 5,3,0 ; 
%               6,1,0 ; 6,2,0 ; 6,3,0 ];
%neumann0   = [ 1,1,0 ; 1,2,0 ; 1,3,0 ; 
%               2,1,0 ; 2,2,0 ; 2,3,0 ; 
%               3,1,0 ; 3,2,0 ; 3,3,0 ; 
%               4,1,0 ; 4,2,0 ; 4,3,0 ; 
%               5,1,0 ; 5,2,0 ; 5,3,0 ; 
%               6,1,0 ; 6,2,0 ; 6,3,0 ];
               
%dirichlet0 = [ 1,1,0 ; 1,2,0 ; 1,3,0 ; 
%               3,1,0 ; 3,2,0 ; 3,3,0 ; 
%               4,1,0 ; 4,2,0 ; 4,3,0 ];
%neumann0   = [ 1,1,0 ; 1,2,0 ; 1,3,0 ; 
%               2,1,0 ; 2,2,0 ; 2,3,0 ; 
%               3,1,0 ; 3,2,0 ; 3,3,0 ; 
%               4,1,0 ; 4,2,0 ; 4,3,0 ; 
%               5,1,0 ; 5,2,0 ; 5,3,0 ; 
%               6,1,0 ; 6,2,0 ; 6,3,0 ];
            
%dirichlet0 = [ 1,1,0 ; 1,2,0 ; 1,3,0 ; 
%               2,1,0 ; 2,2,0 ; 2,3,0 ; 
%               3,1,0 ; 3,2,0 ; 3,3,0 ; 
%               4,1,0 ; 4,2,0 ; 4,3,0 ; 
%               5,1,0 ; 5,2,0 ; 5,3,0 ; 
%               6,1,0 ; 6,2,0 ; 6,3,0 ];
%neumann0   = [ 1,1,0 ; 1,2,0 ; 1,3,0 ; 
%               2,1,0 ; 2,2,0 ; 2,3,0 ; 
%               3,1,0 ; 3,2,0 ; 3,3,0 ; 
%               4,1,0 ; 4,2,0 ; 4,3,0 ; 
%               5,1,0 ; 5,2,0 ; 5,3,0 ; 
%               6,1,0 ; 6,2,0 ; 6,3,0 ];

% Import the mesh
%[ nodes,elements,ntoelem,boundary,order] = readmesh3D( 'meshes/rg3d_crack/plate3d_ur4.msh' );
%[ nodes,elements,ntoelem,boundary,order ] = readmesh3D( 'meshes/rg3d_crack/plate3d_crack0.msh' );
if testkase==1
   [ nodes,elements,ntoelem,boundary,order ] = readmesh3D( 'meshes/rg3d_crack/plate3d_crack1_bis.msh' );
else
   [ nodes,elements,ntoelem,boundary,order ] = readmesh3D( 'meshes/rg3d_crack/plate3d_crack_large.msh' );
end
%[ nodes,elements,ntoelem,boundary,order ] = readmesh3D( 'meshes/rg3d_crack/plate3d_noplane.msh' );
%[ nodes,elements,ntoelem,boundary,order ] = readmesh3D( 'meshes/rg3d_crack/plate3d_smile.msh' );
nnodes = size(nodes,1);

%% Lagrange
%[K,C,nbloq] = Krig3D (nodes,elements,mat,order,boundary,dirichlet);
%Kinter = K(1:3*nnodes, 1:3*nnodes);
%f = zeros(3*nnodes+nbloq,13);
%for i = 1:13
%   f(:,i) = loading3D(nbloq,nodes,boundary,neumann{i});
%end
%
%uin = K\f; uref = uin(1:3*nnodes,:); fref = Kinter*uref;

% % Substitution
% Kinter = Krig3D (nodes,elements,mat,order,boundary,[]);
% C = Cbound2 ( nodes, dirichlet, boundary ); index = find(diag(C*C'));
% 
% indexdex = setdiff(1:3*nnodes,index);
% 
% f = zeros(3*nnodes,13);
% for i = 1:13
% f(:,i) = loading3D(0,nodes,boundary,neumann{i});
% end
% 
% uref = zeros(3*nnodes,13);% relres = zeros(13,1);
%% %tic
% PJag = diag(Kinter); 
% PJac = sparse( 1:3*nnodes, 1:3*nnodes, PJag, 3*nnodes, 3*nnodes ); % Jacobi
% for i=1:13
%    [uref(indexdex,i),flag,relres(i),iter] = pcg( Kinter(indexdex,indexdex), f(indexdex,i), 1e-6, 1000, PJac(indexdex,indexdex) );
% end
% %toc
% relres
%%
% Kinter = .5*(Kinter+Kinter');
% uref(indexdex,:) = Kinter(indexdex,indexdex)\f(indexdex,:);
%%
if testkase==1
   urefl = load('fields/rg3d_Newton_dirichlet.mat');
else
   urefl = load('fields/rg3d_Newton_dirichlet_large.mat');
end
uref = urefl.uref;  fref = urefl.fref;
%%
% fref = Kinter*uref;
%%
% clear Kinter; %clear K; 
% sigma = stress3D(uref,mat,nodes,elements,order,1,ntoelem);
% toplot = { uref(:,1), 'U1' ; uref(:,2), 'U2' ; uref(:,3), 'U3' ;...
%           uref(:,4), 'U4' ; uref(:,5), 'U5' ; uref(:,6), 'U6' ;...
%           uref(:,7), 'U7' ; uref(:,8), 'U8' ; uref(:,9), 'U9' ;...
%           uref(:,10), 'U10' ; uref(:,11), 'U11' ; uref(:,12), 'U12' ;...
%           uref(:,13), 'U13' };
% plotGMSH3D( toplot, elements, nodes, 'output/reference' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ nodes2,elements2,ntoelem2,boundary2,order ] = readmesh3D( 'meshes/rg3d_crack/plate3d_crack1-2.msh' );
%[ nodes2,elements2,ntoelem2,boundary2,order ] = readmesh3D( 'meshes/rg3d_crack/plate3d_crack0.msh' );
%[nodes2,elements2,ntoelem2,boundary2,order] = readmesh3D( 'meshes/rg3d_crack/plate3d_crack2.msh' );
%[ nodes2,elements2,ntoelem2,boundary2,order] = readmesh3D( 'meshes/rg3d_crack/plate3d_u.msh' );

boundary2( find(boundary2(:,1)==7) , : ) = []; % Hack : remove the crack (ie boundary 7)

nnodes2 = size(nodes2,1);
[ ~, nbloq2 ] = Cbound2 ( nodes2, dirichlet, boundary2 );
% [K2,C2,nbloq2,node2c2,c2node2] = Krig3D(nodes2,elements2,mat,order,boundary2,dirichlet);
% Kinter2 = K2( 1:3*nnodes2, 1:3*nnodes2 );
% 
% % Slight contraction of the second mesh for the passmesh
% mean2x = mean(nodes2(:,1)); mean2y = mean(nodes2(:,2)); mean2z = mean(nodes2(:,3));
% nodes2c = nodes2;
% nodes2c(:,1) = (1-1e-6)*(nodes2c(:,1)-mean2x) + mean2x;
% nodes2c(:,2) = (1-1e-6)*(nodes2c(:,2)-mean2y) + mean2y;
% nodes2c(:,3) = (1-1e-6)*(nodes2c(:,3)-mean2z) + mean2z;
% 
% ur = passMesh3D( nodes, elements, nodes2c, elements2, uref, boundary2 );
% fr = Kinter2*ur;
% clear K2; clear Kinter2; clear C2;

ur = urefl.ur; fr = urefl.fr;

% toplot = { ur(:,1), 'U1' ; ur(:,2), 'U2' ; ur(:,3), 'U3' ;...
%            ur(:,4), 'U4' ; ur(:,5), 'U5' ; ur(:,6), 'U6' ;...
%            ur(:,7), 'U7' ; ur(:,8), 'U8' ; ur(:,9), 'U9' ;...
%            ur(:,10), 'U10' ; ur(:,11), 'U11' ; ur(:,12), 'U12' ;...
%            ur(:,13), 'U13' };
% plotGMSH3D(toplot, elements2, nodes2, 'output/bound');

% Recover the force on the top boundary (Dirichlet)
%neumann{3} = [ 1,3,-fscalar ];
for i = 1:13
   fi(:,i) = loading3D(nbloq2,nodes2,boundary2,neumann{i});
end
fi = fr-fi(1:3*nnodes2,:); % Substract the other loads
fe2 = elem_forces(fi,nodes2,boundary2,2,order); % And put it on the edges of bound 2

%try % Debug plot
%figure;
%hold on;
%patch('Faces',boundary2(:,2:4),'Vertices',nodes2,'FaceVertexCData',fi(3:3:end,3),'FaceColor','interp');
%colorbar(); set(colorbar, 'fontsize', 20); axis([0,1,0,1,0,1],'square');
%end

%% Discrete measure point stuff
if min(size(nbDirichlet))>0
   nmin = [];
   for i=1:6
      nbdiscrete = find(nbDirichlet(:,1)==i);
      if min(size(nbdiscrete)) > 0
         nbi = nbDirichlet(nbdiscrete,2);
         if nbi>0
            [~, b2node] = mapBound(i, boundary2, nnodes2);
            xmin = min(nodes2(b2node,1)); xmax = max(nodes2(b2node,1)); Lx1 = xmax-xmin; sx = Lx1/(nbi-1);
            ymin = min(nodes2(b2node,2)); ymax = max(nodes2(b2node,2)); Ly1 = ymax-ymin; sy = Ly1/(nbi-1);
            zmin = min(nodes2(b2node,3)); zmax = max(nodes2(b2node,3)); Lz1 = zmax-zmin; sz = Lz1/(nbi-1);
            if sx==0
               y = ymin:sy:ymax; z = zmin:sz:zmax;
               y1 = y'*ones(1,nbi); z1 = ones(nbi,1)*z;
               xyz = [ xmin*ones(1,nbi^2) ; y1(:)' ; z1(:)' ];
            elseif sy==0
               x = xmin:sx:xmax; z = zmin:sz:zmax;
               x1 = x'*ones(1,nbi); z1 = ones(nbi,1)*z;
               xyz = [ x1(:)' ; ymin*ones(1,nbi^2) ; z1(:)' ];
            elseif sz==0
               x = xmin:sx:xmax; y = ymin:sy:ymax;
               x1 = x'*ones(1,nbi); y1 = ones(nbi,1)*y;
               xyz = [ x1(:)' ; y1(:)' ; zmin*ones(1,nbi^2) ];
            else
               warning('Discrete measurement unimplemented for non parallelepideic domains');
            end
            nm = zeros(nbi^2,1);
            for j=1:nbi^2 % Find the closest nodes
               xyj = xyz(:,j)';
               nodiff = nodes2(b2node,:)-xyj;
               normdiff = nodiff(:,1).^2+nodiff(:,2).^2+nodiff(:,3).^2;
               [~,nm(j)] = min(normdiff);
            end
            nm = unique(nm); nmin = [nmin;b2node(nm)];
         end
      end
   end
   nmin = unique(nmin); % In the global volumic numerotation
else
   nmin = 1:nnodes2;
end

% Add the noise
un = ur; am = zeros(1,13);
% br1 = randn(3*nnodes2,13);
noise = load('noises/cauchyRG3d.mat');
br1 = noise.br1;
for i=1:13
   am(i)   = sqrt(mean(ur(:,i).^2));
   ur(:,i) = ur(:,i) + br*am(i)*br1(:,i);
end

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

% Pass ponctual Dirichlet to boundary
nmin1 = nodeglob2boun(nmin);
nmin1 = nmin1(find(nmin1)); % Remove the zeros (in case interior nodes are in nmin)
nmin1 = [3*nmin1-2;3*nmin1-1;3*nmin1];

nodekD = [];
for i=1:size(knownD,1)
   dof = mod(knownD(i)+2,3)+1; % No of the dof (1=X, 2=Y, 3=Z)
   B = ceil(knownD(i)/3); % No of the boundary element
   N = boun2loc(B,:);
   nodekD = [ nodekD, 3*N-3+dof ];
end
nodekD  = unique(nodekD); knownD = nodekD;  % Erase knownD for name reasons
knownD  = intersect( nmin1', knownD );      % Apply the discrete stuff
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
f_known = zeros(3*nboun2,13); u_known = zeros(3*nnodeboun,13);
for i=1:nboun2
   % Reference force from the BC's
   name = boundary2(i,1);
   fer = zeros(3,13);
   for k=1:13
      nm = neumann{k};
      loads = find(name == nm(:,1));
      if min(size(loads)) > 0
         for j=1:max(size(loads))
            fer(nm(loads(j),2),k) = nm(loads(j),3);
         end
      end
   end
   ind0            = [ 3*i-2, 3*i-1, 3*i ];
   f_known(ind0,:) = fer + fe2(ind0,:); % And add on the boundary2
end
for i=1:nnodeboun
   ind0 = [ 3*i-2, 3*i-1, 3*i ];
   ind1 = [ 3*nodeboun2glob(i)-2, 3*nodeboun2glob(i)-1, 3*nodeboun2glob(i) ];
   u_known(ind0,:) = ur(ind1,:);
end

%%
% Restrict the data on the given part (not necessary, but cleaner)
f_refer = f_known; u_refer = u_known;
f_known(tofindN) = 0; u_known(tofindD) = 0;

%try
%figure; hold on;
%patch('Faces',boundary2(:,2:4),'Vertices',nodes2,'FaceVertexCData',f_refer(3:3:end,3),'FaceColor','flat');
%colorbar; axis equal;
%end

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

nntt = nboun2;
if nntt > nelemax
   times = idivide(uint64(nntt),uint64(nelemax),'ceil');
   nnt   = nelemax*ones(times,1);
   nnt(end) = nntt - (times-1)*nelemax;
else
   nnt   = nntt;
   times = 1;
end

ncumul = 0; % This int cumulates and makes the link between the vectors
Ru = zeros( nftest, 3*nnodeboun );
Rf = zeros( nftest, 3*nboun2);

for cut=1:times
   nntc = nnt(cut);
   %% Build the list of Gauss points, and of construction-functions
   Xxg = zeros( nntc*Ndots,1 ); Yyg = zeros( nntc*Ndots,1 );
   Zzg = zeros( nntc*Ndots,1 ); Wwg = zeros( nntc*Ndots,1 );
   Phi = sparse( 3*nntc*Ndots, 3*nnodeboun );
   Psi = sparse( 3*nntc*Ndots, 3*nboun2 ); % (And not nntc)
   exnor = zeros( nntc*Ndots,3 ); % Normal
   
   for i=1:nntc
      iHi = i+ncumul; % Recover the right index of the element
      bonod = boundary2(iHi,:); exno = extnorm(iHi,:)';
      
      no1 = bonod(2); no2 = bonod(3); no3 = bonod(4); boname = bonod(1);
      x1 = nodes2(no1,1); y1 = nodes2(no1,2); z1 = nodes2(no1,3);
      x2 = nodes2(no2,1); y2 = nodes2(no2,2); z2 = nodes2(no2,3);
      x3 = nodes2(no3,1); y3 = nodes2(no3,2); z3 = nodes2(no3,3);
         
      nsurf = [ (y1-y2)*(z1-z3)-(y1-y3)*(z1-z2);... %
                (x1-x3)*(z1-z2)-(x1-x2)*(z1-z3);... % Vectorial product
                (x1-x2)*(y1-y3)-(x1-x3)*(y1-y2) ];  %
      surfa = .5*norm(nsurf);
        
      % Dofs in the numerotation of the boundary nodes
      N = boun2loc(iHi,:);
      indDtot = [ 3*N(1)-2, 3*N(1)-1, 3*N(1),...
                  3*N(2)-2, 3*N(2)-1, 3*N(2),...
                  3*N(3)-2, 3*N(3)-1, 3*N(3) ];
                     
      % Dofs in the numerotation of the boundary elements
      indNtot = [ 3*iHi-2, 3*iHi-1, 3*iHi ];
                   
      for j=1:Ndots
         xg = Xg(j,:); wg = Wg(j); umxgx = 1-xg(1)-xg(2);
         xgr = umxgx*[x1;y1;z1] + xg(1)*[x2;y2;z2] + xg(2)*[x3;y3;z3] ; % abscissae
   
         X = xgr(1)/Lx; Y = xgr(2)/Lx; Z = xgr(3)/Lx; % Normalize
         %Swg = surfa * wg;
   
         indexg = Ndots*(i-1) + j; % And not iHi
         Xxg( indexg ) = X; Yyg( indexg ) = Y; Zzg( indexg ) = Z; Wwg( indexg ) = surfa * wg;
         exnor( indexg, :) = exno';
   
         Phi( 3*indexg-2, indDtot)  = [ umxgx, 0, 0, xg(1), 0, 0, xg(2), 0, 0 ];
         Phi( 3*indexg-1, indDtot)  = [ 0, umxgx, 0, 0, xg(1), 0, 0, xg(2), 0 ];
         Phi( 3*indexg  , indDtot)  = [ 0, 0, umxgx, 0, 0, xg(1), 0, 0, xg(2) ];
   
         Psi( 3*indexg-2, indNtot) = [1,0,0];
         Psi( 3*indexg-1, indNtot) = [0,1,0];
         Psi( 3*indexg  , indNtot) = [0,0,1];
      end
   end
   
   %% Build the matrix of test-functions
   Vv = zeros( size(coef,1), 3*nntc*Ndots );
   Sv = zeros( size(coef,1), 3*nntc*Ndots );
   exnoX = exnor(:,1); exnoY = exnor(:,2); exnoZ = exnor(:,3);
   
   for ii=0:nmax
      nd1 = nmax-ii;
      Xi = Xxg.^ii;
      if ii>0, XiM1 = Xxg.^(ii-1); end
      for jj=0:nd1
         nd2 = nmax-ii-jj;
         Yj = Yyg.^jj;
         if jj>0, YjM1 = Yyg.^(jj-1); end
         for kk=0:nd2
            Zk = Zzg.^kk; %ZkM1 = Zzg.^(kk-1);   
            % Build the test field's basis
            if ii>0
               %XYZi = ii*Xxg.^(ii-1).*Yyg.^jj.*Zzg.^kk;
               XYZi = ii*XiM1.*Yj.*Zk;
               sloc11a = (lambda+2*mu)*XYZi;
               sloc22a = lambda*XYZi;
               sloc33a = lambda*XYZi;
               sloc12b = mu*XYZi;
               sloc13c = mu*XYZi;
            else
               sloc11a = zeros(nntc*Ndots,1);
               sloc22a = zeros(nntc*Ndots,1);
               sloc33a = zeros(nntc*Ndots,1);
               sloc12b = zeros(nntc*Ndots,1);
               sloc13c = zeros(nntc*Ndots,1);
            end
   
            if jj>0
               %XYZj = jj*Xxg.^ii.*Yyg.^(jj-1).*Zzg.^kk;
               XYZj = jj*Xi.*YjM1.*Zk;
               sloc12a = mu*XYZj; 
               sloc11b = lambda*XYZj; 
               sloc22b = (lambda+2*mu)*XYZj;
               sloc33b = lambda*XYZj;
               sloc23c = mu*XYZj;
            else
               sloc12a = zeros(nntc*Ndots,1);
               sloc11b = zeros(nntc*Ndots,1);
               sloc22b = zeros(nntc*Ndots,1);
               sloc33b = zeros(nntc*Ndots,1);
               sloc23c = zeros(nntc*Ndots,1);
            end
   
            if kk>0
               %XYZk = kk*Xxg.^ii.*Yyg.^jj.*Zzg.^(kk-1);
               XYZk = kk*Xi.*Yj.*Zzg.^(kk-1);
               sloc13a = mu*XYZk; 
               sloc23b = mu*XYZk;
               sloc11c = lambda*XYZk; 
               sloc22c = lambda*XYZk;
               sloc33c = (lambda+2*mu)*XYZk;
            else
               sloc13a = zeros(nntc*Ndots,1);
               sloc23b = zeros(nntc*Ndots,1);
               sloc11c = zeros(nntc*Ndots,1);
               sloc22c = zeros(nntc*Ndots,1);
               sloc33c = zeros(nntc*Ndots,1);
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
   
            %XYZ = Xxg.^ii.*Yyg.^jj.*Zzg.^kk;
            XYZ = Xi.*Yj.*Zk;
            vpaax = XYZ; vpaby = XYZ; vpacz = XYZ;
   
            index = (nmax+1)^2*ii + (nmax+1)*jj + kk + 1;
                     
            Sv( index, 1:3:3*nntc*Ndots-2 )         = Wwg' .* fpaax';
            Sv( nm3 + index, 1:3:3*nntc*Ndots-2 )   = Wwg' .* fpabx';
            Sv( 2*nm3 + index, 1:3:3*nntc*Ndots-2 ) = Wwg' .* fpacx';
   
            Sv( index, 2:3:3*nntc*Ndots-1 )         = Wwg' .* fpaay';
            Sv( nm3 + index, 2:3:3*nntc*Ndots-1 )   = Wwg' .* fpaby';
            Sv( 2*nm3 + index, 2:3:3*nntc*Ndots-1 ) = Wwg' .* fpacy';
   
            Sv( index, 3:3:3*nntc*Ndots )           = Wwg' .* fpaaz';
            Sv( nm3 + index, 3:3:3*nntc*Ndots )     = Wwg' .* fpabz';
            Sv( 2*nm3 + index, 3:3:3*nntc*Ndots )   = Wwg' .* fpacz';
   
            Vv( index, 1:3:3*nntc*Ndots-2 )         = Wwg' .* vpaax';
            Vv( nm3 + index, 2:3:3*nntc*Ndots-1 )   = Wwg' .* vpaby';
            Vv( 2*nm3 + index, 3:3:3*nntc*Ndots )   = Wwg' .* vpacz';
         end
      end
   end
   
   Ruij = Sv*Phi; clear Sv; clear Phi;
   Rfij = Vv*Psi; clear Vv; clear Psi;
   
   Ru = Ru + coef'*Ruij; Rf = Rf + coef'*Rfij;
   clear Ruij; clear Rfij; % We can't afford not to clear

   ncumul = ncumul + nntc;
end

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
Rhs1 = Rur*u_known(knownD,:) - Rfr*f_known(knownN,:); % It's a multiple Rhs
sRhs = size(Rhs1,1);

% Build the matrix that passes f on dofs on the nodes from the bound
Fbton = sparse( 3*nnodeboun, 3*nboun2 ); Fntob = sparse( 3*nboun2, 3*nnodeboun ); 
% And the boundary mass matrix
nodeMass = sparse(3*nnodeboun,3*nnodeboun);  % For u
elemMass = sparse(3*nboun2,3*nboun2);     % For f
% msur = 10000;  % ~= Infinity
for i=1:nboun2
   coefU = [ boun2loc(i,1), boun2loc(i,2), boun2loc(i,3) ];
   
   no1 = boundary2(i,2); no2 = boundary2(i,3); no3 = boundary2(i,4);
   x1 = nodes2(no1,1); y1 = nodes2(no1,2); z1 = nodes2(no1,3);
   x2 = nodes2(no2,1); y2 = nodes2(no2,2); z2 = nodes2(no2,3);
   x3 = nodes2(no3,1); y3 = nodes2(no3,2); z3 = nodes2(no3,3);
   nsurf = [ (y1-y2)*(z1-z3)-(y1-y3)*(z1-z2);... %
             (x1-x3)*(z1-z2)-(x1-x2)*(z1-z3);... % Vectorial product
             (x1-x2)*(y1-y3)-(x1-x3)*(y1-y2) ];  %
   surfa = .5*norm(nsurf);
   
%    if surfa<msur msur = surfa; end % Just a debug asset
   
   Fntob( 3*i-2, 3*coefU-2 ) = 1/(3*surfa) * [1,1,1];
   Fntob( 3*i-1, 3*coefU-1 ) = 1/(3*surfa) * [1,1,1];
   Fntob( 3*i  , 3*coefU   ) = 1/(3*surfa) * [1,1,1];
   Fbton( 3*coefU-2, 3*i-2 ) = surfa/3 * [1;1;1];
   Fbton( 3*coefU-1, 3*i-1 ) = surfa/3 * [1;1;1];
   Fbton( 3*coefU  , 3*i   ) = surfa/3 * [1;1;1];
   
   ico = [ 3*i-2, 3*i-1, 3*i ];
   cco = [ 3*coefU-2, 3*coefU-1, 3*coefU ];
   elemMass(ico,ico) = surfa*eye(3);
   nodeMass(cco,cco) = nodeMass(cco,cco) + ...
                       surfa/3 * [ eye(3), zeros(3), zeros(3) ; ...
                                   zeros(3), eye(3), zeros(3) ; ...
                                   zeros(3), zeros(3), eye(3) ];
end

% Extract interesting parts
Mfm = elemMass(tofindN,tofindN);% Mfr = elemMass(knownN,knownN);
Mum = nodeMass(tofindD,tofindD);% Mur = nodeMass(knownD,knownD);

% Derivative operators
Du  = sparse(6*nboun2,3*nnodeboun);
%Duf = zeros(6*121,3*nnodeboun);
Df  = elemMass; % No derivative for f
for i=1:nboun2    
   no1 = boundary2(i,2); no2 = boundary2(i,3); no3 = boundary2(i,4);
   x1 = nodes2(no1,1); y1 = nodes2(no1,2); z1 = nodes2(no1,3);
   x2 = nodes2(no2,1); y2 = nodes2(no2,2); z2 = nodes2(no2,3);
   x3 = nodes2(no3,1); y3 = nodes2(no3,2); z3 = nodes2(no3,3);

   n = [ (y1-y2)*(z1-z3)-(y1-y3)*(z1-z2);... %
         (x1-x3)*(z1-z2)-(x1-x2)*(z1-z3);... % Normal to the triangle
         (x1-x2)*(y1-y3)-(x1-x3)*(y1-y2) ];  %
   S = .5*norm(n); n = n/norm(n);
   t = [ x1-x3 ; y1-y3 ; z1-z3 ];   % Tangent vector
   t = t/norm(t);
   v = [ n(2)*t(3)-t(2)*n(3) ; t(1)*n(3)-n(1)*t(3) ; n(1)*t(2)-t(1)*n(2) ];
   v = v/norm(v);
   P = [t,v,n];
   
   % Pass the nodes onto the xy plane
   xyz1 = P'*[x1;y1;z1]; xyz2 = P'*[x2;y2;z2]; xyz3 = P'*[x3;y3;z3];
   x1 = xyz1(1); y1 = xyz1(2); z1 = xyz1(3);
   x2 = xyz2(1); y2 = xyz2(2); z2 = xyz2(3);
   x3 = xyz3(1); y3 = xyz3(2); z3 = xyz3(3);
   xg = (x1+x2+x3)/3; yg = (y1+y2+y3)/3; zg = (z1+z2+z3)/3; % Gauss point
      
   Be = 1/(2*S)*[ x3-x2, 0, 0, x1-x3, 0, 0, x2-x1, 0, 0 ;...
                  y2-y3, 0, 0, y3-y1, 0, 0, y1-y2, 0, 0 ;...
                  0, x3-x2, 0, 0, x1-x3, 0, 0, x2-x1, 0 ;...
                  0, y2-y3, 0, 0, y3-y1, 0, 0, y1-y2, 0 ;...
                  0, 0, x3-x2, 0, 0, x1-x3, 0, 0, x2-x1 ;...
                  0, 0, y2-y3, 0, 0, y3-y1, 0, 0, y1-y2 ];
               
   % Product with Fourier functions
   ii = 0:10; jj = 0:10;
   %F  = real(exp(2*I*pi)*(ii'*xg+jj*yg)); F = F(:);
   %Bf = [ F*Be(1,:) ; F*Be(2,:) ; F*Be(3,:) ; F*Be(4,:) ; F*Be(5,:) ; F*Be(6,:) ];

   indelem = [ 6*i-5, 6*i-4, 6*i-3, 6*i-2, 6*i-1, 6*i ];
%   coefU   = [ 3*boundary2(i,2)-2, 3*boundary2(i,2)-1, 3*boundary2(i,2),...
%               3*boundary2(i,3)-2, 3*boundary2(i,3)-1, 3*boundary2(i,3),...
%               3*boundary2(i,4)-2, 3*boundary2(i,4)-1, 3*boundary2(i,4) ];
   N = boun2loc(i,:);
   coefU   = [ 3*N(1)-2, 3*N(1)-1, 3*N(1),...
               3*N(2)-2, 3*N(2)-1, 3*N(2),...
               3*N(3)-2, 3*N(3)-1, 3*N(3) ];
                  
   Du( indelem, coefU ) = Du( indelem, coefU ) + Be*sqrt(S);
   %Duf( :, coefU )      = Duf( :, coefU )      + Bf*sqrt(S);
end
   
Du0 = Du; Df0 = Df*norm(Du,'fro')/norm(Df,'fro'); %Duf0 = Duf;
Du  = Du0(:,tofindD); Df  = Df0(tofindN,tofindN); %Duf  = Duf0(:,tofindD);
Duu = Du'*Du; Dfu = Df'*Df; Duu0 = Du0'*Du0;
if max(size(Duu))>0
   sDu = real(Duu^(1/2));
else
   sDu = Duu;
end
if max(size(Dfu))>0
%  sDf = real(full(Dfu)^(1/2));
   sDf = Df; % The matrix is diagonal
else
   sDf = Dfu;
end

%Duk = Du0(:,knownD);  Dfk = Df0(knownN,knownN);   Dufk = Duf0(:,knownD);
clear Df0; 
% End derivative operator

% Recover the reference parameters of the plane
x1 = nodes(10,1); y1 = nodes(10,2); z1 = nodes(10,3); %
x2 = nodes(11,1); y2 = nodes(11,2); z2 = nodes(11,3); % Three points on the crack
x3 = nodes(12,1); y3 = nodes(12,2); z3 = nodes(12,3); %

M   = [x1,y1,z1;x2,y2,z2;x3,y3,z3];
abc = -M\ones(3,1); abcd = [abc;1]; abcd = abcd/norm(abcd);

if norm(init)==0
   init = abcd;
end
init = init/norm(init);

theta   = init; thetap = theta; thetarec = theta;
phi     = zeros(niter,1);
nori    = zeros(niter,13);
norreg  = zeros(niter,13);
norres  = zeros(niter,13);

if step == 0
   stepstep = .1;
else
   stepstep = step;
end

kB = norm(Ru,'fro')/norm(Rf,'fro');
Z12 = sparse( size(Duu,1) , size(Dfu,1) );
Lt = [ Duu ,Z12 ; Z12', Dfu ];
Lhst = [LhsA,kB*LhsB];
Ahaha = Lhst'*Lhst + mur*Lt; % Yep, I'm getting short of variable names
if schur == 1
   Ahm1 = Ahaha\eye(size(Ahaha));
   ndofix = size(Lt,1);
else
   minL = min(eig(Lt));
   %minL = min( min(eig(Duu)), min(eig(Dfu)) );
end

disp([ 'Pre-computations terminated ', num2str(toc) ]);
tic;

%% Newton algorithm
for iter = 1:niter
   % Build the basis orthogonal to theta (Gramm-Schmidt)
   theta1 = [1;0;0;0];
   theta1 = theta1 - (theta'*theta1)*theta; theta1 = theta1/norm(theta1);
   theta2 = [0;1;0;0];
   theta2 = theta2 - (theta'*theta2)*theta;
   theta2 = theta2 - (theta1'*theta2)*theta1; theta2 = theta2/norm(theta2);
   theta3 = [0;0;1;0];
   theta3 = theta3 - (theta'*theta3)*theta;
   theta3 = theta3 - (theta1'*theta3)*theta1;
   theta3 = theta3 - (theta2'*theta3)*theta2; theta3 = theta3/norm(theta3);
   Q = [theta1,theta2,theta3];
   
   pbmax = 3; % 3 directions
   
   if step == 0 && iter > 1 % Adapt the step
      stepstep = norm(dtheta)/10;
   end
   
   for pb = 0:pbmax

      if pb == 0
         thetac = theta;
      elseif pb == 1
         thetac = theta + stepstep*theta1;
      elseif pb == 2
         thetac = theta + stepstep*theta2;
      else
         thetac = theta + stepstep*theta3;
      end
   
      %% Build the intersection of the plane with the domain
      mm{1} = [ thetac(1),thetac(2),thetac(3) ; 0,1,0 ; 0,0,1 ]; rr{1} = [ -thetac(4) ; 0 ; 0 ];
      mm{2} = [ thetac(1),thetac(2),thetac(3) ; 1,0,0 ; 0,0,1 ]; rr{2} = [ -thetac(4) ; 1 ; 0 ];
      mm{3} = [ thetac(1),thetac(2),thetac(3) ; 0,1,0 ; 0,0,1 ]; rr{3} = [ -thetac(4) ; 1 ; 0 ];
      mm{4} = [ thetac(1),thetac(2),thetac(3) ; 1,0,0 ; 0,0,1 ]; rr{4} = [ -thetac(4) ; 0 ; 0 ];
   
      mm{5} = [ thetac(1),thetac(2),thetac(3) ; 0,1,0 ; 0,0,1 ]; rr{5} = [ -thetac(4) ; 0 ; 1 ];
      mm{6} = [ thetac(1),thetac(2),thetac(3) ; 1,0,0 ; 0,0,1 ]; rr{6} = [ -thetac(4) ; 1 ; 1 ];
      mm{7} = [ thetac(1),thetac(2),thetac(3) ; 0,1,0 ; 0,0,1 ]; rr{7} = [ -thetac(4) ; 1 ; 1 ];
      mm{8} = [ thetac(1),thetac(2),thetac(3) ; 1,0,0 ; 0,0,1 ]; rr{8} = [ -thetac(4) ; 0 ; 1 ];
   
      mm{9}  = [ thetac(1),thetac(2),thetac(3) ; 1,0,0 ; 0,1,0 ]; rr{9}  = [ -thetac(4) ; 0 ; 0 ];
      mm{10} = [ thetac(1),thetac(2),thetac(3) ; 1,0,0 ; 0,1,0 ]; rr{10} = [ -thetac(4) ; 1 ; 0 ];
      mm{11} = [ thetac(1),thetac(2),thetac(3) ; 1,0,0 ; 0,1,0 ]; rr{11} = [ -thetac(4) ; 1 ; 1 ];
      mm{12} = [ thetac(1),thetac(2),thetac(3) ; 1,0,0 ; 0,1,0 ]; rr{12} = [ -thetac(4) ; 0 ; 1 ];
   
      dots = zeros(3,0);
      for i=1:12
         zemat = mm{i};
         if rank(zemat) == 3
            zerhs = rr{i};
            zedot = zemat\zerhs;
            if zedot(1) >= 0 && zedot(1) <= 1 && zedot(2) >= 0 && zedot(2) <= 1 && zedot(3) >= 0 && zedot(3) <= 1
               dots = [dots,zedot];
            end
         end
      end
   
      if min(size(dots)) == 0
         error('identified plane has no intersection with the domain');
      end
   
      %% Find the convex hull (we're sure the intersection is convex)
      vpre = [thetac(1);thetac(2);thetac(3)];
   
      % Test all the vectorial products
      stilltodo = 1;
      curr = dots(:,1);
      doso = curr;
      while size(doso,2) < size(dots,2);
         x1 = curr(1); y1 = curr(2); z1 = curr(3);
         for i=1:size(dots,2) % Find the next one
            x2 = dots(1,i); y2 = dots(2,i); z2 = dots(3,i);
            if x1==x2 && y1==y2 && z1==z2
               continue;
            end
            negative = 0;
            for j=1:size(dots,2)
               x3 = dots(1,j); y3 = dots(2,j); z3 = dots(3,j);
               if (x3==x2 && y3==y2 && z3==z2) || (x3==x1 && y3==y1 && z3==z1)
                  continue;
               end
               vpro = [ (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1) ;...
                        (x3-x1)*(z2-z1) - (x2-x1)*(z3-z1) ;...
                        (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1) ];
               if vpro'*vpre < 0 % There is at least 1 negative vectorial product
                  negative = 1;
                  break;
               end
            end
            if negative == 0 % It's the good one !
               doso = [doso,dots(:,i)];
               curr = dots(:,i);
               break;
            end
         end
      end

     % TODO : if there are 2 dots that are too close : merge them (lc1/2)

      if pb == 0 %% Generate the file for GMSH
         nnso = size(doso,2);
         fmid = fopen(['meshes/rg3d_crack/plane.geo'],'w');
         fprintf(fmid,'%s%i%s\n','lc1 = ',lc1,';');
         for i=1:nnso
            fprintf(fmid,'%s%d%s%f%s%f%s%f%s\n','Point(',i,') = {',doso(1,i),',',doso(2,i),',',doso(3,i),',lc1};');
         end
         for i=1:nnso-1
            fprintf(fmid,'%s%d%s%d%s%d%s\n','Line(',i,') = {',i,',',i+1,'};');
         end
         fprintf(fmid,'%s%d%s%d%s%d%s\n','Line(',nnso,') = {',nnso,',',1,'};');
      
         fprintf(fmid,'%s','Line Loop(11) = {');
         for i=1:nnso-1
            fprintf(fmid,'%d%s',i,',');
         end
         fprintf(fmid,'%d%s\n',nnso,'};');
         fprintf(fmid,'%s\n','Plane Surface(1) = {11};');
         fprintf(fmid,'%s\n','Physical Line(111) = {11};');
         fprintf(fmid,'%s','Physical Surface(11) = {1};');
         fclose(fmid);
   
         % Use GMSH to mesh the surface
         [stat,out] = system('gmsh -2 "meshes/rg3d_crack/plane.geo" -o "meshes/rg3d_crack/plane.msh"');
         [ nodes3,elements3,ntoelem3,boundary3,order3,lin3 ] = readmesh3D( 'meshes/rg3d_crack/plane.msh' );
         boundnodes = lin3(:,2);
      else %% Morph the mesh
         %figure; patch('Faces',boundary3(:,2:4),'Vertices',nodes3,'FaceAlpha',.25);
         nodes3 = morph2D( nodes3, boundary3(:,2:4), boundnodes, doso' );
         %patch('Faces',boundary3(:,2:4),'Vertices',nodes3,'FaceAlpha',0); axis('equal');
      end

      % Write the normal to the elements
      nboun3   = size(boundary3,1); nnodes3 = size(nodes3,1);
      extnorm3 = ones(nboun3,1)*[thetac(1),thetac(2),thetac(3)];
      
      %%%% Build the operator
      %% Build the list of Gauss points, and of construction-functions
      Xxg = zeros( nboun3*Ndots,1 ); Yyg = zeros( nboun3*Ndots,1 );
      Zzg = zeros( nboun3*Ndots,1 ); Wwg = zeros( nboun3*Ndots,1 );
      Phi = sparse( 3*nboun3*Ndots, 3*nnodes3 );
      exnor = zeros( nboun3*Ndots,3); % Normal
     
      for i=1:nboun3
         bonod = boundary3(i,:); exno = extnorm3(i,:)';
        
         no1 = bonod(2); no2 = bonod(3); no3 = bonod(4); boname = bonod(1);
         x1 = nodes3(no1,1); y1 = nodes3(no1,2); z1 = nodes3(no1,3);
         x2 = nodes3(no2,1); y2 = nodes3(no2,2); z2 = nodes3(no2,3);
         x3 = nodes3(no3,1); y3 = nodes3(no3,2); z3 = nodes3(no3,3);
           
         nsurf = [ (y1-y2)*(z1-z3)-(y1-y3)*(z1-z2);... %
                   (x1-x3)*(z1-z2)-(x1-x2)*(z1-z3);... % Vectorial product
                   (x1-x2)*(y1-y3)-(x1-x3)*(y1-y2) ];  %
         surfa = .5*norm(nsurf);
          
         % Dofs in the numerotation of the boundary nodes
         N = boundary3(i,2:4);
         indDtot = [ 3*N(1)-2, 3*N(1)-1, 3*N(1),...
                     3*N(2)-2, 3*N(2)-1, 3*N(2),...
                     3*N(3)-2, 3*N(3)-1, 3*N(3) ];
                       
         % Dofs in the numerotation of the boundary elements
         indNtot = [ 3*i-2, 3*i-1, 3*i ];
                     
         for j=1:Ndots
            xg = Xg(j,:); wg = Wg(j); umxgx = 1-xg(1)-xg(2);
            xgr  = umxgx*[x1;y1;z1] + xg(1)*[x2;y2;z2] + xg(2)*[x3;y3;z3] ; % abscissae
     
            X = xgr(1)/Lx; Y = xgr(2)/Lx; Z = xgr(3)/Lx; % Normalize
            %Swg   = surfa * wg;
     
            indexg = Ndots*(i-1) + j;
            Xxg( indexg ) = X; Yyg( indexg ) = Y; Zzg( indexg ) = Z; Wwg( indexg ) = surfa * wg;
            exnor( indexg, :) = exno';
     
            Phi( 3*indexg-2, indDtot)  = [ umxgx, 0, 0, xg(1), 0, 0, xg(2), 0, 0 ];
            Phi( 3*indexg-1, indDtot)  = [ 0, umxgx, 0, 0, xg(1), 0, 0, xg(2), 0 ];
            Phi( 3*indexg  , indDtot)  = [ 0, 0, umxgx, 0, 0, xg(1), 0, 0, xg(2) ];
         end
      end
     
      %% Build the matrix of test-functions
      Sv = zeros( size(coef,1), 3*nboun3*Ndots );
      exnoX = exnor(:,1); exnoY = exnor(:,2); exnoZ = exnor(:,3);

      for ii=0:nmax
         Xi = Xxg.^ii;
         if ii>0, XiM1 = Xxg.^(ii-1); end
         nd1 = nmax-ii;
         for jj=0:nd1
            nd2 = nmax-ii-jj;
            Yj = Yyg.^jj;
            if jj>0, YjM1 = Yyg.^(jj-1); end
            for kk=0:nd2
               Zk = Zzg.^kk;
               % Build the test field's basis
               if ii>0
                  XYZi = ii*XiM1.*Yj.*Zk;
                  sloc11a = (lambda+2*mu)*XYZi;
                  sloc22a = lambda*XYZi;
                  sloc33a = lambda*XYZi;
                  sloc12b = mu*XYZi;
                  sloc13c = mu*XYZi;
               else
                  sloc11a = zeros(nboun3*Ndots,1);
                  sloc22a = zeros(nboun3*Ndots,1);
                  sloc33a = zeros(nboun3*Ndots,1);
                  sloc12b = zeros(nboun3*Ndots,1);
                  sloc13c = zeros(nboun3*Ndots,1);
               end
     
               if jj>0
                  XYZj = jj*Xi.*YjM1.*Zk;
                  sloc12a = mu*XYZj; 
                  sloc11b = lambda*XYZj; 
                  sloc22b = (lambda+2*mu)*XYZj;
                  sloc33b = lambda*XYZj;
                  sloc23c = mu*XYZj;
               else
                  sloc12a = zeros(nboun3*Ndots,1);
                  sloc11b = zeros(nboun3*Ndots,1);
                  sloc22b = zeros(nboun3*Ndots,1);
                  sloc33b = zeros(nboun3*Ndots,1);
                  sloc23c = zeros(nboun3*Ndots,1);
               end
     
               if kk>0
                  XYZk = kk*Xi.*Yj.*Zzg.^(kk-1);
                  sloc13a = mu*XYZk; 
                  sloc23b = mu*XYZk;
                  sloc11c = lambda*XYZk; 
                  sloc22c = lambda*XYZk;
                  sloc33c = (lambda+2*mu)*XYZk;
               else
                  sloc13a = zeros(nboun3*Ndots,1);
                  sloc23b = zeros(nboun3*Ndots,1);
                  sloc11c = zeros(nboun3*Ndots,1);
                  sloc22c = zeros(nboun3*Ndots,1);
                  sloc33c = zeros(nboun3*Ndots,1);
               end
     
               fpaax = sloc11a.*exnoX + sloc12a.*exnoY + sloc13a.*exnoZ;
               fpaay = sloc12a.*exnoX + sloc22a.*exnoY;
               fpaaz = sloc13a.*exnoX + sloc33a.*exnoZ; 
     
               fpabx = sloc11b.*exnoX + sloc12b.*exnoY;
               fpaby = sloc12b.*exnoX + sloc22b.*exnoY + sloc23b.*exnoZ;
               fpabz = sloc23b.*exnoY + sloc33b.*exnoZ;
     
               fpacx = sloc11c.*exnoX + sloc13c.*exnoZ;
               fpacy = sloc22c.*exnoY + sloc23c.*exnoZ;
               fpacz = sloc13c.*exnoX + sloc23c.*exnoY + sloc33c.*exnoZ;
     
               XYZ = Xi.*Yj.*Zk;
               vpaax = XYZ; vpaby = XYZ; vpacz = XYZ;
     
               index = (nmax+1)^2*ii + (nmax+1)*jj + kk + 1;
                       
               Sv( index, 1:3:3*nboun3*Ndots-2 )         = Wwg' .* fpaax';
               Sv( nm3 + index, 1:3:3*nboun3*Ndots-2 )   = Wwg' .* fpabx';
               Sv( 2*nm3 + index, 1:3:3*nboun3*Ndots-2 ) = Wwg' .* fpacx';
     
               Sv( index, 2:3:3*nboun3*Ndots-1 )         = Wwg' .* fpaay';
               Sv( nm3 + index, 2:3:3*nboun3*Ndots-1 )   = Wwg' .* fpaby';
               Sv( 2*nm3 + index, 2:3:3*nboun3*Ndots-1 ) = Wwg' .* fpacy';
     
               Sv( index, 3:3:3*nboun3*Ndots )           = Wwg' .* fpaaz';
               Sv( nm3 + index, 3:3:3*nboun3*Ndots )     = Wwg' .* fpabz';
               Sv( 2*nm3 + index, 3:3:3*nboun3*Ndots )   = Wwg' .* fpacz';
            end
         end
      end
     
      Rucij = Sv*Phi; clear Sv; clear Phi;
      Ruc = coef'*Rucij; clear Rucij; % don't clear coef as we reuse it
      Ruc(toremove,:) = [];    % Clear zeros lines
      
      %%%% End of the building of the operator

      if regular==0
         %% Surface mass matrix and pass F from elements to nodes
         Fbton3 = zeros( 3*nnodes3, 3*nboun3 ); Fntob3 = zeros( 3*nboun3, 3*nnodes3 ); 
         nodeMass3 = zeros(3*nnodes3);  % For u
         for i=1:nboun3
            coefU = boundary3(i,2:4); % List of the nodes linked to this element
            
            no1 = boundary3(i,2); no2 = boundary3(i,3); no3 = boundary3(i,4);
            x1 = nodes3(no1,1); y1 = nodes3(no1,2); z1 = nodes3(no1,3);
            x2 = nodes3(no2,1); y2 = nodes3(no2,2); z2 = nodes3(no2,3);
            x3 = nodes3(no3,1); y3 = nodes3(no3,2); z3 = nodes3(no3,3);
            nsurf = [ (y1-y2)*(z1-z3)-(y1-y3)*(z1-z2);... %
                      (x1-x3)*(z1-z2)-(x1-x2)*(z1-z3);... % Vectorial product
                      (x1-x2)*(y1-y3)-(x1-x3)*(y1-y2) ];  %
            surfa = .5*norm(nsurf);
         
            Fntob3( 3*i-2, 3*coefU-2 ) = 1/(3*surfa) * [1,1,1];
            Fntob3( 3*i-1, 3*coefU-1 ) = 1/(3*surfa) * [1,1,1];
            Fntob3( 3*i  , 3*coefU   ) = 1/(3*surfa) * [1,1,1];
            Fbton3( 3*coefU-2, 3*i-2 ) = surfa/3 * [1,1,1];
            Fbton3( 3*coefU-1, 3*i-1 ) = surfa/3 * [1,1,1];
            Fbton3( 3*coefU  , 3*i   ) = surfa/3 * [1,1,1];
         
            cco = [ 3*coefU-2, 3*coefU-1, 3*coefU ];
            nodeMass3(cco,cco) = nodeMass3(cco,cco) + ...
                                 surfa/3 * [ eye(3), zeros(3), zeros(3) ; ...
                                             zeros(3), eye(3), zeros(3) ; ...
                                             zeros(3), zeros(3), eye(3) ];
         end
      end
      
      %% Solve the linear system and recover the unknowns
%      if froreg == 1 && min(size(LhsB))>0
%         kB = sqrt(norm(LhsA,'fro')^2+norm(Ruc,'fro')^2)/norm(LhsB,'fro'); % In order to regularize the stuff
%      else
%         kB = 1;
%      end
      
      offset = size(LhsA,2)+size(LhsB,2);
      tormv = [];
      if zerobound == 1 % Remove the boundaries of the surface
         for i=1:nnodes3
            if any(nodes3(i,:)==0) || any(nodes3(i,:)==1)
               tormv = [tormv,3*i-2,3*i-1,3*i];
            end
         end
      end
      tokeep  = setdiff( 1:offset+3*nnodes3 , offset+tormv );
      tokeep3 = setdiff( 1:3*nnodes3 , tormv );
      
%      Lhs = [LhsA,kB*LhsB,Ruc]; sA = size(Lhs,2);
      Lhs = [Lhst,Ruc]; sA = size(Lhs,2);
      
      %% And the derivative operator
      D3  = sparse(6*nboun3,3*nnodes3);
      %D3f = zeros(6*121,3*nnodes3);
      for i=1:nboun3
         no1 = boundary3(i,2); no2 = boundary3(i,3); no3 = boundary3(i,4);
         x1 = nodes3(no1,1); y1 = nodes3(no1,2); z1 = nodes3(no1,3);
         x2 = nodes3(no2,1); y2 = nodes3(no2,2); z2 = nodes3(no2,3);
         x3 = nodes3(no3,1); y3 = nodes3(no3,2); z3 = nodes3(no3,3);
      
         n = [ (y1-y2)*(z1-z3)-(y1-y3)*(z1-z2);... %
               (x1-x3)*(z1-z2)-(x1-x2)*(z1-z3);... % Normal to the triangle
               (x1-x2)*(y1-y3)-(x1-x3)*(y1-y2) ];  %
         S = .5*norm(n); n = n/norm(n);
         t = [ x1-x3 ; y1-y3 ; z1-z3 ];   % Tangent vector
         t = t/norm(t);
         v = [ n(2)*t(3)-t(2)*n(3) ; t(1)*n(3)-n(1)*t(3) ; n(1)*t(2)-t(1)*n(2) ];
         v = v/norm(v);
         P = [t,v,n];
         
         % Pass the nodes onto the xy plane
         xyz1 = P'*[x1;y1;z1]; xyz2 = P'*[x2;y2;z2]; xyz3 = P'*[x3;y3;z3];
         x1 = xyz1(1); y1 = xyz1(2); z1 = xyz1(3);
         x2 = xyz2(1); y2 = xyz2(2); z2 = xyz2(3);
         x3 = xyz3(1); y3 = xyz3(2); z3 = xyz3(3);
         xg = (x1+x2+x3)/3; yg = (y1+y2+y3)/3; zg = (z1+z2+z3)/3; % Gauss point
            
         Be = 1/(2*S)*[ x3-x2, 0, 0, x1-x3, 0, 0, x2-x1, 0, 0 ;...
                        y2-y3, 0, 0, y3-y1, 0, 0, y1-y2, 0, 0 ;...
                        0, x3-x2, 0, 0, x1-x3, 0, 0, x2-x1, 0 ;...
                        0, y2-y3, 0, 0, y3-y1, 0, 0, y1-y2, 0 ;...
                        0, 0, x3-x2, 0, 0, x1-x3, 0, 0, x2-x1 ;...
                        0, 0, y2-y3, 0, 0, y3-y1, 0, 0, y1-y2 ];
      
         % Product with Fourier functions
         ii = 0:10; jj = 0:10;
         %F  = real(exp(2*I*pi)*(ii'*xg+jj*yg)); F = F(:);
         %Bf = [ F*Be(1,:) ; F*Be(2,:) ; F*Be(3,:) ; F*Be(4,:) ; F*Be(5,:) ; F*Be(6,:) ];
      
         indelem = [ 6*i-5, 6*i-4, 6*i-3, 6*i-2, 6*i-1, 6*i ];
         coefU   = [ 3*boundary3(i,2)-2, 3*boundary3(i,2)-1, 3*boundary3(i,2),...
                     3*boundary3(i,3)-2, 3*boundary3(i,3)-1, 3*boundary3(i,3),...
                     3*boundary3(i,4)-2, 3*boundary3(i,4)-1, 3*boundary3(i,4) ];
                        
         D3( indelem, coefU ) = D3( indelem, coefU ) + Be*sqrt(S);
         %D3f( :, coefU )      = D3f( :, coefU )      + Bf*sqrt(S);
      end
      
      if regular == 1
         D3u = D3'*D3; sD3 = real(D3u^(1/2));
         Z12 = sparse( size(Duu,1) , size(Dfu,1) ); Z13 = sparse( size(Duu,1), size(D3u,1) );
         Z23 = sparse( size(Dfu,1), size(D3u,1) );
         Dtot = [ Duu ,Z12, Z13 ; Z12', Dfu, Z23 ; Z13', Z23', D3u ];
         L = Dtot;
         sL = [ sDu ,Z12, Z13 ; Z12', sDf, Z23 ; Z13', Z23', sD3 ];
         %Zuf = zeros( size(Du,1), size(Df,2)); Zfu = zeros( size(Df,1), size(Du,2));
         %Zu3 = zeros( size(Du,1), size(D3,2)); Z3u = zeros( size(D3,1), size(Du,2));
         %Zf3 = zeros( size(Df,1), size(D3,2)); Z3f = zeros( size(D3,1), size(Df,2)); 
         
         L12  = [ Duu0(tofindD,knownD) ; sparse(size(Dfu,2),size(knownD,2)) ; sparse(3*nnodes3,size(knownD,2)) ];
         L12i = L12*u_known(knownD,:);
         L2i  = u_known(knownD,:)'*Duu0(knownD,knownD)*u_known(knownD,:); % Only its diagonal has meaning
      else
         Z12    = sparse(size(Mum,1),size(Mfm,2)); Z13 = sparse(size(Mum,1),size(nodeMass3,2));
         Z23    = sparse(size(Mfm,1),size(nodeMass3,2));
         Mtot   = [ Mum, Z12, Z13 ; Z12', Mfm, Z23 ; Z13', Z23', nodeMass3 ]; % Weighted mass
         L = Mtot; sL = Mtot^(1/2);
         L12i = zeros(ndofs,13); L2i = zeros(13);
      end

      Solu1 = zeros(sA,13);
      AA12 = Lhst'*Ruc; AAA = [ Ahaha, AA12 ; AA12', Ruc'*Ruc ];
      %AAA = Lhs'*Lhs;
      BBB = Lhs'*Rhs1;% CCC = Rhs1'*Rhs1;
      MAT = AAA + mur*L;
      VEC = BBB - mur*L12i;

      %% Impose the inequation U3.n > 0
      C3 = zeros(nnodes3, 3*nnodes3);
%      C(1,2*1-1) = extnorm3(1,1); C(1,2*1) = extnorm3(1,2);
%      C(nnodes3,2*nnodes3-1) = extnorm3(end,1); C(nnodes3,2*nnodes3) = extnorm3(end,2);
      for i=1:nboun3 %ux*nx + uy*ny + uz*nz > 0
         no1 = boundary3(i,2); no2 = boundary3(i,3); no3 = boundary3(i,4);
         C3(no1,3*no1-2) = C3(no1,3*no1-2) + extnorm3(i,1);%ux*nx for i=1,2,3
         C3(no2,3*no2-2) = C3(no2,3*no2-2) + extnorm3(i,1);
         C3(no3,3*no3-2) = C3(no3,3*no3-2) + extnorm3(i,1);
         C3(no1,3*no1-1) = C3(no1,3*no1-1) + extnorm3(i,2);%uy*ny for i=1,2,3
         C3(no2,3*no2-1) = C3(no2,3*no2-1) + extnorm3(i,2);
         C3(no3,3*no3-1) = C3(no3,3*no3-1) + extnorm3(i,2);
         C3(no1,3*no1)   = C3(no1,3*no1)   + extnorm3(i,3);%uz*nz for i=1,2,3
         C3(no2,3*no2)   = C3(no2,3*no2)   + extnorm3(i,3);
         C3(no3,3*no3)   = C3(no3,3*no3)   + extnorm3(i,3);
      end
      
      % Normalize
      for i=1:nnodes3
         C3(i,:) = C3(i,:)/norm(C3(i,:));
      end
      
      % Add the zero rows for the other dofs (and sparsify)
      C = [ sparse(nnodes3,sA-3*nnodes3), C3 ];
      respos = zeros(nuzawa,1); df = zeros(nuzawa,1);

      if schur==0
         if zerobound == 1
            MAT1 = MAT(tokeep,tokeep)\eye(size(MAT(tokeep,tokeep)));
            %[Ll, Uu, Pp] = lu (MAT(tokeep,tokeep));
            %Ll = chol (MAT(tokeep,tokeep),'lower'); % Not a real improvement...
            %kuzawa1 = .999*min(eig(MAT(tokeep,tokeep)));
            kuzawa1 = mur*min( minL, min(eig(D3u(tokeep3,tokeep3))) ); % Quicker sharp lower bound
         else
            MAT1 = MAT\eye(size(MAT));
            %[Ll, Uu, Pp] = lu (MAT);  % P * M = L * U
            kuzawa1 = .999*min(eig(MAT));
         end
      
         % Acually invert the matrix (is is in fact more effecient (at least for the semi-fat case I tried))
         %Uu1 = Uu\eye(size(Uu)); Ll1 = Ll\eye(size(Ll)); MAT1 = Uu1*Ll1*Pp;

         % Just an initialization
         Solu1 = zeros(size(L,1),13);
         f = zeros(nnodes3,13); Ctf = C'*f;
      
         for i=1:nuzawa % Uzawa for the contact problem
            if zerobound == 1
               %Solu1(tokeep,:) = Uu \ ( Ll \ ( Pp * ( VEC(tokeep,:) + Ctf(tokeep,:) ) ) );
               Solu1(tokeep,:) = MAT1 * ( VEC(tokeep,:) + Ctf(tokeep,:) );
               %Solu1(tokeep,:) = (Ll') \ ( Ll \ ( VEC(tokeep,:) + Ctf(tokeep,:) ) );
            else
               %Solu1 = Uu \ ( Ll \ ( Pp * ( VEC + Ctf ) ) );
               Solu1 = MAT1 * ( VEC + Ctf );
            end
            CSolu1 = C*Solu1;
            respos(i) = norm(CSolu1 - abs(CSolu1),'fro');
            fp = f;
         
            f = f - kuzawa1*CSolu1;
            f = .5*(f + abs(f)); Ctf = C'*f;
            df(i) = norm(f-fp,'fro');
         end
         %% End Uzawa

      else
         if zerobound == 1
            Cth = MAT(1:ndofix,ndofix+1:end); Bth = MAT(ndofix+1:end,ndofix+1:end);
            Cth = Cth(:,tokeep3); Bth = Bth(tokeep3,tokeep3);
            b1  = VEC(1:ndofix,:); b3 = VEC(ndofix+1:end,:); b3 = b3(tokeep3,:);
            %kuzawa1 = mur*min( minL, min(eig(D3u(tokeep3,tokeep3))) );
         else
            Cth = MAT(1:ndofix,ndofix+1:end); Bth = MAT(ndofix+1:end,ndofix+1:end);
            b1  = VEC(1:ndofix,:); b3 = VEC(ndofix+1:end,:);
            %Cth = Lhst'*Ruc; Bth = Ruc'*Ruc;
            %kuzawa1 = .999*min(eig(MAT));
         end
         ACth = Ahm1*Cth;
         %ACth = Ahahah\Cth;
         Sth  = (Bth-Cth'*ACth); % Schur complement
         [Ll, Uu, Pp] = lu (Sth);
         %Stm1 = Sth\eye(size(Sth));
         %StAt =  Stm1*ACth';
         %StAt =  Sth\ACth';
         %M11  = Ahm1 + ACth*StAt;
         %M12  = -StAt'; M22 = Stm1;
         %MAT1 = [ M11, M12 ; M12', M22 ];
         SVEC = b3 - ACth'*b1; % Condensed RHS

         kuzawa1 = min(eig(Sth)); % Sharper estimator for the parameter

         % Just an initialization
         x3 = zeros(size(D3u,1),13);
         f = zeros(nnodes3,13); Ctf = C3'*f;

         for i=1:nuzawa % Uzawa only on the D3 part
            if zerobound == 1
               %x3(tokeep3,:) = Sth \ ( SVEC + Ctf(tokeep3,:) );
               x3(tokeep3,:) = Uu \ ( Ll \ ( Pp * ( SVEC + Ctf(tokeep3,:) ) ) );
            else
               %x3 = Sth \ ( SVEC + Ctf );
               x3 = Uu \ ( Ll \ ( Pp * ( SVEC + Ctf ) ) );
            end
            CSolu1 = C3*x3;
            respos(i) = norm(CSolu1 - abs(CSolu1),'fro');
            fp = f;
         
            f = f - kuzawa1*CSolu1;
            f = .5*(f + abs(f)); Ctf = C3'*f;
            df(i) = norm(f-fp,'fro');
         end
         %% End Uzawa, compute the complete solution vector
         x1 = Ahm1 * ( b1 - Cth*x3(tokeep3,:) );
         Solu1 = [x1;x3];
      end      

      % Compute the Gradient and the Hessian
      if pb == 0
         Ax  = Lhs*Solu1; sLx = sL*Solu1;
         L1x = diag(Solu1'*L12i);
         res = Ax(:) - Rhs1(:); rel = sLx(:);
%         for i=1:13
%            norres(iter,i) = Solu1(:,i)'*AAA*Solu1(:,i) - 2*Solu1(:,i)'*BBB(:,i) +...
%                             Rhs1(:,i)'*Rhs1(:,i);
%            norreg(iter,i) = mur*Solu1(:,i)'*L*Solu1(:,i) +...
%                             2*mur*Solu1(:,i)'*L12i(:,i) + mur*L2i(i,i);
%            nori(iter,i)   = norres(iter,i) + norreg(iter,i);
%            phi(iter)    = phi(iter) + nori(iter,i);
%         end
         norres(iter,:) = diag( Solu1'*AAA*Solu1 - 2*Solu1'*BBB + Rhs1'*Rhs1 );
         norreg(iter,:) = diag( mur*Solu1'*L*Solu1 + 2*mur*Solu1'*L12i + mur*L2i );
         nori(iter,:)   = norres(iter,:) + norreg(iter,:); phi(iter) = sum(nori(iter,:));
         
         if phi(iter) == min(phi(1:iter)) % Store this one
            thetap = theta;
            nodes3p = nodes3; Solu1p = Solu1;
            extnorm3p = extnorm3; boundary3p = boundary3;
         end
         nodes3_ref = nodes3; boundary3_ref = boundary3; % Store the Nodes (to passMesh)
         norm_ref   = theta(1:3)/norm(theta(1:3)); theta_ref = theta;
         %disp([ 'iteration ', num2str(iter),' pb 0' ]);
      elseif pb == 1
         Ax1  = Lhs*Solu1; L1x1 = diag(Solu1'*L12i);
         sLx1 = sL*Solu1;% topass0 = sLxt(end-3*size(nodes3,1)+1:end,:); sLx1 = sLxt(1:end-3*size(nodes3,1),:);
         %topass = topass0;
         %topass = passMesh2D3d (nodes3, boundary3(:,2:4), nodes3_ref, [], topass0); % Get it on the reference mesh
         %sLx1   = [ sLx1 ; topass ];
         Dd1    = (Ax1(:)-Ax(:))/stepstep;
         DL1    = (sLx1(:)-sLx(:))/stepstep;
         DL121  = (L1x1-L1x)/stepstep;
         %disp([ 'iteration ', num2str(iter),' pb 1' ]);
      elseif pb == 2
         Ax2  = Lhs*Solu1; L1x2 = diag(Solu1'*L12i);
         sLx2 = sL*Solu1;% topass0 = sLxt(end-3*size(nodes3,1)+1:end,:); sLx2 = sLxt(1:end-3*size(nodes3,1),:);
         %topass = topass0;
         %topass = passMesh2D3d (nodes3, boundary3(:,2:4), nodes3_ref, [], topass0); % Get it on the reference mesh
         %sLx2   = [ sLx2 ; topass ];
         Dd2    = (Ax2(:)-Ax(:))/stepstep;
         DL2    = (sLx2(:)-sLx(:))/stepstep;
         DL122  = (L1x2-L1x)/stepstep;
         %disp([ 'iteration ', num2str(iter),' pb 2' ]);
      else % pb == 3
         Ax3  = Lhs*Solu1; L1x3 = diag(Solu1'*L12i);
         sLx3 = sL*Solu1;% topass0 = sLxt(end-3*size(nodes3,1)+1:end,:); sLx3 = sLxt(1:end-3*size(nodes3,1),:);
         %topass = topass0;
         %topass = passMesh2D3d (nodes3, boundary3(:,2:4), nodes3_ref, [], topass0); % Get it on the reference mesh
         %sLx3   = [ sLx3 ; topass ];
         Dd3    = (Ax3(:)-Ax(:))/stepstep;
         DL3    = (sLx3(:)-sLx(:))/stepstep;
         DL123  = (L1x3-L1x)/stepstep;
         %disp([ 'iteration ', num2str(iter),' pb 3' ]);
      end

   end
   D = [Dd1,Dd2,Dd3]; DL = [DL1,DL2,DL3]; DL12 = [DL121,DL122,DL123];
   dtheta = - ( D'*D + mur*DL'*DL ) \ ( D'*res + mur*DL'*rel + mur*DL12'*ones(size(DL12,1),1) );
   
   % Actualize and normalize
   theta = theta + Q*dtheta;
   theta = theta/norm(theta);
   
   thetarec = [thetarec,theta];
   disp([ 'iteration ', num2str(iter)]);
end

disp([ 'Newton algorithm terminated ', num2str(toc) ]);

Solu1 = Solu1p;

% Post-process : reconstruct u and f from Solu
fsolu1 = zeros(3*nboun2,13); usolu1 = zeros(3*nnodeboun,13);
usolu1(tofindD,:) = Solu1( 1:ndofD, : );
fsolu1(tofindN,:) = kB*Solu1( ndofD+1:ndofD+ndofN, : );
fsoluN = Fbton*fsolu1; f_refeN = Fbton*f_refer;

%try
%figure; hold on;
%patch('Faces',boundary2(:,2:4),'Vertices',nodes2,'FaceVertexCData',fsolu1(3:3:end,3),'FaceColor','flat');
%colorbar; axis equal;
%end

%dofglob = [ 3*nodeboun2glob-2 ;  3*nodeboun2glob-1 ; 3*nodeboun2glob ];
%dofloc = [ 1:3:3*nnodeboun-2 ; 1:3:3*nnodeboun-1 ; 1:3:3*nnodeboun ];
%fsog = zeros(3*nnodes2,13); fsog(dofglob,:) = fsoluN(dofloc,:);
%freg = zeros(3*nnodes2,13); freg(dofglob,:) = f_refeN(dofloc,:);
%plotGMSH3D( {fsog(:,11),'sol';freg(:,11),'ref'}, boundary2, nodes2, 'output/f_output.msh' );

Usol = zeros(3*nnodes2,13); Fsol = zeros(3*nnodes2,13);
Uref = zeros(3*nnodes2,13); Fref = zeros(3*nnodes2,13);
for i=1:nnodeboun
   ind1 = [ 3*nodeboun2glob(i)-2 , 3*nodeboun2glob(i)-1 , 3*nodeboun2glob(i) ];
   ind2 = [ 3*i-2, 3*i-1, 3*i ];
   Usol(ind1,:) = usolu1(ind2,:); Fsol(ind1,:) = fsoluN(ind2,:);
   Uref(ind1,:) = u_refer(ind2,:); Fref(ind1,:) = f_refeN(ind2,:);
end

ucrsol1 = Solu1(end-3*size(nodes3p,1)+1:end,:); nnodes3o = size(nodes3p,1);

% Visulaize the gap in the crack's plane
eno = extnorm3p(1,:); % All the normals are the same
UsN = ucrsol1(1:3:end-2,:)*eno(1) + ucrsol1(2:3:end-1,:)*eno(2) + ucrsol1(3:3:end,:)*eno(3);
try
figure;
colormap(jet(64));
set(gca, 'fontsize', 20);
patch('Faces',boundary3p(:,2:4),'Vertices',nodes3p,'FaceVertexCData',UsN(:,11),'FaceColor','interp');
%caxis( [-1.32e-6, 9.75e-6] );
%caxis( [-0.0000018532, 0.000013567] );
axis([0,1,0,1,0,1],'square'); set(colorbar, 'fontsize', 20); colorbar('SouthOutside');
end

%% (Re-)Export it to gmsh
mesh2GMSH( nodes3p, boundary3p(:,2:4), [], 'output/mesh_usn' );
plotGMSH3D( {UsN(:,11),'UsN'}, boundary3p(:,2:4), nodes3p, 'output/field_usn' );

%%%% And recover the reference
%% Build the intersection of the plane with the domain
mm{1} = [ abcd(1),abcd(2),abcd(3) ; 0,1,0 ; 0,0,1 ]; rr{1} = [ -abcd(4) ; 0 ; 0 ];
mm{2} = [ abcd(1),abcd(2),abcd(3) ; 1,0,0 ; 0,0,1 ]; rr{2} = [ -abcd(4) ; 1 ; 0 ];
mm{3} = [ abcd(1),abcd(2),abcd(3) ; 0,1,0 ; 0,0,1 ]; rr{3} = [ -abcd(4) ; 1 ; 0 ];
mm{4} = [ abcd(1),abcd(2),abcd(3) ; 1,0,0 ; 0,0,1 ]; rr{4} = [ -abcd(4) ; 0 ; 0 ];

mm{5} = [ abcd(1),abcd(2),abcd(3) ; 0,1,0 ; 0,0,1 ]; rr{5} = [ -abcd(4) ; 0 ; 1 ];
mm{6} = [ abcd(1),abcd(2),abcd(3) ; 1,0,0 ; 0,0,1 ]; rr{6} = [ -abcd(4) ; 1 ; 1 ];
mm{7} = [ abcd(1),abcd(2),abcd(3) ; 0,1,0 ; 0,0,1 ]; rr{7} = [ -abcd(4) ; 1 ; 1 ];
mm{8} = [ abcd(1),abcd(2),abcd(3) ; 1,0,0 ; 0,0,1 ]; rr{8} = [ -abcd(4) ; 0 ; 1 ];

mm{9}  = [ abcd(1),abcd(2),abcd(3) ; 1,0,0 ; 0,1,0 ]; rr{9}  = [ -abcd(4) ; 0 ; 0 ];
mm{10} = [ abcd(1),abcd(2),abcd(3) ; 1,0,0 ; 0,1,0 ]; rr{10} = [ -abcd(4) ; 1 ; 0 ];
mm{11} = [ abcd(1),abcd(2),abcd(3) ; 1,0,0 ; 0,1,0 ]; rr{11} = [ -abcd(4) ; 1 ; 1 ];
mm{12} = [ abcd(1),abcd(2),abcd(3) ; 1,0,0 ; 0,1,0 ]; rr{12} = [ -abcd(4) ; 0 ; 1 ];

dots = zeros(3,0);
for i=1:12
   zemat = mm{i};
   if rank(zemat) == 3
      zerhs = rr{i};
      zedot = zemat\zerhs;
      if zedot(1) >= 0 && zedot(1) <= 1 && zedot(2) >= 0 && zedot(2) <= 1 && zedot(3) >= 0 && zedot(3) <= 1
         dots = [dots,zedot];
      end
   end
end

%% Find the convex hull (we're sure the intersection is convex)
% Compute a random vectorial product to get a reference
x1 = dots(1,1); y1 = dots(2,1); z1 = dots(3,1);
x2 = dots(1,2); y2 = dots(2,2); z2 = dots(3,2);
x3 = dots(1,3); y3 = dots(2,3); z3 = dots(3,3);
vpre = [ (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1) ;...
         (x3-x1)*(z2-z1) - (x2-x1)*(z3-z1) ;...
         (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1) ];

% Test all the vectorial products
stilltodo = 1;
curr = dots(:,1);
doso = curr;
while size(doso,2) < size(dots,2);
   x1 = curr(1); y1 = curr(2); z1 = curr(3);
   for i=1:size(dots,2) % Find the next one
      x2 = dots(1,i); y2 = dots(2,i); z2 = dots(3,i);
      if x1==x2 && y1==y2 && z1==z2
         continue;
      end
      negative = 0;
      for j=1:size(dots,2)
         x3 = dots(1,j); y3 = dots(2,j); z3 = dots(3,j);
         if (x3==x2 && y3==y2 && z3==z2) || (x3==x1 && y3==y1 && z3==z1)
            continue;
         end
         vpro = [ (y2-y1)*(z3-z1) - (y3-y1)*(z2-z1) ;...
                  (x3-x1)*(z2-z1) - (x2-x1)*(z3-z1) ;...
                  (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1) ];
         if vpro'*vpre < 0 % There is at least 1 negative vectorial product
            negative = 1;
            break;
         end
      end
      if negative == 0 % It's the good one !
         doso = [doso,dots(:,i)];
         curr = dots(:,i);
         break;
      end
   end
end

%% Generate the file for GMSH
nnso = size(doso,2);
fmid = fopen(['meshes/rg3d_crack/plane.geo'],'w');
fprintf(fmid,'%s%i%s\n','lc1 = ',lc1,';');
for i=1:nnso
   fprintf(fmid,'%s%d%s%f%s%f%s%f%s\n','Point(',i,') = {',doso(1,i),',',doso(2,i),',',doso(3,i),',lc1};');
end
for i=1:nnso-1
   fprintf(fmid,'%s%d%s%d%s%d%s\n','Line(',i,') = {',i,',',i+1,'};');
end
fprintf(fmid,'%s%d%s%d%s%d%s\n','Line(',nnso,') = {',nnso,',',1,'};');

fprintf(fmid,'%s','Line Loop(11) = {');
for i=1:nnso-1
   fprintf(fmid,'%d%s',i,',');
end
fprintf(fmid,'%d%s\n',nnso,'};');
fprintf(fmid,'%s\n','Plane Surface(1) = {11};');
fprintf(fmid,'%s','Physical Surface(11) = {1};');
fclose(fmid);

% Use GMSH to mesh the surface
[stat,out] = system('gmsh -2 "meshes/rg3d_crack/plane.geo" -o "meshes/rg3d_crack/plane.msh"');
[ nodes3,elements3,ntoelem3,boundary3,order3 ] = readmesh3D( 'meshes/rg3d_crack/plane.msh' );
nnodes3 = size(nodes3,1);

% First, make a slight contraction to make sure nothing is outside
xm = mean(nodes3(:,1)); ym = mean(nodes3(:,2)); zm = mean(nodes3(:,3));
nodes3(:,1) = (nodes3(:,1)-xm)*.999 + xm;
nodes3(:,2) = (nodes3(:,2)-ym)*.999 + ym;
nodes3(:,3) = (nodes3(:,3)-zm)*.999 + zm;

normal3 = abcd(1:3)/norm(abcd(1:3));
nodes_up = nodes3 + 1e-4*ones(nnodes3,1)*normal3';
nodes_do = nodes3 - 1e-4*ones(nnodes3,1)*normal3';
U_up = passMesh3D( nodes, elements, nodes_up, [], uref );
U_do = passMesh3D( nodes, elements, nodes_do, [], uref );
U_diff = U_up-U_do;

UrNr = U_diff(1:3:end-2,:)*normal3(1) + U_diff(2:3:end-1,:)*normal3(2) + U_diff(3:3:end,:)*normal3(3);
try
figure;
colormap(jet(64));
set(gca, 'fontsize', 20);
patch('Faces',boundary3(:,2:4),'Vertices',nodes3,'FaceVertexCData',UrNr(:,11),'FaceColor','interp');
%caxis( [-1.32e-6, 9.75e-6] );
%caxis( [-0.0000018532, 0.000013567] );
axis([0,1,0,1,0,1],'square'); set(colorbar, 'fontsize', 20); colorbar('SouthOutside');
end

%% (Re-)Export it to gmsh
mesh2GMSH( nodes3, boundary3(:,2:4), [], 'output/mesh_urnr' );
plotGMSH3D( {UrNr(:,11),'UsN'}, boundary3(:,2:4), nodes3, 'output/field_urnr' );

%%% Plot on a line
%% First : computed stuff : use the projection on the xy plane of nodes3
%nnodes3p = size(nodes3,1);
%x0 = 0; xf = 1; dx = 1/20; Xx = x0+dx:dx:xf-dx; Yy = .5*ones(1,19); % remove the 2 extreme points
%ucx = ucrsol1(1:3:3*nnodes3p-2,:); ucy = ucrsol1(2:3:3*nnodes3p-1,:);
%ucz = ucrsol1(3:3:3*nnodes3p,:);   uca = zeros(nnodes3p,13); % passMesh2D accepts 2 dofs fields
%upass = [ reshape([ucx,ucy]',[13,2*nnodes3])' , reshape([ucz,uca]',[13,2*nnodes3])'];
%U_com = passMesh2D( nodes3p(:,1:2), boundary3p(2:4), [Xx',Yy'], [], upass );
%U_com = [ U_com(1:2:2*19-1,1:13), U_com(2:2:2*19,1:13), U_com(1:2:2*19-1,14:end) ];
%U_com = reshape(U_com',[13,3*19])';
%U_com = [zeros(3,13);U_com;zeros(3,13)]; eno = extnorm3p(1,:);
%UsN2 = U_com(1:3:end-2,:)*eno(1) + U_com(2:3:end-1,:)*eno(2) + U_com(3:3:end,:)*eno(3);
%
%% And the reference
%z0 = -(theta(4)+.5*theta(2))/theta(3); % Pray for having theta(3) \neq 0
%zf = -(theta(1)+theta(4)+.5*theta(2))/theta(3); dz = (zf-z0)/20;
%Zz = z0:dz:zf;
%x0 = 0; xf = 1; dx = 1/20; Xx = x0:dx:xf; Yy = .5*ones(1,21);
%nodes_up2 = [Xx',Yy',Zz'] + 1e-4*ones(21,1)*normal3';
%nodes_do2 = [Xx',Yy',Zz'] - 1e-4*ones(21,1)*normal3';
%U_up2 = passMesh3D( nodes, elements, nodes_up2, [], uref );
%U_do2 = passMesh3D( nodes, elements, nodes_do2, [], uref );
%U_diff2 = U_up2-U_do2;
%eno = extnorm3(1,:);
%UsN2r = U_diff2(1:3:end-2,:)*eno(1) + U_diff2(2:3:end-1,:)*eno(2) + U_diff2(3:3:end,:)*eno(3);
%
%try
%figure
%plot([0,Xx,1],UsN2(:,11),'Color','red');
%plot(Xx,UsN2r(:,11),'Color','black');
%legend('computed gap','reference')
%end
%
%% Extract the identified values on the boundary
% TODO : Neumann is tricky
%if min(size(tofindD)) ~= 0 %11
%   usol0 = Solu1(1:ndofD,:); 
%   uref = zeros(3*nnodeboun,13);
%   uref([1:3:end-2,2:3:end-1,3:3:end],:) = ur([3*nodeboun2glob-2,3*nodeboun2glob-1,3*nodeboun2glob],:); % Pass to local
%   usol = uref; usol(tofindD,:) = usol0;% Replace by the identified field
%   us = zeros(3*nnodes2,13);
%   us([3*nodeboun2glob-2,3*nodeboun2glob-1,3*nodeboun2glob],:) = usol([1:3:end-2,2:3:end-1,3:3:end],:); % Re-pass to global
%   try
%   figure;
%   hold on;
%   patch('Faces',boundary2(:,2:4),'Vertices',nodes2,'FaceVertexCData',ur(1:3:end-2,1),'FaceColor','interp');
%   colorbar(); set(colorbar, 'fontsize', 20); axis([0,1,0,1,0,1],'square');
%   legend('reference');
%   end
%   try
%   figure;
%   hold on;
%   patch('Faces',boundary2(:,2:4),'Vertices',nodes2,'FaceVertexCData',us(1:3:end-2,1),'FaceColor','interp');
%   colorbar(); set(colorbar, 'fontsize', 20); axis([0,1,0,1,0,1],'square');
%   legend('identified');
%   end
%end

% Plot phi
try
   figure;
   plot(log10(phi));
   legend('cost function');
end

% Compare both planes
try
figure;
hold on;
set(gca, 'fontsize', 20);
patch('Faces',boundary3p(:,2:4),'Vertices',nodes3p,'FaceVertexCData',ones(size(nodes3p,1),1),'FaceColor','interp');
patch('Faces',boundary3(:,2:4),'Vertices',nodes3,'FaceVertexCData',zeros(size(nodes3,1),1),'FaceColor','interp');
%patch('Faces',boundary3p(:,2:4),'Vertices',nodes3p,'FaceVertexCData',UsN(:,11),'FaceColor','interp');
%patch('Faces',boundary3(:,2:4),'Vertices',nodes3,'FaceVertexCData',UrNr(:,11),'FaceColor','interp');
colorbar(); set(colorbar, 'fontsize', 20); axis([0,1,0,1,0,1],'square');
end

%% Plot the evolution of the parameters
[~,thisone] = min(phi);
try
figure;
hold on;
plot(thetarec(1,:),'Color','red');
plot(thetarec(2,:),'Color','blue');
plot(thetarec(3,:),'Color','green');
plot(thetarec(4,:),'Color','magenta');
plot(abcd(1)*ones(niter,1),'Color','red');
plot(abcd(2)*ones(niter,1),'Color','blue');
plot(abcd(3)*ones(niter,1),'Color','green');
plot(abcd(4)*ones(niter,1),'Color','magenta');
plot( [thisone,thisone], [-1,1], 'Color', 'black' );
legend('a', 'b','c', 'd');
end
