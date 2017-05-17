% 01/12/2016 : Algo de Schur Dual avec plusieurs sous-domaines

clear all;
close all;

addpath(genpath('./tools'))

% Parameters
E       = 200000; % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 250;    % N.mm-2 : Loading on the plate
mat = [0, E, nu];

nsub    = 3;
niter   = 10;
precond = 1;

[ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/plate.msh' );
nnodes = size(nodes,1);

%dirichlet = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ; 4,1,0 ; 4,2,0 ];
%dirichlet = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ];
dirichlet = [ 1,1,0 ; 1,2,0 ];
%dirichlet = [ 3,1,0 ; 3,2,0 ];
%dirichlet = [ 4,1,0 ; 4,2,0 ];
neumann   = [ 3,2,fscalar ];

%figure
%patch('Faces',elements(:,1:3),'Vertices',nodes,'FaceAlpha',0);
%axis('equal');

%% Submeshes

% Find the dimensions of the plate
L1 = min(nodes(:,1)); L2 = max(nodes(:,1));
H1 = min(nodes(:,2)); H2 = max(nodes(:,2));

newelem  = {};
newnode  = {};
newbouns = {}; % global numerotation of the internal boundaries
bounsloc = {}; % local numerotation of the internal boundaries
map      = {};
boundar  = {}; % restriction of the external boundaries to each subdomain

newelem{1} = elements ;
newnode{1} = nodes ;
boundar{1} = boundary ;
map{1}     = 1:nnodes ;

%% First, split the domain
for i = 1:nsub-1
   elementsc = newelem{i};
   nodesc    = newnode{i};
   boundaryc = boundar{i};

   % Cut along a line
   line = [ L1 ; i*(H1+H2)/nsub ; L2 ; i*(H1+H2)/nsub];
   [ nodes1, elements1, boundary1, map1,...
     nodes2, elements2, boundary2, map2, newbound ] =...
                                   cutMesh (nodesc, elementsc, boundaryc, line);

   % Get the connectivity between the divided domain and the global domain
   if i > 1
      [ b1to2, b2to1 ] = superNodes( nodesc, nodes, 1e-6 );
   else
      b1to2 = (1:nnodes)';
   end
                                   
   newelem{i}   = elements1;
   newelem{i+1} = elements2;
   newnode{i}   = nodes1;
   newnode{i+1} = nodes2;
   
   % Increment the map
   oldmap = map{i};
   map{i} = map1(oldmap);
   map{i+1} = map2(oldmap);

   % get the local boundary
   newbouns{i} = b1to2(newbound);
   
   % The boundaries
   boundar{i} = boundary1; boundar{i+1} = boundary2;
   
   % store b1to2 (always useful)
   b1to2s{i} = b1to2;
end

figure
hold on;
for i = 1:nsub
   elementsc = newelem{i};
   nodesc    = newnode{i};
   ret  = patch('Faces',elementsc(:,1:3),'Vertices',nodesc);
   col1 = rand(); col2 = rand(); col3 = rand(); coto = col1+col2+col3;
   set (ret, 'FaceColor', [col1,col2,col3]);
end
axis('equal');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reference solution
[K,C,nbloq,node2c,c2node] = Krig2 (nodes,elements,mat,order,boundary,dirichlet);
f = loading(nbloq,nodes,boundary,neumann);
u = K\f; u = u(1:2*nnodes,1);
sigma = stress(u,E,nu,nodes,elements,order,1,ntoelem);
plotGMSH({u,'U_vect';sigma,'stress'}, elements, nodes(:,[1,2]), 'output/reference');
uref = u;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin the serious stuff

% Stiffness, rigid modes & Co
K = {}; f = {}; nbloq = {}; G = {}; R = {};
Kp = {}; nbloqp = {}; Ct = {}; % trace operator on the boundary
Gglob = zeros( 2*nnodes, 3*nsub ); % Gglob's size will decrease

% Connectivity tables between the subdomains and the rigid modes
sdtorm = zeros(nsub,1); rmtosd = [];

j = 1; % index for the connectivity tables
for i = 1:nsub
   nodess    = newnode{i};
   elementss = newelem{i};
   boundarys = boundar{i};
   map1      = map{i};
   %map2      = map{i+1};
   
   sdtorm(i) = 0;
   
   addboun = [];
   if i-1 > 0
      boun1       = map1(newbouns{i-1}); % Local boundaries
      bounsloc{2*i-1} = boun1;
      addboun = [addboun ; boun1];
   end
   if i < nsub
      boun2       = map1(newbouns{i});
      bounsloc{2*i} = boun2;
      addboun = [addboun ; boun2];
   end
   
   nno{i} = size(nodess,1);
   [K{i},C,nbloq{i},node2c,c2node] =...
                 Krig2 (nodess,elementss,mat,order,boundarys,dirichlet);
   f{i} = loading(nbloq{i},nodess,boundarys,neumann);
   R{i} = [];
   
   % Primal problem : 
   [Kpt,C,nbloqp{i},node2c,c2node] =...
                 Krig2 (nodess,elementss,mat,order,boundarys,dirichlet);
   Kpinter = Kpt(1:2*nno{i}, 1:2*nno{i});
   nbloqp{i} = nbloqp{i} + 2*size(addboun,1);
   
   % See witch dof has already Dirichlet conditions
   testC = C*ones(size(C,2),1);
   
   % Lagrange multipliers for addboun
   %C1c = zeros( 2*nno{i}, 2*size(addboun) ); % intermediate matrix
   C1c = zeros( 2*nno{i}, 0 ); % intermediate matrix
   index = 0;
   for no=1:size(addboun,1)
      if testC(2*addboun(no)-1) == 0
         index = index+1;
         C1c(2*addboun(no)-1, index) = 1;
      end
      if testC(2*addboun(no)) == 0 
         index = index+1;
         C1c(2*addboun(no), index) = 1;
      end
   end
   Ct{i} = [ C, C1c ];
   Kp{i} = [Kpinter, Ct{i} ; Ct{i}', zeros(size(Ct{i},2))];
   
   % Management of the rigid modes (TOFIX \mu-pp)
   if nbloq{i} == 0 % nbloq < 3 is also problematic, but ...
   
      % Connectivity stuff (subdomains VS rigid modes)
      sdtorm(i) = j;
      rmtosd(j) = i;
   
      nnodess = nno{i};
%      r1 = zeros(2*nnodess,1); r2 = r1; r3 = r1; g1 = r1; g2 = r1; g3 = r1;
%      ind = 2:2:2*nnodess;
%      r1(ind-1,1) = 1; r2(ind,1) = 1;
%      
%      moyx = (max( nodess(:,1) ) - min( nodess(:,1) ))/2;
%      moyy = (max( nodess(:,2) ) - min( nodess(:,2) ))/2;
%      r3(ind-1,1) = -nodess(ind/2,2)+moyy;
%      r3(ind,1) = nodess(ind/2,1)-moyx;
%      g1(2*boun1-1) = r1(2*boun1-1);
%      g2(2*boun1) = r2(2*boun1);
%      g3([2*boun1-1;2*boun1]) = r3([2*boun1-1;2*boun1]);
%   
%      G{i} = [g1,g2,g3]; Gloc = G{i};
%      R{i} = [r1,r2,r3]; Rloc = R{i};
      R{i} = null(K{i}); Rloc = R{i};
      
      nbloq{i} = size(R{i},2);
      
      sign = 1;
      if rem(nsub-i,2) == 0
         sign = -1;
      end
      Gloc = zeros( 2*nnodess, 3 ); % local rigid modes
      if i>1 % Lower boundary
         Gglob( [ 2*newbouns{i-1}-1, 2*newbouns{i-1} ], 3*j-2:3*j ) = ...
%                Gglob( [ 2*newbouns{i-1}-1, 2*newbouns{i-1} ], 3*j-2:3*j ) + ...
                sign*Rloc( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ], : );
         Gloc( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ], : ) = ...
                     sign*Rloc( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ], : );
      end
      if i<nsub
         Gglob( [ 2*newbouns{i}-1, 2*newbouns{i} ], 3*j-2:3*j ) = ...
%                Gglob( [ 2*newbouns{i}-1, 2*newbouns{i} ], 3*j-2:3*j ) + ...
                sign*Rloc( [ 2*bounsloc{2*i}-1, 2*bounsloc{2*i} ], : );
         Gloc( [ 2*bounsloc{2*i}-1, 2*bounsloc{2*i} ], : ) = ...
                     sign*Rloc( [ 2*bounsloc{2*i}-1, 2*bounsloc{2*i} ], : );
      end
      G{i} = Gloc;
      j = j+1;
   end
   
end
Gglob(:,3*j-2:end) = [];  % reshape Gglob

% rework the rigid modes
if j > 1
   % orthonormalize Gglob (no need to do that)
%   g1 = Gglob(:,1); g2 = Gglob(:,2); g3p = Gglob(:,3);
%   g3 = g3p - (g3p'*g1)/(g1'*g1)*g1 - (g3p'*g2)/(g2'*g2)*g2;
%   
%   ng1 = norm(g1); ng2 = norm(g2); ng3 = norm(g3);
%   g1 = g1/ng1; g2 = g2/ng2; g3 = g3/ng3;
   
   eD = zeros( size(Gglob,2), 1 );  % 0 of the same size as Gglob
   
   % Ensure tr(R{i}) = Gglob
   for i = 1:nsub
      if isempty(R{i}) % We do not orthonormalize empty matrices
         continue;
      end
%      Rloc = R{i};
%      
%      nodess = newnode{i}; nnodess = nno{i};
%      [ b1to2, b2to1 ] = superNodes( nodess, nodes, 1e-6 );
%      g1loc = zeros(2*nno{i},1); g2loc = zeros(2*nno{i},1);
%      g1locn = zeros(2*nno{i},1); g2locn = zeros(2*nno{i},1); g3locn = zeros(2*nno{i},1);
%      
%      globind = [ 2*b1to2-1 ; 2*b1to2 ];
%      locind = [ 1:2:2*nnodess-1 , 2:2:2*nnodess ];
%      
%      g1locn( locind ) = g1( globind );
%      g2locn( locind ) = g2( globind );
%      g3locn( locind ) = g3( globind );
%      Gloc = [g1locn, g2locn,g3locn];  % basically Gglob on nodess
%      
%      gloc = G{i}; g1loc = gloc(:,1); g2loc = gloc(:,2); % previous G
%
%      % Modify R such that G = tr(R)
%      r1 = Rloc(:,1); r2 = Rloc(:,2); r3p = Rloc(:,3);
%      r3 = r3p - (r3p'*g1loc)/(g1'*g1)*r1 - (r3p'*g2loc)/(g2'*g2)*r2;
%      r1 = r1/ng1; r2 = r2/ng2; r3 = r3/ng3;
%      R{i} = [r1,r2,r3];
%      K{i} = [ K{i}, R{i} ; R{i}', zeros(size(R{i},2)) ];

      
%      sign = -1;
%      if rem(nsub-i,2) == 0
%         sign = 1;
%      end
      j = sdtorm(i);   % Recover the rigid mode no
      eD(3*j-2:3*j) = eD(3*j-2:3*j) + R{i}'*f{i};

      K{i} = [ K{i}, G{i} ; G{i}', zeros(size(G{i},2)) ];
      f{i} = [ f{i} ; zeros( size(G{i},2), 1 ) ];

   end

   %Gglob = [g1,g2,g3];
   P = eye(2*nnodes) - Gglob*inv(Gglob'*Gglob)*Gglob';  % The projector
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ( constraint - ) CG algorithm
Lambda = zeros(2*nnodes,1);       % Iterate (on the global mesh for convenience)
Res    = zeros(2*nnodes,niter+1);
Rez    = zeros(2*nnodes,niter+1); % Residual without orthogonal projection
Zed    = zeros(2*nnodes,niter+1);
d      = zeros(2*nnodes,niter+1);
Ad     = zeros(2*nnodes,niter+1);
Az     = zeros(2*nnodes,niter+1); % Ad without the projection
resid  = zeros( niter+1,1 );
error  = zeros( niter+1,1 );  % Useless for now (excepted for debug)

% Compute Rhs
b = zeros( 2*nnodes, 1 );
for i = 1:nsub
   ui = K{i}\f{i};
   u = zeros( 2*nnodes, 1 );
   
   sign = -1;
   if rem(nsub-i,2) == 0
      sign = 1;
   end
   
   if i>1 % Lower boundary
      u( [ 2*newbouns{i-1}-1, 2*newbouns{i-1} ]  ) =...
                  sign * ui( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ] );
   end
   if i<nsub
      u( [ 2*newbouns{i}-1, 2*newbouns{i} ] ) =...
                  sign * ui( [ 2*bounsloc{2*i}-1, 2*bounsloc{2*i} ] );
   end
   b = b+u;
end

% initialize
if norm(Gglob) ~= 0
   Lambda = -Gglob*inv(Gglob'*Gglob)*eD;
end

% Compute Ax0
Axz = zeros( 2*nnodes, 1 );
for i = 1:nsub
   floc = zeros( 2*nno{i} + nbloq{i}, 1 );
   if i>1 % Lower boundary
      floc( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ] ) =...
                  Lambda( [ 2*newbouns{i-1}-1, 2*newbouns{i-1} ] );
   end
   if i<nsub
      floc( [ 2*bounsloc{2*i}-1, 2*bounsloc{2*i} ] ) =...
                  Lambda( [ 2*newbouns{i}-1, 2*newbouns{i} ] );
   end
   
   ui = K{i}\floc;
   u = zeros( 2*nnodes, 1 );
   if i>1 % Lower boundary
      u( [ 2*newbouns{i-1}-1, 2*newbouns{i-1} ]  ) =...
                  ui( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ] );
   end
   if i<nsub
      u( [ 2*newbouns{i}-1, 2*newbouns{i} ] ) =...
                  ui( [ 2*bounsloc{2*i}-1, 2*bounsloc{2*i} ] );
   end
   Axz = Axz+u;
end

Rez(:,1) = b - Axz;
if norm(Gglob) ~= 0
   Res(:,1) = P'*(b - Axz);
else
   Res(:,1) = b - Axz;
end

if precond == 1
   Zed(:,1) = zeros( 2*nnodes, 1 );
   for i = 1:nsub

      sign = -1;
      if rem(nsub-i,2) == 0
         sign = 1;
      end
   
      uloc = zeros( 2*nno{i}, 1 );
      if i>1 % Lower boundary
         uloc( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ] ) =...
                     Res( [ 2*newbouns{i-1}-1, 2*newbouns{i-1} ] , 1 );
      end
      if i<nsub
         uloc( [ 2*bounsloc{2*i}-1, 2*bounsloc{2*i} ] ) =...
                     Res( [ 2*newbouns{i}-1, 2*newbouns{i} ] , 1 );
      end
      udir = [ zeros(2*nno{i},1) ; Ct{i}'*uloc ];
      ui = Kp{i}\udir;
      frea = -Ct{i}*ui( 2*nno{i}+1:end );  % Reaction forces on all the bound
      freac = zeros( 2*nnodes, 1 );      % Reaction only on the added boundaries
      if i>1 % Lower boundary
         freac( [ 2*newbouns{i-1}-1, 2*newbouns{i-1} ]  ) =...
                     frea( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ] );
      end
      if i<nsub
         freac( [ 2*newbouns{i}-1, 2*newbouns{i} ] ) =...
                     frea( [ 2*bounsloc{2*i}-1, 2*bounsloc{2*i} ] );
      end
      Zed(:,1) = Zed(:,1) + freac;
   end
   
else
   Zed(:,1) = Res(:,1);
end

if norm(Gglob) ~= 0
   Zed(:,1) = P*Zed(:,1);
else
   Zed(:,1) = Zed(:,1);
end

resid(1) = norm( Res(:,1) );

d(:,1) = Zed(:,1);

%% Loop over the iterations
iter = 0;  % in case niter=0
for iter = 1:niter

   % Compute Ad
   Ad(:,iter) = zeros( 2*nnodes, 1 );
   for i = 1:nsub
      floc = zeros( 2*size(newnode{i},1) + nbloq{i}, 1 );
      % Rem : normally, there should be a sign stuff on Lambda AND u, but by the
      % law of  - * - = +, I put no sign at all.
      if i>1 % Lower boundary
         floc( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ] ) =...
                     d( [ 2*newbouns{i-1}-1, 2*newbouns{i-1} ], iter );
      end
      if i<nsub
         floc( [ 2*bounsloc{2*i}-1, 2*bounsloc{2*i} ] ) =...
                     d( [2*newbouns{i}-1, 2*newbouns{i} ], iter );
      end
      
      ui = K{i}\floc;
      u = zeros( 2*nnodes, 1 );
      if i>1 % Lower boundary
         u( [ 2*newbouns{i-1}-1, 2*newbouns{i-1} ]  ) =...
                     ui( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ] );
      end
      if i<nsub
         u( [ 2*newbouns{i}-1, 2*newbouns{i} ] ) =...
                     ui( [ 2*bounsloc{2*i}-1, 2*bounsloc{2*i} ] );
      end
      Ad(:,iter) = Ad(:,iter) + u;
   end
   Az(:,iter) = Ad(:,iter);
   if norm(Gglob) ~= 0
      Ad(:,iter) = P'*Ad(:,iter);
   end

   den = (d(:,iter)'*Ad(:,iter));
   d(:,iter) = d(:,iter)/sqrt(den); Ad(:,iter) = Ad(:,iter)/sqrt(den);
   Az(:,iter) = Az(:,iter)/sqrt(den);
   num = Res(:,iter)'*d(:,iter);
   Lambda        = Lambda + d(:,iter)*num;
   Res(:,iter+1) = Res(:,iter) - Ad(:,iter)*num;
   Rez(:,iter+1) = Rez(:,iter) - Az(:,iter)*num;
   
   if precond == 1
      Zed(:,iter+1) = zeros( 2*nnodes, 1 );
      for i = 1:nsub
   
         sign = -1;
         if rem(nsub-i,2) == 0
            sign = 1;
         end
      
         uloc = zeros( 2*nno{i}, 1 );
         if i>1 % Lower boundary
            uloc( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ] ) =...
                        Res( [ 2*newbouns{i-1}-1, 2*newbouns{i-1} ] , iter+1 );
         end
         if i<nsub
            uloc( [ 2*bounsloc{2*i}-1, 2*bounsloc{2*i} ] ) =...
                        Res( [ 2*newbouns{i}-1, 2*newbouns{i} ] , iter+1 );
         end
         udir = [ zeros(2*nno{i},1) ; Ct{i}'*uloc ];
         ui = Kp{i}\udir;
         frea = -Ct{i}*ui( 2*nno{i}+1:end );  % Reaction forces on all the bound
         freac = zeros( 2*nnodes, 1 );      % Reaction only on the added boundaries
         if i>1 % Lower boundary
            freac( [ 2*newbouns{i-1}-1, 2*newbouns{i-1} ]  ) =...
                        frea( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ] );
         end
         if i<nsub
            freac( [ 2*newbouns{i}-1, 2*newbouns{i} ] ) =...
                        frea( [ 2*bounsloc{2*i}-1, 2*bounsloc{2*i} ] );
         end
         Zed(:,iter+1) = Zed(:,iter+1) + freac;
      end
   else
      Zed(:,iter+1) = Res(:,iter+1);
   end
   
   if norm(Gglob) ~= 0
      Zed(:,iter+1) = P*Zed(:,iter+1);
      resid(iter+1) = norm( P*Res(:,iter+1) );
   else
      Zed(:,iter+1) = Zed(:,iter+1);
      resid(iter+1) = norm(Res(:,iter+1));
   end
%   resid(iter+1) = norm(Zed(:,iter+1));
   % Orthogonalization
   d(:,iter+1) = Zed(:,iter+1);
   for jter=1:iter
       betaij = ( Zed(:,iter+1)'*Ad(:,jter) );
       d(:,iter+1) = d(:,iter+1) - d(:,jter) * betaij;
   end
   
end

figure; % plot the residual
plot(log10(resid/resid(1)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the final solution and assemble it
usol = zeros( 2*nnodes, 1 );

if norm(Gglob) ~= 0
   alphaD = inv(Gglob'*Gglob)*Gglob'*Rez(:,niter+1);
end

for i = 1:nsub
   nodess    = newnode{i};
%   elementss = newelem{i};
   nnodess   = nno{i};
   
   sign = 1;
   if rem(nsub-i,2) == 0
      sign = -1;
   end

   floc = zeros( 2*size(newnode{i},1) + nbloq{i}, 1 );
   if i>1 % Lower boundary
      floc( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ] ) =...
                  sign * Lambda( [ 2*newbouns{i-1}-1, 2*newbouns{i-1} ] );
   end
   if i<nsub % upper boundary
      floc( [ 2*bounsloc{2*i}-1, 2*bounsloc{2*i} ] ) =...
                  sign * Lambda( [ 2*newbouns{i}-1, 2*newbouns{i} ] );
   end

   ui = K{i}\(floc+f{i});
   if ~isempty(R{i}) % Rigid modes
      j = sdtorm(i);
      ui = ui + [ R{i}*alphaD(3*j-2:3*j) ; zeros(nbloq{i},1) ];
   end
   
   u = zeros( 2*nnodes, 1 );
   [ b1to2, b2to1 ] = superNodes( nodess, nodes, 1e-6 ); % Yeah, I'm a pig
   u( [2*b1to2-1 ; 2*b1to2] ) = ui( [1:2:2*nnodess-1, 2:2:2*nnodess] );
   
   % Special treatment for the interface nodes
   if i>1 % Lower boundary
      u( [ 2*newbouns{i-1}-1, 2*newbouns{i-1} ]  ) =...
                  .5 * ui( [ 2*bounsloc{2*i-1}-1, 2*bounsloc{2*i-1} ] );
   end
   if i<nsub
      u( [ 2*newbouns{i}-1, 2*newbouns{i} ] ) =...
                  .5 * ui( [ 2*bounsloc{2*i}-1, 2*bounsloc{2*i} ] );
   end
   
   usol = usol + u;
end

uti = reshape(usol,2,[])';  ux = uti(:,1);  uy = uti(:,2);
plotGMSH({ux,'U_x';uy,'U_y';usol,'U_vect'}, elements, nodes(:,[1,2]), 'output/solution');
total_error = norm(usol-uref)/norm(uref);
plotGMSH({usol-uref,'Error'}, elements, nodes(:,[1,2]), 'output/error');