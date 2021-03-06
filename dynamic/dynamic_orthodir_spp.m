% 26/01/2018
% Cauchy Dynamique avec ORTHODIR et SPP

close all;
clear all;

% Parameters
E       = 70000;   % MPa : Young modulus
nu      = 0.3;     % Poisson ratio
fscalar = 1;       % N.mm-1 : Loading on the plate
niter   = 800;
br      = 0.;      % noise
dt      = 2e-6;      % s : time discrteization parameter
rho     = 7500e-9; % kg.mm-3 : volumic mass
beta    = .25;     %
gamma   = .5;      % Newmark parameters
precond = 1;       % Use the retro-preconditionner

mat     = [0, E, nu];

% Boundary conditions
% first index  : index of the boundary
% second index : 1=x, 2=y
% third        : value
% [0,1,value] marks a dirichlet regularization therm on x
% dirichlet = [ 1,1,0; 1,2,0 ];
dirichlet = [ 1,1,0; 1,2,0 ; 3,1,0; 3,2,0 ];
%neumann   = [ 2,1,fscalar;
%              4,1,fscalar ];
neumann   = [ 2,1,-fscalar];%,0,-fscalar ];

% Import the mesh
[ nodes,elements,ntoelem,boundary,order ] = readmesh( 'meshes/plate.msh' );
nnodes = size(nodes,1);

% Extract the index of the boundary
[ node2b2, b2node2 ] = mapBound( 2, boundary, nnodes );
index    = 2*b2node2-1;
index    = index(size(index):-1:1);
indto2   = [ 2*b2node2; 2*b2node2-1 ];

% Then, build the stiffness matrix :
[K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
Kinter = K(1:2*nnodes,1:2*nnodes);

% Mass and (zero) damping matrices
M0 = rho * mass_mat(nodes, elements);
M = [ M0 , zeros(2*nnodes,nbloq) ; zeros(2*nnodes,nbloq)' , zeros(nbloq) ];
C = zeros(size(M));
%fatKinter = Kinter + 1/(beta*dt^2)*M0; % fatKinter*(u-u0) = f

% The right hand side :
f = loading(nbloq,nodes,boundary,neumann);
T = 0:1:50;  fa = [ f*T/T(end) , f*(1-T/T(end)), zeros(2*nnodes+nbloq,450) ];
u0 = zeros(2*nnodes+nbloq,1); v0 = zeros(2*nnodes+nbloq,1);
a0 = zeros(2*nnodes+nbloq,1);

% Solve the problem :
[uin,vin,ain] = Newmark (M, C, K, fa, u0, v0, a0, dt, beta, gamma);

% Extract displacement and Lagrange multiplicators :
uref = uin(1:2*nnodes,:); vref = uin(1:2*nnodes,:);
aref = uin(1:2*nnodes,:); ntime = size(uref,2);

lagr = uin(2*nnodes+1:end,:);
urefb = ( 1 + br*randn(2*nnodes,size(uref,2)) ) .* uref;
lagrb = ( 1 + br*randn(nbloq,size(uref,2)) ) .* lagr;
fref = f( 1:2*nnodes,: ); % Reaction forces

% Compute stress :
%sigma = stress(uref,E,nu,nodes,elements,order,1,ntoelem);

N = rigidModes(nodes); P = eye(2*nnodes) - N*((N'*N)\N');
urefo = P*uref;
plotGMSH({uref,'U_vect';urefo,'U_rigide'}, elements, nodes, 'output/reference');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% D problem
% dirichlet1 = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ; 4,1,0 ; 4,2,0];
dirichlet1 = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ; 3,1,0 ; 3,2,0 ; 4,1,0 ; 4,2,0];
% dirichlet1 = [ 2,1,0 ; 2,2,0 ; 4,1,0 ; 4,2,0];
neumann1   = []; % No load
[K1,C1t,nbloq1,node2c1,c2node1] = Krig (nodes,elements,E,nu,order,boundary,dirichlet1);
M1 = [ M0 , zeros(2*nnodes,nbloq1) ; zeros(2*nnodes,nbloq1)' , zeros(nbloq1) ];
C1 = zeros(size(M1));
u01 = zeros(2*nnodes+nbloq1,1); v01 = zeros(2*nnodes+nbloq1,1);
a01 = zeros(2*nnodes+nbloq1,1);
fatK1 = K1 + 1/(beta*dt^2)*M1 + gamma/(beta*dt)*C1;
invK1 = fatK1\eye(size(fatK1));

% N problem
% dirichlet2 = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ];
dirichlet2 = [ 1,1,0 ; 1,2,0 ; 2,1,0 ; 2,2,0 ; 3,1,0 ; 3,2,0 ];
% dirichlet2 = [ 2,1,0 ; 2,2,0 ];
neumann2   = []; % fr = 0
[K2,C2t,nbloq2,node2c2,c2node2] = Krig (nodes,elements,E,nu,order,boundary,dirichlet2);
M2 = [ M0 , zeros(2*nnodes,nbloq2) ; zeros(2*nnodes,nbloq2)' , zeros(nbloq2) ];
C2 = zeros(size(M2));
u02 = zeros(2*nnodes+nbloq2,1); v02 = zeros(2*nnodes+nbloq2,1);
a02 = zeros(2*nnodes+nbloq2,1);
fatK2 = K2 + 1/(beta*dt^2)*M2 + gamma/(beta*dt)*C2;
invK2 = fatK2\eye(size(fatK2));

% Dual N problem
% dirichlet2d = [ 1,1,0 ; 1,2,0 ];
dirichlet2d = [ 1,1,0 ; 1,2,0 ; 3,1,0 ; 3,2,0 ];
neumann2d   = []; % fr = 0
[K2d,C2td,nbloq2d,node2c2d,c2node2d] = Krig (nodes,elements,E,nu,order,boundary,dirichlet2d);
M2d = [ M0 , zeros(2*nnodes,nbloq2d) ; zeros(2*nnodes,nbloq2d)' , zeros(nbloq2d) ];
C2d = zeros(size(M2d));
u02d = zeros(2*nnodes+nbloq2d,1); v02d = zeros(2*nnodes+nbloq2d,1);
a02d = zeros(2*nnodes+nbloq2d,1);
fatK2d = K2d + 1/(beta*dt^2)*M2d + gamma/(beta*dt)*C2d;
invK2d = fatK2d\eye(size(fatK2d));

error    = zeros(niter+1,1);
residual = zeros(niter+1,1);
regulari = zeros(niter+1,1);
Itere    = zeros(2*nnodes,ntime);
Delta    = zeros(niter,1);
qRec     = zeros(2*nnodes*ntime, niter); % Vectorial storing format for q_i
pRec     = zeros(2*nnodes*ntime, niter);
% indtot   = 1:2*nnodes;            % Choose the index for the scalar product
indtot   = indto2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Perform A x0 :
fdir = dirichletRhs2(Itere, 2, c2node1, boundary, nnodes ); % Solve D
f1 = fdir;
[uin1] = Newmark (M1, C1, K1, f1, u01, v01, a01, dt, beta, gamma, invK1);
lagr1 = uin1(2*nnodes+1:end,:);
fri1 = zeros(2*nnodes, ntime);
for i=1:size(lagr1,2)
   fri1(:,i) = lagr2forces2( lagr1(:,i), c2node1, 2, boundary, nnodes );
end
fdir = dirichletRhs2(Itere, 2, c2node2, boundary, nnodes ); % Solve N
f2 = fdir;
[uin2] = Newmark (M2, C2, K2, f2, u02, v02, a02, dt, beta, gamma, invK2);
lagr2 = uin2(2*nnodes+1:end,:);
fri2 = zeros(2*nnodes, ntime);
for i=1:size(lagr1,2)
   fri2(:,i) = lagr2forces2( lagr2(:,i), c2node2, 2, boundary, nnodes );
end
atimesItere = fri1 - fri2;
%% End Ax0

%% Write Rhs :
fdir = dirichletRhs2(urefb, 4, c2node1, boundary, nnodes ); % Solve D
f1 = fdir;
[uin1] = Newmark (M1, C1, K1, f1, u01, v01, a01, dt, beta, gamma, invK1);
lagr1 = uin1(2*nnodes+1:end,:);
fri1 = zeros(2*nnodes, ntime);
for i=1:size(lagr1,2)
   fri1(:,i) = lagr2forces2( lagr1(:,i), c2node1, 2, boundary, nnodes );
end
fri2 = zeros(2*nnodes, ntime); % Solve N : Zero loading
b = fri2 - fri1; % Remember there is a "-"
%% End RHS

Res = b - atimesItere;
if precond == 1
   f2 = zeros( 2*nnodes, ntime ); f2(indtot,:) = Res(indtot,:);
   f2 = fliplr( f2 - f2(:,end)*ones(1,ntime) ); f2 = [ f2 ; zeros( nbloq2d, ntime ) ];
   [uin2, vin2, ain2] = Newmark (M2d, C2d, K2d, f2, u02d, v02d, a02d, -dt, beta, gamma, invK2d);
   ur2 = uin2(1:2*nnodes,:);
   ur2 = ur2 - ur2(:,end)*ones(1,ntime);
   Zed = fliplr( ur2 ); Zed = keepField( Zed, 2, boundary );
else
   Zed = Res;
end
p = Zed;

%plotGMSH({Zed,'Zed';Res,'Res'}, elements, nodes(:,[1,2]), 'output/res');
%bug;
residual(1) = norm(Res(indtot,:),'fro');
error(1)    = norm(Itere(indtot,:)-uref(indtot,:),'fro')/norm(uref(indtot,:),'fro');
regulari(1) = norm(Itere(indtot,:),'fro');

%% Perform Q1 = A P1 :
fdir = dirichletRhs2(p, 2, c2node1, boundary, nnodes ); % Solve D
f1 = fdir;
[uin1] = Newmark (M1, C1, K1, f1, u01, v01, a01, dt, beta, gamma, invK1);
lagr1 = uin1(2*nnodes+1:end,:);
fri1 = zeros(2*nnodes, ntime);
for i=1:size(lagr1,2)
   fri1(:,i) = lagr2forces2( lagr1(:,i), c2node1, 2, boundary, nnodes );
end
fdir = dirichletRhs2(p, 2, c2node2, boundary, nnodes ); % Solve N
f2 = fdir;
[uin2] = Newmark (M2, C2, K2, f2, u02, v02, a02, dt, beta, gamma, invK2);
lagr2 = uin2(2*nnodes+1:end,:);
fri2 = zeros(2*nnodes, ntime);
for i=1:size(lagr1,2)
   fri2(:,i) = lagr2forces2( lagr2(:,i), c2node2, 2, boundary, nnodes );
end
q = fri1 - fri2;
%% End Q1 = A P1
%plotGMSH({uin1(1:2*nnodes,:),'u1';uin2(1:2*nnodes,:),'u2';q,'q';Res,'Res'}, ...
%         elements, nodes(:,[1,2]), 'output/res');
%bug;
pRec(:,1) = p(:);
qRec(:,1) = q(:);

for iter = 1:niter

   Delta(iter)      = norm(q(indtot,:),'fro')^2;  %q(:,iter)'*q(:,iter);
   gammai           = sum(sum(q(indtot,:).*Res(indtot,:))); %q(:,iter)'*Res;
   alphai           = gammai/Delta(iter);

   Itere            = Itere + p*alphai;
   Res              = Res - q*alphai;
   
   if precond == 1
      f2 = zeros( 2*nnodes, ntime ); f2(indtot,:) = Res(indtot,:);
      f2 = fliplr( f2 - f2(:,end)*ones(1,ntime) ); f2 = [ f2 ; zeros( nbloq2d, ntime ) ];
      [uin2, vin2, ain2] = Newmark (M2d, C2d, K2d, f2, u02d, v02d, a02d, -dt, beta, gamma, invK2d);
      ur2 = uin2(1:2*nnodes,:);
      ur2 = ur2 - ur2(:,end)*ones(1,ntime);
      Zed = fliplr( ur2 ); Zed = keepField( Zed, 2, boundary );
   else
      Zed = Res;
   end

   residual(iter+1) = norm(Res(indtot,:),'fro');
   error(iter+1)    = norm(Itere(indtot,:)-uref(indtot,:),'fro') ...
                      / norm(uref(indtot,:),'fro');
   regulari(iter+1) = norm(Itere(indtot,:),'fro');

   %% Perform Ari = A*Zed :
   fdir = dirichletRhs2(Zed, 2, c2node1, boundary, nnodes ); % Solve D
   f1 = fdir;
   [uin1] = Newmark (M1, C1, K1, f1, u01, v01, a01, dt, beta, gamma, invK1);
   lagr1 = uin1(2*nnodes+1:end,:);
   fri1 = zeros(2*nnodes, ntime);
   for i=1:size(lagr1,2)
      fri1(:,i) = lagr2forces2( lagr1(:,i), c2node1, 2, boundary, nnodes );
   end
   fdir = dirichletRhs2(Zed, 2, c2node2, boundary, nnodes ); % Solve N
   f2 = fdir;
   [uin2] = Newmark (M2, C2, K2, f2, u02, v02, a02, dt, beta, gamma, invK2);
   lagr2 = uin2(2*nnodes+1:end,:);
   fri2 = zeros(2*nnodes, ntime);
   for i=1:size(lagr1,2)
      fri2(:,i) = lagr2forces2( lagr2(:,i), c2node2, 2, boundary, nnodes );
   end
   Ari = fri1 - fri2;
   %% End Ari = A*Res
   
   %% Orthogonalization
   p = Zed; % Zed is dead
   q = Ari;

   for jter=1:iter
      pj = reshape( pRec(:,jter), [2*nnodes, ntime] );
      qj = reshape( qRec(:,jter), [2*nnodes, ntime] ); % Recover the space-time matrix
      phiij  = sum(sum(qj(indtot,:).*Ari(indtot,:))); %q(:,jter)'*Ari;
      betaij = phiij/Delta(jter);
      p = p - pj * betaij;
      q = q - qj * betaij;
   end

   pRec(:,iter+1) = p(:);
   qRec(:,iter+1) = q(:);

end

residual = residual/residual(1); % Normalize

figure;
hold on;
set(gca, 'fontsize', 20);
plot(log10(error),'Color','black')
plot(log10(residual),'Color','red')
legend('error (log)', 'residual (log)')

%hold on;
% plot(log10(ferror1),'Color','black')
% plot(log10(ferror2),'Color','blue')
% plot(log10(fresidual),'Color','red')

%%L-curve
%figure
%loglog(residual(2:end),regulari(2:end));

% Plot solution
% hold on;
% set(gca, 'fontsize', 15);
% set(gca,'ylim',[-3e-5 3e-5])
% plot(uref(index));
% plot(u2(index),'Color','red');

% Find the middle node of the bound
xmax = max(nodes(:,1)); ymoy = .5 * ( max(nodes(:,2) + min(nodes(:,2))) );
dmin = xmax^2;
for i=1:nnodes
   x = nodes(i,1); y = nodes(i,2);
   d = (x-xmax)^2 + (y-ymoy)^2;
   if d<dmin
      dmin = d;
      imin = i;
   end
end

%% Compute stress :
%sigma = stress(u2,E,nu,nodes,elements,order,1,ntoelem);
% Output :
Itereo = P*Itere;
plotGMSH({Itere,'U1';Itereo,'U1_rigide'}, elements, nodes(:,[1,2]), 'output/solution');

utref = uref(2*imin-1,:); utsol = Itere(2*imin-1,:);
figure;
hold on;
plot( utref, 'Color', 'blue' );
plot( utsol, 'Color', 'red' );
legend( 'reference', 'solution' );

%% Final resolution
fdir = dirichletRhs2( Itere, 2, c2node1, boundary, nnodes );
f4 = dirichletRhs2( urefb, 4, c2node1, boundary, nnodes );
f1 = fdir + f4;
[uin1] = Newmark (M1, C1, K1, f1, u01, v01, a01, dt, beta, gamma, invK1);
usol = uin1(1:2*nnodes,:);

erroru = norm(usol-uref,'fro')/norm(uref,'fro');