% 24/01/2018
% Cauchy Dynamique avec ORTHODIR

close all;
clear all;

% Parameters
E       = 70000;   % MPa : Young modulus
nu      = 0.3;     % Poisson ratio
fscalar = 1;       % N.mm-1 : Loading on the plate
niter   = 10;
br      = 0.;      % noise
dt      = 2e-6;      % s : time discrteization parameter
rho     = 7500e-9; % kg.mm-3 : volumic mass
beta    = .25;     %
gamma   = .5;      % Newmark parameters
precond = 0;       % Use the retro-preconditionner

mat     = [0, E, nu];

% Boundary conditions
% first index  : index of the boundary
% second index : 1=x, 2=y
% third        : value
% [0,1,value] marks a dirichlet regularization therm on x
%dirichlet = [ 1,1,0; 1,2,0
%              3,1,0; 3,2,0 ];
dirichlet = [];
%neumann   = [ 2,1,fscalar;
%              4,1,fscalar ];
neumann   = [ 2,1,-fscalar ];%,0,-fscalar ];

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

N = null(Kinter); P = eye(2*nnodes) - N*((N'*N)\N');
urefo = P*uref;
plotGMSH({uref,'U_vect';urefo,'U_rigide'}, elements, nodes, 'output/reference');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DN problem
dirichlet1 = [2,1,0;2,2,0];
%dirichlet1 = [ 2,1,0; 2,2,0
%               1,1,0; 1,2,0
%               3,1,0; 3,2,0 ];
neumann1   = []; % Zero load
[K1,C1t,nbloq1,node2c1,c2node1] = Krig (nodes,elements,E,nu,order,boundary,dirichlet1);
M1 = [ M0 , zeros(2*nnodes,nbloq1) ; zeros(2*nnodes,nbloq1)' , zeros(nbloq1) ];
C1 = zeros(size(M1));
u01 = zeros(2*nnodes+nbloq1,1); v01 = zeros(2*nnodes+nbloq1,1);
a01 = zeros(2*nnodes+nbloq1,1);
fatK1 = K1 + 1/(beta*dt^2)*M1 + gamma/(beta*dt)*C1;
invK1 = fatK1\eye(size(fatK1));

% ND problem
%dirichlet2 = [ 4,1,0; 4,2,0
%               1,1,0; 1,2,0
%               3,1,0; 3,2,0 ];
dirichlet2 = [4,1,0;4,2,0];
neumann2   = [];% [2,1,lagr1; 2,2,lagr1]
                % is managed by lagr2forces
[K2,C2t,nbloq2,node2c2,c2node2] = Krig (nodes,elements,E,nu,order,boundary,dirichlet2);
M2 = [ M0 , zeros(2*nnodes,nbloq2) ; zeros(2*nnodes,nbloq2)' , zeros(nbloq2) ];
C2 = zeros(size(M2));
u02 = zeros(2*nnodes+nbloq2,1); v02 = zeros(2*nnodes+nbloq2,1);
a02 = zeros(2*nnodes+nbloq2,1);
fatK2 = K2 + 1/(beta*dt^2)*M2 + gamma/(beta*dt)*C2;
invK2 = fatK2\eye(size(fatK2));

error    = zeros(niter+1,1);
residual = zeros(niter+1,1);
regulari = zeros(niter+1,1);
Itere    = zeros(2*nnodes,ntime);
Delta    = zeros(niter,1);
qRec     = zeros(2*nnodes*ntime, niter); % Vectorial storing format for q_i
pRec     = zeros(2*nnodes*ntime, niter);
indtot   = 1:2*nnodes;            % Choose the index for the scalar product
%indtot   = indto2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Perform A x0 :
fdir = dirichletRhs2(Itere, 2, c2node1, boundary, nnodes ); % Solve DN
f1 = fdir;
[uin1, vin1, ain1] = Newmark (M1, C1, K1, f1, u01, v01, a01, dt, beta, gamma, invK1);
lagr1 = uin1(2*nnodes+1:end,:);
fri = zeros(2*nnodes, ntime); % Solve ND
for i=1:size(lagr1,2)
   fri(:,i) = lagr2forces2( lagr1(:,i), c2node1, 2, boundary, nnodes );
end
fr = [ fri; zeros( nbloq2, ntime ) ];
f2 = fr;
[uin2, vin2, ain2] = Newmark (M2, C2, K2, f2, u02, v02, a02, dt, beta, gamma, invK2);
atimesItere = Itere - uin2(1:2*nnodes,:);
%% End Ax0

%% Write Rhs :
%u1 = zeros(2*nnodes,ntime); % Solve DN (there should be a Newmark, but loading is 0)
%lagr1 = uin1(2*nnodes+1:end,:);
%fri = zeros(size(u1)); % Solve ND
%for i=1:size(lagr1,2)
% fri(:,i) = lagr2forces2( lagr1(:,i), c2node1, 2, boundary, nnodes );
%end
%fr    = [ fri; zeros( size(C2t,2), size(fri,2) ) ];
fdir4 = dirichletRhs2(urefb, 4, c2node2, boundary, nnodes);
f2 = fdir4; % fr = 0 because loading is 0
[uin2, vin2, ain2] = Newmark (M2, C2, K2, f2, u02, v02, a02, dt, beta, gamma, invK2);
b = uin2(1:2*nnodes,:);
%% End RHS

Res = b - atimesItere;
if precond == 1
   udir = Res - Res(:,end)*ones(1,ntime); udir = fliplr(udir);
   fdir = dirichletRhs2( udir, 2, c2node1, boundary, nnodes ); % Solve DN
   f1 = fdir;
   [uin1, vin1, ain1] = Newmark (M1, C1, K1, f1, u01, v01, a01, -dt, beta, gamma, invK1);
   lagr1 = uin1(2*nnodes+1:end,:);
   fri = zeros(2*nnodes, ntime); % Solve ND
   for i=1:size(lagr1,2)
      fri(:,i) = lagr2forces2( lagr1(:,i), c2node1, 2, boundary, nnodes );
   end
   fr = [ fri; zeros( nbloq2, ntime ) ];
   f2 = fr;
   [uin2, vin2, ain2] = Newmark (M2, C2, K2, f2, u02, v02, a02, -dt, beta, gamma, invK2);
   ur2 = uin2(1:2*nnodes,:);
   ur2 = -fliplr(Res) + ur2;
   ur2 = ur2 - ur2(:,end)*ones(1,ntime);
   Zed = fliplr( ur2 );
%for i=1:size(lagr1,2)  %%% I Guess this should be suppressed
%   fri(:,i) = lagr2forces2( lagr1(:,i), c2node1, 2, boundary, nnodes );
%end
%fr = [ fri; zeros( nbloq2, ntime ) ];
%f2 = fr;
%[uin2, vin2, ain2] = Newmark (M2, C2, K2, f2, u02, v02, a02, dt, beta, gamma, invK2);
%atimesItere = Itere - uin2(1:2*nnodes,:);
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
fdir = dirichletRhs2( p, 2, c2node1, boundary, nnodes ); % Solve DN
f1 = fdir;
[uin1, vin1, ain1] = Newmark (M1, C1, K1, f1, u01, v01, a01, dt, beta, gamma, invK1);
lagr1 = uin1(2*nnodes+1:end,:);
fri = zeros(2*nnodes, ntime); % Solve ND
for i=1:size(lagr1,2)
   fri(:,i) = lagr2forces2( lagr1(:,i), c2node1, 2, boundary, nnodes );
end
fr = [ fri; zeros( nbloq2, ntime ) ];
f2 = fr;
[uin2, vin2, ain2] = Newmark (M2, C2, K2, f2, u02, v02, a02, dt, beta, gamma, invK2);
q = p - uin2(1:2*nnodes,:);
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
      udir = Res - Res(:,end)*ones(1,ntime); udir = fliplr(udir);
      fdir = dirichletRhs2( udir, 2, c2node1, boundary, nnodes ); % Solve DN
      f1 = fdir;
      [uin1, vin1, ain1] = Newmark (M1, C1, K1, f1, u01, v01, a01, -dt, beta, gamma, invK1);
      lagr1 = uin1(2*nnodes+1:end,:);
      fri = zeros(2*nnodes, ntime); % Solve ND
      for i=1:size(lagr1,2)
         fri(:,i) = lagr2forces2( lagr1(:,i), c2node1, 2, boundary, nnodes );
      end
      fr = [ fri; zeros( nbloq2, ntime ) ];
      f2 = fr;
      [uin2, vin2, ain2] = Newmark (M2, C2, K2, f2, u02, v02, a02, -dt, beta, gamma, invK2);
      ur2 = uin2(1:2*nnodes,:);
      ur2 = -fliplr(Res) + ur2;
      ur2 = ur2 - ur2(:,end)*ones(1,ntime);
      Zed = fliplr( ur2 );
   else
      Zed = Res;
   end

   residual(iter+1) = norm(Res(indtot,:),'fro');
   error(iter+1)    = norm(Itere(indtot,:)-uref(indtot,:),'fro') ...
                      / norm(uref(indtot,:),'fro');
   regulari(iter+1) = norm(Itere(indtot,:),'fro');

   %% Perform Ari = A*Res :
   fdir = dirichletRhs2( Zed, 2, c2node1, boundary, nnodes ); % Solve DN
   f1 = fdir;
   [uin1, vin1, ain1] = Newmark (M1, C1, K1, f1, u01, v01, a01, dt, beta, gamma, invK1);
   lagr1 = uin1(2*nnodes+1:end,:);
   fri = zeros(2*nnodes, ntime); % Solve ND
   for i=1:size(lagr1,2)
    fri(:,i) = lagr2forces2( lagr1(:,i), c2node1, 2, boundary, nnodes );
   end
   fr = [ fri; zeros( nbloq2, ntime ) ];
   f2 = fr;
   [uin2, vin2, ain2] = Newmark (M2, C2, K2, f2, u02, v02, a02, dt, beta, gamma, invK2);
   Ari = Zed - uin2(1:2*nnodes,:);
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

%% Compute stress :
%sigma = stress(u2,E,nu,nodes,elements,order,1,ntoelem);
% Output :
Itereo = P*Itere;
plotGMSH({Itere,'U1';Itereo,'U1_rigide'}, elements, nodes(:,[1,2]), 'output/solution');