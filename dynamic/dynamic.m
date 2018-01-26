%10/11/2017
%Cauchy Dynamique

close all;
clear all;

addpath(genpath('./tools'))

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
omega   = 1;      % Relaxation parameter

mat     = [0, E, nu];

% Boundary conditions
% first index  : index of the boundary
% second index : 1=x, 2=y
% third        : value
% [0,1,value] marks a dirichlet regularization therm on x
dirichlet = [];%[ 1,1,0; 1,2,0
              %3,1,0; 3,2,0 ];
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
indtot   = [ 2*b2node2; 2*b2node2-1 ]; nindtot = size(indtot,1);

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
T = 0:1:50;  fa = [ f*T/T(end) , f*(1-T/T(end)) ];
u0 = zeros(2*nnodes+nbloq,1); v0 = zeros(2*nnodes+nbloq,1);
a0 = zeros(2*nnodes+nbloq,1);

% Solve the problem :
[uin,vin,ain] = Newmark (M, C, K, fa, u0, v0, a0, dt, beta, gamma);

% Extract displacement and Lagrange multiplicators :
uref = uin(1:2*nnodes,:); vref = uin(1:2*nnodes,:);
aref = uin(1:2*nnodes,:);

lagr = uin(2*nnodes+1:end,:);
urefb = ( 1 + br*randn(2*nnodes,size(uref,2)) ) .* uref;
lagrb = ( 1 + br*randn(nbloq,size(uref,2)) ) .* lagr;
fref = f( 1:2*nnodes,: ); % Reaction forces

% Compute stress :
sigma = stress(uref,E,nu,nodes,elements,order,1,ntoelem);
% Output :
plotGMSH({uref,'U_vect';sigma,'stress'}, elements, nodes, 'output/reference');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init :
u1    = uref-uref;
u2    = u1;
u1p   = u1; v1p   = u1; a1p   = u1;  % Remember, u1 = 0 at this point
u2p   = u1; v2p   = u1; a2p   = u1;
fri   = u1; frip = fri;
v     = u1;
theta = ones(niter+1,1); % First relaxation parameter

% DN problem
dirichlet1 = [2,1,0;2,2,0];
neumann1   = []; % Zero load
[K1,C1t,nbloq1,node2c1,c2node1] = Krig (nodes,elements,E,nu,order,boundary,dirichlet1);
M1 = [ M0 , zeros(2*nnodes,nbloq1) ; zeros(2*nnodes,nbloq1)' , zeros(nbloq1) ];
C1 = zeros(size(M1));
u01 = zeros(2*nnodes+nbloq1,1); v01 = zeros(2*nnodes+nbloq1,1);
a01 = zeros(2*nnodes+nbloq1,1);

%u01 = urefb(:,end); v01 = vref(:,end); a01 = aref(:,end);
%u01 = [ u01 ; zeros(nbloq1,1) ];                                % Reverse time
%v01 = [ v01 ; zeros(nbloq1,1) ];                                %
%a01 = [ a01 ; zeros(nbloq1,1) ];                                %

% ND problem
dirichlet2 = [4,1,0;4,2,0];
neumann2   = [];% [2,1,lagr1; 2,2,lagr1]
                % is managed by lagr2forces
[K2,C2t,nbloq2,node2c2,c2node2] = Krig (nodes,elements,E,nu,order,boundary,dirichlet2);
M2 = [ M0 , zeros(2*nnodes,nbloq2) ; zeros(2*nnodes,nbloq2)' , zeros(nbloq2) ];
C2 = zeros(size(M2));
u02 = zeros(2*nnodes+nbloq2,1); v02 = zeros(2*nnodes+nbloq2,1);
a02 = zeros(2*nnodes+nbloq2,1);

error1   = zeros(niter,1);
error2   = zeros(niter,1);
residual = zeros(niter,1);
regulari = zeros(niter,1);

ntime = 102;
Amat  = zeros( nindtot*ntime ); % ntime = 51

%for iter = 1:niter
for time = 1:ntime
for iter = 1:nindtot
    u2 = zeros(2*nnodes,ntime); u2( indtot(iter), time ) = 1;
%    if iter > 1 %mod(iter,2) == 2 % Reverse Time
%       u01 = [ u2(:,end) ; zeros(nbloq1,1) ];
%       v01 = [ v2(:,end) ; zeros(nbloq1,1) ];
%       a01 = [ a2(:,end) ; zeros(nbloq1,1) ];
%       u2  = fliplr(u2); %u2(:,1) = []; u2 = [u2,zeros(2*nnodes,1)]; % Off-by-one
%    end

    % Solve DN
    fdir = dirichletRhs2(u2, 2, c2node1, boundary, nnodes );
    f1 = fdir; % f = 0 on 4 (and 1,3)

    [uin1, vin1, ain1] = Newmark (M1, C1, K1, f1, u01, v01, a01, dt, beta, gamma);
    u1 = uin1(1:2*nnodes,:); v1 = vin1(1:2*nnodes,:); a1 = ain1(1:2*nnodes,:); 
    lagr1 = uin1(2*nnodes+1:end,:);
    
    fri = zeros(size(u1));
    for i=1:size(lagr1,2)
       fri(:,i) = lagr2forces2( lagr1(:,i), c2node1, 2, boundary, nnodes );
    end
    
%    if iter > 1 %mod(iter,2) == 2
%       fri  = fliplr(fri); % Re-reverse time
%       u1   = fliplr(u1); v1 = fliplr(v1); a1 = fliplr(a1);
%    end
    
    u1e  = omega*u1 + (1-omega)*u1p;             % Relaxation
    v1e  = omega*v1 + (1-omega)*v1p;             %
    a1e  = omega*a1 + (1-omega)*a1p;             %
    frie = omega*fri + (1-omega)*frip;           %
    u1p  = u1; v1p = v1; a1p = a1; frip = fri;   %
    
    error1(iter) = norm(u1-uref,'fro') / norm(uref,'fro');
    
    % Solve ND
    fr = [ frie; zeros( size(C2t,2), size(fri,2) ) ]; % Give to fr the right size
    fdir4 = dirichletRhs2(urefb, 4, c2node2, boundary, nnodes);
    f2 = fr;% + fdir4;
    
    % Reverse time
%    f2 = fliplr(f2);
%    u02 = u1(:,end); v02 = v1(:,end); a02 = a1(:,end);
%    u02 = [ u02 ; zeros(nbloq2,1) ];
%    v02 = [ v02 ; zeros(nbloq2,1) ];
%    a02 = [ a02 ; zeros(nbloq2,1) ];

    [uin2, vin2, ain2] = Newmark (M2, C2, K2, f2, u02, v02, a02, dt, beta, gamma);
    u2 = uin2(1:2*nnodes,:); v2 = vin2(1:2*nnodes,:); a2 = ain2(1:2*nnodes,:); 
    
%    u2 = fliplr(u2); % Re-reverse time
    
    error2(iter) = norm(u2-uref,'fro')/norm(uref,'fro');
        
%    if iter > 1
%       residual(iter) = norm(u1-u2,'fro')/sqrt(norm(u1,'fro')*norm(u2,'fro'));
%    end
                 
    regulari(iter) = 0;%sqrt(u2'*regul(u2, nodes, boundary, 3));
    
    u2 = u2(indtot,:);
    Amat( :, nindtot*(time-1) + iter ) = u2(:);
end
end
bug
%Amat = load('fields/dynamicOperator.mat'); % Recover the full matrix

residual(1) = 1; % tiny hack

figure;
hold on;
set(gca, 'fontsize', 15);
%set(gca,'ylim',[-3 0])
plot(log10(error1),'Color','black')
plot(log10(error2),'Color','blue')
plot(log10(residual),'Color','red')
legend('error1 (log)','error2 (log)','residual (log)')
%hold on;
% plot(log10(ferror1),'Color','black')
% plot(log10(ferror2),'Color','blue')
% plot(log10(fresidual),'Color','red')
% figure
% L-curve
%loglog(residual,regulari);
% plot(regulari.*residual)

% Plot solution
% hold on;
% set(gca, 'fontsize', 15);
% set(gca,'ylim',[-3e-5 3e-5])
% plot(uref(index));
% plot(u2(index),'Color','red');

%% Compute stress :
%sigma = stress(u2,E,nu,nodes,elements,order,1,ntoelem);
% Output :
plotGMSH({u1,'U1';f1(1:2*nnodes,:),'F1';u2,'U2';f2(1:2*nnodes,:),'F2'},...
         elements, nodes(:,[1,2]), 'output/solution');