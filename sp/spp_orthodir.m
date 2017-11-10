% 24/03/2016
% Algo Steklov-Poincaré primal avec Orthodir

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 70000;  % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 1;      % N.mm-1 : Loading on the plate
niter   = 20;
precond = 0;      % Use a dual precond
mu      = 0;    % Regularization parameter
br      = 0.0;     % noise
omega   = 0;%5e9;      % s-1 : pulsation (dynamic case)
rho     = 7500e-9; % kg.mm-3 : volumic mass

% Boundary conditions
% first index  : index of the boundary
% second index : 1=x, 2=y
% third        : value
% [0,1,value] marks a dirichlet regularization therm on x
dirichlet = [4,1,0;
             4,2,0];
neumann   = [1,2,fscalar;
             2,1,fscalar;
             3,2,fscalar];

% Import the mesh
[ nodes,elements,ntoelem,boundary,order ] = readmesh( 'meshes/plate.msh' );
nnodes = size(nodes,1);

% find the nodes in the corners and suppress the element :
xmax = max(nodes(:,1));
xmin = min(nodes(:,1));
ymax = max(nodes(:,2));
ymin = min(nodes(:,2));
no1  = findNode(xmin, ymin, nodes, 1e-5);
no2  = findNode(xmax, ymin, nodes, 1e-5);
no3  = findNode(xmax, ymax, nodes, 1e-5);
no4  = findNode(xmin, ymax, nodes, 1e-5);

boundaryp1 = suppressBound( boundary, [no3;no4], 3 );
%boundaryp2 = suppressBound( boundary, [no4], 3 );
boundaryp2 = boundaryp1;
%boundary = suppressBound( boundary, [no3;no4], 3 );
% Suppress no2 from 1 for a Schur complement computation
boundarym = suppressBound( boundary, [no2], 2 );
boundarym = suppressBound( boundarym, [no1], 1 );

% Then, build the stiffness matrix :
[K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
Kinter = K(1:2*nnodes, 1:2*nnodes);
M0 = mass_mat(nodes, elements);
M0 = rho*M0;
M = [ M0 , zeros(2*nnodes,nbloq) ; zeros(2*nnodes,nbloq)' , zeros(nbloq) ];
[node2b4, b2node4] = mapBound(4, boundaryp1, nnodes);
[node2b3, b2node3] = mapBound(3, boundaryp1, nnodes);
[node2b1, b2node1] = mapBound(1, boundarym, nnodes);
[node2b2, b2node2] = mapBound(2, boundarym, nnodes);
b2node12 = [b2node1;b2node2];
% The right hand side :
f = loading(nbloq,nodes,boundary,neumann);

% Solve the problem :
uin = (K-omega*M)\f;

% Extract displacement and Lagrange multiplicators :
uref = uin(1:2*nnodes,1);
lagr = uin(2*nnodes+1:end,1);

urefb = ( 1 + br*randn(2*nnodes,1) ) .* uref;

sigma = stress(uref,E,nu,nodes,elements,order,1,ntoelem);
plotGMSH({uref,'Vect_U';sigma,'stress'}, elements, nodes, 'output/reference');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition of the operators

% First problem
dirichlet1 = [4,1,0;4,2,0;
              3,1,0;3,2,0;
              2,1,0;2,2,0;
              1,1,0;1,2,0];

[K1,C1,nbloq1,node2c1,c2node1] = Krig (nodes,elements,E,nu,order,boundaryp1,dirichlet1);
M1 = [ M0 , zeros(2*nnodes,nbloq1) ; zeros(2*nnodes,nbloq1)' , zeros(nbloq1) ];
% Second problem
dirichlet2 = [4,1,0;4,2,0;
              3,1,0;3,2,0];
neumann2   = [1,2,fscalar;
              2,1,fscalar];
neumann0   = [];
[K2,C2,nbloq2,node2c2,c2node2] = Krig (nodes,elements,E,nu,order,boundaryp2,dirichlet2);
M2 = [ M0 , zeros(2*nnodes,nbloq2) ; zeros(2*nnodes,nbloq2)' , zeros(nbloq2) ];
error    = zeros(niter+1,1);
residual = zeros(niter+1,1);
regulari = zeros(niter+1,1);

%% Dual problems
% First problem
dirichlet1d = [4,1,0;4,2,0;
               2,1,0;2,2,0;
               1,1,0;1,2,0];
[K1d,C1d,nbloq1d] = Krig (nodes,elements,E,nu,order,boundary,dirichlet1d);
M1d = [ M0 , zeros(2*nnodes,nbloq1d) ; zeros(2*nnodes,nbloq1d)' , zeros(nbloq1d) ];
% Second problem
dirichlet2d = [4,1,0;4,2,0];
[K2d,C2d,nbloq2d] = Krig (nodes,elements,E,nu,order,boundary,dirichlet2d);
M2d = [ M0 , zeros(2*nnodes,nbloq2d) ; zeros(2*nnodes,nbloq2d)' , zeros(nbloq2d) ];
%% Anti-cancellation trick
index = [2*b2node3-1; 2*b2node3];
K1(index,index) = 0;
K2(index,index) = 0;
K1d(index,index) = 0;
K2d(index,index) = 0;

%% Schur operators
% K1s = penalise(  Kinter, [1,2,4], boundary, nnodes, E*1e9 );
% K2s = penalise(  Kinter, [4], boundary, nnodes, E*1e9 );
% %
% f11 = 1e9*E*keepField( uref, 1, boundary ); f12 = 1e9*E*keepField( uref, 2, boundary );
% f1 = assembleDirichlet( [f11,f12] );
% %f1 = f11+f12;
% First problem
% dirichlet1s = [4,1,0;4,2,0;
%               2,1,0;2,2,0;
%               1,1,0;1,2,0];
% [K1s,C1s,nbloq1s,node2c1s,c2node1s] =...
%     Krig (nodes,elements,E,nu,order,boundaryp1,dirichlet1s);
% % Second problem
% dirichlet2s = [4,1,0;4,2,0];
% [K2s,C2s,nbloq2s,node2c2s,c2node2s] =...
%     Krig (nodes,elements,E,nu,order,boundaryp2,dirichlet2s);
% %
% f11 = dirichletRhs2( urefb, 1, c2node1s, boundary, nnodes );
% f12 = dirichletRhs2( urefb, 2, c2node1s, boundary, nnodes );
% f1 = assembleDirichlet( [f11,f12] );
% f2 = loading(nbloq2s,nodes,boundary,neumann2);
% 
% % Schur complement on the redondant boundary :
% % dirichletR = [4,1,0;4,2,0;
% %               3,1,0;3,2,0];
% % [KR,CR,nbloqR,node2cR,c2nodeR] =...
% %     Krig (nodes,elements,E,nu,order,boundaryp1,dirichletR);
% % fR = zeros(2*nnodes+nbloqR,1);
% %
% Krr = Kinter([2*b2node12-1;2*b2node12], [2*b2node12-1;2*b2node12]);
% %
% Kjj = Kinter;
% Kjj([2*b2node12-1;2*b2node12;2*b2node3-1;...
%     2*b2node3;2*b2node4-1;2*b2node4],:) = [];
% Kjj(:,[2*b2node12-1;2*b2node12;2*b2node3-1;...
%     2*b2node3;2*b2node4-1;2*b2node4]) = [];
% %
% Krj = Kinter( [2*b2node12-1;2*b2node12], : );
% Krj( :, [2*b2node12-1;2*b2node12;2*b2node3-1;...
%     2*b2node3;2*b2node4-1;2*b2node4] ) = [];
% %
% Kgj = Kinter( [2*b2node3-1;2*b2node3], : );
% Kgj( :, [2*b2node12-1;2*b2node12;2*b2node3-1;...
%     2*b2node3;2*b2node4-1;2*b2node4] ) = [];
% %
% Ikjj = inv(Kjj);
% %
% SR = full(Krr - Krj*Ikjj*Krj');
% %
% K1k = Kgj*Ikjj*Krj'*inv(SR)*Krj*Ikjj*Kgj';
% %
% [usvdr,svdr] = svd(SR);
% [usvdjj,svdjj] = svd(full(Kjj));
% [usvdrjr,svdrjr] = svd(full(Krj'*SR*Krj));
% [usvdgj,svdgj] = svd(full(Kgj));
% [usvdjjr,svdjjr] = svd(full(Ikjj*Krj'*SR*Krj*Ikjj));
% [usvdrj,svdrj] = svd(full(Krj));
% [usvdk1,svdk1] = svd(K1k);
% %
% % plot(log10(eig(K1k)));
% % figure
% %
% hold on
% plot(log10(diag(svdr)),'Color','red');
% plot(log10(diag(svdjj)),'Color','green');
% plot(log10(diag(svdgj)),'Color','yellow');
% plot(log10(diag(svdrj)),'Color','black');
% plot(log10(diag(svdk1)),'Color','blue');
% legend('svd(Sr)','svd(Kjj)','svd(Kgj)','svd(Krj)','svd(Stot = Sd-Sn)')
% figure
% %
% hold on
% plot(log10(diag(svdr)),'Color','red');
% plot(log10(diag(svdrjr)),'Color','green');
% plot(log10(diag(svdjjr)),'Color','black');
% plot(log10(diag(svdk1)),'Color','blue');
% legend('svd(Sr)','svd(KrjT*SR*Krj)','svd(Ikjj*KrjT*SR*Krj*Ikjj)',...
%     'svd(Stot = Sd-Sn)')
% figure
% %
% for i=11:18
%     u1xy = reshape(usvdk1(:,i),2,[])';
%     hold on
%     plot(u1xy(:,1),'Color','red')
%     plot(u1xy(:,2),'Color','blue')
%     legend(['VP ',num2str(i),' (x)'], ['VP ',num2str(i),' (y)'])
%     figure
% end
%plotGMSH({usvdk1(:,end),'dernier_VP'}, elements, nodes, 'VP');
%%K1k = Kgj*Krj'*SR*Krj*Kgj';
% nz = 0;  % Debug stuff
% for i=1:size(K1k,1)
%     if norm(K1k(:,i))==0
%         nz = nz+1;
%     end
% end
% nz
%
% [ S1, b1, map1 ] = schurComp2( K1s, f1(1:2*nnodes,1), b2node3 );
% [ S2, b2, map2 ] = schurComp2( K2s, f2(1:2*nnodes,1), b2node3 );
% % K1s(index,index) = 0;
% % K2s(index,index) = 0;

% [ S1, b1, map1 ] = schurCompL( K1s, f1, b2node3, nbloq1s, c2node1s );
% [ S2, b2, map2 ] = schurCompL( K2s, f2, b2node3, nbloq2s, c2node2s );
% Stot = S1-S2;
% btot = b1-b2;

% sv1  = real(eig(S1));
% sv2  = real(eig(S2));
% sigv = real(eig(S1-S2));
% sigp = real(eig(S1+S2));
% sv1  = svd(S1);
% sv2  = svd(S2);
% sigv = svd(S1-S2);
% sigp = svd(S1+S2);
% cond(S1-S2)
% hold on;
% plot(eig(S1-S2),'Color','red');
% plot(eig(K1k),'Color','blue');
% % plot(log10(sigp),'Color','blue');
% % plot(log10(sv1),'Color','black');
% % plot(log10(sv2),'Color','black');
% figure

% ucomp = pinv(S1-S2)*(b2-b1);
% hold on;
% urb = reshape(ucomp,2,[])';
% plot(urb(:,1), 'Color', 'red');
% urb = reshape(uref(map1,1),2,[])';
% plot(urb(:,1), 'Color', 'blue');
% figure;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ORTHODIR for the problem : (S10-S20) x = S2-S1
Itere = zeros( 2*nnodes, 1 );
p     = zeros( 2*nnodes, niter+1 );
q     = zeros( 2*nnodes, niter+1 );
Delta = zeros( niter+1, 1 );
Res   = zeros( 2*nnodes, niter+1 );
Zed   = zeros( 2*nnodes, niter+1 );

%% Perform A x0 :
% Solve 1
f1 = dirichletRhs2( Itere, 3, c2node1, boundaryp1, nnodes );
uin1 = (K1-omega*M1)\f1;
lagr1 = uin1(2*nnodes+1:end,1);
lamb1 = lagr2forces2( lagr1, c2node1, 3, boundaryp1, nnodes );
% Solve 2
f2 = dirichletRhs2( Itere, 3, c2node2, boundaryp2, nnodes );
uin2 = (K2-omega*M2)\f2;
lagr2 = uin2(2*nnodes+1:end,1);
lamb2 = lagr2forces2( lagr2, c2node2, 3, boundaryp2, nnodes );
% Regularization term
Nu = regul(Itere, nodes, boundary, 3);
%
Axz = mu*Nu+lamb1-lamb2;
%%%%
%% Compute Rhs :
% Solve 1
f11 = dirichletRhs2( urefb, 1, c2node1, boundary, nnodes );
f12 = dirichletRhs2( urefb, 2, c2node1, boundary, nnodes );
f1 = assembleDirichlet( [f11,f12] );
uin1 = (K1-omega*M1)\f1;
lagr1 = uin1(2*nnodes+1:end,1);
lamb1 = lagr2forces2( lagr1, c2node1, 3, boundaryp1, nnodes );
% Solve 2
f2 = loading(nbloq2,nodes,boundary,neumann2);
uin2 = (K2-omega*M2)\f2;
lagr2 = uin2(2*nnodes+1:end,1);
lamb2 = lagr2forces2( lagr2, c2node2, 3, boundaryp2, nnodes );
%
b = -lamb1+lamb2;
% hold on;
% urb = reshape(b(map1,1),2,[])';
% plot(urb(:,1), 'Color', 'red');
% urb = reshape(btot,2,[])';
% plot(urb(:,1), 'Color', 'blue');
% figure;
%%%%
% %% Compute Residual (alternative way) :
% % Solve 1
% f11 = dirichletRhs(uref, 1, C1, boundary);
% f12 = dirichletRhs(uref, 2, C1, boundary);
% f13 = dirichletRhs(Itere, 3, C1, boundary);
% f1 = f11 + f12 + f13;
% uin1 = K1\f1;
% lagr1 = uin1(2*nnodes+1:end,1);
% lamb1 = lagr2forces( lagr1, C1, 3, boundary );
% % Solve 2
% f21 = loading(nbloq2,nodes,boundary,neumann2);
% f22 = dirichletRhs(Itere, 3, C2, boundary);
% f2 = f21+f22;
% uin2 = K2\f2;
% lagr2 = uin2(2*nnodes+1:end,1);
% lamb2 = lagr2forces( lagr2, C2, 3, boundary );
% %
% Res(:,1) = -lamb1+lamb2;
%%
Res(:,1) = b - Axz;
% urb = reshape(Res(map1,1),2,[])';
% plot(urb(:,1));
% figure;

if precond == 1
    % Solve 1
    f1 = [Res(:,1)/2; zeros(nbloq1d,1)];
    uin1 = (K1d-omega*M1d)\f1;
    u1i = uin1(1:2*nnodes,1);
    u1 = keepField( u1i, 3, boundaryp1 );
%    % Solve 2
%    f2 = [Res(:,1)/2; zeros(nbloq2d,1)];
%    uin2 = (K2d-omega*M2d)\f2;
%    u2i = uin2(1:2*nnodes,1);
%    u2 = keepField( u2i, 3, boundaryp2 );
    %
    Zed(:,1) = u1/2;%-u2/2;
else
    Zed(:,1) = Res(:,1);
end

p(:,1) = Zed(:,1);

%residual(1) = sqrt( myps( Zed(:,1), Zed(:,1), Kinter, boundaryp1, M, nodes ) );
residual(1) = sqrt( myps( Res(:,1), Res(:,1), Kinter, boundary, M, nodes ) );
error(1)    = sqrt( myps( Itere - uref, Itere - uref, Kinter, boundary, M, nodes )...
                 / myps( uref, uref, Kinter, boundary, M, nodes ) );
regulari(1) = sqrt( Itere'*regul(Itere, nodes, boundary, 3) );

%% Perform Q1 = A P1 :
% Solve 1
f1 = dirichletRhs2(p(:,1), 3, c2node1, boundaryp1, nnodes);
uin1 = (K1-omega*M1)\f1;
lagr1 = uin1(2*nnodes+1:end,1);
lamb1 = lagr2forces2( lagr1, c2node1, 3, boundaryp1, nnodes );
% Solve 2
f2 = dirichletRhs2(p(:,1), 3, c2node2, boundaryp2, nnodes);
uin2 = (K2-omega*M2)\f2;
lagr2 = uin2(2*nnodes+1:end,1);
lamb2 = lagr2forces2( lagr2, c2node2, 3, boundaryp2, nnodes );
% Regularization term
Nu = regul(p(:,1), nodes, boundary, 3);
%
q(:,1) = mu*Nu+lamb1-lamb2;
%%%%
%%
for iter = 1:niter
    
    Delta(iter,1) = myps( q(:,iter), q(:,iter), Kinter, boundary, M, nodes ); %q(:,iter)'*q(:,iter);
    gammai        = myps( q(:,iter), Res(:,iter), Kinter, boundary, M, nodes ); %q(:,iter)'*Res;
    alphai        = gammai/Delta(iter,1);
    
    Itere         = Itere + p(:,iter)*alphai;
    Res(:,iter+1) = Res(:,iter) - q(:,iter)*alphai;
    
    if precond == 1
        % Solve 1
        f1 = [Res(:,iter+1)/2; zeros(nbloq1d,1)];
        uin1 = (K1d-omega*M1d)\f1;
        u1i = uin1(1:2*nnodes,1);
        u1 = keepField( u1i, 3, boundaryp1 );
%        % Solve 2
%        f2 = [Res(:,iter+1)/2; zeros(nbloq2d,1)];
%        uin2 = (K2d-omega*M2d)\f2;
%        u2i = uin2(1:2*nnodes,1);
%        u2 = keepField( u2i, 3, boundaryp2 );
        %
        Zed(:,iter+1) = u1/2;%-u2/2;
    else
        Zed(:,iter+1) = Res(:,iter+1);
    end
    
    %residual(iter+1) = sqrt( myps( Zed(:,iter+1), Zed(:,iter+1), Kinter, boundaryp1, M, nodes ) );
    residual(iter+1) = sqrt( myps( Res(:,iter+1), Res(:,iter+1), Kinter, boundary, M, nodes ) );
    error(iter+1)    = sqrt( myps( Itere - uref, Itere - uref, Kinter, boundary, M, nodes )...
                     / myps( uref, uref, Kinter, boundary, M, nodes ) );
    regulari(iter+1) = sqrt( Itere'*regul(Itere, nodes, boundary, 3) );

    %% Perform Ari = A*Res
    % Solve 1
    rhs1 = Zed(:,iter+1);
    f1 = dirichletRhs2(rhs1, 3, c2node1, boundaryp1, nnodes);
    uin1 = (K1-omega*M1)\f1;
    lagr1 = uin1(2*nnodes+1:end,1);
    lamb1 = lagr2forces2( lagr1, c2node1, 3, boundaryp1, nnodes );
    % Solve 2
    rhs2 = Zed(:,iter+1);
    f2 = dirichletRhs2(rhs2, 3, c2node2, boundaryp2, nnodes);
    uin2 = (K2-omega*M2)\f2;
    lagr2 = uin2(2*nnodes+1:end,1);
    lamb2 = lagr2forces2( lagr2, c2node2, 3, boundaryp2, nnodes );
    % Regularization term
    Nu = regul(Zed(:,iter+1), nodes, boundaryp2, 3);
    %
    Ari = mu*Nu+lamb1-lamb2;
    
    %% Orthogonalization
    p(:,iter+1) = Zed(:,iter+1);
    q(:,iter+1) = Ari;
    
    for jter=1:iter
        phiij  = myps( q(:,jter), Ari, Kinter, boundary, M, nodes ); %q(:,jter)'*Ari;
        betaij = phiij/Delta(jter,1);
        p(:,iter+1) = p(:,iter+1) - p(:,jter) * betaij;
        q(:,iter+1) = q(:,iter+1) - q(:,jter) * betaij;
    end
    
end

hold on
% plot(error,'Color','blue')
% plot(residual,'Color','red')
plot(log10(error),'Color','blue')
plot(log10(residual),'Color','red')
legend('error (log)','residual (log)')
%figure;
% L-curve :
%loglog(residual,regulari);
%figure
%%%%
%% Final problem : compute u
% DN problem
dirichlet = [4,1,0;4,2,0;
             3,1,0;3,2,0;
             1,1,0;1,2,0;
             2,1,0;2,2,0];
neumann   = [];
[K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
M = [ M0 , zeros(2*nnodes,nbloq) ; zeros(2*nnodes,nbloq)' , zeros(nbloq) ];
fdir1 = dirichletRhs(urefb, 1, C, boundary);
fdir2 = dirichletRhs(urefb, 2, C, boundary);
fdir3 = dirichletRhs(Itere, 3, C, boundary);
usoli = (K-omega*M) \ assembleDirichlet( [fdir1+fdir3,fdir2] );%( max(fdir1, max(fdir2, fdir3)) );

usol = usoli(1:2*nnodes,1);
fsol = (Kinter-omega*M0)*usol;

%plot(fsol(2*b2node3))

total_error = norm(usol-uref)/norm(uref);
% Compute stress :
sigma = stress(usol,E,nu,nodes,elements,order,1,ntoelem);
plotGMSH({usol,'U_vect';sigma,'stress'}, elements, nodes, 'output/solution');