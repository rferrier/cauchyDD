% 06/06/2016
% Algo Steklov-Poincar� primal avec Gradient Conjugu�

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 70000;  % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 1;      % N.mm-1 : Loading on the plate
niter   = 10;
precond = 1;      % 1 : Use a dual precond, 2 : use regul precond
mu      = 0.;    % Regularization parameter
br      = 0.0;     % noise

% Boundary conditions
% first index  : index of the boundary
% second index : 1=x, 2=y
% third        : value
% [0,1,value] marks a dirichlet regularization therm on x
dirichlet = [1,1,0; 1,2,0 ;
             3,1,0; 3,2,0 ];
neumann   = [2,1,fscalar;
             4,1,fscalar];

% Import the mesh
[ nodes,elements,ntoelem,boundary,order ] = readmesh( 'meshes/platee.msh' );
nnodes = size(nodes,1);

% Import the overkill mesh
[ nodeso,elementso,ntoelemo,boundaryo,ordero ] = readmesh( 'meshes/plateer.msh' );
nnodeso = size(nodeso,1);

% find the nodes in the corners and suppress the element :
xmax = max(nodes(:,1));
xmin = min(nodes(:,1));
ymax = max(nodes(:,2));
ymin = min(nodes(:,2));
no1  = findNode(xmin, ymin, nodes, 1e-5);
no2  = findNode(xmax, ymin, nodes, 1e-5);
no3  = findNode(xmax, ymax, nodes, 1e-5);
no4  = findNode(xmin, ymax, nodes, 1e-5);

boundaryp1 = suppressBound( boundary, no2, 2 );
boundaryp1 = suppressBound( boundaryp1, no3, 2 );
boundaryp1 = suppressBound( boundaryp1, no1, 4 );
boundaryp1 = suppressBound( boundaryp1, no4, 4 );
boundaryp2 = boundaryp1;

[~, b2node2] = mapBound(2, boundaryp1, nnodes);
indexa = [2*b2node2-1; 2*b2node2];

% Then, build the stiffness matrix :
[K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
[Ko,Co,nbloqo] = Krig (nodeso,elementso,E,nu,ordero,boundaryo,dirichlet);
Kinter = K(1:2*nnodes, 1:2*nnodes);
M      = mass_mat(nodes, elements);
[node2b4, b2node4] = mapBound(4, boundaryp1, nnodes);
[node2b3, b2node3] = mapBound(3, boundaryp1, nnodes);
[node2b1, b2node1] = mapBound(1, boundaryp1, nnodes);
[node2b2, b2node2] = mapBound(2, boundaryp1, nnodes);
index = [2*b2node2-1;2*b2node2];
% The right hand side :
f = loading(nbloq,nodes,boundary,neumann);
fo = loading(nbloqo,nodeso,boundaryo,neumann);

% Solve the problem :
uin = K\f;
uino = Ko\fo; urefo = uino(1:2*nnodeso,1);

% Extract displacement and Lagrange multiplicators :
uref = uin(1:2*nnodes,1);
lagr = uin(2*nnodes+1:end,1);

urefb = ( 1 + br*randn(2*nnodes,1) ) .* uref;

sigma = stress(uref,E,nu,nodes,elements,order,1,ntoelem);
plotGMSH({uref,'Vect_U';sigma,'stress'}, elements, nodes, 'reference');

%% Estimate discretization error
%ureftoo = passMesh2D (nodes, elements, nodeso, elementso, uref);
%error_disc = elasticEnergy (ureftoo-urefo,[0,E,nu],nodeso,elementso,ordero); %= 2.0495e-05
%plotGMSH({ureftoo,'U_vect'}, elementso, nodeso, 'output/projected');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition of the operators

% First problem
dirichlet1 = [4,1,0;4,2,0;
              3,1,0;3,2,0;
              2,1,0;2,2,0;
              1,1,0;1,2,0];

[K1,C1,nbloq1,node2c1,c2node1] = Krig (nodes,elements,E,nu,order,boundaryp1,dirichlet1);

% Second problem
dirichlet2 = [1,1,0;1,2,0;
              2,1,0;2,2,0;
              3,1,0;3,2,0];
neumann2   = [4,1,fscalar];
neumann0   = [];
[K2,C2,nbloq2,node2c2,c2node2] = Krig (nodes,elements,E,nu,order,boundaryp2,dirichlet2);

error    = zeros(niter+1,1);
residual = zeros(niter+1,1);
regulari = zeros(niter+1,1);

%% Dual problems
% First problem
dirichlet1d = [4,1,0;4,2,0;
               3,1,0;3,2,0;
               1,1,0;1,2,0];
[K1d,C1d,nbloq1d,node2c1d,c2node1d] = Krig (nodes,elements,E,nu,order,boundary,dirichlet1d);

% Second problem
dirichlet2d = [1,1,0;1,2,0;
               3,1,0;3,2,0];
[K2d,C2d,nbloq2d,node2c1d,c2node1d] = Krig (nodes,elements,E,nu,order,boundary,dirichlet2d);

%% Anti-cancellation trick
% K1(indexa,indexa) = 0;
% K2(indexa,indexa) = 0;
% K1d(indexa,indexa) = 0;
% K2d(indexa,indexa) = 0;

%% Schur operators
% First problem
dirichlet1s = [3,1,0;3,2,0;
               4,1,0;4,2,0;
               1,1,0;1,2,0];
[K1s,C1s,nbloq1s,node2c1s,c2node1s] =...
    Krig (nodes,elements,E,nu,order,boundaryp1,dirichlet1s);
% Second problem
dirichlet2s = [1,1,0;1,2,0;
               3,1,0;3,2,0];
[K2s,C2s,nbloq2s,node2c2s,c2node2s] =...
    Krig (nodes,elements,E,nu,order,boundaryp2,dirichlet2s);
%
f11 = dirichletRhs2( urefb, 1, c2node1s, boundary, nnodes );
f12 = dirichletRhs2( urefb, 2, c2node1s, boundary, nnodes );
f1 = assembleDirichlet( [f11,f12] );
f2 = loading(nbloq2s,nodes,boundary,neumann2);
% 
%Krr = Kinter([2*b2node4-1;2*b2node4], [2*b2node4-1;2*b2node4]);
%%
%Kjj = Kinter;
%Kjj([2*b2node1-1;2*b2node1;2*b2node3-1;2*b2node3;2*b2node2-1;2*b2node2;...
%    2*b2node4-1;2*b2node4],:) = [];
%Kjj(:,[2*b2node1-1;2*b2node1;2*b2node3-1;2*b2node3;2*b2node2-1;2*b2node2;...
%    2*b2node4-1;2*b2node4]) = [];
%%
%Krj = Kinter( [2*b2node4-1;2*b2node4], : );
%Krj( :, [2*b2node1-1;2*b2node1;2*b2node3-1;2*b2node3;2*b2node2-1;2*b2node2;...
%    2*b2node4-1;2*b2node4] ) = [];
%%
%Kgj = Kinter( [2*b2node2-1;2*b2node2], : );
%Kgj( :, [2*b2node1-1;2*b2node1;2*b2node3-1;2*b2node3;2*b2node2-1;2*b2node2;...
%    2*b2node4-1;2*b2node4] ) = [];
%%
%Ikjj = inv(Kjj);
%%
%SR = full(Krr - Krj*Ikjj*Krj');
%%
%Korth = Kgj*Ikjj*Krj';
%%
%K1k = Kgj*Ikjj*Krj'*inv(SR)*Krj*Ikjj*Kgj';
%%
%% solution of the normal equation
%%
%ur     = uref([2*b2node4-1;2*b2node4],:);
%fr     = f([2*b2node4-1;2*b2node4],:);
%btilde = SR*ur + fr;
%um     = Korth\btilde;
%hold on;
%plot(um);
%plot( uref([2*b2node2-1;2*b2node2],:), 'Color', 'red' );
%legend('solution','reference')
%figure;
% [usvdr,svdr] = svd(SR);
% [usvdjj,svdjj] = svd(full(Kjj));
% % [usvdrjr,svdrjr] = svd(full(Krj'*SR*Krj));
% [usvdgj,svdgj] = svd(full(Kgj));
% % [usvdjjr,svdjjr] = svd(full(Ikjj*Krj'*SR*Krj*Ikjj));
% [usvdrj,svdrj] = svd(full(Krj));
% [usvdk1,svdk1] = svd(K1k);
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
%
% [ S1, b1, map1 ] = schurCompL( K1s, f1, b2node2, nbloq1s, c2node1s );
% [ S2, b2, map2 ] = schurCompL( K2s, f2, b2node2, nbloq2s, c2node2s );
% Stot = S1-S2;
% btot = b1-b2;

% hold on;
% plot(log10(svd(S1-S2)),'Color','red');
%plot(log10(svd(K1k)),'Color','blue');
%figure
%plot(log10(svd(Korth)),'Color','blue');
%figure
% cond(K1k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Conjugate Gradient for the problem : (S10-S20) x = S2-S1
Itere = zeros( 2*nnodes, 1 );
usolC = zeros( 2*nnodes, 1 ); % Reconstructed total solution
d     = zeros( 2*nnodes, niter+1 );
Ad    = zeros( 2*nnodes, niter+1 );
Res   = zeros( 2*nnodes, niter+1 );
Zed   = zeros( 2*nnodes, niter+1 );
H12   = H1demi(size(Res,1), nodes, boundary, 2 );
Mgr   = Mgrad(size(Res,1), nodes, boundary, 2 );

%% Error estimation quantities
Eerror    = zeros(niter+1,1); % "real error" with overkill
binferror = zeros(niter+1,1);
lambr = zeros( 2*nnodes, 1 ); % to be incremented
ulamr1 = zeros( 2*nnodes, 1 ); % to be incremented = Z^(-1)*lambr
tr = loading(0,nodes,boundary,neumann2); % loading on the redondant boundary
ur = keepField( uref(1:2*nnodes), 4, boundaryp1 ); % U on the redondant boundary

%% Perform A x0 :
% Solve 1
f1 = dirichletRhs2( Itere, 2, c2node1, boundaryp1, nnodes );
uin1 = K1\f1;
lagr1 = uin1(2*nnodes+1:end,1);
lamb1 = lagr2forces( lagr1, C1, 2, boundaryp1 );
% Solve 2
f2 = dirichletRhs2( Itere, 2, c2node2, boundaryp2, nnodes );
uin2 = K2\f2;
lagr2 = uin2(2*nnodes+1:end,1);
lamb2 = lagr2forces( lagr2, C2, 2, boundaryp2 );
% Regularization term
Nu = regul(Itere, nodes, boundary, 2);
%
Axz = mu*Nu+lamb1-lamb2;

% Error estimation assets
lambr = lambr + lagr2forces2( uin1(2*nnodes+1:end,1), c2node1, 4, boundaryp1, nnodes );
usolC = usolC + uin1(1:2*nnodes);
%ulamr = ulamr + keepField( uin1(1:2*nnodes), 4, boundaryp1 );

%%%%
%% Compute Rhs :
% Solve 1
f1 = dirichletRhs2( urefb, 4, c2node1, boundary, nnodes );
uin1 = K1\f1;
lagr1 = uin1(2*nnodes+1:end,1);
lamb1 = lagr2forces2( lagr1, c2node1, 2, boundaryp1, nnodes );
% Solve 2
f2 = loading(nbloq2,nodes,boundary,neumann2);
uin2 = K2\f2;
lagr2 = uin2(2*nnodes+1:end,1);
lamb2 = lagr2forces2( lagr2, c2node2, 2, boundaryp2, nnodes );
%
b = -lamb1+lamb2;
utr = keepField( uin2(1:2*nnodes), 4, boundaryp2 ); % error estimator assets
lambr = lambr + lagr2forces2( uin1(2*nnodes+1:end,1), c2node1, 4, boundaryp1, nnodes );
usolC = usolC + uin1(1:2*nnodes);
%%
Res(:,1) = b - Axz;

%% Debug : compute Z^(-1)*lambr
%f1 = dirichletRhs2( Itere, 2, c2node1, boundaryp1, nnodes );
%f2 = dirichletRhs2( urefb, 4, c2node1, boundary, nnodes );
%uin1 = K1\(f1+f2);
%lambr1 = lagr2forces2( uin1(2*nnodes+1:end,1), c2node1, 4, boundaryp1, nnodes );
ulamr = K2\[lambr;zeros(nbloq2,1)]; ulamr = ulamr(1:2*nnodes);
utr   = K2\[tr;zeros(nbloq2,1)];    utr = utr(1:2*nnodes);

binferror(1) = (tr-lambr)'*(utr-ulamr); % Error inferior born

if precond == 1
    % Solve 1
    f1 = [Res(:,1)/2; zeros(nbloq1d,1)];
    uin1 = K1d\f1;
    u1i = uin1(1:2*nnodes,1);
    u1 = keepField( u1i, 2, boundaryp1 );
    % Solve 2
    f2 = [Res(:,1)/2; zeros(nbloq2d,1)];
    uin2 = K2d\f2;
    u2i = uin2(1:2*nnodes,1);
    u2 = keepField( u2i, 2, boundaryp2 );
    %
    Zed(:,1) = u1/2-u2/2;
elseif precond == 2
    Zed(index,1) = 1/E*H12(index,index)\Res(index,1);
elseif precond == 3
    Zed(index,1) = 1/E*Mgr(index,index)\Res(index,1);
else
    Zed(:,1) = Res(:,1);
end

d(:,1) = Zed(:,1);

residual(1) = norm(Res( indexa,1));
error(1)    = norm(Itere(indexa) - uref(indexa)) / norm(uref(indexa));
regulari(1) = sqrt( Itere'*regul(Itere, nodes, boundary, 2) );

% Energetic error
%usoltoo = passMesh2D (nodes, elements, nodeso, elementso, usolC);
%Eerror(1) = elasticEnergy (urefo-usoltoo,[0,E,nu],nodes,elements,order);

%%
for iter = 1:niter
    %% Optimal step
    
    % Solve 1
    rhs1 = d(:,iter);
    f1 = dirichletRhs2(rhs1, 2, c2node1, boundaryp1, nnodes);
    uin1 = K1\f1;
    lagr1 = uin1(2*nnodes+1:end,1);
    lamb1 = lagr2forces2( lagr1, c2node1, 2, boundaryp1, nnodes );
    % Solve 2
    rhs2 = d(:,iter);
    f2 = dirichletRhs2(rhs2, 2, c2node2, boundaryp2, nnodes);
    uin2 = K2\f2;
    lagr2 = uin2(2*nnodes+1:end,1);
    lamb2 = lagr2forces2( lagr2, c2node2, 2, boundaryp2, nnodes );
    % Regularization term
    Nu = regul(d(:,iter), nodes, boundaryp2, 2);
    %
    Ad(:,iter) = mu*Nu+lamb1-lamb2;
    
    % Error estimation assets
    lambri = lagr2forces2( uin1(2*nnodes+1:end,1), c2node1, 4, boundaryp1, nnodes );
%    ulamri = keepField( uin1(1:2*nnodes), 4, boundaryp1 );
    
    den = (d(indexa,iter)'*Ad(indexa,iter));
    d(:,iter) = d(:,iter)/sqrt(den); Ad(:,iter) = Ad(:,iter)/sqrt(den);
    num = Res(indexa,iter)'*d(indexa,iter);
    Itere         = Itere + d(:,iter)*num;
    Res(:,iter+1) = Res(:,iter) - Ad(:,iter)*num;
    
    % Increment error estimation
    lambr         = lambr + lambri*num/sqrt(den);
    usolC         = usolC + uin1(1:2*nnodes)*num/sqrt(den); % Total field
%    ulamr         = ulamr + ulamri*num/sqrt(den);

    %% Debug : compute Z^(-1)*lambr
%    f1 = dirichletRhs2( Itere, 2, c2node1, boundaryp1, nnodes );
%    f2 = dirichletRhs2( urefb, 4, c2node1, boundary, nnodes );
%    uin1 = K1\(f1+f2);
%    lambr1 = lagr2forces2( uin1(2*nnodes+1:end,1), c2node1, 4, boundaryp1, nnodes );
    
    ulamr = K2\[lambr;zeros(nbloq2,1)]; ulamr = ulamr(1:2*nnodes);
    binferror(iter+1) = (tr-lambr)'*(utr-ulamr);
%    usoltoo = passMesh2D (nodes, elements, nodeso, elementso, usolC);
%    Eerror(iter+1) = elasticEnergy (urefo-usoltoo,[0,E,nu],nodes,elements,order);
    
    residual(iter+1) = norm(Res(indexa,iter+1));
    error(iter+1)    = norm(Itere(indexa) - uref(indexa)) / norm(uref(indexa));
    regulari(iter+1) = sqrt( Itere'*regul(Itere, nodes, boundary, 2) );
    
    if precond == 1
        % Solve 1
        f1 = [Res(:,iter+1)/2; zeros(nbloq1d,1)];
        uin1 = K1d\f1;
        u1i = uin1(1:2*nnodes,1);
        u1 = keepField( u1i, 2, boundaryp1 );
        % Solve 2
        f2 = [Res(:,iter+1)/2; zeros(nbloq2d,1)];
        uin2 = K2d\f2;
        u2i = uin2(1:2*nnodes,1);
        u2 = keepField( u2i, 2, boundaryp2 );
        %
        Zed(:,iter+1) = u1/2-u2/2;
    elseif precond == 2
        Zed(index,iter+1) = 1/E*H12(index,index)\Res(index,iter+1);
    elseif precond == 3
        Zed(index,iter+1) = 1/E*Mgr(index,index)\Res(index,iter+1);
    else
        Zed(:,iter+1) = Res(:,iter+1);
    end
    
    %% Orthogonalization
    d(:,iter+1) = Zed(:,iter+1);
    
    for jter=1:iter
        betaij = ( Zed(indexa,iter+1)'*Ad(indexa,jter) );%/...
            %( d(indexa,jter)'*Ad(indexa,jter) );
        d(:,iter+1) = d(:,iter+1) - d(:,jter) * betaij;
    end

%     betaij = ( Zed(indexa,iter+1)'*Ad(indexa,iter) )/...
%         ( d(indexa,iter)'*Ad(indexa,iter) );
%     d(:,iter+1) = d(:,iter+1) - d(:,iter) * betaij;
    
end

hold on
plot(log10(error),'Color','blue')
plot(log10(residual),'Color','red')
legend('error (log)','residual (log)')
figure;
% L-curve :
loglog(residual(2:end),regulari(2:end));
figure
%%%%
%% Final problem : compute u
% DN problem
dirichlet = [4,1,0;4,2,0;
             3,1,0;3,2,0;
             1,1,0;1,2,0;
             2,1,0;2,2,0];
neumann   = [];
[K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
fdir2 = dirichletRhs(Itere, 2, C, boundary);
fdir4 = dirichletRhs(urefb, 4, C, boundary);
usoli = K \ (fdir4 + fdir2);

usol = usoli(1:2*nnodes,1);
fsol = Kinter*usol;

plot(fsol(2*b2node2-1),'Color','red');

total_error = norm(usol-uref)/norm(uref);
% Energetic error
error_end = elasticEnergy (uref-usol,[0,E,nu],nodes,elements,order);

%usolo = passMesh2D (nodes, elements, nodeso, elementso, usol);
%error_o = elasticEnergy (usolo-urefo,[0,E,nu],nodeso,elementso,ordero); %= 2.5275e-05
% Compute stress :
sigma = stress(usol,E,nu,nodes,elements,order,1,ntoelem);
plotGMSH({usol,'U_vect';sigma,'stress'}, elements, nodes, 'solution');