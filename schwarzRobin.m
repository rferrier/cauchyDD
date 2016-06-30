%23/03/2016
%Code Schwarz Robin 2D contraintes planes

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 70000;      % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 1;      % N.mm-1 : Loading on the plate
niter   = 15;
robin1  = .6*E;%.6*E; lim stab : 7*E
robin2  = .6*E;%.6*E;            1/7*E
rank    = 9;       % Half of the size of the approximation of Si 
                   % (if 0, don't use it)

dirichlea = [4,1,0;4,2,0;
             1,1,0;1,2,0;
             2,1,0;2,2,0];
dirichleb = [4,1,0;4,2,0;
             3,1,0;3,2,0;
             2,1,0;2,2,0];
neumann   = [];

% First, import the mesh
[ nodea,elementa,ntoelea,boundara,ordea ] = readmesh( 'meshes/plate.msh' );
nodeb = translateMesh( nodea, 0, 10 );
nnodea = size(nodea,1);

% Then, build the stiffness matrix :
[Ka,Ca,nbloqa,node2ca,c2nodea] =...
    Krig (nodea,elementa,E,nu,ordea,boundara,dirichlea);
[Kb,Cb,nbloqb,node2cb,c2nodeb] =...
    Krig (nodeb,elementa,E,nu,ordea,boundara,dirichleb);
%Kab = penalEq( Ka, Kb, E*1e9, nodea, nodeb, 1e-5 );
[ Kab, naj] = lagrEq( Ka, Kb, nodea, nodeb, 1e-5,c2nodea,c2nodeb);

% The right hand side :
fa = volumicLoad( nbloqa, nodea, elementa, 1, fscalar );
fb = volumicLoad( nbloqb, nodeb, elementa, 1, fscalar );

% Solve the problem :
uinab = Kab\[fa;fb;zeros(naj,1)];

% Extract displacement :
u = uinab(1:2*nnodea,1);

% Compute stress :
sigma = stress(u,E,nu,nodea,elementa,ordea,1,ntoelea);

% Output :
% plotNodes(u,elements,nodes); : TODO debug on Matlab r>2013
plotGMSH({u,'U_vect';sigma,'stress'}, elementa, nodea(:,[1,2]), 'ref_field');

%%%%
%% Assembly of the Schur complement

% We work on an half-mesh
[ nodes1,elements,ntoelem,boundary,order ] = readmesh( 'meshes/plate.msh' );
nnodes = size(nodes1,1);
% Second mesh
nodes2 = translateMesh(nodes1, 0, 10); % Move the mesh2
[ b1to2, b2to1 ] = superNodes( nodes1, nodes2, 1e-5 ); % build the connec table

% [ K1s, f1s ] = suppressLine( Kinter, f1in(1:2*nnodes,1), [1,2,4], boundary, nnodes );
% [ K2s, f2s ] = suppressLine( Kinter, f2in(1:2*nnodes,1), [3,2,4], boundary, nnodes );

% K1s = penalise( Kinter, [1,2,4], boundary, nnodes, E*1e9 );
% K2s = penalise( Kinter, [3,2,4], boundary, nnodes, E*1e9 );
% [ S1, b1, map1 ] = schurComp2( K1s, f1in(1:2*nnodes,1), b2node3 );
% [ S2i, b2i, map2 ] = schurComp2( K2s, f2in(1:2*nnodes,1), b2node1 );
% find the nodes in the corners and suppress the elements :
xmax = max(nodes1(:,1));
xmin = min(nodes1(:,1));
ymax = max(nodes1(:,2));
ymin = min(nodes1(:,2));
no1  = findNode(xmin, ymin, nodes1, 1e-5);
no2  = findNode(xmax, ymin, nodes1, 1e-5);
no3  = findNode(xmax, ymax, nodes1, 1e-5);
no4  = findNode(xmin, ymax, nodes1, 1e-5);
boundaryp2 = suppressBound( boundary, [no1;no2], 1 );
boundaryp1 = suppressBound( boundary, [no3;no4], 3 );
% Some useful tables
[node2b1, b2node1] = mapBound(1, boundaryp2, nnodea);
[node2b3, b2node3] = mapBound(3, boundaryp1, nnodea);
%
dirichlet1s = [4,1,0;4,2,0;
               2,1,0;2,2,0;
               1,1,0;1,2,0];
[K1s,C1s,nbloq1s,node2c1s,c2node1s] =...
    Krig (nodes1,elements,E,nu,order,boundaryp1,dirichlet1s);
Kinter = K1s(1:2*nnodes, 1:2*nnodes);
f1ins = volumicLoad( nbloq1s, nodes1, elements, 1, fscalar );
% Second problem
dirichlet2s = [3,1,0;3,2,0;
               2,1,0;2,2,0;
               4,1,0;4,2,0];
[K2s,C2s,nbloq2s,node2c2s,c2node2s] =...
    Krig (nodes2,elements,E,nu,order,boundaryp2,dirichlet2s);
f2ins = volumicLoad( nbloq2s, nodes2, elements, 1, fscalar );
%
[ S1, b1, map1 ]   = schurCompL( K1s, f1ins, b2node3, nbloq1s, c2node1s );
[ S2i, b2i, map2 ] = schurCompL( K2s, f2ins, b2node1, nbloq2s, c2node2s );

ndof = size(S1,1);
% Re-index S2 in order to match S1
% S1 on b1 -> b2node1 -> S1 on n2 -> b2ot1 -> S1 on n1 -> node2b3 -> S1 on 3
i0 = 1:1:ndof/2;   % Indices on b1
i1 = b2node1(i0,1)';% indices on n2
i2 = b2to1(i1,1)';  % indices on n1
i3 = node2b3(i2,1)';% indices on b2
m0 = 1:1:ndof;           % The same with dofs (1node = 2dofs)
m3 = zeros(1,ndof);
m3(1, 2.*i0) = 2.*i3;
m3(1, 2.*i0-1) = 2.*i3-1;
% Actual re-indexation
S2(m0, m0) = S2i(m3, m3);
b2(m0, 1) = b2i(m3, 1);

% Solution of the boundary problem (to validate the schur complement)
% urs = (S1+S2)\(b1+b2);
% ub1sb = reshape(urs,2,[])';
% plot(ub1sb(:,1));
% legend('Solution sur le bord avec complément de Schur x')
% figure;
% plot(ub1sb(:,2));
% legend('Solution sur le bord avec complément de Schur y')
% figure;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Obtain an approxiamtion of the Schur complements

if rank >= 1
    
    % First problem
    dirichlet1s = [ 4,1,0;4,2,0 ; 2,1,0;2,2,0 ; 1,1,0;1,2,0 ; 3,1,0;3,2,0 ];
    [K1s,C1s,nbloq1s,node2c1s,c2node1s] =...
        Krig (nodes1,elements,E,nu,order,boundaryp1,dirichlet1s);
    % Second problem
    dirichlet2s = [ 3,1,0;3,2,0 ; 2,1,0;2,2,0 ; 4,1,0;4,2,0 ; 1,1,0;1,2,0 ];
    [K2s,C2s,nbloq2s,node2c2s,c2node2s] =...
        Krig (nodes2,elements,E,nu,order,boundaryp2,dirichlet2s);
    
    bsize3 = size(b2node3,1);
    index1 = zeros( 2*bsize3 ,1 );
    index1(1:2:2*bsize3-1) = 2*b2node3-1;
    index1(2:2:2*bsize3)   = 2*b2node3;
    
    bsize1 = size(b2node1,1);
    index2 = zeros( 2*bsize1 ,1 );
    index2(1:2:2*bsize1-1) = 2*b2node1-1;
    index2(2:2:2*bsize1)   = 2*b2node1;
    
    xmaxs = max(nodes1(b2node1,1));
    xmins = min(nodes1(b2node1,1));
    
    L     = xmax-xmin;
    Ls    = xmaxs-xmins;
    xmoys = (xmaxs+xmins)/2;
    P1w = zeros( 2*bsize3, 2*rank ); % Projection matrices
    P2w = zeros( 2*bsize1, 2*rank );
    P1f = zeros( 2*bsize3, 2*rank ); % Projection matrices
    P2f = zeros( 2*bsize1, 2*rank );
    
    for i=1:rank
        for j=1:bsize3 % Loop over nodes of the boundary
            x      = nodes1( b2node3(j), 1 ) - xmin;
            P1w(2*j-1,2*i-1) = sin(i*pi*x/L);
            P1w(2*j,2*i)     = sin(i*pi*x/L);
            P1f(2*j-1,2*i-1) = x^(i-1);
            P1f(2*j,2*i)     = x^(i-1);
        end
        for j=1:bsize1 % Loop over nodes of the boundary
            x = nodes2( b2node1(j), 1 ) - xmin;
            P2w(2*j-1,2*i-1) = sin(i*pi*x/L);
            P2w(2*j,2*i)     = sin(i*pi*x/L);
            P2f(2*j-1,2*i-1) = x^(i-1);
            P2f(2*j,2*i)     = x^(i-1);
        end
        
        % normalize
        P1w(:,2*i-1) = P1w(:,2*i-1)/norm(P1w(:,2*i-1));
        P1w(:,2*i) = P1w(:,2*i)/norm(P1w(:,2*i));
        P2w(:,2*i-1) = P2w(:,2*i-1)/norm(P2w(:,2*i-1));
        P2w(:,2*i) = P2w(:,2*i)/norm(P2w(:,2*i));
        
        % Orthonormalize
        for k=1:i-1
            P1f(:,2*i-1) = P1f(:,2*i-1) -...
                P1f(:,2*k-1) * (P1f(:,2*i-1)'*P1f(:,2*k-1))/norm(P1f(:,2*k-1))^2;
            P1f(:,2*i) = P1f(:,2*i) -...
                P1f(:,2*k) * (P1f(:,2*i)'*P1f(:,2*k))/norm(P1f(:,2*k))^2;
            P2f(:,2*i-1) = P2f(:,2*i-1) -...
                P2f(:,2*k-1) * (P2f(:,2*i-1)'*P2f(:,2*k-1))/norm(P2f(:,2*k-1))^2;
            P2f(:,2*i) = P2f(:,2*i) -...
                P2f(:,2*k) * (P2f(:,2*i)'*P2f(:,2*k))/norm(P2f(:,2*k))^2;
        end
        P1f(:,2*i-1) = P1f(:,2*i-1)/norm(P1f(:,2*i-1));
        P1f(:,2*i) = P1f(:,2*i)/norm(P1f(:,2*i));
        P2f(:,2*i-1) = P2f(:,2*i-1)/norm(P2f(:,2*i-1));
        P2f(:,2*i) = P2f(:,2*i)/norm(P2f(:,2*i));
        
    end
    
%     P1f = P1w;
%     P2f = P2w;
    
    %% Solve the 2*2*rank problems
    f1h = zeros( 2*bsize3, 2*rank );
    f2h = zeros( 2*bsize3, 2*rank );
    
    for i=1:rank
        % First problem, x
        udir = zeros(2*nnodes, 1);
        udir(index1) = P1w(:,2*i-1);
        fdir = dirichletRhs2( udir, 3, c2node1s, boundaryp1, nnodes );
        us = K1s\fdir;
        f = Kinter*us(1:2*nnodes,1);
        f1h(:,2*i-1) = f(index1);
        
        % First problem, y
        udir = zeros(2*nnodes, 1);
        udir(index1) = P1w(:,2*i);
        fdir = dirichletRhs2( udir, 3, c2node1s, boundaryp1, nnodes );
        us = K1s\fdir;
        f = Kinter*us(1:2*nnodes,1);
        f1h(:,2*i) = f(index1);
        
        % Second problem, x
        udir = zeros(2*nnodes, 1);
        udir(index2) = P2w(:,2*i-1);
        fdir = dirichletRhs2( udir, 1, c2node2s, boundaryp2, nnodes );
        us = K2s\fdir;
        f = Kinter*us(1:2*nnodes,1);
        f2h(:,2*i-1) = f(index2);
        
        % Second problem, y
        udir = zeros(2*nnodes, 1);
        udir(index2) = P2w(:,2*i);
        fdir = dirichletRhs2( udir, 1, c2node2s, boundaryp2, nnodes );
        us = K2s\fdir;
        f = Kinter*us(1:2*nnodes,1);
        f2h(:,2*i) = f(index2);
        
    end
    
    % The macro forces
    F1 = P1f'*f1h;
    F2 = P2f'*f2h;
    
    % The approximation of the Schur complement : Sa = F*W-1, with W = Id
    Sa1 = P1f*F1*P1w';
    %Sa1 = 0.9*Sa1 + 0.1*norm(Sa1)*eye(size(Sa1));
    Sa2i = P2f*F2*P2w';
    Sa2(m0, m0) = Sa2i(m3, m3); % Re-indexation (to match S2)
    %Sa2 = 0.9*Sa2 + 0.1*norm(Sa2)*eye(size(Sa2));
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Schwarz Robin algorithm

% R1 problem
dirichlet1 = [4,1,0;4,2,0;
              1,1,0;1,2,0;
              2,1,0;2,2,0];
neumann1   = [];
[K1,C1,nbloq1] = Krig (nodes1,elements,E,nu,order,boundary,dirichlet1);
f1in = volumicLoad( nbloq1, nodes1, elements, 1, fscalar );

% Computation of a fictive C matrix (to extract u on the boundaries)
[~,C1f,~] = Krig (nodes1,elements,E,nu,order,boundary, [3,1,0;3,2,0]);
[~,C2f,~] = Krig (nodes1,elements,E,nu,order,boundary, [1,1,0;1,2,0]);

% R2 problem
dirichlet2 = [4,1,0;4,2,0;
              3,1,0;3,2,0;
              2,1,0;2,2,0];
neumann2   = [];

[K2,C2,nbloq2] = Krig (nodes2,elements,E,nu,order,boundary,dirichlet2);
f2in = volumicLoad( nbloq2, nodes2, elements, 1, fscalar );

error1    = zeros(niter,1);
error2    = zeros(niter,1);
residual  = zeros(niter,1);
error1b   = zeros(niter,1);
error2b   = zeros(niter,1);
residualb = zeros(niter,1);
u1        = zeros(2*nnodes,1);
u2        = zeros(2*nnodes,1);

% Directions of research :
robmat1 = robin1 * eye(ndof);
robmat2 = robin2 * eye(ndof);
ub1 = zeros(ndof,1);
ub2 = zeros(ndof,1);
lam1 = zeros(ndof,1);
lam2 = zeros(ndof,1);
%%%%
%%
%urni = projectField(u, nodea, nodes1, 1e-5);
urn = dirichletRhs(u, 3, C1f, boundary);
uref1 = u(map1,1); % uref on the boundary

% uref12 = reshape(uref1,2,[])';
% plot(uref12(:,1));
% figure;

% utest11 = utest1(map1,1); % on the boundary
% uref12 = reshape(utest11,2,[])';
% plot(uref12(:,1));
% figure;
% ub112 = reshape(ub1t,2,[])';
% plot(ub112(:,1));
% figure;
% plot(uref12(:,2));
% figure;
%%
f2 = zeros(2*nnodes,1);
for iter = 1:niter
    %% Solve R1
    fri = projectField( f2, nodes2, nodes1, 1e-5 );
    ui2 = projectField( u2, nodes2, nodes1, 1e-5 );
    kuplusf = ui2 + fri/robin1;

    [Kp1, fro1] = robinRHS( nbloq1, nodes1, boundary, kuplusf, robin1, 3 );
    f1 = f1in + fro1;
    K1e = K1 + Kp1;

    uin1 = K1e\f1;
    u1 = uin1(1:2*nnodes,1);
    f1 = -Kinter*u1 + f1in(1:2*nnodes,1); % Get the forces at the boundary
    
%     hold on;
%     ub1sb = reshape(u1(map1,1),2,[])';
%     plot(ub1sb(:,1),'color','blue');
    %legend('u1 ss assembler')
    %figure;
    
    % Computation of the error :
    u1n = dirichletRhs(u1, 3, C1f, boundary);
    error1(iter) = sqrt( (u1n-urn)'*(u1n-urn)/(urn'*urn) );
    
    % Schur conjugate version
    %ub1 = (S1+robmat1)\(b1+robmat1*ub2-lam2);
    %ub1 = (S1+Sa2+robmat1)\(b1+Sa2*ub2+robmat1*ub2-lam2);
    ub1 = (S1+Sa2)\(b1+Sa2*ub2-lam2);
    lam1 = S1*ub1 - b1;
    u1no = zeros(2*nnodes,1); % Put back ub1 on the mesh
    u1no(map1,1) = ub1;
    
%     ub1sb = reshape(ub1,2,[])';
%     plot(ub1sb(:,1),'color','red');
    %legend('u1 en assemblant')
    %figure;

    % Computation of the error :
    error1b(iter) = sqrt( (ub1-uref1)'*(ub1-uref1)/(uref1'*uref1) );
%     u1nb = dirichletRhs(u1no, 3, C1f, boundary);
%     error1b(iter) = (u1nb-urn)'*(u1nb-urn)/(urn'*urn);
    
    %% Solve R2
    fri = projectField( f1, nodes1, nodes2, 1e-5 );
    ui1 = projectField( u1, nodes1, nodes2, 1e-5 );
    kuplusf = ui1 + fri/robin2;

    [Kp2, fro2] = robinRHS( nbloq2, nodes2, boundary, kuplusf, robin2, 1 );
    f2 = f2in + fro2;
    K2e = K2 + Kp2;

    uin2 = K2e\f2;
    u2 = uin2(1:2*nnodes,1);
    f2 = -Kinter*u2 + f2in(1:2*nnodes,1); % Get the forces at the boundary
    
    % Computation of the error :
    u2ni = projectField(u2, nodes2, nodes1, 1e-5);
    u2n = dirichletRhs(u2ni, 3, C1f, boundary);
    %urni = projectField(u, nodes, nodes1, .5);
    %urn = dirichletRhs(urni, 3, C1f, boundary);
    error2(iter) = sqrt( (u2n-urn)'*(u2n-urn)/(urn'*urn) );
    residual(iter) = sqrt( (u2n-u1n)'*(u2n-u1n)/sqrt((u2n'*u2n)*(u1n'*u1n)) );
    
    % Schur conjugate version
    %ub2 = (S2+robmat2)\(b2+robmat2*ub1-lam1);
    %ub2 = (S2+Sa1+robmat2)\(b2+Sa1*ub1+robmat2*ub1-lam1);
    ub2 = (S2+Sa1)\(b2+Sa1*ub1-lam1);
    lam2 = S2*ub2 - b2;
    u2no = zeros(2*nnodes,1); % Put back ub1 on the mesh
    u2no(map2,1) = ub2;
    % Computation of the error :
    error2b(iter) = sqrt( (ub2-uref1)'*(ub2-uref1)/(uref1'*uref1) );
    residualb(iter) = sqrt( (ub2-ub1)'*(ub2-ub1)/sqrt((ub2'*ub2)*(ub1'*ub1)) );
    
    % Update Sa1 and Sa2 (because rank)
    if iter<2
        Sa1 = norm(Sa1)*eye(size(Sa1));
        Sa2 = norm(Sa2)*eye(size(Sa2));
    end
    
end
%%
sigma1 = stress(u1,E,nu,nodes1,elements,order,1,ntoelem);
sigma2 = stress(u2,E,nu,nodes2,elements,order,1,ntoelem);

plotGMSH({u1,'U_vect';sigma1,'stress'}, elements, nodes1(:,[1,2]), 'field1');
plotGMSH({u2,'U_vect';sigma2,'stress'}, elements, nodes1(:,[1,2]), 'field2');

% ub112 = reshape(ub1,2,[])';
% plot(ub112(:,1));
% figure;
% ub212 = reshape(ub2,2,[])';
% plot(ub212(:,1));
% figure;

hold on;
% plot(error1,'Color','black')
% plot(error2,'Color','blue')
% plot(residual,'Color','red')
set(gca, 'fontsize', 15);
set(gca,'ylim',[-16 0])
plot(log10(error1b),'Color','black')
plot(log10(error2b),'Color','blue')
plot(log10(residualb),'Color','red')
legend('error1 (log)','error2 (log)','residual (log)')
title('Avec les compléments de Schur')
figure;
hold on;
set(gca, 'fontsize', 15);
set(gca,'ylim',[-11 0])
plot(log10(error1),'Color','black')
plot(log10(error2),'Color','blue')
plot(log10(residual),'Color','red')
legend('error1 (log)','error2 (log)','residual (log)')
%title('Avec les résolutions non assemblées')