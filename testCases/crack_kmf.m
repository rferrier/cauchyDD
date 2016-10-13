% 12/04/2016
% Identification de fissures en 2D

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 70000;  % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 1;      % N.mm-1 : Loading on the plate
br      = 0.0;      % noise

% Methods : 1=KMF, 2=KMF Orthodir, 3=KMF Robin, 4=SPP, 5=SPD,
% 6 = Evanescent regul, 7 = CG + Ritz filter
methods = [7];

% Boundary conditions
% first index  : index of the boundary
% second index : 1=x, 2=y
% third        : value
% [0,1,value] marks a dirichlet regularization therm on x
dirichlet_up   = [5,1,0; 5,2,0;
                  1,2,0; 3,2,0 ];
dirichlet_down = [1,1,0; 1,2,0];

% Import the meshes
[ nodes_up,elements_up,ntoelem_up,boundary_up,order ] =...
    readmesh( 'meshes/plate_up.msh' );
nnodes_up = size(nodes_up,1);
[ nodes_down,elements_down,ntoelem_down,boundary_down,order ] =...
    readmesh( 'meshes/plate_down.msh' );
nnodes_down = size(nodes_down,1);

% Then, build the stiffness matrix :
[K_up,C_up,nbloq_up,node2c_up,c2node_up] =...
    Krig (nodes_up,elements_up,E,nu,order,boundary_up,dirichlet_up);
Kinter_up = K_up(1:2*nnodes_up, 1:2*nnodes_up);
M_up      = mass_mat(nodes_up, elements_up);

% find the nodes in the corners and suppress the element :
xmax = max(nodes_up(:,1));
xmin = min(nodes_up(:,1));
ymax = max(nodes_up(:,2));
ymin = min(nodes_up(:,2));
no1  = findNode(xmin, ymin, nodes_up, 1e-5);
no4  = findNode(xmax, ymin, nodes_up, 1e-5);
no5  = findNode(xmax, ymax, nodes_up, 1e-5);
no6  = findNode(xmin, ymax, nodes_up, 1e-5);

% Suppress some nodes from the boundaries
boundaryp1 = suppressBound( boundary_up, [no1;no4], 9 );
boundaryp1 = suppressBound( boundaryp1, [no5], 4 );
boundaryp1 = suppressBound( boundaryp1, [no6], 6 );

[ node2b1, b2node1 ] = mapBound( 1, boundaryp1, nnodes_up );
[ node2b2, b2node2 ] = mapBound( 2, boundaryp1, nnodes_up );
[ node2b3, b2node3 ] = mapBound( 3, boundaryp1, nnodes_up );
[ node2b4, b2node4 ] = mapBound( 4, boundaryp1, nnodes_up );
[ node2b5, b2node5 ] = mapBound( 5, boundaryp1, nnodes_up );
[ node2b6, b2node6 ] = mapBound( 6, boundaryp1, nnodes_up );
[ node2b9, b2node9 ] = mapBound( 9, boundaryp1, nnodes_up );

[K_down,C_down,nbloq_down] =...
    Krig (nodes_down,elements_down,E,nu,order,boundary_down,dirichlet_down);
Kinter_down = K_down(1:2*nnodes_down, 1:2*nnodes_down);
M_down      = mass_mat(nodes_down, elements_down);

% The right hand side :
udi = zeros( nnodes_up, 1 );
ind = 2:2:2*nnodes_up;
udi(ind,1) = 1;
f_up = dirichletRhs2( udi, 5, c2node_up, boundary_up, nnodes_up );

% Solve the problem :
uin = K_up\f_up;

% Extract displacement and Lagrange multiplicators :
uref = uin(1:2*nnodes_up,1);
lagr = uin(2*nnodes_up+1:end,1);
fref = lagr2forces2( lagr, c2node_up, 9, boundary_up, nnodes_up );
urefb = ( 1 + br*randn(2*nnodes_up,1) ) .* uref;
lagrb = ( 1 + br*randn(nbloq_up,1) ) .* lagr;

% Compute stress :
sigma = stress(uref,E,nu,nodes_up,elements_up,order,1,ntoelem_up);
% Output :
plotGMSH({uref,'U_vect';sigma,'stress'}, elements_up, nodes_up, 'reference');

% Plot displacement on the interface :
%index = 2*[b2node1;b2node2;b2node3];
index = 2*b2node9;
plot(uref(index,1));
figure
plot(fref(index,1),'Color','red');

indexxy = [index;index+1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KMF prequisites
if size(find(methods==1),1) == 1 || size(find(methods==2),1) == 1
    % DN problem
    dirichlet1 = [5,1,0;5,2,0;
                  1,1,0;1,2,0;
                  2,1,0;2,2,0;
                  3,1,0;3,2,0;];

    [K1,C1,nbloq1,node2c1,c2node1] =...
        Krig (nodes_up,elements_up,E,nu,order,boundary_up,dirichlet1);
    % Dirichlet upper loading :
    fdir51 = dirichletRhs2( udi, 5, c2node1, boundary_up, nnodes_up );

    % ND problem
    dirichlet2 = [5,1,0;5,2,0;
                  4,1,0;4,2,0;
                  6,1,0;6,2,0];

    [K2,C2,nbloq2,node2c2,c2node2] =...
        Krig (nodes_up,elements_up,E,nu,order,boundary_up,dirichlet2);
    % Dirichlet upper loading :
    fdir52 = dirichletRhs2( udi, 5, c2node2, boundary_up, nnodes_up );
    fdir4  = dirichletRhs2( urefb, 4, c2node2, boundary_up, nnodes_up );
    fdir6  = dirichletRhs2( urefb, 6, c2node2, boundary_up, nnodes_up );

    f2i = fdir4 + fdir6;
    f2 = assembleDirichlet( [f2i,fdir52] );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KMF algo : 10000 iterations
if find(methods==1)

    niter = 1000;
    relax = 0;  % Computation of the best relaxtion parameter
    
    % init :
    u1    = uref-uref;
    u2    = u1;
    v     = u1;
    fri   = u1;

    error1   = zeros(niter,1);
    error2   = zeros(niter,1);
    residual = zeros(niter,1);
    regulari = zeros(niter,1);
    theta    = ones(niter,1); % Relxation parameter

    for iter = 1:niter
        % Solve DN
        u1o = u1;
        f1 = dirichletRhs2(v, 9, c2node1, boundary_up, nnodes_up) + fdir51;
        uin1 = K1\f1;
        u1 = uin1(1:2*nnodes_up,1);
        error1(iter) = norm(u1(indexxy)-uref(indexxy)) / norm(uref(indexxy));
        
        % Keep only the forces at the bottom boundary
        fri = Kinter_up*u1;% - keepField( Kinter_up*u1, 5, boundary_up );
        
        % Solve ND
        u2o = u2;
        fr = [ fri ; zeros(size(C2,2),1) ]; % Give to fr the right size
        uin2 = K2\(fr+f2);
        u2 = uin2(1:2*nnodes_up,1);
        
        vo = v;

        v = theta(iter)*u2 + (1-theta(iter))*vo;
        if relax == 1
            e1 = u1-u1o;
            e2 = u2-u2o;
            theta(iter+1) = e1(indexxy)'*(e1(indexxy)-e2(indexxy)) /...
                norm(e1(indexxy)-e2(indexxy))^2;
        end
        
        error2(iter) = norm(u2(indexxy)-uref(indexxy)) / norm(uref(indexxy));
        residual(iter) = norm(u1(indexxy)-u2(indexxy)) /...
                         sqrt ( norm(u1(indexxy))*norm(u2(indexxy)) );             
        regulari(iter) = sqrt(u2'*( regul(u2, nodes_up, boundary_up, 1) +...
                            regul(u2, nodes_up, boundary_up, 2) +...
                            regul(u2, nodes_up, boundary_up, 3) ));
    end

    figure
    hold on;
    plot(log10(error1),'Color','black')
    plot(log10(error2),'Color','blue')
    plot(log10(residual),'Color','red')
    legend('error1 (log)','error2 (log)','residual (log)')
    % L-curve
    figure
    loglog(residual,regulari);
    % Compute stress :
    sigma = stress(u2,E,nu,nodes_up,elements_up,order,1,ntoelem_up);
    % Output :
    plotGMSH({u2,'U_vect';sigma,'stress'}, elements_up, nodes_up, 'field2');
    % Plot displacement on the interface :
    figure
    plot(u2(index,1));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KMF Orthodir 20 iterations
if find(methods==2)
    niter   = 20;
    mu      = 0;      % Regularization parameter
    
    % Init
    Itere    = zeros( 2*nnodes_up, 1 );
    p        = zeros( 2*nnodes_up, niter+1 );
    q        = zeros( 2*nnodes_up, niter+1 );
    Delta    = zeros( niter+1, 1 );
    Res      = zeros( 2*nnodes_up, niter+1 );
    error    = zeros(niter+1,1);
    residual = zeros(niter+1,1);
    regulari = zeros(niter+1,1);
    
    %% Compute residual
    % Ax0 : Solve DN
%     fdir1 = dirichletRhs2(Itere, 1, c2node1, boundary_up, nnodes_up );
%     fdir2 = dirichletRhs2(Itere, 2, c2node1, boundary_up, nnodes_up );
%     fdir3 = dirichletRhs2(Itere, 3, c2node1, boundary_up, nnodes_up );
%     f1 = assembleDirichlet( [fdir1,fdir2,fdir3] );
    f1 = dirichletRhs2(Itere, 9, c2node1, boundary_up, nnodes_up );
    uin1 = K1\f1;
    u1 = uin1(1:2*nnodes_up,1);
    % Keep only the forces at the bottom boundary
    fri = Kinter_up*u1;
    % Solve ND
    fr = [ fri ; zeros(size(C2,2),1) ]; % Give to fr the right size
    uin2 = K2\fr;
    %
    Nu = regul(Itere, nodes_up, boundary_up, 9);
    atimesItere = mu*Nu + Itere - uin2(1:2*nnodes_up,1);
    %
    % RHS : Solve DN
    f1 = fdir51;
    uin1 = K1\f1;
    u1 = uin1(1:2*nnodes_up,1);
    % Keep only the forces at the bottom boundary
    fri = Kinter_up*u1;
    % Solve ND
    fr = [ fri ; zeros(size(C2,2),1) ]; % Give to fr the right size
    uin2 = K2\(fr+f2);
    %
    b = uin2(1:2*nnodes_up,1);
    %
    Res(:,1) = b - atimesItere;
    p(:,1) = Res(:,1);

    residual(1) = norm(Res(indexxy,1));
    error(1)    = norm( Itere(indexxy) - uref(indexxy )) / norm(uref(indexxy));
    regulari(1) = sqrt(Itere'*( regul(Itere, nodes_up, boundary_up, 9)));

    
    %% Compute Q1 = AP1
    % Solve DN
    f1 = dirichletRhs2(p(:,1), 9, c2node1, boundary_up, nnodes_up );
    uin1 = K1\f1;
    u1 = uin1(1:2*nnodes_up,1);
    % Keep only the forces at the bottom boundary
    fri = Kinter_up*u1;
    % Solve ND
    fr = [ fri ; zeros(size(C2,2),1) ]; % Give to fr the right size
    uin2 = K2\fr;
    %
    Nu = regul(p(:,1), nodes_up, boundary_up, 9);
    q(:,1) = mu*Nu + p(:,1) - uin2(1:2*nnodes_up,1);
    
    for iter = 1:niter

        Delta(iter,1) = norm(q(indexxy,iter))^2;
        gammai        = q(indexxy,iter)'*Res(indexxy,iter);
        alphai        = gammai/Delta(iter,1);

        Itere         = Itere + p(:,iter)*alphai;
        Res(:,iter+1) = Res(:,iter) - q(:,iter)*alphai;

        residual(iter+1) = norm(Res(indexxy,iter+1));
        error(iter+1)    = norm( Itere(indexxy) - uref(indexxy)) / norm(uref(indexxy));
        regulari(iter+1) = sqrt(Itere'*( regul(Itere, nodes_up, boundary_up, 9)));

        %% Perform Ari = A*Res
        % Solve DN
        f1 = dirichletRhs2(Res(:,iter+1), 9, c2node1, boundary_up, nnodes_up );
        uin1 = K1\f1;
        u1 = uin1(1:2*nnodes_up,1);
        % Keep only the forces at the bottom boundary
        fri = Kinter_up*u1;
        % Solve ND
        fr = [ fri ; zeros(size(C2,2),1) ]; % Give to fr the right size
        uin2 = K2\fr;
        %
        Nu = regul(Res(:,iter+1), nodes_up, boundary_up, 9);
        Ari = mu*Nu + Res(:,iter+1) - uin2(1:2*nnodes_up,1);

        %% Orthogonalization
        p(:,iter+1) = Res(:,iter+1);
        q(:,iter+1) = Ari;

        for jter=1:iter
            phiij  = q(indexxy,jter)'*Ari(indexxy);
            betaij = phiij/Delta(jter,1);
            p(:,iter+1) = p(:,iter+1) - p(:,jter) * betaij;
            q(:,iter+1) = q(:,iter+1) - q(:,jter) * betaij;
        end

    end
    figure
    hold on
    plot(log10(error(2:end)),'Color','blue')
    plot(log10(residual(2:end)),'Color','red')
    legend('error (log)','residual (log)')
    % L-curve
    figure
    loglog(residual,regulari);
    % Output
    figure
    plot(Itere(index));
    
    % Compute stress :
    sigma = stress(Itere,E,nu,nodes_up,elements_up,order,1,ntoelem_up);
    % Output :
    plotGMSH({Itere,'U_vect';sigma,'stress'}, elements_up, nodes_up, 'field2');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KMF Robin
if find(methods==3)

    niter = 10;
    robin1  = 1e5*E;%-7*E;%-1e1*E;      % Robin parameters
    %robin2  = -1/7*E;%-1/7*E;%-1e-1*E;
    
    %% Computation of the best stiffness
    dirichlet2 = [9,1,0;9,2,0;
                  %4,1,0;4,2,0;
                  5,1,0;5,2,0;
                  %6,1,0;6,2,0;
                  ];
    [K2,C2,nbloq2,node2c2,c2node2] = Krig (nodes_up,elements_up,E,nu,order,boundary_up,dirichlet2);
    
    fd4 = dirichletRhs2( urefb, 4, c2node2, boundary_up, nnodes_up );
    fd5 = dirichletRhs2( udi, 5, c2node2, boundary_up, nnodes_up );
    fd6 = dirichletRhs2( urefb, 6, c2node2, boundary_up, nnodes_up );
    u2 = K2\fd5;%assembleDirichlet( [fd4+fd6,fd5] );
    f2 = Kinter_up*u2(1:2*nnodes_up);
    b2 = f2([2*b2node9-1 ; 2*b2node9]);
    
    dirichlet1 = [5,1,0;5,2,0;];
    [K1,C1,nbloq1,node2c1,c2node1] = Krig (nodes_up,elements_up,E,nu,order,boundary_up,dirichlet1);
    f1 = zeros(size(K1,1),1); f1([2*b2node9-1 ; 2*b2node9]) = b2;
    u1 = K1\f1;
    D1b2 = u1([2*b2node9-1 ; 2*b2node9]);

    robin2 = -norm(b2)/norm(D1b2);
%     norm(b2)/norm(D1b2)/E
    %% Definition of the problems
    % First problem
    dirichlet1 = [5,1,0;5,2,0];

    [K1,C1,nbloq1,node2c1,c2node1] =...
        Krig (nodes_up,elements_up,E,nu,order,boundary_up,dirichlet1);
    % Dirichlet upper loading :
    fdir51 = dirichletRhs2( udi, 5, c2node1, boundary_up, nnodes_up );

    % Second problem
    dirichlet2 = [5,1,0;5,2,0;
                  4,1,0;4,2,0;
                  6,1,0;6,2,0];

    [K2,C2,nbloq2,node2c2,c2node2] =...
        Krig (nodes_up,elements_up,E,nu,order,boundary_up,dirichlet2);
    % Dirichlet upper loading :
    fdir52 = dirichletRhs2( udi, 5, c2node2, boundary_up, nnodes_up );
    fdir4  = dirichletRhs2( urefb, 4, c2node2, boundary_up, nnodes_up );
    fdir6  = dirichletRhs2( urefb, 6, c2node2, boundary_up, nnodes_up );

    f2i = fdir4 + fdir6;
    f2c = assembleDirichlet( [f2i,fdir52] );
    
    % init :
    u1    = uref-uref;
    u2    = u1;
    f2    = u1;

    error1   = zeros(niter,1);
    error2   = zeros(niter,1);
    residual = zeros(niter,1);
    regulari = zeros(niter,1);

    for iter = 1:niter
    
        % Solve 1
        kuplusf = u2 + f2/robin1;
        [Kp1, fro1] =...
            robinRHS( nbloq1, nodes_up, boundary_up, kuplusf, robin1, 9 );
        f1 = fro1 + fdir51;
        uin1 = (K1+Kp1)\f1;
        u1 = uin1(1:2*nnodes_up,1);
        error1(iter) = norm(u1(indexxy)-uref(indexxy)) / norm(uref(indexxy));
        % Keep only the forces at the bottom boundary
        f1 = Kinter_up*u1;%

        % Solve 2
        kuplusf = u1 + f1/robin2;
        [Kp2, fro2] =...
            robinRHS( nbloq2, nodes_up, boundary_up, kuplusf, robin2, 9 );
        uin2 = (K2+Kp2)\(fro2+f2c);
        u2 = uin2(1:2*nnodes_up,1);
        f2 = Kinter_up*u2;

        error2(iter) = norm(u2(indexxy)-uref(indexxy)) / norm(uref(indexxy));
        residual(iter) = norm(u1(indexxy)-u2(indexxy)) /...
                         sqrt ( norm(u1(indexxy))*norm(u2(indexxy)) );             
        regulari(iter) = sqrt(u2'*( regul(u2, nodes_up, boundary_up, 1) +...
                            regul(u2, nodes_up, boundary_up, 2) +...
                            regul(u2, nodes_up, boundary_up, 3) ));
    end

    figure
    hold on;
    plot(log10(error1),'Color','black')
    plot(log10(error2),'Color','blue')
    plot(log10(residual),'Color','red')
    legend('error1 (log)','error2 (log)','residual (log)')
    % L-curve
    figure
    loglog(residual,regulari);
    % Compute stress :
    sigma = stress(u2,E,nu,nodes_up,elements_up,order,1,ntoelem_up);
    % Output :
    plotGMSH({u2,'U_vect';sigma,'stress'}, elements_up, nodes_up, 'field2');
    % Plot displacement on the interface :
    figure
    plot(u2(index,1));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SP prequisites
if size(find(methods==4),1) == 1 || size(find(methods==5),1) == 1 ||...
   size(find(methods==7),1) == 1
    % find the nodes in the corners and suppress the element :
    xmax = max(nodes_up(:,1));
    xmin = min(nodes_up(:,1));
    ymax = max(nodes_up(:,2));
    ymin = min(nodes_up(:,2));
    no1  = findNode(xmin, ymin, nodes_up, 1e-5);
    no4  = findNode(xmax, ymin, nodes_up, 1e-5);
    boundaryp = suppressBound( boundary_up, [no1;no4], 1 );
    boundaryp = suppressBound( boundaryp, [no1;no4], 3 );
    boundaryp = suppressBound( boundaryp, [no1;no4], 9 );
    %% Definition of the operators
    % First problem
    dirichlet1 = [4,1,0;4,2,0;
                  5,1,0;5,2,0;
                  6,1,0;6,2,0;
                  9,1,0;9,2,0];
    [K1,C1,nbloq1,node2c1,c2node1] =...
        Krig (nodes_up,elements_up,E,nu,order,boundaryp,dirichlet1);

    % Second problem
    dirichlet2 = [5,1,0;5,2,0;
                  9,1,0;9,2,0];
    [K2,C2,nbloq2,node2c2,c2node2] =...
        Krig (nodes_up,elements_up,E,nu,order,boundaryp,dirichlet2);

    %% Dual problems
    % First problem
    dirichlet1d = [4,1,0;4,2,0;
                   5,1,0;5,2,0;
                   6,1,0;6,2,0];
    [K1d,C1d,nbloq1d,node2c1s,c2node1s] =...
        Krig (nodes_up,elements_up,E,nu,order,boundaryp,dirichlet1d);

    % Second problem
    dirichlet2d = [5,1,0;5,2,0];
    [K2d,C2d,nbloq2d,node2c2s,c2node2s] =...
        Krig (nodes_up,elements_up,E,nu,order,boundaryp,dirichlet2d);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Primal SP with Orthodir (niter = 10)
if find(methods==4)
    niter   = 6;
    mu      = 0.01;      % Regularization parameter
    precond = 1;      % use a dual precond ?
    
    % Init
    Itere    = zeros( 2*nnodes_up, 1 );
    p        = zeros( 2*nnodes_up, niter+1 );
    q        = zeros( 2*nnodes_up, niter+1 );
    Delta    = zeros( niter+1, 1 );
    Res      = zeros( 2*nnodes_up, niter+1 );
    Zed      = zeros( 2*nnodes_up, niter+1 );
    error    = zeros(niter+1,1);
    residual = zeros(niter+1,1);
    regulari = zeros(niter+1,1);
    
    %% Perform A x0 :
    % Solve 1
    f1 = dirichletRhs2( Itere, 9, c2node1, boundaryp, nnodes_up );
    uin1 = K1\f1;
    lagr1 = uin1(2*nnodes_up+1:end,1);
    lamb1 = lagr2forces( lagr1, C1, 9, boundaryp );
    % Solve 2
    f2 = dirichletRhs2( Itere, 9, c2node2, boundaryp, nnodes_up );
    uin2 = K2\f2;
    lagr2 = uin2(2*nnodes_up+1:end,1);
    lamb2 = lagr2forces( lagr2, C2, 9, boundaryp );
    % Regularization term
    Nu = regul(Itere, nodes_up, boundaryp, 9);
    %
    Axz = mu*Nu+lamb1-lamb2;
    %%%%
    %% Compute Rhs :
    % Solve 1
    f11 = dirichletRhs2( urefb, 4, c2node1, boundaryp, nnodes_up );
    f12 = dirichletRhs2( urefb, 6, c2node1, boundaryp, nnodes_up );
    f_up = dirichletRhs2( udi, 5, c2node1, boundaryp, nnodes_up );
    f1 = assembleDirichlet( [f11+f12,f_up] );
    uin1 = K1\f1;
    lagr1 = uin1(2*nnodes_up+1:end,1);
    lamb1 = lagr2forces2( lagr1, c2node1, 9, boundaryp, nnodes_up );
    % Solve 2
    f_up = dirichletRhs2( udi, 5, c2node2, boundaryp, nnodes_up );
    uin2 = K2\f_up;
    lagr2 = uin2(2*nnodes_up+1:end,1);
    lamb2 = lagr2forces2( lagr2, c2node2, 9, boundaryp, nnodes_up );
    b = lamb2-lamb1;
    %
    Res(:,1) = b - Axz;
    %% Preconditionning
    if precond == 1
        % Solve 1
        f1 = [Res(:,1)/2; zeros(nbloq1d,1)];
        uin1 = K1d\f1;
        u1i = uin1(1:2*nnodes_up,1);
        u1 = keepField( u1i, 9, boundaryp );
        % Solve 2
        f2 = [Res(:,1)/2; zeros(nbloq2d,1)];
        uin2 = K2d\f2;
        u2i = uin2(1:2*nnodes_up,1);
        u2 = keepField( u2i, 9, boundaryp );
        %
        Zed(:,1) = u1/2-u2/2;
    else
        Zed(:,1) = Res(:,1);
    end
    p(:,1) = Zed(:,1);
    
    residual(1) = norm(Res(indexxy,1));
    error(1)    = norm( Itere(indexxy) - uref(indexxy )) / norm(uref(indexxy));
    regulari(1) = sqrt(Itere'*( regul(Itere, nodes_up, boundary_up, 9)));

    %% Perform Q1 = A P1 :
    % Solve 1
    f1 = dirichletRhs(p(:,1), 9, C1, boundaryp);
    uin1 = K1\f1;
    lagr1 = uin1(2*nnodes_up+1:end,1);
    lamb1 = lagr2forces( lagr1, C1, 9, boundaryp );
    % Solve 2
    f2 = dirichletRhs(p(:,1), 9, C2, boundaryp);
    uin2 = K2\f2;
    lagr2 = uin2(2*nnodes_up+1:end,1);
    lamb2 = lagr2forces( lagr2, C2, 9, boundaryp );
    % Regularization term
    Nu = regul(p(:,1), nodes_up, boundaryp, 9);
    %
    q(:,1) = mu*Nu+lamb1-lamb2;
    %%%%
    
    for iter = 1:niter

        Delta(iter,1) = norm(q(indexxy,iter))^2;
        gammai        = q(indexxy,iter)'*Res(indexxy,iter);
        alphai        = gammai/Delta(iter,1);

        Itere         = Itere + p(:,iter)*alphai;
        Res(:,iter+1) = Res(:,iter) - q(:,iter)*alphai;

        if precond == 1
            % Solve 1
            f1 = [Res(:,iter+1)/2; zeros(nbloq1d,1)];
            uin1 = K1d\f1;
            u1i = uin1(1:2*nnodes_up,1);
            u1 = keepField( u1i, 9, boundaryp );
            % Solve 2
            f2 = [Res(:,iter+1)/2; zeros(nbloq2d,1)];
            uin2 = K2d\f2;
            u2i = uin2(1:2*nnodes_up,1);
            u2 = keepField( u2i, 9, boundaryp );
            %
            Zed(:,iter+1) = u1/2-u2/2;
        else
            Zed(:,iter+1) = Res(:,iter+1);
        end

        residual(iter+1) = norm(Res(indexxy,iter+1));
        error(iter+1)    = norm( Itere(indexxy) - uref(indexxy)) / norm(uref(indexxy));
        regulari(iter+1) = sqrt(Itere'*( regul(Itere, nodes_up, boundary_up, 9)));

        %% Perform Ari = A*Res
        % Solve 1
        rhs1 = Zed(:,iter+1);
        f1 = dirichletRhs(rhs1, 9, C1, boundaryp);
        uin1 = K1\f1;
        lagr1 = uin1(2*nnodes_up+1:end,1);
        lamb1 = lagr2forces( lagr1, C1, 9, boundaryp );
        % Solve 2
        rhs2 = Zed(:,iter+1);
        f2 = dirichletRhs(rhs2, 9, C2, boundaryp);
        uin2 = K2\f2;
        lagr2 = uin2(2*nnodes_up+1:end,1);
        lamb2 = lagr2forces( lagr2, C2, 9, boundaryp );
        % Regularization term
        Nu = regul(Zed(:,iter+1), nodes_up, boundaryp, 9);
        %
        Ari = mu*Nu+lamb1-lamb2;

        %% Orthogonalization
        p(:,iter+1) = Zed(:,iter+1);
        q(:,iter+1) = Ari;

        for jter=1:iter
            phiij  = q(indexxy,jter)'*Ari(indexxy);
            betaij = phiij/Delta(jter,1);
            p(:,iter+1) = p(:,iter+1) - p(:,jter) * betaij;
            q(:,iter+1) = q(:,iter+1) - q(:,jter) * betaij;
        end
    end
    
    figure
    hold on
    plot(log10(error),'Color','blue')
    plot(log10(residual),'Color','red')
    legend('error (log)','residual (log)')
    % L-curve
    figure
    loglog(residual,regulari);
    
    % Output
    figure
    plot(Itere(index));
    
    %% Final problem
    dirichlet = [9,1,0;9,2,0;
                 4,1,0;4,2,0;
                 5,1,0;5,2,0;
                 6,1,0;6,2,0];
    [K,C,nbloq] = Krig (nodes_up,elements_up,E,nu,order,boundary_up,dirichlet);
    fdir6 = dirichletRhs(urefb, 6, C, boundary_up);
    fdir4 = dirichletRhs(urefb, 4, C, boundary_up);
    fdir5 = dirichletRhs(udi, 5, C, boundary_up);
    fdir9 = dirichletRhs(Itere, 9, C, boundary_up);
    usoli = K \ assembleDirichlet( [fdir9+fdir5,fdir4+fdir6] );
    usol = usoli(1:2*nnodes_up,1);

    total_error = norm(usol-uref)/norm(uref);
    % Compute stress :
    sigma = stress(usol,E,nu,nodes_up,elements_up,order,1,ntoelem_up);
    plotGMSH({usol,'U_vect';sigma,'stress'}, elements_up, nodes_up, 'solution');
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dual SP with Orthodir (niter = 8)
if find(methods==5)
    niter   = 8;
    mu      = 0.0/E;      % Regularization parameter
    precond = 1;      % use a primal precond ?
    
    % Init
    Itere    = zeros( 2*nnodes_up, 1 );
    p        = zeros( 2*nnodes_up, niter+1 );
    q        = zeros( 2*nnodes_up, niter+1 );
    Delta    = zeros( niter+1, 1 );
    Res      = zeros( 2*nnodes_up, niter+1 );
    Zed      = zeros( 2*nnodes_up, niter+1 );
    error    = zeros(niter+1,1);
    residual = zeros(niter+1,1);
    regulari = zeros(niter+1,1);
    
    %% Perform A x0 :
    % Solve 1
    f1 = [Itere; zeros(nbloq1d,1)];
    uin1 = K1d\f1;
    u1i = uin1(1:2*nnodes_up,1);
    u1 = keepField( u1i, 9, boundaryp );
    % Solve 2
    f2 = [Itere; zeros(nbloq2d,1)];
    uin2 = K2d\f2;
    u2i = uin2(1:2*nnodes_up,1);
    u2 = keepField( u2i, 9, boundaryp );
    % Regularization term
    Nu = regul(Itere, nodes_up, boundaryp, 9);
    %
    Axz = mu*Nu+u1-u2;
    %%%%
    
    %% Compute Rhs :
    % Solve 1
    f11 = dirichletRhs2( urefb, 4, c2node1s, boundaryp, nnodes_up );
    f12 = dirichletRhs2( urefb, 6, c2node1s, boundaryp, nnodes_up );
    f_up = dirichletRhs2( udi, 5, c2node1s, boundaryp, nnodes_up );
    f1 = assembleDirichlet( [f11+f12,f_up] );
    uin1 = K1d\f1;
    u1i = uin1(1:2*nnodes_up,1);
    u1 = keepField( u1i, 9, boundaryp );
    % Solve 2
    f_up = dirichletRhs2( udi, 5, c2node2s, boundaryp, nnodes_up );
    uin2 = K2d\f_up;
    u2i = uin2(1:2*nnodes_up,1);
    u2 = keepField( u2i, 9, boundaryp );
    b = u2-u1;
    %
    Res(:,1) = b - Axz;
    %% Preconditionning
    if precond == 1
        % Solve 1
        f1 = dirichletRhs(Res(:,1)/2, 9, C1, boundaryp);
        uin1 = K1\f1;
        lagr1 = uin1(2*nnodes_up+1:end,1);
        lamb1 = lagr2forces( lagr1, C1, 9, boundaryp );
        % Solve 2
        f2 = dirichletRhs(Res(:,1)/2, 9, C2, boundaryp);
        uin2 = K2\f2;
        lagr2 = uin2(2*nnodes_up+1:end,1);
        lamb2 = lagr2forces( lagr2, C2, 9, boundaryp );
        %
        Zed(:,1) = lamb1/2-lamb2/2;        
    else
        Zed(:,1) = Res(:,1);
    end
    p(:,1) = Zed(:,1);
    
    residual(1) = norm(Res(indexxy,1));
    error(1)    = norm( Itere(indexxy) - fref(indexxy )) / norm(fref(indexxy));
    regulari(1) = sqrt(Itere'*( regul(Itere, nodes_up, boundary_up, 9)));

    %% Perform Q1 = A P1 :
    % Solve 1
    f1 = [p(:,1); zeros(nbloq1d,1)];
    uin1 = K1d\f1;
    u1i = uin1(1:2*nnodes_up,1);
    u1 = keepField( u1i, 9, boundaryp );
    % Solve 2
    f2 = [p(:,1); zeros(nbloq2d,1)];
    uin2 = K2d\f2;
    u2i = uin2(1:2*nnodes_up,1);
    u2 = keepField( u2i, 9, boundaryp );
    % Regularization term
    Nu = regul(p(:,1), nodes_up, boundaryp, 9);
    %
    q(:,1) = mu*Nu+u1-u2;
    %%%%
    
    for iter = 1:niter

        Delta(iter,1) = norm(q(indexxy,iter))^2;
        gammai        = q(indexxy,iter)'*Res(indexxy,iter);
        alphai        = gammai/Delta(iter,1);

        Itere         = Itere + p(:,iter)*alphai;
        Res(:,iter+1) = Res(:,iter) - q(:,iter)*alphai;

        if precond == 1
            % Solve 1
            f1 = dirichletRhs(Res(:,iter+1)/2, 9, C1, boundaryp);
            uin1 = K1\f1;
            lagr1 = uin1(2*nnodes_up+1:end,1);
            lamb1 = lagr2forces( lagr1, C1, 9, boundaryp );
            % Solve 2
            f2 = dirichletRhs(Res(:,iter+1)/2, 9, C2, boundaryp);
            uin2 = K2\f2;
            lagr2 = uin2(2*nnodes_up+1:end,1);
            lamb2 = lagr2forces( lagr2, C2, 9, boundaryp );
            %
            Zed(:,iter+1) = lamb1/2-lamb2/2; 
        else
            Zed(:,iter+1) = Res(:,iter+1);
        end

        residual(iter+1) = norm(Res(indexxy,iter+1));
        error(iter+1)    = norm( Itere(indexxy) - fref(indexxy)) / norm(fref(indexxy));
        regulari(iter+1) = sqrt(Itere'*( regul(Itere, nodes_up, boundary_up, 9)));

        %% Perform Ari = A*Res
        % Solve 1
        f1 = [Zed(:,iter+1); zeros(nbloq1d,1)];
        uin1 = K1d\f1;
        u1i = uin1(1:2*nnodes_up,1);
        u1 = keepField( u1i, 9, boundaryp );
        % Solve 2
        f2 = [Zed(:,iter+1); zeros(nbloq2d,1)];
        uin2 = K2d\f2;
        u2i = uin2(1:2*nnodes_up,1);
        u2 = keepField( u2i, 9, boundaryp );
        % Regularization term
        Nu = regul(Zed(:,iter+1), nodes_up, boundaryp, 9);
        %
        Ari = mu*Nu+u1-u2;
        %%%%

        %% Orthogonalization
        p(:,iter+1) = Zed(:,iter+1);
        q(:,iter+1) = Ari;

        for jter=1:iter
            phiij  = q(indexxy,jter)'*Ari(indexxy);
            betaij = phiij/Delta(jter,1);
            p(:,iter+1) = p(:,iter+1) - p(:,jter) * betaij;
            q(:,iter+1) = q(:,iter+1) - q(:,jter) * betaij;
        end
    end
    
    figure
    hold on
    plot(log10(error),'Color','blue')
    plot(log10(residual),'Color','red')
    legend('error (log)','residual (log)')
    % L-curve
    figure
    loglog(residual,regulari);
    
    % Output
    figure
    plot(Itere(index),'Color','red');
    
    %% Final problem
    dirichlet = [4,1,0;4,2,0;
                 5,1,0;5,2,0;
                 6,1,0;6,2,0];
    [K,C,nbloq] = Krig (nodes_up,elements_up,E,nu,order,boundary_up,dirichlet);
    fdir6 = dirichletRhs(urefb, 6, C, boundary_up);
    fdir4 = dirichletRhs(urefb, 4, C, boundary_up);
    fdir5 = dirichletRhs(udi, 5, C, boundary_up);
    f1 = [Itere; zeros(nbloq,1)];
    usoli = K \ ( assembleDirichlet( [fdir5,fdir4+fdir6] ) + f1 );
    usol = usoli(1:2*nnodes_up,1);

    total_error = norm(usol-uref)/norm(uref);
    % Compute stress :
    sigma = stress(usol,E,nu,nodes_up,elements_up,order,1,ntoelem_up);
    plotGMSH({usol,'U_vect';sigma,'stress'}, elements_up, nodes_up, 'solution');
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Evanescent regularization method
if find(methods==6)
   niter = 10;
   mu = .001;
   
   %% Computation of the inner stiffness
   dirichlet1 = [5,1,0;5,2,0];
   [K1,C1,nbloq1,node2c1,c2node1] = Krig (nodes_up,elements_up,E,nu,order,boundary_up,dirichlet1);
   neumann1   = []; % There is no alone Neumann
   f1 = loading( nbloq_up, nodes_up, boundary_up, neumann1 );

   % map of the nodes
   b2node469 = [b2node4;b2node6;b2node9];
   bbound = [2*b2node469-1;2*b2node469];
   bbzero = [2*b2node5-1;2*b2node5];
   nbound = size(bbound,1);
   
   %% Schur operator
   [ S, b, map ] = schurComp2( Kinter_up, f1(1:2*nnodes_up), bbound, bbzero, udi );

   error    = zeros(niter,1);
   residual = zeros(niter,1);
   regulari = zeros(niter,1);

   %% Mass matrices
   Mr  = bMass_mat(nodes_up, boundary_up, [9]);
   Mrt = 1/E*Mr;  % Ideally 1/EL
   M   = bMass_mat(nodes_up, boundary_up, [4;6;9]);
   Mt  = 1/E*M;
   Mm  = bMass_mat(nodes_up, boundary_up, 5);

   % Extract coords
   Mr  = Mr(bbound, bbound);
   Mrt = Mrt(bbound, bbound);
   M   = M(bbound, bbound);
   Mt  = Mt(bbound, bbound);
   Mm  = Mm(bbound, bbound);

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% Evanescent regularization method
   Itere  = zeros( 2*nnodes_up, 1 );
   Iteref = zeros( 2*nnodes_up, 1 );

   % Compute errors
   %error(1)    = (Itere(bbound)-uref(bbound))'*Mm*...
    %  (Itere(bbound)-uref(bbound)) / (uref(bbound)'*Mm*uref(bbound));
   error(1) = norm(Itere(indexxy)-uref(indexxy))/norm(uref(indexxy));
   residual(1) = (Itere(bbound)-urefb(bbound))'*Mr*...
      (Itere(bbound)-urefb(bbound)) / (uref(bbound)'*Mr*uref(bbound));
   regulari(1) = 0;

   % Build the fat matrix
   Atot = [Mr+mu*M, zeros(size(M)), S'
           zeros(size(M)), Mrt+mu*Mt, -eye(size(M,1),size(S,1))
           S, -eye(size(S,1), size(M,1)), zeros(size(M))];
        
   %disp( [ 'Log of the cond of the problem :', num2str(log10(cond(Atot))) ] )

   %figure;
   %plot(log10(abs(eig(Atot))));
   % plot(eig(Atot));

   for i = 2:niter

      % Rhs
      btot = [Mr*urefb(bbound) + mu*M*Itere(bbound)
              Mrt*fref(bbound) + mu*Mt*Iteref(bbound)
              b]; % ud is contained in b
              %b + S*udi(bbound)];

      % Solve and extract the relevant parts
      Iterep = Itere;
      xtot   = Atot\btot;
      Itere(bbound)  = xtot(1:nbound);
      Iteref(bbound) = xtot(nbound+1:2*nbound);
   
      % Compute errors
      error(i) = norm(Itere(indexxy)-uref(indexxy))/norm(uref(indexxy));
      %error(i)    = (Itere(bbound)-uref(bbound))'*Mm*...
       %  (Itere(bbound)-uref(bbound)) / (uref(bbound)'*Mm*uref(bbound));
      residual(i) = (Itere(bbound)-urefb(bbound))'*Mr*...
         (Itere(bbound)-urefb(bbound)) / (uref(bbound)'*Mr*uref(bbound));
      regulari(i) = Itere(bbound)'*M*Itere(bbound) / (uref(bbound)'*M*uref(bbound));
   
   end

   figure;
   hold on
   plot(log10(error),'Color','blue')
   plot(log10(residual),'Color','red')
   legend('error (log)','residual (log)')

   %%%%
   %% Final problem : compute u
   % DN problem
   dirichlet = [4,1,0;4,2,0;
                5,1,0;5,2,0;
                6,1,0;6,2,0;
                9,1,0;9,2,0];
   neumann   = [];
   [K,C,nbloq] = Krig (nodes_up,elements_up,E,nu,order,boundary_up,dirichlet_up);
   fdir4 = dirichletRhs(urefb, 4, C, boundary_up);
   fdir6 = dirichletRhs(urefb, 6, C, boundary_up);
   fdir9 = dirichletRhs(Itere, 9, C, boundary_up);
   fdir5 = dirichletRhs(udi  , 5, C, boundary_up);
   usoli = K \ assembleDirichlet( [fdir4+fdir6,fdir5+fdir9] );

   usol = usoli(1:2*nnodes_up,1);
   fsol = Kinter_up*usol;

   % Plot displacement on the interface :
   efe = Kinter_up*usol;
   figure
   hold on
   set(gca, 'fontsize', 15);
   plot(usol(index,1));
   plot(usol(index-1,1), 'Color', 'red');
   figure
   hold on
   set(gca, 'fontsize', 15);
   plot(efe(index,1));
   plot(efe(index-1,1), 'Color', 'red');
    
   total_error = norm(uref-usol)/norm(uref);
   total_errorf = norm(fref-efe)/norm(fref);

   % Compute stress :
   sigma = stress(usol,E,nu,nodes_up,elements_up,order,1,ntoelem_up);
   plotGMSH({usol,'U_vect';sigma,'stress'}, elements_up, nodes_up, 'solution');
end

if find(methods == 7)

   niter   = 15;
   precond = 0;      % 1 : Use a dual precond
   ratio   = 1e-200;    % Maximal ratio (for eigenfilter)
   epsilon = 1e-1;   % Convergence criterion for ritz value

   %% Conjugate Gradient for the problem : (S10-S20) x = S2-S1
   Itere = zeros( 2*nnodes_up, 1 );
   d     = zeros( 2*nnodes_up, niter+1 );
   Ad    = zeros( 2*nnodes_up, niter+1 );
   Res   = zeros( 2*nnodes_up, niter+1 );
   Zed   = zeros( 2*nnodes_up, niter+1 );
   alpha = zeros( niter+1, 1 );
   beta  = zeros( niter+1, 1 );
   ntrunc = 10;  % In case the algo finishes at niter
   
   %% Perform A x0 :
   % Solve 1
   f1 = dirichletRhs2( Itere, 9, c2node1, boundaryp, nnodes_up );
   uin1 = K1\f1;
   lagr1 = uin1(2*nnodes_up+1:end,1);
   lamb1 = lagr2forces( lagr1, C1, 9, boundaryp );
   % Solve 2
   f2 = dirichletRhs2( Itere, 9, c2node2, boundaryp, nnodes_up );
   uin2 = K2\f2;
   lagr2 = uin2(2*nnodes_up+1:end,1);
   lamb2 = lagr2forces( lagr2, C2, 9, boundaryp );
   %
    Axz = lamb1-lamb2;
   %%%%
   %% Compute RHS
   % Solve 1
   f11 = dirichletRhs2( urefb, 4, c2node1, boundaryp, nnodes_up );
   f12 = dirichletRhs2( urefb, 6, c2node1, boundaryp, nnodes_up );
   f_up = dirichletRhs2( udi, 5, c2node1, boundaryp, nnodes_up );
   f1 = assembleDirichlet( [f11+f12,f_up] );
   uin1 = K1\f1;
   lagr1 = uin1(2*nnodes_up+1:end,1);
   lamb1 = lagr2forces2( lagr1, c2node1, 9, boundaryp, nnodes_up );
   % Solve 2
   f_up = dirichletRhs2( udi, 5, c2node2, boundaryp, nnodes_up );
   uin2 = K2\f_up;
   lagr2 = uin2(2*nnodes_up+1:end,1);
   lamb2 = lagr2forces2( lagr2, c2node2, 9, boundaryp, nnodes_up );
   b = lamb2-lamb1;
   %
   Res(:,1) = b - Axz;
   
   if precond == 1
       % Solve 1
       f1 = [Res(:,1); zeros(nbloq1d,1)];
       uin1 = K1d\f1;
       u1i = uin1(1:2*nnodes_up,1);
       u1 = keepField( u1i, 9, boundaryp );
       Zed(:,1) = u1;
   else
       Zed(:,1) = Res(:,1);
   end
   
   d(:,1) = Zed(:,1);
   
   residual(1) = sqrt( norm(Res( indexxy,1)));
   error(1)    = sqrt( norm(Itere(indexxy) - uref(indexxy)) / norm(uref(indexxy)));
   regulari(1) = sqrt( Itere'*regul(Itere, nodes_up, boundary_up, 2) );
   
   ritzval  = 0; % Last ritz value that converged
   oldtheta = 0;
   eta      = 0;
   %%
   for iter = 1:niter
       %% Optimal step
       
       % Solve 1
       f1 = dirichletRhs2( d(:,iter), 9, c2node1, boundaryp, nnodes_up );
       uin1 = K1\f1;
       lagr1 = uin1(2*nnodes_up+1:end,1);
       lamb1 = lagr2forces( lagr1, C1, 9, boundaryp );
       % Solve 2
       f2 = dirichletRhs2( d(:,iter), 9, c2node2, boundaryp, nnodes_up );
       uin2 = K2\f2;
       lagr2 = uin2(2*nnodes_up+1:end,1);
       lamb2 = lagr2forces( lagr2, C2, 9, boundaryp );
       %
       Ad(:,iter) = lamb1-lamb2;
       
       den = (d(indexxy,iter)'*Ad(indexxy,iter));
       d(:,iter) = d(:,iter)/sqrt(den); Ad(:,iter) = Ad(:,iter)/sqrt(den);
       num = Res(indexxy,iter)'*d(indexxy,iter);
       
       Itere         = Itere + d(:,iter)*num;%/den;
       Res(:,iter+1) = Res(:,iter) - Ad(:,iter)*num;%/den;
       
       residual(iter+1) = sqrt( norm(Res(indexxy,iter+1)));
       error(iter+1)    = sqrt( norm(Itere(indexxy) - uref(indexxy)) / norm(uref(indexxy)));
       regulari(iter+1) = sqrt( Itere'*regul(Itere, nodes_up, boundary_up, 2) );
       
       if precond == 1
           % Solve 1
           f1 = [Res(:,iter+1); zeros(nbloq1d,1)];
           uin1 = K1d\f1;
           u1i = uin1(1:2*nnodes_up,1);
           u1 = keepField( u1i, 9, boundaryp );
           Zed(:,iter+1) = u1;
       else
           Zed(:,iter+1) = Res(:,iter+1);
       end
       
       % Needed values for the Ritz stuff
       alpha(iter) = Res(indexxy,iter)'*Res(indexxy,iter) / den;
       beta(iter)  = Zed(indexxy,iter+1)'*Res(indexxy,iter+1) /... 
                                   (Zed(indexxy,iter)'*Res(indexxy,iter));
       
       % First Reorthogonalize the residual (as we use it next), in sense of M
       for jter=1:iter-1
           betac = Zed(indexxy,iter+1)'*Res(indexxy,jter) / (Zed(indexxy,jter)'*Res(indexxy,jter));
           Zed(:,iter+1) = Zed(:,iter+1) - Zed(:,jter) * betac;
           Res(:,iter+1) = Res(:,iter+1) - Res(:,jter) * betac;
       end
       
       %% Orthogonalization
       d(:,iter+1) = Zed(:,iter+1);
       
       for jter=iter:iter % No need to reorthogonalize (see above)
           betaij = ( Zed(indexxy,iter+1)'*Ad(indexxy,jter) );%/...
               %( d(indexxy,jter)'*Ad(indexxy,jter) );
           d(:,iter+1) = d(:,iter+1) - d(:,jter) * betaij;
   
       end
       
       %% Ritz algo : find the Ritz elements
       % Build the matrices
       V(:,iter) = zeros(2*nnodes_up,1);
       V(indexxy,iter) = (-1)^(iter-1)*Zed(indexxy,iter)/(sqrt(Res(indexxy,iter)'*Zed(indexxy,iter)));
       etap   = eta;
       delta  = 1/alpha(iter);
       if iter > 1
          delta  = delta + beta(iter-1)/alpha(iter-1);
       end
       eta    = sqrt(beta(iter))/alpha(iter);
       
       if iter > 1
          H(iter,[iter-1,iter]) = [etap, delta];
          H(iter-1,iter)        = etap;
       else
          H(iter,iter) = delta;
       end
       
       % Compute eigenelems of the Hessenberg :
       [Q,Theta1] = eig(H);
       theta = diag(Theta1);
       % Sort it
       [theta,Ind] = sort(theta,'descend');
       Q = Q(:,Ind);
       Theta1 = Theta1(Ind,Ind);
       Y = V*Q;
       
       % See if the current one converged
       if abs(theta(ritzval+1) - oldtheta) < epsilon*oldtheta
          % increment oldtheta
          ritzval = ritzval+1;
          if size(theta,1) > ritzval
             oldtheta = theta(ritzval+1);
          else
             oldtheta = 0;
          end
          
          % Check small value / hold value
          if ritzval > 1
             if theta(ritzval) < ratio*theta(1)
                ntrunc = ritzval
                getmeout = 1;
                break;
             end
          end
       else
          oldtheta = theta(ritzval+1);
       end
   end
   
   % Compute the solution
   chi = inv(Theta1)*Y'*b;
   if ntrunc > 0
      chi(ntrunc:end) = 0;
   end
   ItereR = Y*chi;
   
   figure;
   hold on;
   plot(log10(theta),'Color','blue')
   plot(log10(abs(Y'*b)),'Color','red')
   plot(log10(abs(chi)),'Color','black')
   legend('Ritz Values','RHS values','solution coefficients')
   
   figure;
   hold on;
   plot(Itere(index),'Color','red')
   plot(ItereR(index),'Color','blue')
   plot(uref(index),'Color','green')
   legend('brutal solution','filtred solution', 'reference')
   
   %% Final problem : compute u
   dirichlet = [4,1,0;4,2,0;
                3,1,0;3,2,0;
                1,1,0;1,2,0;
                2,1,0;2,2,0];
   neumann   = [];
   [K,C,nbloq] = Krig (nodes_up,elements_up,E,nu,order,boundary_up,dirichlet_up);
   fdir2 = dirichletRhs(ItereR, 2, C, boundary_up);
   fdir4 = dirichletRhs(urefb, 4, C, boundary_up);
   usoli = K \ (fdir4 + fdir2);
   
   usol = usoli(1:2*nnodes_up,1);
   fsol = Kinter_up*usol;
   
   total_error = norm(usol-uref)/norm(uref);
   % Compute stress :
   sigma = stress(usol,E,nu,nodes_up,elements_up,order,1,ntoelem_up);
   plotGMSH({usol,'U_vect';sigma,'stress'}, elements_up, nodes_up, 'solution');
end