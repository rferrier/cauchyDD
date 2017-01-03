% 20/05/2016
% Identification de fissures en 2D (domaine mince)

close all;
clear all;

addpath(genpath('./tools'))
addpath(genpath('./regu'))

% Parameters
E       = 70000;  % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 1;      % N.mm-1 : Loading on the plate
br      = 0.;      % noise

% Methods : 1=KMF, 2=KMF Orthodir, 3=KMF Robin, 4=SPP, 5=SPD rigid modes
% 6=SPD with 4=0
methods = [6];

% Boundary conditions
% first index  : index of the boundary
% second index : 1=x, 2=y
% third        : value
% [0,1,value] marks a dirichlet regularization therm on x
dirichlet_up   = [0,1,0; 1,2,0; 3,2,0 ];
neumann        = [5,2,fscalar];
dirichlet_down = [1,1,0; 1,2,0];

% Import the meshes
[ nodes_up,elements_up,ntoelem_up,boundary_up,order ] =...
    readmesh( 'meshes/plate_up2.msh' );
nnodes_up = size(nodes_up,1);
[ nodes_down,elements_down,ntoelem_down,boundary_down,order ] =...
    readmesh( 'meshes/plate_down.msh' );
nnodes_down = size(nodes_down,1);

% Then, build the stiffness matrix :
[K_up,C_up,nbloq_up,node2c_up,c2node_up] =...
    Krig (nodes_up,elements_up,E,nu,order,boundary_up,dirichlet_up);
Kinter_up = K_up(1:2*nnodes_up, 1:2*nnodes_up);
M_up      = mass_mat(nodes_up, elements_up);
[ node2b1, b2node1 ] = mapBound( 1, boundary_up, nnodes_up );
[ node2b2, b2node2 ] = mapBound( 2, boundary_up, nnodes_up );
[ node2b3, b2node3 ] = mapBound( 3, boundary_up, nnodes_up );
[ node2b9, b2node9 ] = mapBound( 9, boundary_up, nnodes_up );
[K_down,C_down,nbloq_down] =...
    Krig (nodes_down,elements_down,E,nu,order,boundary_down,dirichlet_down);
Kinter_down = K_down(1:2*nnodes_down, 1:2*nnodes_down);
M_down      = mass_mat(nodes_down, elements_down);

% The right hand side :
f_up = loading(nbloq_up,nodes_up,boundary_up,neumann);

% Solve the problem :
uin = K_up\f_up;

% Extract displacement and Lagrange multiplicators :
uref = uin(1:2*nnodes_up,1);
lagr = uin(2*nnodes_up+1:end,1);
fref = Kinter_up*uref - f_up(1:2*nnodes_up,1);
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
    dirichlet1 = [1,1,0;1,2,0;
                  2,1,0;2,2,0;
                  3,1,0;3,2,0;];
    neumann1 = [5,2,fscalar];
              
    [K1,C1,nbloq1,node2c1,c2node1] =...
        Krig (nodes_up,elements_up,E,nu,order,boundary_up,dirichlet1);

    fn1 = loading(nbloq1,nodes_up,boundary_up,neumann1);
    
    % ND problem
    dirichlet2 = [5,1,0;5,2,0;
                  4,1,0;4,2,0;
                  6,1,0;6,2,0];

    [K2,C2,nbloq2,node2c2,c2node2] =...
        Krig (nodes_up,elements_up,E,nu,order,boundary_up,dirichlet2);
    % Dirichlet upper loading :
    fdir5 = dirichletRhs2( urefb, 5, c2node2, boundary_up, nnodes_up );
    fdir4  = dirichletRhs2( urefb, 4, c2node2, boundary_up, nnodes_up );
    fdir6  = dirichletRhs2( urefb, 6, c2node2, boundary_up, nnodes_up );

    f2i = fdir4 + fdir6;
    f2 = assembleDirichlet( [f2i,fdir5] );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KMF algo : 200 iterations
if find(methods==1)

    niter = 10;
    relax = 0;  % Computation of the best relaxtion parameter
    
    % init :
    u1    = uref-uref;
    u2    = zeros(2*nnodes_up,niter);
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
        f1 = dirichletRhs2(v, 9, c2node1, boundary_up, nnodes_up);
        uin1 = K1\(f1+fn1);
        u1 = uin1(1:2*nnodes_up,1);
        error1(iter) = norm(u1(indexxy)-uref(indexxy)) / norm(uref(indexxy));
        
        % Keep only the forces at the bottom boundary
        fri = Kinter_up*u1;% - keepField( Kinter_up*u1, 5, boundary_up );
        
        % Solve ND
        u2o = u2(:,iter);
        fr = [ fri ; zeros(size(C2,2),1) ]; % Give to fr the right size
        uin2 = K2\(fr+f2);
        u2(:,iter) = uin2(1:2*nnodes_up,1);
        
        vo = v;

        v = theta(iter)*u2(:,iter) + (1-theta(iter))*vo;
        if relax == 1
            e1 = u1-u1o;
            e2 = u2(:,iter)-u2o;
            theta(iter+1) = e1(indexxy)'*(e1(indexxy)-e2(indexxy)) /...
                norm(e1(indexxy)-e2(indexxy))^2;
        end
        
        error2(iter) = norm(u2(indexxy,iter)-uref(indexxy)) / norm(uref(indexxy));
        residual(iter) = norm(u1(indexxy)-u2(indexxy,iter)) /...
                         sqrt ( norm(u1(indexxy))*norm(u2(indexxy,iter)) );             
        regulari(iter) = sqrt(u2(:,iter)'*( regul(u2(:,iter), nodes_up, boundary_up, 1) +...
                            regul(u2(:,iter), nodes_up, boundary_up, 2) +...
                            regul(u2(:,iter), nodes_up, boundary_up, 3) ));
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
    % Automatic stop (regularization tools) :
    [stopHere,~,~] = l_corner(residual,regulari,1:1:niter);
    % Compute stress :
    sigma = stress(u2(:,stopHere),E,nu,nodes_up,elements_up,order,1,ntoelem_up);
    % Output :
    plotGMSH({u2(:,stopHere),'U_vect';sigma,'stress'}, elements_up, nodes_up, 'field2');
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
    Itere    = zeros( 2*nnodes_up, niter+1 );
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
    f1 = dirichletRhs2(Itere(:,1), 9, c2node1, boundary_up, nnodes_up );
    uin1 = K1\f1;
    u1 = uin1(1:2*nnodes_up,1);
    % Keep only the forces at the bottom boundary
    fri = Kinter_up*u1;
    % Solve ND
    fr = [ fri ; zeros(size(C2,2),1) ]; % Give to fr the right size
    uin2 = K2\fr;
    %
    Nu = regul(Itere(:,1), nodes_up, boundary_up, 9);
    atimesItere = mu*Nu + Itere(:,1) - uin2(1:2*nnodes_up,1);
    %
    % RHS : Solve DN
    f1 = fn1;
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
    error(1)    = norm( Itere(indexxy,1) - uref(indexxy )) / norm(uref(indexxy));
    regulari(1) = sqrt(Itere(:,1)'*( regul(Itere(:,1), nodes_up, boundary_up, 9)));

    
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

        Itere(:,iter+1) = Itere(:,iter) + p(:,iter)*alphai;
        Res(:,iter+1)   = Res(:,iter) - q(:,iter)*alphai;

        residual(iter+1) = norm(Res(indexxy,iter+1));
        error(iter+1)    = norm( Itere(indexxy,iter+1) - uref(indexxy)) / norm(uref(indexxy));
        regulari(iter+1) = sqrt(Itere(:,iter+1)'*( regul(Itere(:,iter+1), nodes_up, boundary_up, 9)));

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
    plot(log10(error),'Color','blue')
    plot(log10(residual),'Color','red')
    legend('error (log)','residual (log)')
    % L-curve
    figure
    loglog(residual,regulari);
    % Automatic stop (regularization tools) :
    [stopHere,~,~] = l_corner(residual,regulari,1:1:niter+1);
    % Output
    figure
    plot(Itere(index,stopHere));
    
    % Compute stress :
    sigma = stress(Itere(:,stopHere),E,nu,nodes_up,elements_up,order,1,ntoelem_up);
    % Output :
    plotGMSH({Itere(:,stopHere),'U_vect';sigma,'stress'}, elements_up, nodes_up, 'field2');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KMF Robin
if find(methods==3)

    niter = 20;
    robin1  = 1e5*E;    % Robin parameters
    robin2  = 1e-5*E;
    
    %% Definition of the problems
    % First problem
    dirichlet1 = [];
    neumann1 = [5,2,fscalar];

    [K1,C1,nbloq1,node2c1,c2node1] =...
        Krig (nodes_up,elements_up,E,nu,order,boundary_up,dirichlet1);
    % upper loading :
    fn1 = loading(nbloq1,nodes_up,boundary_up,neumann1);

    % Second problem
    dirichlet2 = [5,1,0;5,2,0;
                  4,1,0;4,2,0;
                  6,1,0;6,2,0];

    [K2,C2,nbloq2,node2c2,c2node2] =...
        Krig (nodes_up,elements_up,E,nu,order,boundary_up,dirichlet2);
    % Dirichlet upper loading :
    fdir5 = dirichletRhs2( urefb, 5, c2node2, boundary_up, nnodes_up );
    fdir4  = dirichletRhs2( urefb, 4, c2node2, boundary_up, nnodes_up );
    fdir6  = dirichletRhs2( urefb, 6, c2node2, boundary_up, nnodes_up );

    f2i = fdir4 + fdir6;
    f2c = assembleDirichlet( [f2i,fdir5] );
    
    % init :
    u1    = uref-uref;
    u2 = zeros(2*nnodes_up,niter+1);
    f2    = u1;

    error1   = zeros(niter,1);
    error2   = zeros(niter,1);
    residual = zeros(niter,1);
    regulari = zeros(niter,1);

    for iter = 1:niter
    
        % Solve 1
        kuplusf = u2(:,iter) + f2/robin1;
        [Kp1, fro1] =...
            robinRHS( nbloq1, nodes_up, boundary_up, kuplusf, robin1, 9 );
        f1 = fro1 + fn1;
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
        u2(:,iter+1) = uin2(1:2*nnodes_up,1);
        f2 = Kinter_up*u2(:,iter+1);

        error2(iter) = norm(u2(indexxy,iter+1)-uref(indexxy)) / norm(uref(indexxy));
        residual(iter) = norm(u1(indexxy)-u2(indexxy,iter+1)) /...
                         sqrt ( norm(u1(indexxy))*norm(u2(indexxy+1,iter)) );             
        regulari(iter) = sqrt(u2(:,iter+1)'*( regul(u2(:,iter+1), nodes_up, boundary_up, 1) +...
                            regul(u2(:,iter+1), nodes_up, boundary_up, 2) +...
                            regul(u2(:,iter+1), nodes_up, boundary_up, 3) ));
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
    % Automatic stop (regularization tools) :
    [stopHere,~,~] = l_corner(residual,regulari,1:1:niter);
    % Compute stress :
    sigma = stress(u2(:,stopHere+1),E,nu,nodes_up,elements_up,order,1,ntoelem_up);
    % Output :
    plotGMSH({u2(:,stopHere+1),'U_vect';sigma,'stress'}, elements_up, nodes_up, 'field2');
    % Plot displacement on the interface :
    figure
    plot(u2(index,stopHere+1));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SP prequisites
if size(find(methods==4),1) == 1 || size(find(methods==5),1) == 1 || size(find(methods==6),1) == 1
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
    dirichlet2 = [9,1,0;9,2,0];
    [K2,C2,nbloq2,node2c2,c2node2] =...
        Krig (nodes_up,elements_up,E,nu,order,boundaryp,dirichlet2);

    %% Dual problems
    % First problem
    dirichlet1d = [4,1,0;4,2,0;
                   5,1,0;5,2,0;
                   6,1,0;6,2,0];
    [K1d,C1d,nbloq1d,node2c1s,c2node1s] =...
        Krig (nodes_up,elements_up,E,nu,order,boundaryp,dirichlet1d);

    neumann2 = [5,2,fscalar];
    
    % Second problem
    dirichlet2d = [];
    [K2d,C2d,~,node2c2s,c2node2s] =...
        Krig (nodes_up,elements_up,E,nu,order,boundaryp,dirichlet2d);
    
    % Anti-cancellation stuff :
    %K1d(indexxy,indexxy) = 0;
    %K2d(indexxy,indexxy) = 0;
    
    % Management of the rigid modes
    % G = null(full(K2d));
    % More efficient method provided dim(ker) = 3
    r1 = zeros(2*nnodes_up,1); r2 = r1; r3p = r1;
    ind = 2:2:2*nnodes_up;
    r1(ind-1,1) = 1; r2(ind,1) = 1;
    g1 = keepField( r1, 9, boundary_up );%boundaryp ?
    g2 = keepField( r2, 9, boundary_up );
    r3p(ind,1) = nodes_up(ind/2,1);
    r3p(ind-1,1) = -nodes_up(ind/2,2);
    g3p = keepField( r3p, 9, boundary_up );
    % orthogonalize G
    g3 = g3p - (g3p'*g1)/(g1'*g1)*g1 - (g3p'*g2)/(g2'*g2)*g2;
    r3 = r3p - (r3p'*g1)/(g1'*g1)*r1 - (r3p'*g2)/(g2'*g2)*r2;
    % Normalize G
    ng1 = norm(g1); ng2 = norm(g2); ng3 = norm(g3);
    g1 = g1/ng1; g2 = g2/ng2; g3 = g3/ng3;
    r1 = r1/ng1; r2 = r2/ng2; r3 = r3/ng3;
    
    G = [g1,g2,g3];
    R = [r1,r2,r3];
    K2d = [ K2d, R ; R', zeros(size(R,2)) ];
    nbloq2d = size(R,2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Primal SP with Orthodir (niter = 20)
if find(methods==4)
    niter   = 20;
    mu      = 0.;      % Regularization parameter
    precond = 0;      % use a dual precond ?
    
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
    f_up = dirichletRhs2( urefb, 5, c2node1, boundaryp, nnodes_up );
    f1 = assembleDirichlet( [f11+f12,f_up] );
    uin1 = K1\f1;
    lagr1 = uin1(2*nnodes_up+1:end,1);
    lamb1 = lagr2forces2( lagr1, c2node1, 9, boundaryp, nnodes_up );
    % Solve 2
    f_up = loading(nbloq2,nodes_up,boundary_up,neumann2);
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
    fdir5 = dirichletRhs(urefb, 5, C, boundary_up);
    fdir9 = dirichletRhs(Itere, 9, C, boundary_up);
    usoli = K \ assembleDirichlet( [fdir9+fdir5,fdir4+fdir6] );
    usol = usoli(1:2*nnodes_up,1);

    total_error = norm(usol-uref)/norm(uref);
    % Compute stress :
    sigma = stress(usol,E,nu,nodes_up,elements_up,order,1,ntoelem_up);
    plotGMSH({usol,'U_vect';sigma,'stress'}, elements_up, nodes_up, 'solution');
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dual SP with Orthodir (niter = 5)
if find(methods==5)
    niter   = 5;
    mu      = 0.0/E;      % Regularization parameter
    precond = 0;      % use a primal precond ?
    
    P = eye(2*nnodes_up) - G*G';  % Projector
    
    % Init
    p        = zeros( 2*nnodes_up, niter+1 );
    q        = zeros( 2*nnodes_up, niter+1 );
    Delta    = zeros( niter+1, 1 );
    Res      = zeros( 2*nnodes_up, niter+1 );
    Zed      = zeros( 2*nnodes_up, niter+1 );
    error    = zeros(niter+1,1);
    residual = zeros(niter+1,1);
    regulari = zeros(niter+1,1);

    %% Compute Rhs :
    % Solve 1
    f11 = dirichletRhs2( urefb, 4, c2node1s, boundaryp, nnodes_up );
    f12 = dirichletRhs2( urefb, 6, c2node1s, boundaryp, nnodes_up );
    f_up = dirichletRhs2( urefb, 5, c2node1s, boundaryp, nnodes_up );
    f1 = assembleDirichlet( [f11+f12,f_up] );
    uin1 = K1d\f1;
    u1i = uin1(1:2*nnodes_up,1);
    u1 = keepField( u1i, 9, boundaryp );
    % Solve 2
    f_up = loading(nbloq2d,nodes_up,boundary_up,neumann2);
    uin2 = K2d\f_up;
    u2i = uin2(1:2*nnodes_up,1);
    u2 = keepField( u2i, 9, boundaryp );
    b = u2-u1 ;
    
    Itere = -G*R'*f_up(1:2*nnodes_up);

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
    %
    Res(:,1) = P'*(b - Axz) ;
    %%%%

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
        Zed(:,1) = P*(lamb1/2-lamb2/2);        
    else
        Zed(:,1) = P*Res(:,1);
    end
    p(:,1) = Zed(:,1);
    
    residual(1) = 1; % norm(Res(indexxy,1));
    error(1)    = norm( Itere(indexxy) - fref(indexxy )) / norm(fref(indexxy));
    regulari(1) = sqrt( abs( Itere'*( regul(Itere, nodes_up, boundary_up, 9))));

    %% Perform Q1 = A P1 :
    % Solve 1
    f1 = [p(1:2*nnodes_up,1); zeros(nbloq1d,1)];
    uin1 = K1d\f1;
    u1i = uin1(1:2*nnodes_up,1);
    u1 = keepField( u1i, 9, boundaryp );
    % Solve 2
    f2 = [p(1:2*nnodes_up,1); zeros(nbloq2d,1)];
    uin2 = K2d\f2;
    u2i = uin2(1:2*nnodes_up,1);
    u2 = keepField( u2i, 9, boundaryp );
    % Regularization term
    Nu = regul(p(1:2*nnodes_up,1), nodes_up, boundaryp, 9);
    %
    q(:,1) = P'*(mu*Nu+u1-u2);
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
            Zed(:,iter+1) = P*(lamb1/2-lamb2/2); 
        else
            Zed(:,iter+1) = P*(Res(:,iter+1));
        end

        residual(iter+1) = norm(Res(indexxy,iter+1)) / norm(Res(indexxy,1));
        error(iter+1)    = norm( Itere(indexxy) - fref(indexxy)) / norm(fref(indexxy));
        regulari(iter+1) = sqrt( abs(Itere'*( regul(Itere, nodes_up, boundary_up, 9))));

        %% Perform Ari = A*Res
        % Solve 1
        f1 = [Zed(1:2*nnodes_up,iter+1); zeros(nbloq1d,1)];
        uin1 = K1d\f1;
        u1i = uin1(1:2*nnodes_up,1);
        u1 = keepField( u1i, 9, boundaryp );
        % Solve 2
        f2 = [Zed(1:2*nnodes_up,iter+1); zeros(nbloq2d,1)];
        uin2 = K2d\f2;
        u2i = uin2(1:2*nnodes_up,1);
        u2 = keepField( u2i, 9, boundaryp );
        % Regularization term
        Nu = regul(Zed(1:2*nnodes_up,iter+1), nodes_up, boundaryp, 9);
        %
        Ari = P'*(mu*Nu+u1-u2) ;
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
    fdir5 = dirichletRhs(urefb, 5, C, boundary_up);
    f1 = [Itere(1:2*nnodes_up); zeros(nbloq,1)];
    usoli = K \ ( assembleDirichlet( [fdir5,fdir4+fdir6] ) + f1 );
    usol = usoli(1:2*nnodes_up,1);

    figure
    plot(usol(index));
    total_error = norm(usol-uref)/norm(uref);
    % Compute stress :
    sigma = stress(usol,E,nu,nodes_up,elements_up,order,1,ntoelem_up);
    plotGMSH({usol,'U_vect';sigma,'stress'}, elements_up, nodes_up, 'solution');
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dual SP with Orthodir (no rigid mode) (niter = 8)
if find(methods==6)
    
    % Non singular second problem
    dirichlet2d = [4,1,0;4,2,0];
    [K2d,C2d,nbloq2d,node2c2s,c2node2s] =...
        Krig (nodes_up,elements_up,E,nu,order,boundaryp,dirichlet2d);
    
    niter   = 8;
    mu      = 0.0/E;      % Regularization parameter
    precond = 0;      % use a primal precond ?
    
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
    f_up = dirichletRhs2( urefb, 5, c2node1s, boundaryp, nnodes_up );
    f1 = assembleDirichlet( [f11+f12,f_up] );
    uin1 = K1d\f1;
    u1i = uin1(1:2*nnodes_up,1);
    u1 = keepField( u1i, 9, boundaryp );
    % Solve 2
    f_up = loading(nbloq2d,nodes_up,boundary_up,neumann2);
    uin2 = K2d\f_up;
    u2i = uin2(1:2*nnodes_up,1);
    u2 = keepField( u2i, 9, boundaryp );
    b = u2-u1 ;
    %
    Res(:,1) = b - Axz ;
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
    Nu = regul(p(1:2*nnodes_up,1), nodes_up, boundaryp, 9);
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
        Nu = regul(Zed(1:2*nnodes_up,iter+1), nodes_up, boundaryp, 9);
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
    fdir5 = dirichletRhs(urefb, 5, C, boundary_up);
    f1 = [Itere(1:2*nnodes_up); zeros(nbloq,1)];
    usoli = K \ ( assembleDirichlet( [fdir5,fdir4+fdir6] ) + f1 );
    usol = usoli(1:2*nnodes_up,1);

    total_error = norm(usol-uref)/norm(uref);
    % Compute stress :
    sigma = stress(usol,E,nu,nodes_up,elements_up,order,1,ntoelem_up);
    plotGMSH({usol,'U_vect';sigma,'stress'}, elements_up, nodes_up, 'solution');
    
end