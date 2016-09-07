% 13/05/2016
% Probl�me de tuyau

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 70000;  % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 1;      % N.mm-1 : Loading on the plate
br      = 0.0;      % noise

% Methods : 1=KMF, 2=KMF Orthodir, 3=KMF Robin, 4=SPP, 5=SPD,
% 6=SPD flottant, 7=SPD flottant constraint, 8=evanescent regu

methods = [8];

% Boundary conditions
% first index  : index of the boundary
% second index : 1=x, 2=y
% third        : value
% [0,1,value] marks a dirichlet regularization therm on x
dirichlet   = [2,1,0; 2,2,0];

% Import the meshes
[ nodes,elements,ntoelem,boundary,order ] =...
    readmesh( 'meshes/tube.msh' );
nnodes = size(nodes,1);

% patch('Faces',elements,'Vertices',nodes,'FaceAlpha',0);
% figure

% Then, build the stiffness matrix :
[K,C,nbloq,node2c,c2node] =...
    Krig (nodes,elements,E,nu,order,boundary,dirichlet,1);
Kinter = K(1:2*nnodes, 1:2*nnodes);
M      = mass_mat(nodes, elements);
[ node2b1, b2node1 ] = mapBound( 1, boundary, nnodes );
[ node2b2, b2node2 ] = mapBound( 2, boundary, nnodes );
[ node2b3, b2node3 ] = mapBound( 3, boundary, nnodes );
[ node2b9, b2node9 ] = mapBound( 9, boundary, nnodes );

% The right hand side :
f = pressureLoad( nbloq, nodes, boundary, fscalar*[10*sin(pi/6),0;-1,0], 8 );

% Solve the problem :
uin = K\f;

% Extract displacement and Lagrange multiplicators :
uref = uin(1:2*nnodes,1);
fref = Kinter*uref;
lagr = uin(2*nnodes+1:end,1);
urefb = ( 1 + br*randn(2*nnodes,1) ) .* uref;
frefb = ( 1 + br*randn(2*nnodes,1) ) .* fref;

% Compute stress :
sigma = stress(uref,E,nu,nodes,elements,order,1,ntoelem,1);
% Output :
plotGMSH({uref,'U_vect';sigma,'stress'}, elements, nodes, 'reference');

% Plot displacement on the interface :
%index = 2*[b2node1;b2node2;b2node3];
index = 2*b2node3;
thetax = 0:2*pi/size(index,1):2*pi*(1-1/size(index,1));
hold on
set(gca, 'fontsize', 15);
plot(thetax,uref(index,1));
plot(thetax,uref(index-1,1),'Color','red');
legend('uy','ux')
xlabel('angle(rad)')
figure
hold on
set(gca, 'fontsize', 15);
plot(thetax,fref(index,1));
plot(thetax,fref(index-1,1),'Color','red');
legend('fy','fx')
xlabel('angle(rad)')

indexxy = [index;index+1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KMF prequisites
if size(find(methods==1),1) == 1 || size(find(methods==2),1) == 1
    % DN problem
    dirichlet1 = [2,1,0;2,2,0;
                  3,1,0;3,2,0];

    [K1,C1,nbloq1,node2c1,c2node1] =...
        Krig (nodes,elements,E,nu,order,boundary,dirichlet1,1);

    % ND problem
    dirichlet2 = [2,1,0;2,2,0;
                  1,1,0;1,2,0];

    [K2,C2,nbloq2,node2c2,c2node2] =...
        Krig (nodes,elements,E,nu,order,boundary,dirichlet2,1);
    % Dirichlet loading :
    f2  = dirichletRhs2( urefb, 1, c2node2, boundary, nnodes );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KMF algo : 10 iterations
if find(methods==1)

    niter = 100;
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
        f1 = dirichletRhs2(v, 3, c2node1, boundary, nnodes);
        uin1 = K1\f1;
        u1 = uin1(1:2*nnodes,1);
        error1(iter) = norm(u1(indexxy)-uref(indexxy)) / norm(uref(indexxy));
        
        % Keep only the forces at the bottom boundary
        fri = Kinter*u1;% - keepField( Kinter_up*u1, 5, boundary_up );
        
        % Solve ND
        u2o = u2;
        fr = [ fri ; zeros(size(C2,2),1) ]; % Give to fr the right size
        uin2 = K2\(fr+f2);
        u2 = uin2(1:2*nnodes,1);
        
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
        regulari(iter) = sqrt(u2'*( regul(u2, nodes, boundary, 3) ));
    end
    
    figure
    hold on;
    set(gca, 'fontsize', 15);
    plot(log10(error1),'Color','black')
    plot(log10(error2),'Color','blue')
    plot(log10(residual),'Color','red')
    legend('error1 (log)','error2 (log)','residual (log)')
    % L-curve
    figure
    loglog(residual,regulari);
    % Compute stress :
    sigma = stress(u2,E,nu,nodes,elements,order,1,ntoelem,1);
    % Output :
    plotGMSH({u2,'U_vect';sigma,'stress'}, elements, nodes, 'field2');
    % Plot displacement on the interface :
    efe = Kinter*u2;
    figure
    hold on
    set(gca, 'fontsize', 15);
    plot(thetax,u2(index,1));
    plot(thetax,u2(index-1,1), 'Color', 'red');
    legend('uy','ux')
    xlabel('angle(rad)')
    figure
    hold on
    set(gca, 'fontsize', 15);
    plot(thetax,efe(index,1));
    plot(thetax,efe(index-1,1), 'Color', 'red');
    legend('fy','fx')
    xlabel('angle(rad)')

    total_error = norm(uref-u2)/norm(uref);
    total_errorf = norm(fref-efe)/norm(fref);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KMF Orthodir 10 iterations
if find(methods==2)
    niter   = 10;
    mu      = 0;      % Regularization parameter
    
    % Init
    Itere    = zeros( 2*nnodes, 1 );
    p        = zeros( 2*nnodes, niter+1 );
    q        = zeros( 2*nnodes, niter+1 );
    Delta    = zeros( niter+1, 1 );
    Res      = zeros( 2*nnodes, niter+1 );
    error    = zeros(niter+1,1);
    residual = zeros(niter+1,1);
    regulari = zeros(niter+1,1);
    
    %% Compute residual
    % Ax0 : Solve DN
    f1 = dirichletRhs2(Itere, 3, c2node1, boundary, nnodes );
    uin1 = K1\f1;
    u1 = uin1(1:2*nnodes,1);
    % Keep only the forces at the bottom boundary
    fri = Kinter*u1;
    % Solve ND
    fr = [ fri ; zeros(size(C2,2),1) ]; % Give to fr the right size
    uin2 = K2\fr;
    %
    Nu = regul(Itere, nodes, boundary, 3);
    atimesItere = mu*Nu + Itere - uin2(1:2*nnodes,1);
    %
    % RHS : Solve DN (this one is useless because 0)
    f1 = zeros(2*nnodes+nbloq1,1);
    uin1 = K1\f1;
    u1 = uin1(1:2*nnodes,1);
    % Keep only the forces at the bottom boundary
    fri = Kinter*u1;
    % Solve ND
    fr = [ fri ; zeros(size(C2,2),1) ]; % Give to fr the right size
    uin2 = K2\(fr+f2);
    %
    b = uin2(1:2*nnodes,1);
    %
    Res(:,1) = b - atimesItere;
    p(:,1) = Res(:,1);

    residual(1) = 1; %norm(Res(indexxy,1));
    error(1)    = norm( Itere(indexxy) - uref(indexxy )) / norm(uref(indexxy));
    regulari(1) = sqrt(Itere'*( regul(Itere, nodes, boundary, 3)));

    
    %% Compute Q1 = AP1
    % Solve DN
    f1 = dirichletRhs2(p(:,1), 3, c2node1, boundary, nnodes );
    uin1 = K1\f1;
    u1 = uin1(1:2*nnodes,1);
    % Keep only the forces at the bottom boundary
    fri = Kinter*u1;
    % Solve ND
    fr = [ fri ; zeros(size(C2,2),1) ]; % Give to fr the right size
    uin2 = K2\fr;
    %
    Nu = regul(p(:,1), nodes, boundary, 3);
    q(:,1) = mu*Nu + p(:,1) - uin2(1:2*nnodes,1);
    
    for iter = 1:niter

        Delta(iter,1) = norm(q(indexxy,iter))^2;
        gammai        = q(indexxy,iter)'*Res(indexxy,iter);
        alphai        = gammai/Delta(iter,1);

        Itere         = Itere + p(:,iter)*alphai;
        Res(:,iter+1) = Res(:,iter) - q(:,iter)*alphai;

        residual(iter+1) = norm(Res(indexxy,iter+1)) / norm(Res(indexxy,1));
        error(iter+1)    = norm( Itere(indexxy) - uref(indexxy)) / norm(uref(indexxy));
        regulari(iter+1) = sqrt(Itere'*( regul(Itere, nodes, boundary, 3)));

        %% Perform Ari = A*Res
        % Solve DN
        f1 = dirichletRhs2(Res(:,iter+1), 3, c2node1, boundary, nnodes );
        uin1 = K1\f1;
        u1 = uin1(1:2*nnodes,1);
        % Keep only the forces at the bottom boundary
        fri = Kinter*u1;
        % Solve ND
        fr = [ fri ; zeros(size(C2,2),1) ]; % Give to fr the right size
        uin2 = K2\fr;
        %
        Nu = regul(Res(:,iter+1), nodes, boundary, 3);
        Ari = mu*Nu + Res(:,iter+1) - uin2(1:2*nnodes,1);

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
    set(gca, 'fontsize', 15);
    plot(log10(error),'Color','blue')
    plot(log10(residual),'Color','red')
    legend('error (log)','residual (log)')
    % L-curve
    figure
    loglog(residual,regulari);
    % Output : terminal computation
    dirichlet = [1,1,0;1,2,0;
                 2,1,0;2,2,0;
                 3,1,0;3,2,0];
    [K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet,1);
    fdir1 = dirichletRhs(urefb, 1, C, boundary);
    fdir3 = dirichletRhs(Itere, 3, C, boundary);
    usoli = K \ ( fdir1 + fdir3 );
    usol = usoli(1:2*nnodes,1);
    
    efe = Kinter*usol;
    figure
    hold on
    set(gca, 'fontsize', 15);
    plot(thetax,usol(index,1));
    plot(thetax,usol(index-1,1), 'Color', 'red');
    legend('uy','ux')
    xlabel('angle(rad)')
    figure
    hold on
    set(gca, 'fontsize', 15);
    plot(thetax,efe(index,1));
    plot(thetax,efe(index-1,1), 'Color', 'red');
    legend('fy','fx')
    xlabel('angle(rad)')
    
    total_error = norm(uref-usol)/norm(uref);
    total_errorf = norm(fref-efe)/norm(fref);
    
    % Compute stress :
    sigma = stress(usol,E,nu,nodes,elements,order,1,ntoelem,1);
    % Output :
    plotGMSH({usol,'U_vect';sigma,'stress'}, elements, nodes, 'field2');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KMF Robin
if find(methods==3)

    niter = 10;
    robin1  = 1e5*E; % Robin parameters
    %robin2  = 1e-5*E;
    
    %% Computation of the best stiffness
    dirichlet2 = [3,1,0;3,2,0;
                  1,1,0;1,2,0;
                  2,1,0;2,2,0];
    [K2,C2,nbloq2,node2c2,c2node2] = Krig (nodes,elements,E,nu,order,boundary,dirichlet2);
    
    fd1 = dirichletRhs2( urefb, 1, c2node2, boundary, nnodes );
    fd2 = dirichletRhs2( urefb, 2, c2node2, boundary, nnodes );
    u2 = K2\assembleDirichlet( [fd1,fd2] );
    f2 = Kinter*u2(1:2*nnodes);
    b2 = f2([2*b2node3-1 ; 2*b2node3]);
    bug
    dirichlet1 = [2,1,0;2,2,0];
    [K1,C1,nbloq1,node2c1,c2node1] = Krig (nodes,elements,E,nu,order,boundary,dirichlet1);
    f1 = zeros(size(K1,1),1); f1([2*b2node3-1 ; 2*b2node3]) = b2;
    u1 = K1\f1;
    D1b2 = u1([2*b2node3-1 ; 2*b2node3]);

    robin2 = -norm(b2)/norm(D1b2);
    %norm(b2)/norm(D1b2)/E
    %% Definition of the problems
    % First problem
    dirichlet1 = [2,1,0;2,2,0];

    [K1,C1,nbloq1,node2c1,c2node1] =...
        Krig (nodes,elements,E,nu,order,boundary,dirichlet1,1);
  
    % Second problem
    dirichlet2 = [2,1,0;2,2,0;
                  1,1,0;1,2,0];

    [K2,C2,nbloq2,node2c2,c2node2] =...
        Krig (nodes,elements,E,nu,order,boundary,dirichlet2,1);
    % Dirichlet loading :
    f2c = dirichletRhs2( urefb, 1, c2node2, boundary, nnodes );
    
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
            robinRHS( nbloq1, nodes, boundary, kuplusf, robin1, 3 );
        f1 = fro1;
        uin1 = (K1+Kp1)\f1;
        u1 = uin1(1:2*nnodes,1);
        error1(iter) = norm(u1(indexxy)-uref(indexxy)) / norm(uref(indexxy));
        % Keep only the forces at the bottom boundary
        f1 = Kinter*u1;%

        % Solve 2
        kuplusf = u1 + f1/robin2;
        [Kp2, fro2] =...
            robinRHS( nbloq2, nodes, boundary, kuplusf, robin2, 3 );
        uin2 = (K2+Kp2)\(fro2+f2c);
        u2 = uin2(1:2*nnodes,1);
        f2 = Kinter*u2;

        error2(iter) = norm(u2(indexxy)-uref(indexxy)) / norm(uref(indexxy));
        residual(iter) = norm(u1(indexxy)-u2(indexxy)) /...
                         sqrt ( norm(u1(indexxy))*norm(u2(indexxy)) );             
        regulari(iter) = sqrt(u2'*( regul(u2, nodes, boundary, 3) ));
    end

    figure
    hold on;
    set(gca, 'fontsize', 15);
    plot(log10(error1),'Color','black')
    plot(log10(error2),'Color','blue')
    plot(log10(residual),'Color','red')
    legend('error1 (log)','error2 (log)','residual (log)')
    % L-curve
    figure
    loglog(residual,regulari);
    % Compute stress :
    sigma = stress(u2,E,nu,nodes,elements,order,1,ntoelem,1);
    % Output :
    plotGMSH({u2,'U_vect';sigma,'stress'}, elements, nodes, 'field2');
    % Plot displacement on the interface :
    efe = Kinter*u2;
    figure
    hold on
    set(gca, 'fontsize', 15);
    plot(thetax,u2(index,1));
    plot(thetax,u2(index-1,1), 'Color', 'red');
    xlabel('angle(rad)')
    figure
    hold on
    set(gca, 'fontsize', 15);
    plot(thetax,efe(index,1));
    plot(thetax,efe(index-1,1), 'Color', 'red');
    xlabel('angle(rad)')

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SP prequisites
if size(find(methods==4),1) == 1 || size(find(methods==5),1) == 1
    % no node to suppress
    boundaryp = boundary;
    %% Definition of the operators
    % First problem
    dirichlet1 = [2,1,0;2,2,0;
                  1,1,0;1,2,0
                  3,1,0;3,2,0];
    [K1,C1,nbloq1,node2c1,c2node1] =...
        Krig (nodes,elements,E,nu,order,boundaryp,dirichlet1,1);

    % Second problem
    dirichlet2 = [2,1,0;2,2,0;
                  3,1,0;3,2,0];
    [K2,C2,nbloq2,node2c2,c2node2] =...
        Krig (nodes,elements,E,nu,order,boundaryp,dirichlet2,1);

    %% Dual problems
    % First problem
    dirichlet1d = [2,1,0;2,2,0;
                   1,1,0;1,2,0];
    [K1d,C1d,nbloq1d,node2c1s,c2node1s] =...
        Krig (nodes,elements,E,nu,order,boundaryp,dirichlet1d,1);

    % Second problem
    dirichlet2d = [2,1,0;2,2,0];
    [K2d,C2d,nbloq2d,node2c2s,c2node2s] =...
        Krig (nodes,elements,E,nu,order,boundaryp,dirichlet2d,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Primal SP with Orthodir (niter = 10)
if find(methods==4)
    niter   = 100;
    mu      = 0;      % Regularization parameter
    precond = 0;      % use a dual precond ?
    
    % Init
    Itere    = zeros( 2*nnodes, 1 );
    p        = zeros( 2*nnodes, niter+1 );
    q        = zeros( 2*nnodes, niter+1 );
    Delta    = zeros( niter+1, 1 );
    Res      = zeros( 2*nnodes, niter+1 );
    Zed      = zeros( 2*nnodes, niter+1 );
    error    = zeros(niter+1,1);
    residual = zeros(niter+1,1);
    regulari = zeros(niter+1,1);
    
    %% Perform A x0 :
    % Solve 1
    f1 = dirichletRhs2( Itere, 3, c2node1, boundaryp, nnodes );
    uin1 = K1\f1;
    lagr1 = uin1(2*nnodes+1:end,1);
    lamb1 = lagr2forces( lagr1, C1, 3, boundaryp );
    % Solve 2
    f2 = dirichletRhs2( Itere, 3, c2node2, boundaryp, nnodes );
    uin2 = K2\f2;
    lagr2 = uin2(2*nnodes+1:end,1);
    lamb2 = lagr2forces( lagr2, C2, 3, boundaryp );
    % Regularization term
    Nu = regul(Itere, nodes, boundaryp, 3);
    %
    Axz = mu*Nu+lamb1-lamb2;
    %%%%
    %% Compute Rhs :
    % Solve 1
    f = dirichletRhs2( urefb, 1, c2node1, boundaryp, nnodes );
    uin1 = K1\f;
    lagr1 = uin1(2*nnodes+1:end,1);
    lamb1 = lagr2forces2( lagr1, c2node1, 3, boundaryp, nnodes );
    % Solve 2 (zero)
    f = zeros(2*nnodes+nbloq2);
    uin2 = K2\f;
    lagr2 = uin2(2*nnodes+1:end,1);
    lamb2 = lagr2forces2( lagr2, c2node2, 3, boundaryp, nnodes );
    b = lamb2-lamb1;
    %
    Res(:,1) = b - Axz;
    %% Preconditionning
    if precond == 1
        % Solve 1
        f1 = [Res(:,1)/2; zeros(nbloq1d,1)];
        uin1 = K1d\f1;
        u1i = uin1(1:2*nnodes,1);
        u1 = keepField( u1i, 3, boundaryp );
        % Solve 2
        f2 = [Res(:,1)/2; zeros(nbloq2d,1)];
        uin2 = K2d\f2;
        u2i = uin2(1:2*nnodes,1);
        u2 = keepField( u2i, 3, boundaryp );
        %
        Zed(:,1) = u1/2-u2/2;
    else
        Zed(:,1) = Res(:,1);
    end
    p(:,1) = Zed(:,1);
    
    residual(1) = 1;%norm(Res(indexxy,1));
    error(1)    = norm( Itere(indexxy) - uref(indexxy )) / norm(uref(indexxy));
    regulari(1) = sqrt(Itere'*( regul(Itere, nodes, boundary, 9)));

    %% Perform Q1 = A P1 :
    % Solve 1
    f1 = dirichletRhs(p(:,1), 3, C1, boundaryp);
    uin1 = K1\f1;
    lagr1 = uin1(2*nnodes+1:end,1);
    lamb1 = lagr2forces( lagr1, C1, 3, boundaryp );
    % Solve 2
    f2 = dirichletRhs(p(:,1), 3, C2, boundaryp);
    uin2 = K2\f2;
    lagr2 = uin2(2*nnodes+1:end,1);
    lamb2 = lagr2forces( lagr2, C2, 3, boundaryp );
    % Regularization term
    Nu = regul(p(:,1), nodes, boundaryp, 3);
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
            u1i = uin1(1:2*nnodes,1);
            u1 = keepField( u1i, 3, boundaryp );
            % Solve 2
            f2 = [Res(:,iter+1)/2; zeros(nbloq2d,1)];
            uin2 = K2d\f2;
            u2i = uin2(1:2*nnodes,1);
            u2 = keepField( u2i, 3, boundaryp );
            %
            Zed(:,iter+1) = +u1/2-u2/2;
        else
            Zed(:,iter+1) = Res(:,iter+1);
        end

        residual(iter+1) = norm(Res(indexxy,iter+1)) / norm(Res(indexxy,1));
        error(iter+1)    = norm( Itere(indexxy) - uref(indexxy)) / norm(uref(indexxy));
        regulari(iter+1) = sqrt(Itere'*( regul(Itere, nodes, boundary, 3)));

        %% Perform Ari = A*Res
        % Solve 1
        rhs1 = Zed(:,iter+1);
        f1 = dirichletRhs(rhs1, 3, C1, boundaryp);
        uin1 = K1\f1;
        lagr1 = uin1(2*nnodes+1:end,1);
        lamb1 = lagr2forces( lagr1, C1, 3, boundaryp );
        % Solve 2
        rhs2 = Zed(:,iter+1);
        f2 = dirichletRhs(rhs2, 3, C2, boundaryp);
        uin2 = K2\f2;
        lagr2 = uin2(2*nnodes+1:end,1);
        lamb2 = lagr2forces( lagr2, C2, 3, boundaryp );
        % Regularization term
        Nu = regul(Zed(:,iter+1), nodes, boundaryp, 3);
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
    set(gca, 'fontsize', 15);
    plot(error,'Color','blue')
    plot(residual,'Color','red')
    legend('error','residual')
    % L-curve
    figure
    loglog(residual,regulari);
    
    %% Final problem
    dirichlet = [2,1,0;2,2,0;
                 3,1,0;3,2,0];
    [K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet,1);
    %fdir1 = dirichletRhs(urefb, 1, C, boundary);
    fdir3 = dirichletRhs(Itere, 3, C, boundary);
    usoli = K \ fdir3 ;
    usol = usoli(1:2*nnodes,1);

    % Compute stress :
    sigma = stress(usol,E,nu,nodes,elements,order,1,ntoelem,1);
    plotGMSH({usol,'U_vect';sigma,'stress'}, elements, nodes, 'solution');
    
    % Plot displacement on the interface :
    efe = Kinter*usol;
    figure
    hold on
    set(gca, 'fontsize', 15);
    plot(thetax,usol(index,1));
    plot(thetax,usol(index-1,1), 'Color', 'red');
    xlabel('angle(rad)')
    figure
    hold on
    set(gca, 'fontsize', 15);
    plot(thetax,efe(index,1));
    plot(thetax,efe(index-1,1), 'Color', 'red');
    xlabel('angle(rad)')
    
    total_error = norm(uref-usol)/norm(uref);
    total_errorf = norm(fref-efe)/norm(fref);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dual SP with Orthodir (niter = 70)
if find(methods==5)
    niter   = 100;
    mu      = 0.0/E;      % Regularization parameter
    precond = 0;      % use a primal precond ?
    
    % Init
    Itere    = zeros( 2*nnodes, 1 );
    p        = zeros( 2*nnodes, niter+1 );
    q        = zeros( 2*nnodes, niter+1 );
    Delta    = zeros( niter+1, 1 );
    Res      = zeros( 2*nnodes, niter+1 );
    Zed      = zeros( 2*nnodes, niter+1 );
    error    = zeros(niter+1,1);
    residual = zeros(niter+1,1);
    regulari = zeros(niter+1,1);
    
    %% Perform A x0 :
    % Solve 1
    f1 = [Itere; zeros(nbloq1d,1)];
    uin1 = K1d\f1;
    u1i = uin1(1:2*nnodes,1);
    u1 = keepField( u1i, 3, boundaryp );
    % Solve 2
    f2 = [Itere; zeros(nbloq2d,1)];
    uin2 = K2d\f2;
    u2i = uin2(1:2*nnodes,1);
    u2 = keepField( u2i, 3, boundaryp );
    % Regularization term
    Nu = regul(Itere, nodes, boundaryp, 3);
    %
    Axz = mu*Nu+u1-u2;
    %%%%
    
    %% Compute Rhs :
    % Solve 1
    f1 = dirichletRhs2( urefb, 1, c2node1s, boundaryp, nnodes );
    uin1 = K1d\f1;
    u1i = uin1(1:2*nnodes,1);
    u1 = keepField( u1i, 3, boundaryp );
    % Solve 2 (ZERO)
    f = zeros(2*nnodes+nbloq2d,1);
    uin2 = K2d\f;
    u2i = uin2(1:2*nnodes,1);
    u2 = keepField( u2i, 3, boundaryp );
    b = u2-u1;
    %
    Res(:,1) = b - Axz;
    %% Preconditionning
    if precond == 1
        % Solve 1
        f1 = dirichletRhs(Res(:,1)/2, 3, C1, boundaryp);
        uin1 = K1\f1;
        lagr1 = uin1(2*nnodes+1:end,1);
        lamb1 = lagr2forces( lagr1, C1, 3, boundaryp );
        % Solve 2
        f2 = dirichletRhs(Res(:,1)/2, 3, C2, boundaryp);
        uin2 = K2\f2;
        lagr2 = uin2(2*nnodes+1:end,1);
        lamb2 = lagr2forces( lagr2, C2, 3, boundaryp );
        %
        Zed(:,1) = lamb1/2-lamb2/2;        
    else
        Zed(:,1) = Res(:,1);
    end
    p(:,1) = Zed(:,1);
    
    residual(1) = 1;%norm(Res(indexxy,1));
    error(1)    = norm( Itere(indexxy) - fref(indexxy )) / norm(fref(indexxy));
    regulari(1) = sqrt(Itere'*( regul(Itere, nodes, boundary, 3)));

    %% Perform Q1 = A P1 :
    % Solve 1
    f1 = [p(:,1); zeros(nbloq1d,1)];
    uin1 = K1d\f1;
    u1i = uin1(1:2*nnodes,1);
    u1 = keepField( u1i, 3, boundaryp );
    % Solve 2
    f2 = [p(:,1); zeros(nbloq2d,1)];
    uin2 = K2d\f2;
    u2i = uin2(1:2*nnodes,1);
    u2 = keepField( u2i, 3, boundaryp );
    % Regularization term
    Nu = regul(p(:,1), nodes, boundaryp, 3);
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
            f1 = dirichletRhs(Res(:,iter+1)/2, 3, C1, boundaryp);
            uin1 = K1\f1;
            lagr1 = uin1(2*nnodes+1:end,1);
            lamb1 = lagr2forces( lagr1, C1, 3, boundaryp );
            % Solve 2
            f2 = dirichletRhs(Res(:,iter+1)/2, 3, C2, boundaryp);
            uin2 = K2\f2;
            lagr2 = uin2(2*nnodes+1:end,1);
            lamb2 = lagr2forces( lagr2, C2, 3, boundaryp );
            %
            Zed(:,iter+1) = lamb1/2-lamb2/2; 
        else
            Zed(:,iter+1) = Res(:,iter+1);
        end

        residual(iter+1) = norm(Res(indexxy,iter+1)) / norm(Res(indexxy,1));
        error(iter+1)    = norm( Itere(indexxy) - fref(indexxy)) / norm(fref(indexxy));
        regulari(iter+1) = sqrt(Itere'*( regul(Itere, nodes, boundary, 3)));

        %% Perform Ari = A*Res
        % Solve 1
        f1 = [Zed(:,iter+1); zeros(nbloq1d,1)];
        uin1 = K1d\f1;
        u1i = uin1(1:2*nnodes,1);
        u1 = keepField( u1i, 3, boundaryp );
        % Solve 2
        f2 = [Zed(:,iter+1); zeros(nbloq2d,1)];
        uin2 = K2d\f2;
        u2i = uin2(1:2*nnodes,1);
        u2 = keepField( u2i, 3, boundaryp );
        % Regularization term
        Nu = regul(Zed(:,iter+1), nodes, boundaryp, 3);
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
    set(gca, 'fontsize', 15);
    plot(log10(error),'Color','blue')
    plot(log10(residual),'Color','red')
    legend('error (log)','residual (log)')
    % L-curve
    figure
    loglog(residual,regulari);
    
    %% Final problem
    dirichlet = [2,1,0;2,2,0];
    [K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet,1);
    %fdir1 = dirichletRhs(urefb, 1, C, boundary);
    fdir2 = dirichletRhs(urefb, 2, C, boundary);
    f1 = [Itere; zeros(nbloq,1)];
    usoli = K \ ( fdir2 + f1 );
    usol = usoli(1:2*nnodes,1);

    % Compute stress :
    sigma = stress(usol,E,nu,nodes,elements,order,1,ntoelem,1);
    plotGMSH({usol,'U_vect';sigma,'stress'}, elements, nodes, 'solution');
    
    % Plot displacement on the interface :
    efe = Kinter*usol;
    figure
    hold on
    set(gca, 'fontsize', 15);
    plot(thetax,usol(index,1));
    plot(thetax,usol(index-1,1), 'Color', 'red');
    xlabel('angle(rad)')
    figure
    hold on
    set(gca, 'fontsize', 15);
    plot(thetax,efe(index,1));
    plot(thetax,efe(index-1,1), 'Color', 'red');
    xlabel('angle(rad)')
    
    total_error = norm(uref-usol)/norm(uref);
    total_errorf = norm(fref-efe)/norm(fref);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SP floating prequisites
if size(find(methods==6),1) == 1 || size(find(methods==7),1) == 1
    %% Definition of the operators (quite particular)
    boundaryp = boundary;
    % First problem
    dirichlet1d = [1,1,0;1,2,0;
                   2,1,0;2,2,0];
    [K1d,C1d,nbloq1d,node2c1s,c2node1s] =...
        Krig (nodes,elements,E,nu,order,boundary,dirichlet1d,1);
     
    % Second problem
    K2d = Kinter;
    
    % Management of the rigid modes
    r1 = zeros(2*nnodes,1); r2 = r1; r3p = r1;
    ind = 2:2:2*nnodes;
    r1(ind-1,1) = 1; r2(ind,1) = 1;
    g1 = keepField( r1, 3, boundary );
    g2 = keepField( r2, 3, boundary );
    r3p(ind-1,1) = -nodes(ind/2,2);
    r3p(ind,1) = nodes(ind/2,1);
    g3p = keepField( r3p, 3, boundary );
    % orthonormalize G and ensure G = trace(R)
    g3 = g3p - (g3p'*g1)/(g1'*g1)*g1 - (g3p'*g2)/(g2'*g2)*g2;
    r3 = r3p - (r3p'*g1)/(g1'*g1)*r1 - (r3p'*g2)/(g2'*g2)*r2;
    
    ng1 = norm(g1); ng2 = norm(g2); ng3 = norm(g3);
    g1 = g1/ng1; g2 = g2/ng2; g3 = g3/ng3;
    r1 = r1/ng1; r2 = r2/ng2; r3 = r3/ng3;

    G = [g1,g2,g3];
    R = [r1,r2,r3];
    K2d = [ K2d, R ; R', zeros(size(R,2)) ];
    nbloq2d = size(R,2);
    
    %% Primal precond problems
    % First problem
    dirichlet1 = [2,1,0;2,2,0;
                  1,1,0;1,2,0
                  3,1,0;3,2,0];
    [K1,C1,nbloq1,node2c1,c2node1] =...
        Krig (nodes,elements,E,nu,order,boundaryp,dirichlet1,1);
    % Second problem
    dirichlet2 = [3,1,0;3,2,0];
    [K2,C2,nbloq2,node2c2,c2node2] =...
        Krig (nodes,elements,E,nu,order,boundaryp,dirichlet2,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dual SP with Orthodir and rigid modes (niter = )
if find(methods==6)
    niter   = 10;
    mu      = 0.0/E;      % Regularization parameter
    precond = 0;      % use a primal precond ?
    indexpp = [indexxy ; 2*nnodes+1 ; 2*nnodes+2 ; 2*nnodes+3]; % indices of the vector to choose
    
    % Init
    Itere    = zeros( 2*nnodes+nbloq2d, 1 );
    p        = zeros( 2*nnodes+nbloq2d, niter+1 );
    q        = zeros( 2*nnodes+nbloq2d, niter+1 );
    Delta    = zeros( niter+1, 1 );
    Res      = zeros( 2*nnodes+nbloq2d, niter+1 );
    Zed      = zeros( 2*nnodes+nbloq2d, niter+1 );
    error    = zeros(niter+1,1);
    residual = zeros(niter+1,1);
    regulari = zeros(niter+1,1);
    
    %% Compute Rhs :
    % Solve 1
    f1 = dirichletRhs2( urefb, 1, c2node1s, boundaryp, nnodes );
    uin1 = K1d\f1;
    u1i = uin1(1:2*nnodes,1);
    u1 = keepField( u1i, 3, boundaryp );
    % Solve 2
    f = [keepField( frefb, 2, boundaryp ) ; zeros(nbloq2d,1) ];
    uin2 = K2d\f;
    u2i = uin2(1:2*nnodes,1);
    u2 = keepField( u2i, 3, boundaryp );
    b = [u2-u1 ; R'*f(1:2*nnodes,1)] ;
    
    %% Init the Iterate
    Itere(1:2*nnodes,1) = -G*R'*f(1:2*nnodes,1);
    
    %% Perform A x0 :
    % Solve 1
    f1 = [Itere(1:2*nnodes); zeros(nbloq1d,1)];
    uin1 = K1d\f1;
    u1i = uin1(1:2*nnodes,1);
    u1 = keepField( u1i, 3, boundaryp );
    % Solve 2
    f2 = [Itere(1:2*nnodes); zeros(nbloq2d,1)];
    uin2 = K2d\f2;
    u2i = uin2(1:2*nnodes,1);
    u2 = keepField( u2i, 3, boundaryp );
    alpha = Itere(end-nbloq2d+1:end,1);
    % Regularization term
    Nu = regul(Itere(1:2*nnodes), nodes, boundaryp, 3);
    %
    Axz = [mu*Nu+u1-u2 - G*alpha ;...
           -G'*Itere(1:2*nnodes)];
    %
    Res(:,1) = b - Axz;
    %%%%
    %% Preconditionning
    if precond == 1
        % Solve 1
        f1 = dirichletRhs(Res(:,1)/2, 3, C1, boundaryp);
        uin1 = K1\f1;
        lagr1 = uin1(2*nnodes+1:end,1);
        lamb1 = lagr2forces( lagr1, C1, 3, boundaryp );
        % Solve 2
        f2 = dirichletRhs(Res(:,1)/2, 3, C2, boundaryp);
        uin2 = K2\f2;
        lagr2 = uin2(2*nnodes+1:end,1);
        lamb2 = lagr2forces( lagr2, C2, 3, boundaryp );
        %
        Zed(:,1) = lamb1/2-lamb2/2;        
    else
        Zed(:,1) = Res(:,1);
    end
    p(:,1) = Zed(:,1);
    
    residual(1) = 1;%norm(Res(indexxy,1));
    error(1)    = norm( Itere(indexxy) - fref(indexxy )) / norm(fref(indexxy));
    regulari(1) = sqrt(Itere'*( regul(Itere, nodes, boundary, 3)));

    %% Perform Q1 = A P1 :
    % Solve 1
    f1 = [p(1:2*nnodes,1); zeros(nbloq1d,1)];
    uin1 = K1d\f1;
    u1i = uin1(1:2*nnodes,1);
    u1 = keepField( u1i, 3, boundaryp );
    % Solve 2
    f2 = [p(1:2*nnodes,1); zeros(nbloq2d,1)];
    uin2 = K2d\f2;
    u2i = uin2(1:2*nnodes,1);
    u2 = keepField( u2i, 3, boundaryp );
    alpha = p(end-nbloq2d+1:end,1);
    % Regularization term
    Nu = regul(p(1:2*nnodes,1), nodes, boundaryp, 3);
    %
    q(:,1) = [mu*Nu+u1-u2 - G*alpha ;...
              -G'*p(1:2*nnodes,1)];
    %%%%
    
    for iter = 1:niter

        Delta(iter,1) = norm(q(indexpp,iter))^2;
        gammai        = q(indexpp,iter)'*Res(indexpp,iter);
        alphai        = gammai/Delta(iter,1);

        Itere         = Itere + p(:,iter)*alphai;
        Res(:,iter+1) = Res(:,iter) - q(:,iter)*alphai;

        if precond == 1
            % Solve 1
            f1 = dirichletRhs(Res(:,iter+1)/2, 3, C1, boundaryp);
            uin1 = K1\f1;
            lagr1 = uin1(2*nnodes+1:end,1);
            lamb1 = lagr2forces( lagr1, C1, 3, boundaryp );
            % Solve 2
            f2 = dirichletRhs(Res(:,iter+1)/2, 3, C2, boundaryp);
            uin2 = K2\f2;
            lagr2 = uin2(2*nnodes+1:end,1);
            lamb2 = lagr2forces( lagr2, C2, 3, boundaryp );
            %
            Zed(:,iter+1) = lamb1/2-lamb2/2; 
        else
            Zed(:,iter+1) = Res(:,iter+1);
        end

        residual(iter+1) = norm(Res(indexpp,iter+1)) / norm(Res(indexpp,1));
        error(iter+1)    = norm( Itere(indexxy) - fref(indexxy)) / norm(fref(indexxy));
        regulari(iter+1) = sqrt(Itere'*( regul(Itere, nodes, boundary, 3)));

        %% Perform Ari = A*Res
        % Solve 1
        f1 = [Zed(1:2*nnodes,iter+1); zeros(nbloq1d,1)];
        uin1 = K1d\f1;
        u1i = uin1(1:2*nnodes,1);
        u1 = keepField( u1i, 3, boundaryp );
        % Solve 2
        f2 = [Zed(1:2*nnodes,iter+1); zeros(nbloq2d,1)];
        uin2 = K2d\f2;
        u2i = uin2(1:2*nnodes,1);
        u2 = keepField( u2i, 3, boundaryp );
        alpha = Zed(end-nbloq2d+1:end,iter+1);
        % Regularization term
        Nu = regul(Zed(1:2*nnodes,iter+1), nodes, boundaryp, 3);
        %
        Ari = [mu*Nu+u1-u2 - G*alpha ;...
               -G'*p(1:2*nnodes,1)];
        %%%%

        %% Orthogonalization
        p(:,iter+1) = Zed(:,iter+1);
        q(:,iter+1) = Ari;

        for jter=1:iter
            phiij  = q(indexpp,jter)'*Ari(indexpp);
            betaij = phiij/Delta(jter,1);
            p(:,iter+1) = p(:,iter+1) - p(:,jter) * betaij;
            q(:,iter+1) = q(:,iter+1) - q(:,jter) * betaij;
        end
    end
    
    figure
    hold on
    set(gca, 'fontsize', 15);
    plot(log10(error),'Color','blue')
    plot(log10(residual),'Color','red')
    legend('error (log)','residual (log)')
    % L-curve
    figure
    loglog(residual,regulari);
    
    %% Final problem
    K = K2d;
    usoli = K \ Itere;%([Itere(1:2*nnodes);Itere(2*nnodes+1:end)]);
    usol = usoli(1:2*nnodes,1);

    % Compute stress :
    sigma = stress(usol,E,nu,nodes,elements,order,1,ntoelem,1);
    plotGMSH({usol,'U_vect';sigma,'stress'}, elements, nodes, 'solution');
    
    % Plot displacement on the interface :
    efe = Kinter*usol;
    figure
    hold on
    set(gca, 'fontsize', 15);
    plot(usol(index,1));
    plot(usol(index-1,1), 'Color', 'red');
    xlabel('angle(rad)')
    figure
    hold on
    set(gca, 'fontsize', 15);
    plot(efe(index,1));
    plot(efe(index-1,1), 'Color', 'red');
    xlabel('angle(rad)')
    
    total_error = norm(uref-usol)/norm(uref);
    total_errorf = norm(fref-efe)/norm(fref);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dual SP with Orthodir (niter = 80), floating domain and constraint Krylov
if find(methods==7)
    niter   = 90;
    mu      = 0.0/E;      % Regularization parameter
    precond = 0;      % use a primal precond ?
    
    P = eye(2*nnodes) - G*G';  % Projector
    
    % Init
    p        = zeros( 2*nnodes, niter+1 );
    q        = zeros( 2*nnodes, niter+1 );
    Delta    = zeros( niter+1, 1 );
    Res      = zeros( 2*nnodes, niter+1 );
    Zed      = zeros( 2*nnodes, niter+1 );
    error    = zeros(niter+1,1);
    residual = zeros(niter+1,1);
    regulari = zeros(niter+1,1);
    
    %% Compute Rhs :
    % Solve 1
    f1 = dirichletRhs2( urefb, 1, c2node1s, boundaryp, nnodes );
    uin1 = K1d\f1;
    u1i = uin1(1:2*nnodes,1);
    u1 = keepField( u1i, 3, boundaryp );
    % Solve 2 (ZERO)
    f = [keepField( frefb, 2, boundaryp ) ; zeros(nbloq2d,1)];
    uin2 = K2d\f;
    u2i = uin2(1:2*nnodes,1);
    u2 = keepField( u2i, 3, boundaryp );
    b = u2-u1;
    %%%%
    
    Itere(1:2*nnodes,1) = -G*R'*f(1:2*nnodes);
    
    %% Perform A x0 :
    % Solve 1
    f1 = [Itere; zeros(nbloq1d,1)];
    uin1 = K1d\f1;
    u1i = uin1(1:2*nnodes,1);
    u1 = keepField( u1i, 3, boundaryp );
    % Solve 2
    f2 = [Itere; zeros(nbloq2d,1)];
    uin2 = K2d\f2;
    u2i = uin2(1:2*nnodes,1);
    u2 = keepField( u2i, 3, boundaryp );
    % Regularization term
    Nu = regul(Itere, nodes, boundaryp, 3);
    %
    Axz = mu*Nu+u1-u2;
    %
    Res(:,1) = P'*(b - Axz);
    
    %% Preconditionning
    if precond == 1
        % Solve 1
        f1 = dirichletRhs(Res(:,1)/2, 3, C1, boundaryp);
        uin1 = K1\f1;
        lagr1 = uin1(2*nnodes+1:end,1);
        lamb1 = lagr2forces( lagr1, C1, 3, boundaryp );
        % Solve 2
        f2 = dirichletRhs(Res(:,1)/2, 3, C2, boundaryp);
        uin2 = K2\f2;
        lagr2 = uin2(2*nnodes+1:end,1);
        lamb2 = lagr2forces( lagr2, C2, 3, boundaryp );
        %
        Zed(:,1) = P*(lamb1/2-lamb2/2);        
    else
        Zed(:,1) = P*Res(:,1);
    end
    p(:,1) = Zed(:,1);
    
    residual(1) = 1;%norm(Res(indexxy,1));
    error(1)    = norm( Itere(indexxy) - fref(indexxy )) / norm(fref(indexxy));
    regulari(1) = sqrt(Itere'*( regul(Itere, nodes, boundary, 3)));

    %% Perform Q1 = P A P1 :
    % Solve 1
    f1 = [p(:,1); zeros(nbloq1d,1)];
    uin1 = K1d\f1;
    u1i = uin1(1:2*nnodes,1);
    u1 = keepField( u1i, 3, boundaryp );
    % Solve 2
    f2 = [p(:,1); zeros(nbloq2d,1)];
    uin2 = K2d\f2;
    u2i = uin2(1:2*nnodes,1);
    u2 = keepField( u2i, 3, boundaryp );
    % Regularization term
    Nu = regul(p(:,1), nodes, boundaryp, 3);
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
            f1 = dirichletRhs(Res(:,iter+1)/2, 3, C1, boundaryp);
            uin1 = K1\f1;
            lagr1 = uin1(2*nnodes+1:end,1);
            lamb1 = lagr2forces( lagr1, C1, 3, boundaryp );
            % Solve 2
            f2 = dirichletRhs(Res(:,iter+1)/2, 3, C2, boundaryp);
            uin2 = K2\f2;
            lagr2 = uin2(2*nnodes+1:end,1);
            lamb2 = lagr2forces( lagr2, C2, 3, boundaryp );
            %
            Zed(:,iter+1) = P*(lamb1/2-lamb2/2); 
        else
            Zed(:,iter+1) = P*Res(:,iter+1);
        end

        residual(iter+1) = norm(Res(indexxy,iter+1)) / norm(Res(indexxy,1));
        error(iter+1)    = norm( Itere(indexxy) - fref(indexxy)) / norm(fref(indexxy));
        regulari(iter+1) = sqrt(Itere'*( regul(Itere, nodes, boundary, 3)));

        %% Perform Ari = P*A*Res
        % Solve 1
        f1 = [Zed(:,iter+1); zeros(nbloq1d,1)];
        uin1 = K1d\f1;
        u1i = uin1(1:2*nnodes,1);
        u1 = keepField( u1i, 3, boundaryp );
        % Solve 2
        f2 = [Zed(:,iter+1); zeros(nbloq2d,1)];
        uin2 = K2d\f2;
        u2i = uin2(1:2*nnodes,1);
        u2 = keepField( u2i, 3, boundaryp );
        % Regularization term
        Nu = regul(Zed(:,iter+1), nodes, boundaryp, 3);
        %
        Ari = P'*(mu*Nu+u1-u2);
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
    set(gca, 'fontsize', 15);
    plot(log10(error),'Color','blue')
    plot(log10(residual),'Color','red')
    legend('error (log)','residual (log)')
    % L-curve
    figure
    loglog(residual,regulari);
    
    %% Final problem
    dirichlet = [2,1,0;2,2,0];
    [K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet,1);
    %fdir1 = dirichletRhs(urefb, 1, C, boundary);
    fdir2 = dirichletRhs(urefb, 2, C, boundary);
    f1 = [Itere; zeros(nbloq,1)];
    usoli = K \ ( fdir2 + f1 );
    usol = usoli(1:2*nnodes,1);

    % Compute stress :
    sigma = stress(usol,E,nu,nodes,elements,order,1,ntoelem,1);
    plotGMSH({usol,'U_vect';sigma,'stress'}, elements, nodes, 'solution');
    
    % Plot displacement on the interface :
    efe = Kinter*usol;
    figure
    hold on
    set(gca, 'fontsize', 15);
    plot(thetax,usol(index,1));
    plot(thetax,usol(index-1,1), 'Color', 'red');
    xlabel('angle(rad)')
    figure
    hold on
    set(gca, 'fontsize', 15);
    plot(thetax,efe(index,1));
    plot(thetax,efe(index-1,1), 'Color', 'red');
    xlabel('angle(rad)')
    
    total_error = norm(uref-usol)/norm(uref);
    total_errorf = norm(fref-efe)/norm(fref);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Evanescent regularization method
if find(methods==8)
   niter = 100;
   mu = .1;
   
   %% Computation of the inner stiffness
   dirichlet1 = [2,1,0;2,2,0];
   [K1,C1,nbloq1,node2c1,c2node1] = Krig (nodes,elements,E,nu,order,boundary,dirichlet1,1);
   neumann1   = []; % There is no alone Neumann
   f1 = loading( nbloq, nodes, boundary, neumann1 );

   % map of the nodes
   b2node13 = [b2node1;b2node3];
   bbound = [2*b2node13-1;2*b2node13];
   bbzero = [2*b2node2-1;2*b2node2];
   nbound = size(bbound,1);
   
   %% Schur operator
   [ S, b, map ] = schurComp2( Kinter, f1(1:2*nnodes), bbound, bbzero );

   error    = zeros(niter,1);
   residual = zeros(niter,1);
   regulari = zeros(niter,1);

   %% Mass matrices
   Mr  = bMass_mat(nodes, boundary, [1]);
   Mrt = 1/E*Mr;  % Ideally 1/EL
   M   = bMass_mat(nodes, boundary, [3;1]);
   Mt  = 1/E*M;
   Mm  = bMass_mat(nodes, boundary, 3);

   % Extract coords
   Mr  = Mr(bbound, bbound);
   Mrt = Mrt(bbound, bbound);
   M   = M(bbound, bbound);
   Mt  = Mt(bbound, bbound);
   Mm  = Mm(bbound, bbound);

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% Evanescent regularization method
   Itere  = zeros( 2*nnodes, 1 );
   Iteref = zeros( 2*nnodes, 1 );

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
              Mrt*frefb(bbound) + mu*Mt*Iteref(bbound)
              b]; % Don't forget to add Kinter*uimp if needed

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
   dirichlet = [3,1,0;3,2,0;
                1,1,0;1,2,0;
                2,1,0;2,2,0];
   neumann   = [];
   [K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet,1);
   fdir1 = dirichletRhs(urefb, 1, C, boundary);
   fdir2 = dirichletRhs(urefb, 2, C, boundary);
   fdir3 = dirichletRhs(Itere, 3, C, boundary);
   usoli = K \ assembleDirichlet( [fdir1+fdir3,fdir2] );

   usol = usoli(1:2*nnodes,1);
   fsol = Kinter*usol;

   % Plot displacement on the interface :
   efe = Kinter*usol;
   figure
   hold on
   set(gca, 'fontsize', 15);
   plot(thetax,usol(index,1));
   plot(thetax,usol(index-1,1), 'Color', 'red');
   xlabel('angle(rad)')
   figure
   hold on
   set(gca, 'fontsize', 15);
   plot(thetax,efe(index,1));
   plot(thetax,efe(index-1,1), 'Color', 'red');
   xlabel('angle(rad)')
    
   total_error = norm(uref-usol)/norm(uref);
   total_errorf = norm(fref-efe)/norm(fref);

   % Compute stress :
   sigma = stress(usol,E,nu,nodes,elements,order,1,ntoelem);
   plotGMSH({usol,'U_vect';sigma,'stress'}, elements, nodes, 'solution');
end