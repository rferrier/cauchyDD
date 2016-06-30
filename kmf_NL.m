%20/04/2016
%Algo KMF pour un problème non-linéire (hyperélasticité)

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 70000;  % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
alpha   = 1e10;   % Compressibility nonliner parameter
fscalar = 1;      % N.mm-1 : Loading on the plate
niter   = 10;
br      = 0.;     % noise
kn      = 1;      % use the KMF(Newton) method
nk      = 0;      % use the Newton(KMF) method
relax   = 0;

noises = load('./noises/noise2.mat'); % Particular noise vector
noise  = noises.bruit1;
% noise = randn(2*nnodes,1);

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

% Then, build the stiffness matrix (we'll need it later) :
[K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
Kinter = K(1:2*nnodes, 1:2*nnodes);
M      = mass_mat(nodes, elements);
% 
% The right hand side :
f = loading(nbloq,nodes,boundary,neumann);

% Solve the non-linear problem
uin = solveHyp(nodes,elements,E,nu,order,boundary,dirichlet,...
    alpha,1e-5,50,ntoelem,0,f);
disp('reference problem converged')

% Extract displacement and Lagrange multiplicators :
uref = uin(1:2*nnodes,1);
lagr = uin(2*nnodes+1:end,1);
urefb = ( 1 + br*noise ) .* uref;
%lagrb = ( 1 + br*randn(nbloq,1) ) .* lagr;
fref = f( 1:2*nnodes,1 ); % Reaction forces

% Compute stress :
sigma = stressHyp(uref,E,nu,nodes,elements,order,1,ntoelem,alpha,1);

% Output :
plotGMSH({uref,'U_vect';sigma,'stress'}, elements, nodes(:,[1,2]), 'reference');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DN problem
dirichlet1 = [4,1,0;4,2,0;
              3,1,0;3,2,0];
neumann1   = [1,2,fscalar;
              2,1,fscalar];
[K1,C1,nbloq1,node2c1,c2node1] = Krig (nodes,elements,E,nu,order,boundary,dirichlet1,1);

% ND problem
dirichlet2 = [4,1,0;4,2,0;
              1,1,0;1,2,0;
              2,1,0;2,2,0];

[K2,C2,nbloq2,node2c2,c2node2] = Krig (nodes,elements,E,nu,order,boundary,dirichlet2,1);

if kn == 1
    tic;
    % init :
    u2    = uref-uref;
    u1    = u2;
    v     = u2;
    theta = ones(niter+1,1); % First relaxation parameter

    error1   = zeros(niter,1);
    error2   = zeros(niter,1);
    residual = zeros(niter,1);
    regulari = zeros(niter,1);

    for iter = 1:niter
        
        u1o = u1;
        
        % Solve DN
        f1in = loading(nbloq1,nodes,boundary,neumann1);
        fdir = dirichletRhs2(v, 3, c2node1, boundary, nnodes );

        % Solve the non-linear problem
        uin1 = solveHyp(nodes,elements,E,nu,order,boundary,dirichlet1,...
            alpha,1e-5,50,ntoelem,fdir,f1in);

        u1 = uin1(1:2*nnodes,1);
        lagr1 = uin1(2*nnodes+1:end,1);
        fri = lagr2forces2( lagr1, c2node1, 3, boundary, nnodes );

        error1(iter) = sqrt( myps(u1-uref,u1-uref,Kinter,boundary,M,nodes)/...
            myps(uref,uref,Kinter,boundary,M,nodes) );

        disp( [ 'DN problem at iteration ', num2str(iter), ' converged' ] )

        u2o = u2;
        
        % Solve ND
        fr = [ fri; zeros(size(C2,2),1) ]; % Give to fr the right size

        fdir1 = dirichletRhs2(urefb, 1, c2node2, boundary, nnodes);
        fdir2 = dirichletRhs2(urefb, 2, c2node2, boundary, nnodes);
        fdir = assembleDirichlet( [fdir1,fdir2] );

        % Solve the non-linear problem
        uin2 = solveHyp(nodes,elements,E,nu,order,boundary,dirichlet2,...
            alpha,1e-5,50,ntoelem,fdir,fr);

        disp( [ 'ND problem at iteration ', num2str(iter), ' converged' ] )

        u2 = uin2(1:2*nnodes,1);

        vo = v;
        v = theta(iter)*u2 + (1-theta(iter))*vo;
    
        if relax == 1 && iter > 1
            e1 = u1-u1o;
            e2 = u2-u2o;
            theta(iter+1) = myps(e1,e1-e2,Kinter,boundary,M,nodes) /...
                myps(e1-e2,e1-e2,Kinter,boundary,M,nodes);
        end

        error2(iter) = sqrt( myps(u2-uref,u2-uref,Kinter,boundary,M,nodes)/...
            myps(uref,uref,Kinter,boundary,M,nodes) );

        residual(iter) = sqrt( myps(u1-u2,u1-u2,Kinter,boundary,M,nodes)/...
                         sqrt(myps(u1,u1,Kinter,boundary,M,nodes)*...
                         myps(u2,u2,Kinter,boundary,M,nodes)) );

        regulari(iter) = sqrt(u2'*regul(u2, nodes, boundary, 3));
    %     Nu = regul(u2, nodes, boundary, 3);
    end

    hold on;
    plot(error1,'Color','black')
    plot(error2,'Color','blue')
    plot(residual,'Color','red')
    legend('error1','error2','residual')
    figure
    % L-curve
    loglog(residual,regulari);
    % figure
    % plot(regulari.*residual)
    error_kn = norm(u2-uref)/norm(uref);
    toc
end

if nk == 1
    tic
    u2 = uref-uref;
    nnewt = 10;
    residual = zeros(nnewt+1,1);
    error    = ones(nnewt+1,1);
    
    % compute the residual
    f1in = loading(nbloq1,nodes,boundary,neumann1);
    fdir = dirichletRhs2(u2, 3, c2node1, boundary, nnodes );
    f1 = f1in + fdir;
    uin1 = solveHyp(nodes,elements,E,nu,order,boundary,dirichlet1,...
        alpha,1e-5,50,ntoelem,fdir,f1in);
    u1 = uin1(1:2*nnodes,1);
    lagr1 = uin1(2*nnodes+1:end,1);
    fri = lagr2forces2( lagr1, c2node1, 3, boundary, nnodes );
    % Solve ND
    fr = [ fri; zeros(size(C2,2),1) ]; % Give to fr the right size
    fdir1 = dirichletRhs2(urefb, 1, c2node2, boundary, nnodes);
    fdir2 = dirichletRhs2(urefb, 2, c2node2, boundary, nnodes);
    f2 = assembleDirichlet( [fdir1,fdir2] );
    uin2 = solveHyp(nodes,elements,E,nu,order,boundary,dirichlet2,...
        alpha,1e-5,50,ntoelem,f2,fr);
    %
    Res = uin2(1:2*nnodes,1) - u2;
    residual(1) = sqrt( myps(Res,Res,Kinter,boundary,M,nodes) );
    disp('residual computation converged')
    
    for i=1:nnewt % Newton scheme
        
        % Richardson for the problem : x = (D0t o S0t - 1) Res
        it = Res-Res;
        % compute the tangent matrixes
        [K1t,~,~] = KrigHyp (nodes,elements,E,nu,order,boundary,...
            dirichlet1,alpha, u2, 1);
        [K2t,~,~] = KrigHyp (nodes,elements,E,nu,order,boundary,...
            dirichlet2,alpha, u2, 1);
        
        for iter=1:niter
            % Solve DN
            fdir = dirichletRhs2( it, 3, c2node1, boundary, nnodes );
            f1 = fdir;
            uin1 = K1t\f1;
            lagr1 = uin1(2*nnodes+1:end,1);
            fri = lagr2forces2( lagr1, c2node1, 3, boundary, nnodes );
            % Solve ND
            fr = [ fri; zeros(size(C2,2),1) ]; % Give to fr the right size
            f2 = fr;
            uin2 = K2t\f2;
            %
            it = uin2(1:2*nnodes,1) + Res;
        end
        % Actualize the estimation
        u2 = u2 + it;
        
        % compute the residual
        f1in = loading(nbloq1,nodes,boundary,neumann1);
        fdir = dirichletRhs2(u2, 3, c2node1, boundary, nnodes );
        f1 = f1in + fdir;
        uin1 = solveHyp(nodes,elements,E,nu,order,boundary,dirichlet1,...
            alpha,1e-5,50,ntoelem,fdir,f1in);
        u1 = uin1(1:2*nnodes,1);
        lagr1 = uin1(2*nnodes+1:end,1);
        fri = lagr2forces2( lagr1, c2node1, 3, boundary, nnodes );
        % Solve ND
        fr = [ fri; zeros(size(C2,2),1) ]; % Give to fr the right size
        fdir1 = dirichletRhs2(urefb, 1, c2node2, boundary, nnodes);
        fdir2 = dirichletRhs2(urefb, 2, c2node2, boundary, nnodes);
        f2 = assembleDirichlet( [fdir1,fdir2] );
        uin2 = solveHyp(nodes,elements,E,nu,order,boundary,dirichlet2,...
            alpha,1e-5,50,ntoelem,f2,fr);
        %
        Res = uin2(1:2*nnodes,1) - u2;
        residual(i+1) = sqrt( myps(Res,Res,Kinter,boundary,M,nodes) );
        error(i+1) = sqrt( myps(u1-uref,u1-uref,Kinter,boundary,M,nodes) ) /...
            sqrt( myps(uref,uref,Kinter,boundary,M,nodes) );
        disp([ 'residual computation at Newton iteration ', ...
            num2str(i), ' converged' ])
        
    end
    hold on
    plot(residual/residual(1),'Color','red')
    plot(error)
    error_nk = norm(u2-uref)/norm(uref);
    toc
end

% Compute stress :
sigma = stress(u2,E,nu,nodes,elements,order,1,ntoelem);
% Output :
plotGMSH({u2,'U_vect';sigma,'stress'}, elements, nodes(:,[1,2]), 'field2');