%10/05/2016
%Algo KMF orthodir pour un problème non-linéire (hyperélasticité)

close all;
clear all;

addpath(genpath('./tools'))
addpath(genpath('./regu'))

% Parameters
E       = 70000;  % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
alpha   = 1e10;   % Compressibility nonliner parameter
fscalar = 1;      % N.mm-1 : Loading on the plate
niter   = 6;      % nb of iterations of the Orthodir algo
nnewt   = 10;      % nb of iterations of the Newton algo
br      = 0.;     % noise
mu      = 0.;      % Regularization parameter (use the second reg method)
actual  = 1;      % Wether to actualize the tangent matrix
% Note : if actual==0, a relaxation is applied in order to ensure
% convergence

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
urefb = ( 1 + br*randn(2*nnodes,1) ) .* uref;
lagrb = ( 1 + br*randn(nbloq,1) ) .* lagr;
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
[K1,C1,nbloq1,node2c1,c2node1] = Krig (nodes,elements,E,nu,order,boundary,dirichlet1);

% ND problem
dirichlet2 = [4,1,0;4,2,0;
              1,1,0;1,2,0;
              2,1,0;2,2,0];

[K2,C2,nbloq2,node2c2,c2node2] = Krig (nodes,elements,E,nu,order,boundary,dirichlet2);

tic
u2 = uref-uref;

residual = ones(nnewt+1,1);
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

for i=1:nnewt % quasi-Newton scheme

    if actual == 1
        [K1t,~,~] = KrigHyp (nodes,elements,E,nu,order,boundary,...
            dirichlet1,alpha, u2, 1);
        [K2t,~,~] = KrigHyp (nodes,elements,E,nu,order,boundary,...
            dirichlet2,alpha, u2, 1);
    else
        K1t = K1;
        K2t = K2;
    end
    
    %% Orthodir for the problem : Res = (1 - D0oS0) x
    
    % initialization
    Itere     = zeros( 2*nnodes, niter+1 );
    p         = zeros( 2*nnodes, niter+1 );
    q         = zeros( 2*nnodes, niter+1 );
    Delta     = zeros( niter+1, 1 );
    Resi      = zeros( 2*nnodes, niter+1 );
    residuali = zeros(niter+1,1);
    regulari  = zeros(niter+1,1);

    %% Perform A x0 :
    % Solve DN
    fdir = dirichletRhs(Itere(:,1), 3, C1, boundary);
    f1 = fdir;
    uin1 = K1t\f1;
    lagr1 = uin1(2*nnodes+1:end,1);
    % Solve ND
    fr1 = lagr2forces( lagr1, C1, 3, boundary );
    fr = [ fr1; zeros(size(C2,2),1) ];
    f2 = fr;
    uin2 = K2t\f2;
    %
    Nu = regul(Itere(:,1), nodes, boundary, 3);
    atimesItere = mu*Nu + uin2(1:2*nnodes,1) - Itere(:,1);
    %%%%

    %% Write Rhs (for once, it's explicit) :
    b = Res;%keepField( Res, 3, boundary);
    %%%%
    %%
    Resi(:,1) = b - atimesItere;
    p(:,1) = Resi(:,1);

    residuali(1) = sqrt( myps( Resi(:,1), Resi(:,1), Kinter, boundary, M, nodes ) );
    regulari(1)  = sqrt( Itere(:,1)'*regul(Itere(:,1), nodes, boundary, 3) );

    %% Perform Q1 = A P1 :
    % Solve DN
    fdir = dirichletRhs(p(:,1), 3, C1, boundary);
    f1 = fdir;
    uin1 = K1t\f1;
    lagr1 = uin1(2*nnodes+1:end,1);
    % Solve ND
    fr1 = lagr2forces( lagr1, C1, 3, boundary );
    fr = [ fr1; zeros(size(C2,2),1) ]; % Give to fr the right size
    f2 = fr;
    uin2 = K2t\f2;
    %
    Nu = regul(p(:,1), nodes, boundary, 3);
    q(:,1) = mu*Nu - uin2(1:2*nnodes,1) + p(:,1) ;
    %%%%
    %%
    for iter = 1:niter

        Delta(iter,1)   = myps( q(:,iter), q(:,iter), Kinter, boundary, M, nodes ); %q(:,iter)'*q(:,iter);
        gammai          = myps( q(:,iter), Resi(:,iter), Kinter, boundary, M, nodes ); %q(:,iter)'*Res;
        alphai          = gammai/Delta(iter,1);

        Itere(:,iter+1)  = Itere(:,iter) + p(:,iter)*alphai;
        Resi(:,iter+1)   = Resi(:,iter) - q(:,iter)*alphai;

        residuali(iter+1) = sqrt( myps( Resi(:,iter+1), Resi(:,iter+1), Kinter, boundary, M, nodes ) );
        regulari(iter+1)  = sqrt( Itere(:,iter+1)'*regul(Itere(:,iter+1), nodes, boundary, 3) );

        %% Perform Ari = A*Res
        % Solve DN
        fdir = dirichletRhs(Resi(:,iter+1), 3, C1, boundary);
        f1 = fdir;
        uin1 = K1t\f1;
        lagr1 = uin1(2*nnodes+1:end,1);
        % Solve ND
        fr1 = lagr2forces( lagr1, C1, 3, boundary );
        fr = [ fr1; zeros(size(C2,2),1) ]; % Give to fr the right size
        f2 = fr;
        uin2 = K2t\f2;
        %
        Nu = regul(Resi(:,iter+1), nodes, boundary, 3);
        Ari = mu*Nu - uin2(1:2*nnodes,1) + Resi(:,iter+1) ;

        %% Orthogonalization
        p(:,iter+1) = Resi(:,iter+1);
        q(:,iter+1) = Ari;

        for jter=1:iter
            phiij  = myps( q(:,jter), Ari, Kinter, boundary, M, nodes );
            betaij = phiij/Delta(jter,1);
            p(:,iter+1) = p(:,iter+1) - p(:,jter) * betaij;
            q(:,iter+1) = q(:,iter+1) - q(:,jter) * betaij;
        end
    end
    %% Automatic stop (regularization tools) :
    [stopHere,~,~] = l_corner(residuali,regulari,1:1:niter+1);
%     hold on;
%     plot(log10(residuali),'Color','red')
%     plot(log10(regulari))
%     figure
%     loglog(residuali,regulari)
%     figure

    %% Actualize the estimation
    if actual == 1
        u2 = u2 + Itere(:,stopHere);
    else
        u2 = u2 + 0.5*Itere(:,stopHere); % relaxation in order to ensure convergence
    end

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
plot(log10(residual/residual(1)),'Color','red')
plot(log10(error))
error_nk = norm(u2-uref)/norm(uref);
toc

% Compute stress :
sigma = stress(u2,E,nu,nodes,elements,order,1,ntoelem);
% Output :
plotGMSH({u2,'U_vect';sigma,'stress'}, elements, nodes(:,[1,2]), 'field2');