%21/03/2016
%Algo KMF avec Orthodir et r�gularisation

close all;
clear all;

addpath(genpath('./tools'))
addpath(genpath('./regu'))

% Parameters
E        = 70000;      % MPa : Young modulus
nu       = 0.3;    % Poisson ratio
fscalar  = 1;      % N.mm-1 : Loading on the plate
niter    = 20;
nrej     = 1;       % number of iterations of the noise rejection algo
br       = .01;      % noise
brt      = .0;     % "Translation" noise : constant noise
mu       = 0;%30;      % Regularization parameter
ev0      = 0;      % Wether or not to use the evanescent regularization (Cimetière)
bestiter = 1;      % Wether to choose the best iteration

noises = load('./noises/noise0r.mat'); % Particular noise vector
noise  = noises.bruit1;

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
[ nodes,elements,ntoelem,boundary,order ] = readmesh( 'meshes/plate_reg.msh' );
nnodes = size(nodes,1);

%noise = randn(2*nnodes,1);

[node2t, t2node] = mapBound(3, boundary, nnodes);

% Then, build the stiffness matrix :
[K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
Kinter = K(1:2*nnodes, 1:2*nnodes);
M      = mass_mat(nodes, elements);
% The right hand side :
f = loading(nbloq,nodes,boundary,neumann);

% Solve the problem :
uin = K\f;

% Extract the index of the boundary
[ node2b3, b2node3 ] = mapBound( 3, boundary, nnodes );
index    = 2*b2node3-1;
indextot = [2*b2node3-1 ; 2*b2node3];  % the same, with x and y

% Extract on the suraboundant surfaces
[ node2b1, b2node1 ] = mapBound( 1, boundary, nnodes );
[ node2b2, b2node2 ] = mapBound( 2, boundary, nnodes );
index12 = [2*b2node1; 2*b2node1-1; 2*b2node2; 2*b2node2-1];
index1 = [2*b2node1; 2*b2node1-1];
index2 = [2*b2node2; 2*b2node2-1];

% Extract displacement and Lagrange multiplicators :
uref = uin(1:2*nnodes,1);
lagr = uin(2*nnodes+1:end,1);

urefb = ( 1 + br*noise ) .* uref + brt*uref;
%urefb = uref;
%lagrb = ( 1 + br*randn(nbloq,1) ) .* lagr;

% Post-pro :
sigma = stress(uref,E,nu,nodes,elements,order,1,ntoelem);
ui = reshape(uref,2,[])';  ux = ui(:,1);  uy = ui(:,2);
plotGMSH({ux,'U_x';uy,'U_y';uref,'U_vect';sigma,'stress'}, elements, nodes, 'reference');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition of the operators

% DN problem
dirichlet1 = [4,1,0;
              4,2,0;
              3,1,0;
              3,2,0];
neumann1   = [1,2,fscalar;
              2,1,fscalar];
neumann0   = [];
[K1,C1,nbloq1] = Krig (nodes,elements,E,nu,order,boundary,dirichlet1);

% ND problem
dirichlet2 = [4,1,0;
              4,2,0;
              1,1,0;
              1,2,0;
              2,1,0;
              2,2,0];
neumann2   = [];% [3,1,lagr1; 3,2,lagr1]
                % is managed by lagr2forces
[K2,C2,nbloq2] = Krig (nodes,elements,E,nu,order,boundary,dirichlet2);

error   = zeros(niter+1,1);
residual = zeros(niter+1,1);
regulari = zeros(niter+1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computation of an explicit form of the operator (1 - D0 o S0 + mu*N)

Ax = zeros(2*nnodes,size(indextot,1));

for i=1:size(indextot,1)
    x = zeros( 2*nnodes, 1 );
    x(indextot(i)) = 1;
    % Solve DN
    f1in = loading(nbloq1,nodes,boundary,neumann0);
    fdir = dirichletRhs(x, 3, C1, boundary);
    f1 = f1in + fdir;
    uin1 = K1\f1;
    u1 = uin1(1:2*nnodes,1);
    lagr1 = uin1(2*nnodes+1:end,1);
    % Solve ND
    fr = lagr2forces( lagr1, C1, 3, boundary );
    fr = [ fr; zeros(size(C2,2),1) ]; 
    f2 = fr;
    uin2 = K2\f2;

    % Regularization term
    Nu = regul(x, nodes, boundary, 3);
    Ax(:,i) = mu*Nu + x - uin2(1:2*nnodes,1);
end

A = Ax(indextot,:); % extract the useful components of A
[U,s,V] = csvd(A);
% cond(A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iterations of the noise rejection algo
total_error = zeros(nrej,1);
for brit = 1:nrej
    %% ORTHODIR for the problem : (1 - D0 o S0 + mu*N) x = DoS(0), with N, a 
    % regulalization term
    Itere = zeros( 2*nnodes, niter+1 );
    p     = zeros( 2*nnodes, niter+1 );
    q     = zeros( 2*nnodes, niter+1 );
    Delta = zeros( niter+1, 1 );
    Res   = zeros( 2*nnodes, niter+1 );

    %% Perform A x0 :
    % Solve DN
    f1in = loading(nbloq1,nodes,boundary,neumann0);
    fdir = dirichletRhs(Itere(:,1), 3, C1, boundary);
    f1 = f1in + fdir;

    uin1 = K1\f1;
    u1 = uin1(1:2*nnodes,1);
    lagr1 = uin1(2*nnodes+1:end,1);

    % Solve ND
    fr = lagr2forces( lagr1, C1, 3, boundary );
    fr = [ fr; zeros(size(C2,2),1) ]; % Give to fr the right size (because Lagrange)
    %fdir1 = dirichletRhs(uref, 1, C2, boundary);
    %fdir2 = dirichletRhs(uref, 2, C2, boundary);
    f2 = fr;% + fdir1 + fdir2;

    uin2 = K2\f2;

    % Regularization term
    Nu = regul(Itere(:,1), nodes, boundary, 3);

    atimesItere = mu*Nu + Itere(:,1) - uin2(1:2*nnodes,1);
    %%%%

    %% Write Rhs :
    % Solve DN
    f1in = loading(nbloq1,nodes,boundary,neumann1);
    fdir = dirichletRhs(zeros( 2*nnodes,1 ), 3, C1, boundary);
    f1 = f1in + fdir;

    uin1 = K1\f1;
    u1 = uin1(1:2*nnodes,1);
    lagr1 = uin1(2*nnodes+1:end,1);

    % Solve ND
    fr = lagr2forces( lagr1, C1, 3, boundary );
    fr = [ fr; zeros(size(C2,2),1) ]; % Give to fr the right size
    fdir1 = dirichletRhs(urefb, 1, C2, boundary);
    fdir2 = dirichletRhs(urefb, 2, C2, boundary);
    f2 = fr + assembleDirichlet( [fdir1,fdir2] );

    uin2 = K2\f2;
    b = uin2(1:2*nnodes,1);
    %%%%

    %% Picard Stuff
    bu = b(indextot);  % extract components
    set(gca, 'fontsize', 15);
    picard(U,s,bu);
    figure
    
    %%
    Res(:,1) = b - atimesItere;
    p(:,1) = Res(:,1);

    residual(1) = sqrt( myps( Res(:,1), Res(:,1), Kinter, boundary, M, nodes ) );
    error(1)    = sqrt( myps( Itere(:,1) - uref, Itere(:,1) - uref, Kinter, boundary, M, nodes )...
                     / myps( uref, uref, Kinter, boundary, M, nodes ) );
    regulari(1) = sqrt( Itere(:,1)'*regul(Itere(:,1), nodes, boundary, 3) );

    %% Perform Q1 = A P1 :
    % Solve DN
    f1in = loading(nbloq1,nodes,boundary,neumann0);
    fdir = dirichletRhs(p(:,1), 3, C1, boundary);
    f1 = f1in + fdir;

    uin1 = K1\f1;
    u1 = uin1(1:2*nnodes,1);
    lagr1 = uin1(2*nnodes+1:end,1);

    % Solve ND
    fr1 = lagr2forces( lagr1, C1, 3, boundary );
    fr = [ fr1; zeros(size(C2,2),1) ]; % Give to fr the right size
    %fdir1 = dirichletRhs(uref, 1, C2, boundary); Those guys have nothing to do here : we are in a linear problem
    %fdir2 = dirichletRhs(uref, 2, C2, boundary);
    f2 = fr;% + fdir1 + fdir2;

    uin2 = K2\f2;

    % Regularization term
    Nu = regul(p(:,1), nodes, boundary, 3);

    q(:,1) = mu*Nu + p(:,1) - uin2(1:2*nnodes,1);
    %%%%
    %%
    for iter = 1:niter

        Delta(iter,1)   = myps( q(:,iter), q(:,iter), Kinter, boundary, M, nodes ); %q(:,iter)'*q(:,iter);
        gammai          = myps( q(:,iter), Res(:,iter), Kinter, boundary, M, nodes ); %q(:,iter)'*Res;
        alphai          = gammai/Delta(iter,1);

        Itere(:,iter+1) = Itere(:,iter) + p(:,iter)*alphai;
        Res(:,iter+1)   = Res(:,iter) - q(:,iter)*alphai;

        residual(iter+1) = sqrt( myps( Res(:,iter+1), Res(:,iter+1), Kinter, boundary, M, nodes ) );
        error(iter+1)    = sqrt( myps( Itere(:,iter+1) - uref, Itere(:,iter+1) - uref, Kinter, boundary, M, nodes )...
                         / myps( uref, uref, Kinter, boundary, M, nodes ) );
        regulari(iter+1) = sqrt( Itere(:,iter+1)'*regul(Itere(:,iter+1), nodes, boundary, 3) );

        %% Perform Ari = A*Res
        % Solve DN
        f1in = loading(nbloq1,nodes,boundary,neumann0);
        fdir = dirichletRhs(Res(:,iter+1), 3, C1, boundary);
        f1 = f1in + fdir;

        uin1 = K1\f1;
        u1 = uin1(1:2*nnodes,1);
        lagr1 = uin1(2*nnodes+1:end,1);

        % Solve ND
        fr1 = lagr2forces( lagr1, C1, 3, boundary );
        fr = [ fr1; zeros(size(C2,2),1) ]; % Give to fr the right size
        f2 = fr;

        uin2 = K2\f2;

        % Regularization term
        if ev0 == 1
           Nu = regul(Res(:,iter+1) - Res(:,iter), nodes, boundary, 3);
        else
           Nu = regul(Res(:,iter+1), nodes, boundary, 3);
        end

        Ari = mu*Nu + Res(:,iter+1) - uin2(1:2*nnodes,1);

        %% Orthogonalization
        p(:,iter+1) = Res(:,iter+1);
        q(:,iter+1) = Ari;

        for jter=1:iter
            phiij  = myps( q(:,jter), Ari, Kinter, boundary, M, nodes ); %q(:,jter)'*Ari;
            betaij = phiij/Delta(jter,1);
            p(:,iter+1) = p(:,iter+1) - p(:,jter) * betaij;
            q(:,iter+1) = q(:,iter+1) - q(:,jter) * betaij;
            %myps( q(:,iter+1), q(:,jter), Kinter, boundary, M, nodes )
        end

    end

    hold on
    set(gca, 'fontsize', 15);
    set(gca,'ylim',[-11 1])
    plot(log10(error(2:end)),'Color','blue')
    plot(log10(residual(2:end)/residual(1)),'Color','red')
    legend('error (log)','residual (log)')
    figure
    
    % L-curve
    loglog(residual,regulari);
    figure

    %% Automatic stop (regularization tools, broken on Octave) :
    %[stopHere,~,~] = l_corner(residual,regulari,1:1:niter+1);
    if bestiter == 1
       stopHere = findCorner(residual(2:end),regulari(2:end)) + 1;
    else
       stopHere = niter+1;
    end
    
    % Plot chosen Itere
    hold on;
    set(gca, 'fontsize', 15);
    set(gca,'ylim',[-3e-5 3e-5])
    plot(uref(index));
    plot(Itere(index,stopHere),'Color','red');
    figure
    %%%%
    %% Final problem : compute u
    % DN problem
    dirichlet = [4,1,0;4,2,0;
                 3,1,0;3,2,0];
    neumann   = [];
    [K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
    f1in = loading(nbloq,nodes,boundary,neumann1);
    fdir3 = dirichletRhs(Itere(:,stopHere), 3, C, boundary);
    usoli = K \ ( assembleDirichlet(fdir3) + f1in );
    usol = usoli(1:2*nnodes,1);

    %% Boundary error => mesh adaptation

%     hold on
%     plot( abs( urefb(index12)-usol(index12) ) / norm( urefb(index12) ) )
%     plot( abs( urefb(index12)-uref(index12) ) / norm( uref(index12) ), 'Color', 'red' )
%     plot( ( uref(index12)-usol(index12) ) / norm( uref(index12) ) )
%     plot( ( uref(index12)-urefb(index12) ) / norm( uref(index12) ), 'Color', 'red' )
%     legend('�cart entre la solution et le champ de r�f�rence','bruit')
%     figure
    
    urefb(index12) = usol(index12); % remplace urefb by usol, that is a smoothed version
    total_error(brit) = norm(usol-uref)/norm(uref);
end

% plot(total_error)

% Post-pro :
sigma = stress(usol,E,nu,nodes,elements,order,1,ntoelem);
ui = reshape(usol,2,[])';  ux = ui(:,1);  uy = ui(:,2);
plotGMSH({ux,'U_x';uy,'U_y';usol,'U_vect';sigma,'stress'}, elements, nodes, 'solution');

plotGMSH({(usol-uref)/norm(uref),'Error'}, elements, nodes, 'error');