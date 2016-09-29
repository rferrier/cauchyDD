%18/04/2016
%Algo KMF avec Orthodir et PGD

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 70000;  % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 1;      % N.mm-1 : Loading on the plate
niter   = 20;     % Nb of Orthodir iterations
br      = 0.;      % noise
mu      = 0.;      % Regularization parameter
nt      = 10;     % Nb of time steps
dt      = 1;      % 
nmodes  = 1;      % Nb of PGD modes
fxpoer  = 1e-3;   % Fixed point error
jlim    = 5;      % Max iterations of the fixed point

Time    = 1:dt:nt*dt;

% Boundary conditions
% first index  : index of the boundary
% second index : 1=x, 2=y
% third        : value
% [0,1,value] marks a dirichlet regularization therm on x
dirichlet = [4,1,0;
             4,2,0];
neumann   = [1,2,fscalar;
             2,1,fscalar];
neumannt  = [3,2,fscalar]; % Time-dependant Neumann loading

% Import the mesh
[ nodes,elements,ntoelem,boundary,order ] = readmesh( 'meshes/plate.msh' );
nnodes = size(nodes,1);

[node2t, t2node] = mapBound(3, boundary, nnodes);

% Then, build the stiffness matrix :
[K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
Kinter = K(1:2*nnodes, 1:2*nnodes);
M      = mass_mat(nodes, elements);
MI     = tempIntegral(Time);

% The right hand side :
fnt = loading(nbloq,nodes,boundary,neumann);
ft = loading(nbloq,nodes,boundary,neumannt);

f = zeros(size(fnt,1),nt);
uin = zeros(size(f,1),nt);
uref = zeros(2*nnodes,nt);
lagr = zeros(nbloq,nt);
for i=1:nt
    f(:,i) = dt*i*fnt + dt^2*(nt-i)*i*ft; % *dt*i, dt*(nt-i)*, sin(2*pi*i/nt)
    % Solve the problem :
    uin(:,i) = K\f(:,i);

    % Extract displacement and Lagrange multiplicators :
    uref(:,i) = uin(1:2*nnodes,i);
    lagr(:,i) = uin(2*nnodes+1:end,i);
end

urefb = ( 1 + br*randn(2*nnodes,nt) ) .* uref;
lagrb = ( 1 + br*randn(nbloq,nt) ) .* lagr;

% Post-pro :
sigma = stress(uref(:,1),E,nu,nodes,elements,order,1,ntoelem);
ui = reshape(uref(:,1),2,[])';  ux = ui(:,1);  uy = ui(:,2);
plotGMSH({ux,'U_x';uy,'U_y';uref(:,1),'U_vect';sigma,'stress'}, elements, nodes, 'reference');

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PGD Process
uPGD     = zeros(size(uref));

%% Write Rhs : (todo, maybe, an other PGD here)
b = zeros(size(uref));
for k=1:nt
    % Solve DN
    f1in = loading(nbloq1,nodes,boundary,neumann1)*dt*k;%*dt*(nt-k);
    fdir = dirichletRhs(zeros( 2*nnodes,1 ), 3, C1, boundary);
    f1 = f1in + fdir;
    uin1 = K1\f1;
    u1 = uin1(1:2*nnodes,1);
    lagr1 = uin1(2*nnodes+1:end,1);
    % Solve ND
    fr = lagr2forces( lagr1, C1, 3, boundary );
    fr = [ fr; zeros(size(C2,2),1) ]; % Give to fr the right size
    fdir1 = dirichletRhs(urefb(:,k), 1, C2, boundary);
    fdir2 = dirichletRhs(urefb(:,k), 2, C2, boundary);
    f2 = fr + assembleDirichlet( [fdir1,fdir2] );
    uin2 = K2\f2;
    b(:,k) = uin2(1:2*nnodes,1);
end

for i=1:nmodes
    % Fixed point
    Lambda   = ones(size(uref,1),1);
    j        = 0; % Nb iterations of fixed point
    errorPGD = fxpoer+1;
    
    %% Perform A*uPGD : (todo, maybe, an other PGD here)
    auPGD = zeros(size(uref));
    for k=1:nt
        % Solve DN
        f1in = loading(nbloq1,nodes,boundary,neumann0);
        fdir = dirichletRhs(uPGD(:,k), 3, C1, boundary);
        f1 = f1in + fdir;
        uin1 = K1\f1;
        u1 = uin1(1:2*nnodes,1);
        lagr1 = uin1(2*nnodes+1:end,1);
        % Solve ND
        fr = lagr2forces( lagr1, C1, 3, boundary );
        fr = [ fr; zeros(size(C2,2),1) ]; % Give to fr the right size (because Lagrange)
        f2 = fr;
        uin2 = K2\f2;
        % Regularization term
        Nu = regul(uPGD(:,k), nodes, boundary, 3);
        auPGD(:,k) = mu*Nu + uPGD(:,k) - uin2(1:2*nnodes,1);
    end
    
    %%
    while errorPGD > fxpoer
        j = j+1;
        if j>jlim  % Emergency quit
            warning('Fixed-point excited without convergence')
            break
        end
        
        %% Perform A Lambda :
        % Solve DN
        f1in = loading(nbloq1,nodes,boundary,neumann0);
        fdir = dirichletRhs(Lambda, 3, C1, boundary);
        f1 = f1in + fdir;
        uin1 = K1\f1;
        u1 = uin1(1:2*nnodes,1);
        lagr1 = uin1(2*nnodes+1:end,1);
        % Solve ND
        fr = lagr2forces( lagr1, C1, 3, boundary );
        fr = [ fr; zeros(size(C2,2),1) ]; % Give to fr the right size (because Lagrange)
        f2 = fr;
        uin2 = K2\f2;
        % Regularization term
        Nu = regul(Lambda, nodes, boundary, 3);
        aLambda = mu*Nu + Lambda - uin2(1:2*nnodes,1);
        
        %% Compute lambda
        mumaupgd = zeros(size(b));
        for k=1:nt
            mumaupgd(:,k) = norm_bound(b(:,k)-auPGD(:,k), nodes, boundary, 3);
        end
        lambda = Lambda'*mumaupgd /(Lambda'*norm_bound(aLambda, nodes, boundary, 3));
        lambda = lambda/norm(lambda);

        Lambdao = Lambda;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% ORTHODIR for the problem : 
        %(1 - D0 o S0 + mu*N) Lambda =
        % int,t lambda'(DoS(0) - (1 - D0 o S0 + mu*N)uPGD) / int,t (lambda'lambda),
        % with N, a regulalization term
        Itere = zeros( 2*nnodes, 1 );
        p     = zeros( 2*nnodes, niter+1 );
        q     = zeros( 2*nnodes, niter+1 );
        Delta = zeros( niter+1, 1 );
        Res   = zeros( 2*nnodes, niter+1 );
        
        RHS = (b-auPGD)*MI*lambda' / (lambda*MI*lambda');

        %% Perform A x0 :
        % Solve DN
        f1in = loading(nbloq1,nodes,boundary,neumann0);
        fdir = dirichletRhs(Itere, 3, C1, boundary);
        f1 = f1in + fdir;
        uin1 = K1\f1;
        u1 = uin1(1:2*nnodes,1);
        lagr1 = uin1(2*nnodes+1:end,1);
        % Solve ND
        fr = lagr2forces( lagr1, C1, 3, boundary );
        fr = [ fr; zeros(size(C2,2),1) ]; % Give to fr the right size (because Lagrange)
        f2 = fr;
        uin2 = K2\f2;
        % Regularization term
        Nu = regul(Itere, nodes, boundary, 3);
        atimesItere = mu*Nu + Itere - uin2(1:2*nnodes,1);
        %%%%

        %%
        Res(:,1) = RHS - atimesItere;
        p(:,1) = Res(:,1);

%         residual(1) = sqrt( myps( Res(:,1), Res(:,1), Kinter, boundary, M, nodes ) );
%         error(1)    = sqrt( myps( Itere - uref, Itere - uref, Kinter, boundary, M, nodes )...
%                          / myps( uref, uref, Kinter, boundary, M, nodes ) );
%         regulari(1) = sqrt( Itere'*regul(Itere, nodes, boundary, 3) );

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
        f2 = fr;% + fdir1 + fdir2;
        uin2 = K2\f2;
        % Regularization term
        Nu = regul(p(:,1), nodes, boundary, 3);
        q(:,1) = mu*Nu + p(:,1) - uin2(1:2*nnodes,1);
        %%%%
        %%
        for iter = 1:niter

            Delta(iter,1) = myps( q(:,iter), q(:,iter), Kinter, boundary, M, nodes ); %q(:,iter)'*q(:,iter);
            gammai        = myps( q(:,iter), Res(:,iter), Kinter, boundary, M, nodes ); %q(:,iter)'*Res;
            alphai        = gammai/Delta(iter,1);

            Itere         = Itere + p(:,iter)*alphai;
            Res(:,iter+1) = Res(:,iter) - q(:,iter)*alphai;

%             residual(iter+1) = sqrt( myps( Res(:,iter+1), Res(:,iter+1), Kinter, boundary, M, nodes ) );
%             error(iter+1)    = sqrt( myps( Itere - uref, Itere - uref, Kinter, boundary, M, nodes )...
%                              / myps( uref, uref, Kinter, boundary, M, nodes ) );
%             regulari(iter+1) = sqrt( Itere'*regul(Itere, nodes, boundary, 3) );

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
            Nu = regul(Res(:,iter+1), nodes, boundary, 3);
            Ari = mu*Nu + Res(:,iter+1) - uin2(1:2*nnodes,1);

            %% Orthogonalization
            p(:,iter+1) = Res(:,iter+1);
            q(:,iter+1) = Ari;

            for jter=1:iter
                phiij  = myps( q(:,jter), Ari, Kinter, boundary, M, nodes ); %q(:,jter)'*Ari;
                betaij = phiij/Delta(jter,1);
                p(:,iter+1) = p(:,iter+1) - p(:,jter) * betaij;
                q(:,iter+1) = q(:,iter+1) - q(:,jter) * betaij;
            end

        end

    %     hold on
    %     plot(log10(error),'Color','blue')
    %     plot(log10(residual),'Color','red')
    %     legend('error (log)','residual (log)')
    %     figure
    %     % L-curve
    %     loglog(regulari, residual);
    
        % Upgrade the error
        Lambda = Itere;
        errorPGD = sqrt( myps( Lambdao-Lambda, Lambdao-Lambda, Kinter, boundary, M, nodes ) /...
            myps( Lambda, Lambda, Kinter, boundary, M, nodes ));
    end
    % Upgrade the solution
    uPGD = uPGD + Lambda*lambda;
end

surf(uref(2*t2node,:))
shading interp;
figure
surf(uPGD(2*t2node,:))
shading interp;
figure
surf(uPGD(2*t2node,:)-uref(2*t2node,:))
shading interp;

% figure
% hold on
% plot(uref(2*t2node(1),:),'Color','red')
% affn = uref(2*t2node(1),1):(uref(2*t2node(1),nt)-uref(2*t2node(1),1))/(nt-1)...
%     :uref(2*t2node(1),nt);
% plot(affn)
% delt = affn-uref(2*t2node(1),:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Final problem : compute u
% DN problem
dirichlet = [4,1,0;
             4,2,0;
             3,1,0;
             3,2,0;
             1,1,0;
             1,2,0;
             2,1,0;
             2,2,0];
neumann   = [];
[K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);

usol = zeros(2*nnodes, nt);

for i=1:nt
    fdir1 = dirichletRhs(urefb(:,1), 1, C, boundary);
    fdir2 = dirichletRhs(urefb(:,1), 2, C, boundary);
    fdir3 = dirichletRhs(uPGD(:,1), 3, C, boundary);
    usoli = K \ assembleDirichlet( [fdir1,fdir2,fdir3] );
    usol(:,i) = usoli(1:2*nnodes,1);
end

total_error = norm(usol-uref)/norm(uref);

% Post-pro :
sigma = stress(usol(:,1),E,nu,nodes,elements,order,1,ntoelem);
ui = reshape(usol(:,1),2,[])';  ux = ui(:,1);  uy = ui(:,2);
plotGMSH({ux,'U_x';uy,'U_y';usol(:,1),'U_vect';sigma,'stress'}, elements, nodes, 'solution');

plotGMSH({(usol(:,1)-uref(:,1))/norm(uref(:,1)),'Error'}, elements, nodes, 'error');