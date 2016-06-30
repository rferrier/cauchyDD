%30/03/2016
%Algo KMF avec CL de Robin (Latin)

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E         = 70000;      % MPa : Young modulus
nu        = 0.3;    % Poisson ratio
fscalar   = 1;      % N.mm-1 : Loading on the plate
niter     = 15;
br        = 0.;     % noise
robin1    = 1e5*E;%1e5*E;%0.38*E;%1e5*E;%-7*E;%-1e1*E;  % Robin parameters
robin2    = -0.097*E;%-0.097*E;%-0.16*E;%1e-5*E;
rank      = 0;       % Half of the size of the approximation of Si 
                     % (if 0, don't use it)

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

% Then, build the stiffness matrix :
[K,C,nbloq] = Krig (nodes,elements,E,nu,order,boundary,dirichlet);
Kinter = K(1:2*nnodes, 1:2*nnodes);
M      = mass_mat(nodes, elements);

% The right hand side :
f = loading(nbloq,nodes,boundary,neumann);

% Solve the problem :
uin = K\f;

% Extract displacement and Lagrange multiplicators :
uref = uin(1:2*nnodes,1);
lagr = uin(2*nnodes+1:end,1);
urefb = ( 1 + br*randn(2*nnodes,1) ) .* uref;
lagrb = ( 1 + br*randn(nbloq,1) ) .* lagr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Obtain an approxiamtion of the Schur complements

xmax = max(nodes(:,1));
xmin = min(nodes(:,1));
ymax = max(nodes(:,2));
no3  = findNode(xmax, ymax, nodes, 1e-5);
no4  = findNode(xmin, ymax, nodes, 1e-5);
boundaryd = suppressBound( boundary, [no3;no4], 3 );

if rank > 0

    %% Preliminary stuff
    [ node2b3, b2node3 ] = mapBound( 3, boundaryd, nnodes );
    bsize3 = size(b2node3,1);
    index = zeros( 2*bsize3 ,1 );
    index(1:2:2*bsize3-1) = 2*b2node3-1;
    index(2:2:2*bsize3)   = 2*b2node3;
    
    %% First compute b1 and b2
    
    % problem 1
    dirichlet1 = [4,1,0;4,2,0 ; 3,1,0;3,2,0];
    neumann1   = [1,2,fscalar;
                  2,1,fscalar];
    [K1,C1,nbloq1,node2c1,c2node1] = Krig (nodes,elements,E,nu,order,boundary,dirichlet1);
    f1in = loading(nbloq1,nodes,boundary,neumann1);
    % problem 2
    dirichlet2 = [4,1,0;4,2,0;
                  3,1,0;3,2,0;
                  1,1,0;1,2,0;
                  2,1,0;2,2,0];
    [K2,C2,nbloq2,node2c2,c2node2] = Krig (nodes,elements,E,nu,order,boundary,dirichlet2);
    
    u1 = K1\f1in;
    f1 = Kinter*u1(1:2*nnodes);
    b1 = f1([2*b2node3-1 ; 2*b2node3]);
    
    fd1 = dirichletRhs2( urefb, 1, c2node2, boundary, nnodes );
    fd2 = dirichletRhs2( urefb, 2, c2node2, boundary, nnodes );
    u2 = K2\assembleDirichlet( [fd1,fd2] );
    f2 = Kinter*u2(1:2*nnodes);
    b2 = f2([2*b2node3-1 ; 2*b2node3]);
    
    %% Compute D1*b2 and D2*b1
    % problem 1
    dirichlet1 = [4,1,0;4,2,0];
    [K1,C1,nbloq1,node2c1,c2node1] = Krig (nodes,elements,E,nu,order,boundary,dirichlet1);
    % problem 2
    dirichlet2 = [4,1,0;4,2,0;
                  1,1,0;1,2,0;
                  2,1,0;2,2,0];
    [K2,C2,nbloq2,node2c2,c2node2] = Krig (nodes,elements,E,nu,order,boundary,dirichlet2);
    
    f1 = zeros(size(K1,1),1); f1([2*b2node3-1 ; 2*b2node3]) = b2;
    u1 = K1\f1;
    S1b2 = u1([2*b2node3-1 ; 2*b2node3]);
    
    f2 = zeros(size(K2,1),1); f2([2*b2node3-1 ; 2*b2node3]) = b1;
    u2 = K2\f2;
    S2b1 = u2([2*b2node3-1 ; 2*b2node3]);
    
    %% Values
    kkedevrai1 = norm(b1)/norm(S2b1) /E;
    kkedevrai2 = norm(b2)/norm(S1b2) /E;
    
    %% Then, estimate S1 and S2
    L = xmax-xmin;
    P = zeros( 2*bsize3, 2*rank ); % Projection matrix
    
    for i=1:rank
        for j=1:bsize3 % Loop over nodes of the boundary
            x = nodes( b2node3(j), 1 ) - xmin;
            P(2*j-1,2*i-1) = sin(i*pi*x/L);
            P(2*j,2*i) = sin(i*pi*x/L);
        end
        % normalize
        P(:,2*i-1) = P(:,2*i-1)/norm(P(:,2*i-1));
        P(:,2*i) = P(:,2*i)/norm(P(:,2*i));
    end
    
    %% Define the Dirichlet problems
    % problem 1
    dirichlet1 = [ 4,1,0;4,2,0 ; 3,1,0;3,2,0 ];
    [K1,~,~,~,c2node1] = Krig (nodes,elements,E,nu,order,boundaryd,dirichlet1);
    % problem 2
    dirichlet2 = [ 4,1,0;4,2,0 ; 3,1,0;3,2,0 ; 1,1,0;1,2,0 ; 2,1,0;2,2,0 ];
    [K2,~,~,~,c2node2] = Krig (nodes,elements,E,nu,order,boundaryd,dirichlet2);
    
    % Solve the 2*2*rank problems
    f1h = zeros( 2*bsize3, 2*rank );
    f2h = zeros( 2*bsize3, 2*rank );
    
    for i=1:rank
        % First problem, x
        udir = zeros(2*nnodes, 1);
        udir(index) = P(:,2*i-1);
        fdir = dirichletRhs2( udir, 3, c2node1, boundaryd, nnodes );
        u = K1\fdir;
        f = Kinter*u(1:2*nnodes,1);
        f1h(:,2*i-1) = f(index);
        
        % First problem, y
        udir = zeros(2*nnodes, 1);
        udir(index) = P(:,2*i);
        fdir = dirichletRhs2( udir, 3, c2node1, boundaryd, nnodes );
        u = K1\fdir;
        f = Kinter*u(1:2*nnodes,1);
        f1h(:,2*i) = f(index);
        
        % Second problem, x
        udir = zeros(2*nnodes, 1);
        udir(index) = P(:,2*i-1);
        fdir = dirichletRhs2( udir, 3, c2node2, boundaryd, nnodes );
        u = K2\fdir;
        f = Kinter*u(1:2*nnodes,1);
        f2h(:,2*i-1) = f(index);
        
        % Second problem, y
        udir = zeros(2*nnodes, 1);
        udir(index) = P(:,2*i);
        fdir = dirichletRhs2( udir, 3, c2node2, boundaryd, nnodes );
        u = K2\fdir;
        f = Kinter*u(1:2*nnodes,1);
        f2h(:,2*i) = f(index);
        
    end

    % The macro forces
    F1 = P'*f1h;
    F2 = P'*f2h;

    % The approximation of the Schur complement : Sa = F*W-1, with W = Id
    Sa1 = P*F1*P';
    Sa2 = P*F2*P'; 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init : get null vectors
u1 = uref-uref;
u2 = u1;
f2 = u2;

% problem 1
dirichlet1 = [4,1,0;4,2,0];

neumann1   = [1,2,fscalar;
              2,1,fscalar];
[K1,C1,nbloq1,node2c1,c2node1] = Krig (nodes,elements,E,nu,order,boundary,dirichlet1);

f1in = loading(nbloq1,nodes,boundary,neumann1);

% problem 2
dirichlet2 = [4,1,0;4,2,0;
              1,1,0;1,2,0;
              2,1,0;2,2,0];

[K2,C2,nbloq2,node2c2,c2node2] = Krig (nodes,elements,E,nu,order,boundary,dirichlet2);

if rank > 0
    % full size matrix to add
    Kaj1 = K1-K1;
    Kaj2 = K2-K2;
    Kaj1(index,index) = robin1*eye(size(index,1));%-Sa2;
    Kaj2(index,index) = -Sa1;%robin2*eye(size(index,1));

%     kkedevrai = norm(Sa1*b2)/norm(b2) /E
%     kkedevrai = norm(Sa1*b1)/norm(b1) /E
%     kkedevrai = norm(Sa2*b1)/norm(b1) /E
%     kkedevrai = norm(Sa2*b2)/norm(b2) /E
    %kkedevrai = sqrt( - 0.5* ( norm(Sa1)^2 + 1 - norm(Sa1+eye(size(Sa1)))^2 ) ) /E
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fixed point algorithm

error1   = zeros(niter,1);
error2   = zeros(niter,1);
residual = zeros(niter,1);
regulari = zeros(niter,1);
for iter = 1:niter
    % Solve 1
    
    if rank == 0
        kuplusf = u2 + f2/robin1;
        [Kp1, fro1] = robinRHS( nbloq1, nodes, boundary, kuplusf, robin1, 3 );
        f1 = f1in + fro1;
        K1e = K1 + Kp1;
    else
        faj1 = Kaj1*[ u2 ; zeros(nbloq1,1) ] +...
            [ keepField( f2, 3, boundaryd ) ; zeros(nbloq1,1) ];
        f1 = f1in + faj1;
        K1e = K1 + Kaj1;
    end

    uin1 = K1e\f1;
    u1 = uin1(1:2*nnodes,1);
    f1 = Kinter*u1 - f1in(1:2*nnodes,1); % Get the forces at the boundary

    error1(iter) = sqrt( myps(u1-uref,u1-uref,Kinter,boundary,M,nodes)/...
        myps(uref,uref,Kinter,boundary,M,nodes) );
    
    % Solve 2
    
    fdir1 = dirichletRhs2(urefb, 1, c2node2, boundary, nnodes);
    fdir2 = dirichletRhs2(urefb, 2, c2node2, boundary, nnodes);
    
    if rank == 0
        kuplusf = u1 + f1/robin2;
        [Kp2, fro2] = robinRHS( nbloq2, nodes, boundary, kuplusf, robin2, 3 );
        f2 = fro2 + assembleDirichlet( [fdir1,fdir2] );
        K2e = K2 + Kp2;
    else
        faj2 = Kaj2*[ u1 ; zeros(nbloq2,1) ] +...
            [ keepField( f1, 3, boundaryd ) ; zeros(nbloq2,1) ];
        f2 = assembleDirichlet( [fdir1,fdir2] ) + faj2;
        K2e = K2 + Kaj2;
    end
    
    uin2 = K2e\f2;
    u2 = uin2(1:2*nnodes,1);
    f2 = Kinter*u2; % Get the forces at the boundary

    error2(iter) = sqrt( myps(u2-uref,u2-uref,Kinter,boundary,M,nodes)/...
        myps(uref,uref,Kinter,boundary,M,nodes) );
    
    residual(iter) = sqrt( myps(u1-u2,u1-u2,Kinter,boundary,M,nodes)/...
                     sqrt(myps(u1,u1,Kinter,boundary,M,nodes)*...
                     myps(u2,u2,Kinter,boundary,M,nodes)) );
                 
    regulari(iter) = sqrt(u2'*regul(u2, nodes, boundary, 3));
end

residual(1) = 1; % trick

hold on;
set(gca, 'fontsize', 15);
set(gca,'ylim',[-3 0])
plot(log10(error1),'Color','black')
plot(log10(error2),'Color','blue')
plot(log10(residual),'Color','red')
legend('error1 (log)','error2 (log)','residual (log)')
figure;
% L-curve
loglog(residual,regulari);