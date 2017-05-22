% 19/05/2017
% Détection de fissure 3D plane par SP

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E        = 210000; % MPa : Young modulus
nu       = 0.3;    % Poisson ratio
fscalar  = 250;    % N.mm-2 : Loading
mat      = [0, E, nu];
regmu    = 0;      % Regularization parameter
br       = .0;
niter_up = 4;      % Number of iterations
niter_do = 4;

cracked_mesh   = 'meshes/rg3dm/platemSPc.msh';
uncracked_mesh = 'meshes/rg3dm/platemSPu.msh';

centCrack = [4;3;1]; % Point on the crack (for reference)
alpha     = pi/15;
tic
% Boundary conditions
% first index  : index of the boundary
% second index : 1=x, 2=y
% third        : value
% [0,1,value] marks a dirichlet regularization therm on x
dirichlet  = [0,1,0 ; 0,2,0 ; 0,3,0 ; 0,4,0 ; 0,5,0 ; 0,6,0];
neumann1   = [2,3,fscalar ; 1,3,-fscalar];
neumann2   = [4,1,fscalar ; 6,1,-fscalar];

% First, import the mesh
[ nodes,elements,ntoelem,boundary,order] = readmesh3D( cracked_mesh );
nnodes = size(nodes,1);

% mapBounds
[ node2b1, b2node1 ]   = mapBound3D( 1, boundary, nnodes );
[ node2b2, b2node2 ]   = mapBound3D( 2, boundary, nnodes );
[ node2b3, b2node3 ]   = mapBound3D( 3, boundary, nnodes );
[ node2b4, b2node4 ]   = mapBound3D( 4, boundary, nnodes );
[ node2b5, b2node5 ]   = mapBound3D( 5, boundary, nnodes );
[ node2b6, b2node6 ]   = mapBound3D( 6, boundary, nnodes );

indexbound  = [3*b2node1-2 ; 3*b2node1-1 ; 3*b2node1 ;...
               3*b2node2-2 ; 3*b2node2-1 ; 3*b2node2 ;...
               3*b2node3-2 ; 3*b2node3-1 ; 3*b2node3 ;...
               3*b2node4-2 ; 3*b2node4-1 ; 3*b2node4 ;...
               3*b2node5-2 ; 3*b2node5-1 ; 3*b2node5 ;...
               3*b2node6-2 ; 3*b2node6-1 ; 3*b2node6 ];

% Then, build the stiffness matrix :
[K,C,nbloq,node2c,c2node] = Krig3D (nodes,elements,mat,order,boundary,dirichlet);
Kinter = K( 1:3*nnodes, 1:3*nnodes );

f1  = loading3D(nbloq,nodes,boundary,neumann1);
uin = K\f1;
resno = norm(K*uin-f1)/norm(f1);
u1 = uin(1:3*nnodes,1);
f1 = Kinter*u1;

%f2  = loading3D(nbloq,nodes,boundary,neumann2);
%uin = K\f2;
%u2 = uin(1:3*nnodes,1);
%f2 = Kinter*u2;

ui = reshape(u1,3,[])'; ux = ui(:,1); uy = ui(:,2); uz = ui(:,3);

% Compute stress :
sigma = stress3D(u1,mat,nodes,elements,order,1,ntoelem);

% Output :
plotGMSH3D({ux,'U_x';uy,'U_y';uz,'U_z';u1,'U_vect';sigma,'stress'}, elements, nodes, 'output/reference');
disp([ 'Direct problem solved ', num2str(toc) ]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build the data : Adapt mesh
tic
% Import the uncracked mesh
%[ nodes,elements,ntoelem,boundary,order] = readmesh3D( uncracked_mesh );
%nnodes = size(nodes,1);

nodes_up = nodes; elements_up = elements; ntoelem_up = ntoelem; boundary_up = boundary;
nodes_do = nodes; elements_do = elements; ntoelem_do = ntoelem; boundary_do = boundary;

u1_up = u1; f1_up = f1; u1_do = u1; f1_do = f1;

crit = -centCrack(1)*sin(alpha) + centCrack(3)*cos(alpha);
epsi = 1e-6;

toremoveup = []; soustracup = zeros(nnodes,1);
toremovedo = []; soustracdo = zeros(nnodes,1);

jup = 0; jdo = 0;
for i=1:nnodes
   zed = -nodes(i,1)*sin(alpha) + nodes(i,3)*cos(alpha);
   if zed+epsi < crit
      jup = jup+1;
      toremoveup(end+1) = i;
   elseif zed-epsi > crit
      jdo = jdo+1;
      toremovedo(end+1) = i;
   end
   soustracup(i,1) = jup; % This list gives the way to renumber the nodes in elements_up
   soustracdo(i,1) = jdo;
end
nodes_up(toremoveup,:) = []; ntoelem_up(toremoveup) = [];
nodes_do(toremovedo,:) = []; ntoelem_do(toremovedo) = [];

u1_up([3*toremoveup-2;3*toremoveup-1;3*toremoveup]) = [];
f1_up([3*toremoveup-2;3*toremoveup-1;3*toremoveup]) = [];
u1_do([3*toremovedo-2;3*toremovedo-1;3*toremovedo]) = [];
f1_do([3*toremovedo-2;3*toremovedo-1;3*toremovedo]) = [];

torup = []; tordo = [];
for i=1:size(toremoveup,2)
   tr = find(elements == toremoveup(i));
   tr = mod(tr-1,size(elements,1))+1; % Because there are several raws in elements
   torup = [torup;tr];
end
for i=1:size(toremovedo,2)
   tr = find(elements == toremovedo(i));
   tr = mod(tr-1,size(elements,1))+1; % Because there are several raws in elements
   tordo = [tordo;tr];
end
elements_up(torup,:) = [];
elements_do(tordo,:) = [];

torup = []; tordo = [];
for i=1:size(toremoveup,2)
   tr = find(boundary(:,2:end) == toremoveup(i));
   tr = mod(tr-1,size(boundary,1))+1; % Because there are several raws in elements
   torup = [torup;tr];
end
for i=1:size(toremovedo,2)
   tr = find(boundary(:,2:end) == toremovedo(i));
   tr = mod(tr-1,size(boundary,1))+1; % Because there are several raws in elements
   tordo = [tordo;tr];
end
boundary_up(torup,:) = [];
boundary_do(tordo,:) = [];

% Renumber elements and boundary
elements_up = elements_up - soustracup(elements_up);
boundary_up(:,2:end) = boundary_up(:,2:end) - soustracup(boundary_up(:,2:end));
elements_do = elements_do - soustracdo(elements_do);
boundary_do(:,2:end) = boundary_do(:,2:end) - soustracdo(boundary_do(:,2:end));

nnodes_up = size(nodes_up,1); nnodes_do = size(nodes_do,1);

%% Second pass : seek and destroy orphan nodes (= doubled nodes on the crack)
toremoveup = []; soustracup = zeros(nnodes,1);
toremovedo = []; soustracdo = zeros(nnodes,1);

jup = 0; jdo = 0;
for i=1:nnodes_up
   if size( find(elements_up(:,2:end)==i),1 ) == 0
      jup = jup+1;
      toremoveup(end+1) = i;
   end
   soustracup(i,1) = jup;
end
for i=1:nnodes_do
   if size( find(elements_do(:,2:end)==i),1 ) == 0
      jdo = jdo+1;
      toremovedo(end+1) = i;
   end
   soustracdo(i,1) = jdo;
end
nodes_up(toremoveup,:) = []; ntoelem_up(toremoveup) = [];
nodes_do(toremovedo,:) = []; ntoelem_do(toremovedo) = [];

u1_up([3*toremoveup-2;3*toremoveup-1;3*toremoveup]) = [];
f1_up([3*toremoveup-2;3*toremoveup-1;3*toremoveup]) = [];
u1_do([3*toremovedo-2;3*toremovedo-1;3*toremovedo]) = [];
f1_do([3*toremovedo-2;3*toremovedo-1;3*toremovedo]) = [];
% And of course, no need to suppress elements (cause we were chasing orphan nodes)

torup = []; tordo = [];
for i=1:size(toremoveup,2)
   tr = find(boundary_up(:,2:end) == toremoveup(i));
   tr = mod(tr-1,size(boundary_up,1))+1; % Because there are several raws in elements
   torup = [torup;tr];
end
for i=1:size(toremovedo,2)
   tr = find(boundary_do(:,2:end) == toremovedo(i));
   tr = mod(tr-1,size(boundary_do,1))+1;
   tordo = [tordo;tr];
end
boundary_up(torup,:) = [];
boundary_do(tordo,:) = [];

% Renumber elements and boundary
elements_up = elements_up - soustracup(elements_up);
boundary_up(:,2:end) = boundary_up(:,2:end) - soustracup(boundary_up(:,2:end));
elements_do = elements_do - soustracdo(elements_do);
boundary_do(:,2:end) = boundary_do(:,2:end) - soustracdo(boundary_do(:,2:end));

%% Work on the boundary 8 : first, puppress it, then rebuild it
boundary_up( find( boundary_up(:,1)==8 ), : ) = [];
boundary_do( find( boundary_do(:,1)==8 ), : ) = [];

for i=1:size(elements_up,1)
   no1 = elements_up(i,1); no2 = elements_up(i,2);
   no3 = elements_up(i,3); no4 = elements_up(i,4);
   
   zed1 = -nodes_up(no1,1)*sin(alpha) + nodes_up(no1,3)*cos(alpha);
   zed2 = -nodes_up(no2,1)*sin(alpha) + nodes_up(no2,3)*cos(alpha);
   zed3 = -nodes_up(no3,1)*sin(alpha) + nodes_up(no3,3)*cos(alpha);
   zed4 = -nodes_up(no4,1)*sin(alpha) + nodes_up(no4,3)*cos(alpha);
   
   score = [];
   if zed1-epsi < crit %zed1+epsi > crit && 
      score(end+1) = no1;
   end
   if zed2-epsi < crit
      score(end+1) = no2;
   end
   if zed3-epsi < crit
      score(end+1) = no3;
   end
   if zed4-epsi < crit
      score(end+1) = no4;
   end
   
   if size(score,2) == 3
      boundary_up(end+1,:) = [8,score];
   end
end

for i=1:size(elements_do,1)
   no1 = elements_do(i,1); no2 = elements_do(i,2);
   no3 = elements_do(i,3); no4 = elements_do(i,4);
   
   zed1 = -nodes_do(no1,1)*sin(alpha) + nodes_do(no1,3)*cos(alpha);
   zed2 = -nodes_do(no2,1)*sin(alpha) + nodes_do(no2,3)*cos(alpha);
   zed3 = -nodes_do(no3,1)*sin(alpha) + nodes_do(no3,3)*cos(alpha);
   zed4 = -nodes_do(no4,1)*sin(alpha) + nodes_do(no4,3)*cos(alpha);
   
   score = [];
   if zed1+epsi > crit
      score(end+1) = no1;
   end
   if zed2+epsi > crit
      score(end+1) = no2;
   end
   if zed3+epsi > crit
      score(end+1) = no3;
   end
   if zed4+epsi > crit
      score(end+1) = no4;
   end
   
   if size(score,2) == 3
      boundary_do(end+1,:) = [8,score];
   end
end

%% Extract the boundary 8 (missing boundary)
bounindex_up = find( boundary_up(:,1)==8 );
patch2_up    = boundary_up(bounindex_up,2:end);
bounindex_do = find( boundary_do(:,1)==8 );
patch2_do    = boundary_do(bounindex_do,2:end);

boundaryp_up = suppressBound( boundary_up, 'extreme', 8, nodes_up, 1e-5 );
boundaryp_do = suppressBound( boundary_do, 'extreme', 8, nodes_do, 1e-5 );

nnodes_up = size(nodes_up,1); nnodes_do = size(nodes_do,1);

disp([ 'Mesh adaptated ', num2str(toc) ]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build the operators & such
tic

dirichlet1d_up = [6,1,0;6,2,0;6,3,0;
                  5,1,0;5,2,0;5,3,0;
                  4,1,0;4,2,0;4,3,0;
                  3,1,0;3,2,0;3,3,0;
                  2,1,0;2,2,0;2,3,0];
[K1d_up,C1d_up,nbloq1d_up] = Krig3D(nodes_up,elements_up,mat,order,boundary_up,dirichlet1d_up);

% Second problem
dirichlet2d_up = [6,1,0;6,2,0;6,3,0;
                  5,1,0;5,2,0;5,3,0;
                  4,1,0;4,2,0;4,3,0;
                  3,1,0;3,2,0;3,3,0];
[K2d_up,C2d_up,nbloq2d_up] = Krig3D(nodes_up,elements_up,mat,order,boundary_up,dirichlet2d_up);

f2d = loading3D(nbloq2d_up,nodes_up,boundary_up,neumann1);

Cm = Cbound ( nodes_up, [8,1,0;8,2,0;8,3,0], boundaryp_up );
Cr = Cbound ( nodes_up, [2,1,0;2,2,0;2,3,0], boundary );

% find the index of the nodes on the crack plane
indexC = sum(Cm'); indexa = find(indexC);
C8 = Cbound ( nodes, [8,1,0;8,2,0;8,3,0], boundary );
indexC8 = sum(C8'); indexa8 = find(indexC8);

%%%%%%%%%%%%%%%%%%%
%% Down subdomain

dirichlet1d_do = [6,1,0;6,2,0;6,3,0;
                  5,1,0;5,2,0;5,3,0;
                  4,1,0;4,2,0;4,3,0;
                  3,1,0;3,2,0;3,3,0;
                  1,1,0;1,2,0;1,3,0];
[K1d_do,C1d_do,nbloq1d_do] = Krig3D(nodes_do,elements_do,mat,order,boundary_do,dirichlet1d_do);

% Second problem
dirichlet2d_do = [6,1,0;6,2,0;6,3,0;
                  5,1,0;5,2,0;5,3,0;
                  4,1,0;4,2,0;4,3,0;
                  3,1,0;3,2,0;3,3,0];
[K2d_do,C2d_do,nbloq2d_do] = Krig3D(nodes_do,elements_do,mat,order,boundary_do,dirichlet2d_do);

f2ddo = loading3D(nbloq2d_do,nodes_do,boundary_do,neumann1);

Cmdo = Cbound ( nodes_do, [8,1,0;8,2,0;8,3,0], boundaryp_do );
Crdo = Cbound ( nodes_do, [1,1,0;1,2,0;1,3,0], boundary );

% find the index of the nodes on the crack plane
indexCdo = sum(Cmdo'); indexado = find(indexCdo);
%C8 = Cbound ( nodes, [8,1,0;8,2,0;8,3,0], boundary );
%indexC8 = sum(C8'); indexa8 = find(indexC8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CG resolution

%% First, plot the reference
%ret = patch('Faces',patch2_up,'Vertices',nodes_up, ... 
%            'FaceVertexCData',u1_up(3:3:end),'FaceColor','interp');
%colorbar;

niter = niter_up;
error    = zeros(niter+1,1);
residual = zeros(niter+1,1);
regulari = zeros(niter+1,1);

Itere = zeros( 3*nnodes_up, 1 );
d     = zeros( 3*nnodes_up, niter+1 );
Ad    = zeros( 3*nnodes_up, niter+1 );
Res   = zeros( 3*nnodes_up, niter+1 );
Zed   = zeros( 3*nnodes_up, niter+1 );
%alpha = zeros( niter+1, 1 );
%beta  = zeros( niter+1, 1 );
%alpha2 = zeros( niter+1, 1 );

%% Perform A x0 :
% Solve 1
f1 = [Itere; zeros(nbloq1d_up,1)];
uin1 = K1d_up\f1;
u1i = uin1(1:3*nnodes_up,1);
u1 = Cm*Cm'*u1i;
% Solve 2
f2 = [Itere; zeros(nbloq2d_up,1)];
uin2 = K2d_up\f2;
u2i = uin2(1:3*nnodes_up,1);
u2 = Cm*Cm'*u2i;
%
Axz = u2-u1;

%%%%
%% Compute Rhs :
% Solve 1
f1 = [ zeros(3*nnodes_up,1) ; C1d_up'*Cr*Cr'*u1_up ];
uin1 = K1d_up\f1; resno = norm(K1d_up*uin1-f1)/norm(f1);
u1i = uin1(1:3*nnodes_up);
u1 = Cm*Cm'*u1i;
% Solve 2
f2 = f2d;
uin2 = K2d_up\f2; resno = norm(K1d_up*uin1-f1)/norm(f1);
u2i = uin2(1:3*nnodes_up);
u2 = Cm*Cm'*u2i;
%
b = u1-u2;
%%
Res(:,1) = b - Axz;
Zed(:,1) = Res(:,1);
d(:,1) = Zed(:,1);

residual(1) = norm(Res( :,1));
error(1)    = norm(Itere(indexa) - f1_up(indexa)) / norm(f1_up(indexa));
regulari(1) = norm(Itere);

for iter = 1:niter
    %% Optimal step
    % Solve 1
    f1 = [d(:,iter); zeros(nbloq1d_up,1)];
    uin1 = K1d_up\f1; norm(K1d_up*uin1-f1)/norm(f1);
    u1i = uin1(1:3*nnodes_up,1);
    u1 = Cm*Cm'*u1i;
    % Solve 2
    f2 = [d(:,iter); zeros(nbloq2d_up,1)];
    uin2 = K2d_up\f2; norm(K2d_up*uin2-f2)/norm(f2);
    u2i = uin2(1:3*nnodes_up,1);
    u2 = Cm*Cm'*u2i;
    %
    Ad(:,iter) = u2-u1;
    
    den = (d(:,iter)'*Ad(:,iter));
    d(:,iter) = d(:,iter)/sqrt(den); Ad(:,iter) = Ad(:,iter)/sqrt(den);
    num = Res(:,iter)'*d(:,iter);
    
    Itere         = Itere + d(:,iter)*num;%/den;
    Res(:,iter+1) = Res(:,iter) - Ad(:,iter)*num;%/den;
    
    Zed(:,iter+1) = Res(:,iter+1);
    
    residual(iter+1) = norm(Res(:,iter+1));
    error(iter+1)    = norm(Itere(indexa) - f1_up(indexa)) / norm(f1_up(indexa));
%    regulari(iter+1) = sqrt( Itere'*regul(Itere, nodes, boundary, 2) );
    regulari(iter+1) = norm(Itere);

%    % First Reorthogonalize the residual (as we use it next), in sense of M
%    for jter=1:iter-1
%        betac = Zed(:,iter+1)'*Res(:,jter) / (Zed(:,jter)'*Res(:,jter));
%        Zed(:,iter+1) = Zed(:,iter+1) - Zed(:,jter) * betac;
%    end
    
    %% Orthogonalization
    d(:,iter+1) = Zed(:,iter+1);
    
    for jter=1:iter % No need to reorthogonalize (see above)
        betaij = ( Zed(:,iter+1)'*Ad(:,jter) );%/...
            %( d(:,jter)'*Ad(:,jter) );
        d(:,iter+1) = d(:,iter+1) - d(:,jter) * betaij;
    end
    
end

disp([ 'Inverse problem solved ', num2str(toc) ]);

%% Convergence curve
%figure;
%hold on;
%%plot(error,'Color','blue'); % The error is wrong
%plot(residual,'Color','black');
%legend('residual');

% L-curve
figure;
loglog( residual(2:end) , regulari(2:end) );

% Final problem
f_upD = [ zeros(3*nnodes_up,1) ; C1d_up'*C1d_up*C1d_up'*u1_up ];
f_upN = [Itere; zeros(nbloq1d_up,1)];
u_up  = K1d_up\(f_upN+f_upD);

error_up = norm( u_up(indexa) - u1_up(indexa) ) / norm( u1_up(indexa) );

%figure;
%ret = patch('Faces',patch2_up,'Vertices',nodes_up, ... 
%            'FaceVertexCData',f_upN(3:3:end),'FaceColor','interp');
%colorbar;

%figure;
%ret = patch('Faces',patch2_up,'Vertices',nodes_up, ... 
%            'FaceVertexCData',u_up(3:3:end),'FaceColor','interp');
%colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CG resolution (down)

%% First, plot the reference
%ret = patch('Faces',patch2_do,'Vertices',nodes_do, ... 
%            'FaceVertexCData',u1_do(3:3:end),'FaceColor','interp');
%colorbar;

niter      = niter_do;
errordo    = zeros(niter+1,1);
residualdo = zeros(niter+1,1);
regularido = zeros(niter+1,1);

Iteredo = zeros( 3*nnodes_do, 1 );
d       = zeros( 3*nnodes_do, niter+1 );
Ad      = zeros( 3*nnodes_do, niter+1 );
Res     = zeros( 3*nnodes_do, niter+1 );
Zed     = zeros( 3*nnodes_do, niter+1 );
%alpha   = zeros( niter+1, 1 );
%beta    = zeros( niter+1, 1 );
%alpha2   = zeros( niter+1, 1 );

%% Perform A x0 :
% Solve 1
f1 = [Iteredo; zeros(nbloq1d_do,1)];
uin1 = K1d_do\f1;
u1i = uin1(1:3*nnodes_do,1);
u1 = Cmdo*Cmdo'*u1i;
% Solve 2
f2 = [Iteredo; zeros(nbloq2d_do,1)];
uin2 = K2d_do\f2;
u2i = uin2(1:3*nnodes_do,1);
u2 = Cmdo*Cmdo'*u2i;
%
Axz = u2-u1;

%%%%
%% Compute Rhs :
% Solve 1
f1 = [ zeros(3*nnodes_do,1) ; C1d_do'*Crdo*Crdo'*u1_do ];
uin1 = K1d_do\f1; resno = norm(K1d_do*uin1-f1)/norm(f1);
u1i = uin1(1:3*nnodes_do);
u1 = Cmdo*Cmdo'*u1i;
% Solve 2
f2 = f2ddo;
uin2 = K2d_do\f2; resno = norm(K1d_do*uin1-f1)/norm(f1);
u2i = uin2(1:3*nnodes_do);
u2 = Cmdo*Cmdo'*u2i;
%
b = u1-u2;
%%
Res(:,1) = b - Axz;
Zed(:,1) = Res(:,1);
d(:,1) = Zed(:,1);

residualdo(1) = norm(Res( :,1));
%errordo(1)    = norm(Itere(indexado) - f1_do(indexado)) / norm(f1_do(indexado));
regularido(1) = norm(Iteredo);

for iter = 1:niter
    %% Optimal step
    % Solve 1
    f1 = [d(:,iter); zeros(nbloq1d_do,1)];
    uin1 = K1d_do\f1;
    u1i = uin1(1:3*nnodes_do,1);
    u1 = Cmdo*Cmdo'*u1i;
    % Solve 2
    f2 = [d(:,iter); zeros(nbloq2d_do,1)];
    uin2 = K2d_do\f2;
    u2i = uin2(1:3*nnodes_do,1);
    u2 = Cmdo*Cmdo'*u2i;
    %
    Ad(:,iter) = u2-u1;
    
    den = (d(:,iter)'*Ad(:,iter));
    d(:,iter) = d(:,iter)/sqrt(den); Ad(:,iter) = Ad(:,iter)/sqrt(den);
    num = Res(:,iter)'*d(:,iter);
    
    Iteredo       = Iteredo + d(:,iter)*num;%/den;
    Res(:,iter+1) = Res(:,iter) - Ad(:,iter)*num;%/den;
    
    Zed(:,iter+1) = Res(:,iter+1);
    
    residualdo(iter+1) = norm(Res(:,iter+1));
    %errordo(iter+1)    = norm(Itere(indexa) - f1_do(indexa)) / norm(f1_do(indexa));
%    regulari(iter+1) = sqrt( Itere'*regul(Itere, nodes, boundary, 2) );
    regularido(iter+1) = norm(Iteredo);

%    % First Reorthogonalize the residual (as we use it next), in sense of M
%    for jter=1:iter-1
%        betac = Zed(:,iter+1)'*Res(:,jter) / (Zed(:,jter)'*Res(:,jter));
%        Zed(:,iter+1) = Zed(:,iter+1) - Zed(:,jter) * betac;
%    end
    
    %% Orthogonalization
    d(:,iter+1) = Zed(:,iter+1);
    
    for jter=1:iter % No need to reorthogonalize (see above)
        betaij = ( Zed(:,iter+1)'*Ad(:,jter) );%/...
            %( d(:,jter)'*Ad(:,jter) );
        d(:,iter+1) = d(:,iter+1) - d(:,jter) * betaij;
    end
    
end

disp([ 'Inverse problem solved ', num2str(toc) ]);

%% Convergence curve
%figure;
%hold on;
%%plot(error,'Color','blue'); % The error is wrong
%plot(residualdo,'Color','black');
%legend('residual');

% L-curve
figure;
loglog( residualdo(2:end) , regularido(2:end) );

% Final problem
f_doD = [ zeros(3*nnodes_do,1) ; C1d_do'*C1d_do*C1d_do'*u1_do ];
f_doN = [Iteredo; zeros(nbloq1d_do,1)];
u_do  = K1d_do\(f_doN+f_doD);

error_do = norm( u_do(indexa) - u1_do(indexa) ) / norm( u1_do(indexa) );

%figure;
%ret = patch('Faces',patch2_do,'Vertices',nodes_do, ... 
%            'FaceVertexCData',f_doN(3:3:end),'FaceColor','interp');
%colorbar;

%figure;
%ret = patch('Faces',patch2_do,'Vertices',nodes_do, ... 
%            'FaceVertexCData',u_do(3:3:end),'FaceColor','interp');
%colorbar;

%%%%%%%%%%%%%%%%%%
%% Change the coordinates
Q = [ cos(alpha), 0, sin(alpha)
      0, 1, 0
      -sin(alpha), 0, cos(alpha) ];
      
nodes2_up = transpose(Q*nodes_up');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Soustract fields

u1_doonup = passMesh3D ( nodes_do, elements_do,...
                        nodes_up, elements_up, u1_do(1:3*nnodes_do) );
u1_diff = u1_up(1:3*nnodes_up)-u1_doonup;

u_doonup = passMesh3D ( nodes_do, elements_do,...
                        nodes_up, elements_up, u_do(1:3*nnodes_do) );
u_diff = u_up(1:3*nnodes_up)-u_doonup;

ui = reshape(u_diff,3,[])'; ui = transpose(Q*ui');
u_diff = u_diff-u_diff;
u_diff(1:3:end-2) = ui(:,1); u_diff(2:3:end-1) = ui(:,2); u_diff(3:3:end) = ui(:,3);

ui = reshape(u1_diff,3,[])'; ui = transpose(Q*ui');
u1_diff = u1_diff-u1_diff;
u1_diff(1:3:end-2) = ui(:,1); u1_diff(2:3:end-1) = ui(:,2); u1_diff(3:3:end) = ui(:,3);

figure;
ret = patch('Faces',patch2_up,'Vertices',nodes2_up, ... 
            'FaceVertexCData',u1_diff(3:3:end),'FaceColor','interp');
colorbar;

figure;
ret = patch('Faces',patch2_up,'Vertices',nodes2_up, ... 
            'FaceVertexCData',u_diff(3:3:end),'FaceColor','interp');
colorbar;

Y = 0:.1:10; Y = Y';
X = 4*ones(101,1); Z = (1+1e-8)*ones(101,1);

u_lin = passMesh3D ( nodes_up, elements_up,...
                     [X,Y,Z], [], u_diff );

u1_lin = passMesh3D ( nodes_up, elements_up,...
                     [X,Y,Z], [], u1_diff );
figure;
hold on;
plot(Y,u_lin(3:3:end),'Color','blue');
plot(Y,u1_lin(3:3:end),'Color','black');
legend('identified gap','reference');

error_tot = .5*(error_up+error_do);

csvwrite('fields/sp3d.csv',[Y,u_lin(3:3:end),u1_lin(3:3:end)]);