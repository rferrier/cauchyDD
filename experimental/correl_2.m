% Transformations de données de corrélation
clear all;
close all;

E     = 70000; Ey = E;
nu    = 0.3;
alpha = 100/E^2;%1000/E^2; % Regularization parameter
mat   = [0,E,nu];
niter = 30;

load('e3_r0_000-015-Mesh.mat'); %15

nodes_c = Mesh.Znode;
eleme_c = Mesh.TRI;
u_c = U;

% Read the mesh
nnodes = size(nodes_c,1); nelem = size(eleme_c,1);
nodes = zeros(nnodes,2);
nodes(:,1) = real(nodes_c); nodes(:,2) = imag(nodes_c);
elements = eleme_c;

 ntoelem = zeros(size(nodes,1), 1);
 for i=1:size(elements,1)
     for j=1:3
         ntoelem(elements(i,j),1) = ntoelem(elements(i,j),1) + 1;
     end
 end
 
 boundary = []; nbound = size(boundary,1);
 order = 1;

% Read u
u = u_c;

% Export the mesh on gmsh
fmid = fopen('meshes/correlation/correli_mesh.msh','w');
fprintf(fmid,'%s\n%s\n%s\n','$MeshFormat','2.2 0 8','$EndMeshFormat');

% Nodes
fprintf(fmid,'%s\n','$Nodes');
fprintf(fmid,'%d\n',nnodes);
for n=1:nnodes
    fprintf(fmid,'%d %d %d %d \n',n,nodes(n,1),nodes(n,2),0);   
end
fprintf(fmid,'%s\n','$EndNodes');

% Elements
fprintf(fmid,'%s\n','$Elements');
fprintf(fmid,'%d\n',nelem+nbound);
for n=1:nbound
    fprintf(fmid,'%d %d %d %d %d %d %d \n',n,1,2,boundary(n,1),...
        boundary(n,1),boundary(n,2),boundary(n,3));   
end
for n=1:nelem
    fprintf(fmid,'%d %d %d %d %d %d %d %d \n',n+nbound,2,2,5,1,...
        elements(n,1),elements(n,2),elements(n,3));   
end
fprintf(fmid,'%s\n','$EndElements');

% Compute the stiffness matrix & such
[K,C,nbloq] = Krig2 (nodes,elements,mat,order,boundary,[]);

% Transform u
P = null(full(K));
u  = u - P*inv(P'*P)*P'*u;

% toremove = [254:264,212:222,233:243];
toremove = [500:521,558:578,532:547];
toremove = [2*toremove-1,2*toremove];

K1 = K; K1(toremove,:) = [];

% Equilibrium gap
egg = u'*(K1'*K1)*u;

% for i=1:100
    Op = eye(2*nnodes) + alpha* (K1'*K1);
    v = Op\u;
    residual = norm(v-u)/norm(u);
    regul    = v'*K*v;
% end
% loglog(residual,regul);
% bug
sigma  = stress(u,Ey,nu,nodes,elements,order,1,ntoelem);
sigmav  = stress(v,Ey,nu,nodes,elements,order,1,ntoelem);

% Export u
plotGMSH({u,'Vect_U';sigma,'stress'},elements, nodes, 'output/correlation');
plotGMSH({v,'Vect_U';sigmav,'stress'},elements, nodes, 'output/correlation_reg');

% Build the bound
% b2 = [269:297,207];
% b31 = [212:222]; b32 = [254:264];
% b4 = [233:243];
b2 = [585:629,494];
b31 = [500:512]; b32 = [565:578];
b4 = [532:547];

bo2 = zeros(size(b2,2)-1,2);
bo31 = zeros(size(b31,2)-1,2); bo32 = zeros(size(b32,2)-1,2);
bo4 = zeros(size(b4,2)-1,2);

bo2(:,1) = b2(1:end-1);   bo2(:,2) = b2(2:end);
bo31(:,1) = b31(1:end-1); bo31(:,2) = b31(2:end);
bo32(:,1) = b32(1:end-1); bo32(:,2) = b32(2:end);
bo4(:,1) = b4(1:end-1);   bo4(:,2) = b4(2:end);
bo3 = [bo31;bo32];

boundary = zeros(size(bo2,1)+size(bo3,1)+size(bo4,1),3);
i = 0;
boundary(i+1:i+size(bo2,1),1) = 2;
boundary(i+1:i+size(bo2,1),2:3) = bo2;
i = i+size(bo2,1);
boundary(i+1:i+size(bo3,1),1) = 3;
boundary(i+1:i+size(bo3,1),2:3) = bo3;
i = i+size(bo3,1);
boundary(i+1:i+size(bo4,1),1) = 4;
boundary(i+1:i+size(bo4,1),2:3) = bo4;

uref = u; urefu = u;

% Compare regularized and un-regularized stuff
fu = K*u(1:2*nnodes); fv = K*v;
% figure;
% hold on;
% plot(fu(2*b32-1),'Color','red');
% plot(fv(2*b32-1),'Color','blue');

% figure;
% hold on;
% plot(fu(2*b2-1),'Color','red');
% plot(fv(2*b2-1),'Color','blue');

% Solve the validation problem to compute residual for DIC
dirichlet = [4,1,0 ; 4,2,0 ; 3,1,0 ; 3,2,0];
[K1,C,nbloq,node2c1,c2node1] = Krig2 (nodes,elements,mat,order,boundary,dirichlet);
f4 = dirichletRhs2(urefu, 4, c2node1, boundary, nnodes );
f3 = dirichletRhs2(urefu, 3, c2node1, boundary, nnodes );

f = (f4+f3);
ucomp = K1\f; ucomp = ucomp(1:2*nnodes);
DICresidual = norm(ucomp-urefu)/norm(urefu);
plotGMSH({(ucomp-urefu)/max(abs(urefu)),'Vect_U'},elements, nodes, 'output/correlation_residual');

%% KMF Stuff
E      = Ey; % because it got ecrased
Kinter = K;
M      = mass_mat(nodes, elements);
Mb     = Mgrad( 2*nnodes, nodes, boundary, 3 );

% init :
u1    = uref-uref;
u2    = u1;
fri   = u1;
v     = u1;
theta = ones(niter+1,1); % First relaxation parameter

% DN problem
dirichlet1 = [4,1,0;4,2,0;
              3,1,0;3,2,0];
neumann1   = [];
[K1,C1,nbloq1,node2c1,c2node1] = Krig (nodes,elements,E,nu,order,boundary,dirichlet1);
f1in = dirichletRhs2(uref, 4, c2node1, boundary, nnodes );

% ND problem
dirichlet2 = [4,1,0;4,2,0;
              1,1,0;1,2,0;
              2,1,0;2,2,0];
neumann2   = [];
[K2,C2,nbloq2,node2c2,c2node2] = Krig (nodes,elements,E,nu,order,boundary,dirichlet2);
f2in = dirichletRhs2(uref, 4, c2node2, boundary, nnodes );

error1   = zeros(niter,1);
error2   = zeros(niter,1);
residual = zeros(niter,1);
regulari = zeros(niter,1);
serror1   = zeros(niter,1); % Error for sigma
serror2   = zeros(niter,1);
sresidual = zeros(niter,1);
ferror1   = zeros(niter,1); % Error for reaction force
ferror2   = zeros(niter,1);
fresidual = zeros(niter,1);

for iter = 1:niter
    % Solve DN
    u1o = u1;
    fdir = dirichletRhs2(v, 3, c2node1, boundary, nnodes );
    f1 = f1in + fdir;

    uin1 = K1\f1;
    u1 = uin1(1:2*nnodes,1);
    lagr1 = uin1(2*nnodes+1:end,1);
    frio = fri; % Store fri for the residual computation
    fri = lagr2forces2( lagr1, c2node1, 3, boundary, nnodes );
    
    error1(iter) = sqrt( myps(u1-urefu,u1-urefu,Kinter,boundary,M,nodes)/...
        myps(urefu,urefu,Kinter,boundary,M,nodes) );
    
    % Solve ND
    u2o = u2;
    fr = [ fri; zeros(size(C2,2),1) ]; % Give to fr the right size
    fdir2 = dirichletRhs2(urefu, 2, c2node2, boundary, nnodes);
    f2 = fr + fdir2 + f2in;

    uin2 = K2\f2;
    u2 = uin2(1:2*nnodes,1);
    
    vo = v;
    v = theta(iter)*u2 + (1-theta(iter))*vo;
    
    error2(iter) = sqrt( myps(u2-urefu,u2-urefu,Kinter,boundary,M,nodes)/...
        myps(urefu,urefu,Kinter,boundary,M,nodes) );

    if iter > 1
       residual(iter) = sqrt( myps(u1-u2,u1-u2,Kinter,boundary,M,nodes)/...
                        sqrt(myps(u1,u1,Kinter,boundary,M,nodes)*...
                        myps(u2,u2,Kinter,boundary,M,nodes)) );
    end
                 
    regulari(iter) = u2(1:2*nnodes)'*Mb*u2(1:2*nnodes);
end

sigma2  = stress(u2,Ey,nu,nodes,elements,order,1,ntoelem);
plotGMSH({u2,'Vect_U';sigma2,'stress'},elements, nodes, 'output/correlation_id');

figure
plot(error2);

% figure;
% hold on;
% plot(u2(2*b32),'Color','red');
% plot(u(2*b32),'Color','blue');

% figure;
% hold on;
% plot(u2(2*b32-1),'Color','red');
% plot(urefu(2*b32-1),'Color','blue');

f2 = Kinter*u2(1:2*nnodes); fref = Kinter*urefu; fru = Kinter*u;

figure;
hold on;
plot(f2(2*b32-1),'Color','red');
plot(fref(2*b32-1),'Color','blue');

figure;
hold on;
plot(f2(2*b31-1),'Color','red');
plot(fref(2*b31-1),'Color','blue');

index = [2*b31,2*b32,2*b31-1,2*b32-1];

% figure;
% hold on;
% plot(u2(index),'Color','red');
% plot(urefu(index),'Color','blue');
% plot(u(index),'Color','green');
% legend('Cauchy id','reg_correl','correlation')
% 
% figure;
% hold on;
% plot(f2(index),'Color','red');
% plot(fref(index),'Color','blue');
% plot(fru(index),'Color','green');
% legend('Cauchy id','reg_correl','correlation')
