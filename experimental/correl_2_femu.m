% Transformations de données de corrélation et FEMU
clear all;
close all;

E     = 70000; Ey = E;
nu    = 0.3; mat   = [0,E,nu];
Ex0   = 70000; % MPa reference young
Ey0   = 70000; % MPa reference young
nuxy0 = .3;
Gxy0  = Ex0/(2*(1+nuxy0));    % reference Shear
dt    = .01;   % Differential step
maxit = 170;
crit  = 1e-3;  % Sensibility kernel criterion

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

%% Export the mesh on gmsh
%fmid = fopen('meshes/correlation/correli_mesh.msh','w');
%fprintf(fmid,'%s\n%s\n%s\n','$MeshFormat','2.2 0 8','$EndMeshFormat');
%
%% Nodes
%fprintf(fmid,'%s\n','$Nodes');
%fprintf(fmid,'%d\n',nnodes);
%for n=1:nnodes
%    fprintf(fmid,'%d %d %d %d \n',n,nodes(n,1),nodes(n,2),0);   
%end
%fprintf(fmid,'%s\n','$EndNodes');
%
%% Elements
%fprintf(fmid,'%s\n','$Elements');
%fprintf(fmid,'%d\n',nelem+nbound);
%for n=1:nbound
%    fprintf(fmid,'%d %d %d %d %d %d %d \n',n,1,2,boundary(n,1),...
%        boundary(n,1),boundary(n,2),boundary(n,3));   
%end
%for n=1:nelem%1000:1120%1:nelem
%    fprintf(fmid,'%d %d %d %d %d %d %d %d \n',n+nbound,2,2,5,1,...
%        elements(n,1),elements(n,2),elements(n,3));   
%end
%fprintf(fmid,'%s\n','$EndElements');
%%for n=1:nelem%1000:1120%1:nelem
%%    fprintf(fmid,'%d %d %d %d %d %d %d %d \n',n+nbound,2,2,5,1,...
%%        elements(n,1),elements(n,2),elements(n,3));   
%%end

% Compute the stiffness matrix & such
[K,C,nbloq] = Krig2 (nodes,elements,mat,order,boundary,[]);

% Transform u
n1 = zeros(2*nnodes,1); n1(1:2:2*nnodes-1) = 1;
n2 = zeros(2*nnodes,1); n2(2:2:2*nnodes) = 1;
n3 = zeros(2*nnodes,1);
n3(1:2:2*nnodes-1) = -nodes(:,2); n3(2:2:2*nnodes) = nodes(:,1);

%P = null(full(K));
P = [n1,n2,n3];
u = u - P*inv(P'*P)*P'*u;

% toremove = [254:264,212:222,233:243];
toremove = [500:521,558:578,532:547];
toremove = [2*toremove-1,2*toremove];

K1 = K; K1(toremove,:) = [];

% Equilibrium gap
egg = u'*(K1'*K1)*u;

v = u;

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
ball = [494:629];

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

umes = uref;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FEMU

[K,C,nbloq,node2c,c2node] = Krig2 (nodes,elements,mat,order,boundary,dirichlet);
f = [ zeros(2*nnodes,1) ; (C'*C)\C'*umes ];


S0 = [1/Ex0, -nuxy0/Ex0, 0;...
      -nuxy0/Ex0, 1/Ey0, 0;...
      0, 0, 1/Gxy0];
Sm10 = inv(S0);

k110 = Sm10(1,1); % Give the initial parameters
k120 = Sm10(1,2);
k220 = Sm10(2,2);
k330 = Sm10(3,3);

Kdofs0 = [k110;k120;k220;k330];
Kdofs  = [k110;k120;k220;k330];

dtd    = dt*Kdofs0;

%% Inverse identification
ndof  = size(Kdofs,1);
res   = zeros(maxit,1);
err   = zeros(maxit,1);
delta = zeros(maxit,1);
ndK   = zeros(maxit,1);

S1 = zeros(2*nnodes,ndof);
u0 = zeros(2*nnodes,1);

for i = 1:maxit
   % Compute the Gradient and Hessian
   up = u0;
   Ktest = Kdofs;
   mat = [1.5, Ktest(1), Ktest(2), Ktest(3), Ktest(4)];
   [K,C,nbloq,node2c,c2node] = Krig2 (nodes,elements,mat,order,boundary,dirichlet);
   u0 = K\f; u0 = u0(1:2*nnodes,:); u0 = u0(:);
   phi1 = norm(u0(:,1)-umes(:,1))^2;
   res(i)   = norm(u0(:,1)-umes(:,1))/norm(umes);
   delta(i) = norm(up-u0);
   
   for j=1:ndof
      Ktest    = Kdofs;
      Ktest(j) = Ktest(j) + dtd(j);
      mat = [1.5, Ktest(1), Ktest(2), Ktest(3), Ktest(4)];
      [K,C,nbloq,node2c,c2node] = Krig2 (nodes,elements,mat,order,boundary,dirichlet);
      uGj = K\f; uGj = uGj(1:2*nnodes,:); uGj = uGj(:);
      S1(:,j) = 1/dtd(j) .* (umes(:,1)-uGj(:,1));
   end

   H = ( S1'*S1 ); % TODO if needed : put the inverse of correlation matrix here (from correli)
   b = ( S1'*(umes(:,1)-u0(:,1)) ); % TODO : add Tychonov
   
   [ U,Theta,V ] = svd(H);  % At this point, H = U*Theta*U'
   theta = diag(Theta);
   toremove = find( theta/theta(1) < crit );
   N = U(:,toremove);
   P = eye(ndof) - N*((N'*N)\N'); % Orthogonal projector
   
   % Solve
   dKdof  = - H\b; % TODO : Variable relaxation
   Kdofs  = Kdofs + P*dKdof;
   ndK(i) = norm(P*dKdof);
   
%   err(i) = norm(Kdofs-Kref)/norm(Kref);
end

%res = res/res(1);

figure;
plot(res,'Color','red');
legend('residual');

%% Recovery of the parameters
Sm1f = [ Kdofs(1),Kdofs(2),0 ; Kdofs(2),Kdofs(3),0 ; 0,0,Kdofs(4) ];
Sf = inv(Sm1f);

Exf   = 1/Sf(1,1);
Eyf   = 1/Sf(2,2);
Gxyf  = Kdofs(4);
nuxyf = - Exf*Sf(1,2);