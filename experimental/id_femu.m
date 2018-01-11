clear all;
close all;

% Identification de matrice homogène orthotrope FEMU

Ex0   = 70000; % MPa reference young
Ey0   = 70000; % MPa reference young
nuxy0 = .3;
Gxy0  = Ex0/2*(1+nuxy0);    % reference Shear
dt    = .01;   % Differential step
maxit = 100;
br    = 0.0;   % Noise

S0 = [1/Ex0, -nuxy0/Ex0, 0;...
      -nuxy0/Ex0, 1/Ey0, 0;...
      0, 0, 1/Gxy0];
Sm10 = inv(S0);

%% reference problem
fscalar = 250;
Ex = 90000; Ey = 50000; nuxy = .33; Gxy = 30000; mat = [1, Ex, Ey, nuxy, Gxy];

S = [1/Ex, -nuxy/Ex, 0;...
     -nuxy/Ex, 1/Ey, 0;...
     0, 0, 1/Gxy];
Sm1 = inv(S);
Kref = [Sm1(1,1);Sm1(1,2);Sm1(2,2);Sm1(3,3)];

dirichlet = [ 1,1,0 ; 1,2,0 ];
neumann   = [ 3,2,fscalar ];
neumann1  = [ 3,2,fscalar ];
neumann2  = [ 3,1,fscalar ];
neumann3  = [ 2,1,fscalar ; 4,1,-fscalar ];
%neumann   = [ 3,1,fscalar ; 3,2,fscalar ; 2,1,fscalar ];
[ nodes,elements,ntoelem,boundary,order] = readmesh( 'meshes/plate.msh' );
nnodes = size(nodes,1);
[K,C,nbloq,node2c,c2node] = Krig2 (nodes,elements,mat,order,boundary,dirichlet);
f  = loading(nbloq,nodes,boundary,neumann);
f1 = loading(nbloq,nodes,boundary,neumann1);
f2 = loading(nbloq,nodes,boundary,neumann2);
f3 = loading(nbloq,nodes,boundary,neumann3);

umes = K\[f1,f2,f3];
umes = umes(1:2*nnodes,:);

amp1 = sqrt(mean(umes(:,1).^2)); % Amplitudes (for noise)
amp2 = sqrt(mean(umes(:,2).^2));
amp3 = sqrt(mean(umes(:,3).^2));

umes(:,1) = umes(:,1) + amp1*br*rand(2*nnodes,1);
umes(:,2) = umes(:,2) + amp2*br*rand(2*nnodes,1);
umes(:,3) = umes(:,3) + amp3*br*rand(2*nnodes,1);

%umes(:,1) = umes(:,1) .* (1+br*rand(2*nnodes,1));
%umes(:,2) = umes(:,2) .* (1+br*rand(2*nnodes,1));
%umes(:,3) = umes(:,3) .* (1+br*rand(2*nnodes,1));

plotGMSH({umes(:,1),'Ux';umes(:,2),'Uxy';umes(:,3),'Uy'}, elements, nodes, 'output/reference');
umes = umes(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The dofs in K : [1,2,0 ; 2,3,0 ; 0,0,4]

k110 = Sm10(1,1); % Give the initial parameters
k120 = Sm10(1,2);
k220 = Sm10(2,2);
k330 = Sm10(3,3);

%k110 = 1.1*Kref(1); % Give the initial parameters
%k120 = 1.1*Kref(2);
%k220 = 1.1*Kref(3);
%k330 = 1.1*Kref(4);

Kdofs0 = [k110;k120;k220;k330];
Kdofs  = [k110;k120;k220;k330];

dtd    = dt*Kdofs0;

%% Inverse identification
ndof = size(Kdofs,1);
uG   = cell(ndof,1);
res  = zeros(maxit,1);
err  = zeros(maxit,1);

S = zeros(6*nnodes,ndof);

Grad = zeros(ndof,1);
Hess = zeros(ndof,1);

for i = 1:maxit
   % Compute the dots
   Ktest = Kdofs;
   mat = [1.5, Ktest(1), Ktest(2), Ktest(3), Ktest(4)];
   [K,C,nbloq,node2c,c2node] = Krig2 (nodes,elements,mat,order,boundary,dirichlet);
   u0 = K\[f1,f2,f3]; u0 = u0(1:2*nnodes,:); u0 = u0(:);
   phi1 = norm(u0(:,1)-umes(:,1))^2;
%   phi2 = norm(u0(:,2)-umes(:,2))^2;
%   phi3 = norm(u0(:,2)-umes(:,3))^2;
%   res(i) = phi1+phi2+phi3;
   res(i) = phi1;
   
   for j=1:ndof
      Ktest    = Kdofs;
      Ktest(j) = Ktest(j) + dtd(j);
      mat = [1.5, Ktest(1), Ktest(2), Ktest(3), Ktest(4)];
      [K,C,nbloq,node2c,c2node] = Krig2 (nodes,elements,mat,order,boundary,dirichlet);
      uGj = K\[f1,f2,f3]; uGj = uGj(1:2*nnodes,:); uGj = uGj(:);
      S1(:,j) = 1/dtd(j) .* (umes(:,1)-uGj(:,1));
%      S2(:,j) = 1/dtd(j) .* (umes(:,2)-uGj(:,2));
%      S3(:,j) = 1/dtd(j) .* (umes(:,3)-uGj(:,3));
   end
      
%   H = (S1'*S1 + S2'*S2 + S3'*S3);
%   b = (S1'*(umes(:,1)-u0(:,1)) + S2'*(umes(:,2)-u0(:,2)) + S3'*(umes(:,3)-u0(:,3)));
   H = ( S1'*S1 );
   b = ( S1'*(umes(:,1)-u0(:,1)) );
   
   % Solve
   dKdof = - H\b;
   Kdofs = Kdofs + dKdof;
   
   err(i) = norm(Kdofs-Kref)/norm(Kref);
end

res = res/res(1);

figure;
hold on;
plot(res,'Color','red');
plot(err);
legend('residual','error');

%% Recovery of the parameters
Sm1f = [ Kdofs(1),Kdofs(2),0 ; Kdofs(2),Kdofs(3),0 ; 0,0,Kdofs(4) ];
Sf = inv(Sm1f);

Exf   = 1/Sf(1,1);
Eyf   = 1/Sf(2,2);
Gxyf  = Kdofs(4);
nuxyf = - Exf*Sf(1,2);
