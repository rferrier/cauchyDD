% 30/03/2017
% Calculs de spectres

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 70000;  % MPa : Young modulus
nu      = 0.3;    % Poisson ratio
fscalar = 1;      % N.mm-1 : Loading on the plate

mat = [2, E, nu, 1, 1];
%mat = [0, E, nu];

dirichlet = [1,1,0; 1,2,0 ;
             3,1,0; 3,2,0 ];
neumann   = [2,1,fscalar,0,fscalar;
             4,1,fscalar,0,-fscalar];

% Import the mesh
[ nodes,elements,ntoelem,boundary,order ] = readmesh( 'meshes/platee.msh' );
nnodes = size(nodes,1);

% find the nodes in the corners and suppress the element :
xmax = max(nodes(:,1));
xmin = min(nodes(:,1));
ymax = max(nodes(:,2));
ymin = min(nodes(:,2));
no1  = findNode(xmin, ymin, nodes, 1e-5);
no2  = findNode(xmax, ymin, nodes, 1e-5);
no3  = findNode(xmax, ymax, nodes, 1e-5);
no4  = findNode(xmin, ymax, nodes, 1e-5);

boundaryp1 = suppressBound( boundary, no2, 2 );
boundaryp1 = suppressBound( boundaryp1, no3, 2 );
boundaryp1 = suppressBound( boundaryp1, no1, 4 );
boundaryp1 = suppressBound( boundaryp1, no4, 4 );
boundaryp2 = boundaryp1;

[~, b2node2] = mapBound(2, boundaryp1, nnodes);
indexa = [2*b2node2-1; 2*b2node2];

% Then, build the stiffness matrix :
[K,C,nbloq] = Krig2 (nodes,elements,mat,order,boundary,dirichlet);


Kinter = K(1:2*nnodes, 1:2*nnodes);
M      = mass_mat(nodes, elements);
[node2b4, b2node4] = mapBound(4, boundaryp1, nnodes);
[node2b3, b2node3] = mapBound(3, boundaryp1, nnodes);
[node2b1, b2node1] = mapBound(1, boundaryp1, nnodes);
[node2b2, b2node2] = mapBound(2, boundaryp1, nnodes);
index = [2*b2node2-1;2*b2node2];

% Solve the problem
f = loading(nbloq,nodes,boundary,neumann);
uin = K\f;

%uin = uin.*(1+100*randn(2*nnodes+nbloq,1));
%uin = randn(2*nnodes+nbloq,1);
%f = randn(2*nnodes+nbloq,1);

% Extract data
ur = uin([2*b2node4-1;2*b2node4]);
tr = f([2*b2node4-1;2*b2node4]);
um = uin([2*b2node2-1;2*b2node2]);
tm = f([2*b2node2-1;2*b2node2]);

%%
Kbb1 = Kinter([2*b2node2-1;2*b2node2], [2*b2node2-1;2*b2node2]);
Kbb2 = Kinter([2*b2node2-1;2*b2node2], [2*b2node2-1;2*b2node2]);

Kii1 = Kinter; 
Kii1(:,[2*b2node1-1;2*b2node1;2*b2node3-1;2*b2node3;2*b2node4-1;2*b2node4; ...
        2*b2node2-1;2*b2node2]) = []; 
Kii1([2*b2node1-1;2*b2node1;2*b2node3-1;2*b2node3;2*b2node4-1;2*b2node4; ...
      2*b2node2-1;2*b2node2],:) = [];
Kii2 = Kinter; 
Kii2(:,[2*b2node1-1;2*b2node1;2*b2node3-1;2*b2node3; ...
        2*b2node2-1;2*b2node2]) = []; 
Kii2([2*b2node1-1;2*b2node1;2*b2node3-1;2*b2node3; ...
      2*b2node2-1;2*b2node2],:) = [];

Kib1 = Kinter(:, [2*b2node2-1;2*b2node2]);
Kib1([2*b2node1-1;2*b2node1;2*b2node3-1;2*b2node3;2*b2node4-1;2*b2node4; ...
      2*b2node2-1;2*b2node2],:) = [];
Kib2 = Kinter(:, [2*b2node2-1;2*b2node2]); 
Kib2([2*b2node1-1;2*b2node1;2*b2node3-1;2*b2node3; ...
      2*b2node2-1;2*b2node2],:) = [];

Kbi1 = Kinter([2*b2node2-1;2*b2node2],:); 
Kbi1(:,[2*b2node1-1;2*b2node1;2*b2node3-1;2*b2node3;2*b2node4-1;2*b2node4; ...
        2*b2node2-1;2*b2node2]) = [];
Kbi2 = Kinter([2*b2node2-1;2*b2node2],:); 
Kbi2(:,[2*b2node1-1;2*b2node1;2*b2node3-1;2*b2node3; ...
        2*b2node2-1;2*b2node2]) = [];
%%
Iii1 = inv(Kii1); Iii2 = inv(Kii2);
S1 = full(Kbb1 - Kbi1*Iii1*Kib1);
S2 = full(Kbb2 - Kbi2*Iii2*Kib2);
D1 = inv(S1); D2 = inv(S2);

%%
Krr = Kinter([2*b2node4-1;2*b2node4], [2*b2node4-1;2*b2node4]);
Kjj = Kinter;
Kjj([2*b2node1-1;2*b2node1;2*b2node3-1;2*b2node3;2*b2node2-1;2*b2node2;...
    2*b2node4-1;2*b2node4],:) = [];
Kjj(:,[2*b2node1-1;2*b2node1;2*b2node3-1;2*b2node3;2*b2node2-1;2*b2node2;...
    2*b2node4-1;2*b2node4]) = [];

Krj = Kinter( [2*b2node4-1;2*b2node4], : );
Krj( :, [2*b2node1-1;2*b2node1;2*b2node3-1;2*b2node3;2*b2node2-1;2*b2node2;...
    2*b2node4-1;2*b2node4] ) = [];

Kgj = Kinter( [2*b2node2-1;2*b2node2], : );
Kgj( :, [2*b2node1-1;2*b2node1;2*b2node3-1;2*b2node3;2*b2node2-1;2*b2node2;...
    2*b2node4-1;2*b2node4] ) = [];

Ikjj = inv(Kjj);
SR = full(Krr - Krj*Ikjj*Krj');
K1k = Kgj*Ikjj*Krj'*inv(SR)*Krj*Ikjj*Kgj';
%%
Kgg = Kinter( [2*b2node2-1;2*b2node2], [2*b2node2-1;2*b2node2] );

Kjg = Kinter( :,[2*b2node2-1;2*b2node2] );
Kjg( [2*b2node1-1;2*b2node1;2*b2node3-1;2*b2node3;2*b2node2-1;2*b2node2;...
    2*b2node4-1;2*b2node4], : ) = [];
    
Kjr = Kinter( :, [2*b2node4-1;2*b2node4] );
Kjr( [2*b2node1-1;2*b2node1;2*b2node3-1;2*b2node3;2*b2node2-1;2*b2node2;...
    2*b2node4-1;2*b2node4], : ) = [];
    
Ikgg = inv(Kgg);
L = Kjj-Kjg*Ikgg*Kgj; Il = inv(L);
SR2 = Krr-Krj*Il*Kjr;
K2k = Ikgg*Kgj*Il*Kjr*inv(SR2)*Krj*Il*Kjg*Ikgg;
%%
Stot = S1-S2; Dtot = D2-D1;
%cond(Stot)
%cond(Dtot)
Korth = Kgj*Ikjj*Krj';
Korth2 = Ikgg*Kgj*Il*Kjr;

%[Vp,vp] = eig(Korth);
%[Vd,vd] = eig(Korth2);
%vp = diag(abs(vp)); vd = diag(abs(vd));
%[vp,Ind] = sort(vp,'descend');
%Vp = Vp(:,Ind);
%[vd,Ind] = sort(vd,'descend');
%Vd = Vd(:,Ind);
%
%[Vpp,vpp] = eig(D1^(1/2)*Korth);
%[Vdp,vdp] = eig(S2^(1/2)*Korth2);
%vpp = diag(abs(vpp)); vdp = diag(abs(vdp));
%[vpp,Ind] = sort(vpp,'descend');
%Vpp = Vpp(:,Ind);
%[vdp,Ind] = sort(vdp,'descend');
%Vdp = Vdp(:,Ind);

[Vp,vp] = eig(K1k);
[Vd,vd] = eig(full(K2k));
vp = diag(abs(vp)); vd = diag(abs(vd));
[vp,Ind] = sort(vp,'descend');
Vp = Vp(:,Ind);
[vd,Ind] = sort(vd,'descend');
Vd = Vd(:,Ind);

%[Vpp,vpp] = eig(D1*K1k);
%[Vdp,vdp] = eig(S2*K2k);
%vpp = diag(abs(vpp)); vdp = diag(abs(vdp));
%[vpp,Ind] = sort(vpp,'descend');
%Vpp = Vpp(:,Ind);
%[vdp,Ind] = sort(vdp,'descend');
%Vdp = Vdp(:,Ind);

[Vpp,vpp] = eig(K1k,S1); Vpp = Vpp*(Vpp'*S1*Vpp)^(-1/2);
[Vdp,vdp] = eig(full(K2k),D2); Vdp = Vdp*(Vdp'*D2*Vdp)^(-1/2);
vpp = diag(abs(vpp)); vdp = diag(abs(vdp));
[vpp,Ind] = sort(vpp,'descend');
Vpp = Vpp(:,Ind);
[vdp,Ind] = sort(vdp,'descend');
Vdp = Vdp(:,Ind);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RHS
%borth = SR*ur - tr;
%eib = abs(real(Vp)'*borth);
%eibp = abs(real(Vpp)'*D1^(1/2)*borth);
%
%dorth = SR2*ur - tr;
%eid = abs(real(Vd)'*dorth);
%eidp = abs(real(Vdp)'*S2^(1/2)*dorth);

b = Kgj*Ikjj*Krj'*(ur - SR\tr);
eib = real(Vp)'*b;
eibp = real(Vpp)'*(b);

d = Ikgg*Kgj*Il*Kjr*(ur - SR2\tr);
eid = real(Vd)'*d;
eidp = real(Vdp)'*(d);

%hold on;
%plot(log10(svd(K1k)),'Color','red');
%plot(log10(svd(Stot)),'Color','magenta');
%plot(log10(svd(K2k)),'Color','blue');
%plot(log10(svd(Dtot)),'Color','cyan');
%plot(log10(svd(D1*K1k)),'Color','green');
%plot(log10(svd(S2*K2k)),'Color','black');
%plot(log10(svd(Korth)),'Color','yellow');
%legend('Sd-Sn','Sd-Sn cancel','Dn-Dd','Dn-Dd cancel','Dd(Sd-Sn)','Sn(Dn-Dd)','Korth')

%hold on;
%plot(log10(svd(K1k)),'Color','red');
%plot(log10(svd(K2k)),'Color','blue');
%plot(log10(svd(D1*K1k)),'Color','green');
%plot(log10(svd(S2*K2k)),'Color','black');
%legend('Sd-Sn','Dn-Dd','Dd(Sd-Sn)','Sn(Dn-Dd)')

%hold on;
%plot(log10(svd(Stot)),'Color','red');
%plot(log10(svd(K1k)),'Color','blue');
%legend('Sd-Sn','Sd-Sn formule')

%norm(D1^(1/2)*Korth - S1^(1/2)*Korth2)
%norm(D1^(1/2)*Korth)

%hold on;
%plot(log10(vp),'Color','red');
%plot(log10(vpp),'Color','magenta');
%plot(log10(vd),'Color','blue');
%plot(log10(vdp),'Color','cyan');
%legend('Korth primal','Korth precond','Korth dual','Korth precond');

% %% Iterative resolutions
% [x,flag,relres,iter] = pcg(K1k,b,1e-15,20);
% figure; plot(x); relres
% [x,flag,relres,iter] = symmlq(K1k,b,1e-15,20);
% figure; plot(x); relres
% [x,flag,relres,iter] = gmres(K1k,b,[],1e-15,20);
% figure; plot(x); relres
% x = K1k\b;
% figure; plot(x); norm(K1k*x-b)
% bug
%% Primal noprec
chib = eib./vp; chiD = chib; niter = size(chiD,1);
%[indm2,pol] = findPicard2(log10(abs(chiD)), 3, 1);
%n = size(pol,1);
%t = 1:.05:niter; tt = zeros(n,20*(niter-1)+1);
%for j=1:n
%   tt(j,:) = t.^(n-j);
%end
%px = pol'*tt;

%figure;
%hold on;
%plot(log10(vp),'Color','red');
%plot(log10(abs(eib)),'Color','blue');
%plot(log10(abs(eib)./vp),'Color','green');
%%plot(t,px,'Color','cyan')
%legend('primal','rhs','solution');%,'interpolation');
%figure;
%plot(log10(abs(eib)./vp),'Color','black');

%% Primal prec
chibp = eibp./vpp; chiD = chibp; niter = size(chiD,1);
%[indm2,pol] = findPicard2(log10(abs(chiD)), 3, 1);
%n = size(pol,1);
%t = 1:.05:niter; tt = zeros(n,20*(niter-1)+1);
%for j=1:n
%   tt(j,:) = t.^(n-j);
%end
%px = pol'*tt;

%figure;
%hold on;
%plot(log10(vpp),'Color','red');
%plot(log10(abs(eibp)),'Color','blue');
%plot(log10(abs(eibp)./vpp),'Color','green');
%%plot(t,px,'Color','cyan')
%legend('primal prec','rhs','solution');%,'interpolation');

figure;
hold on;
plot(log10(abs(eib)./vp),'Color','blue');
plot(log10(abs(eibp)./vpp),'Color','black');
legend('no prec','prec');

%% Dual noprec
chid = eid./vd; chiD = chid; niter = size(chiD,1);
%[indm2,pol] = findPicard2(log10(abs(chiD)), 3, 1);
%n = size(pol,1);
%t = 1:.05:niter; tt = zeros(n,20*(niter-1)+1);
%for j=1:n
%   tt(j,:) = t.^(n-j);
%end
%px = pol'*tt;

%figure;
%hold on;
%plot(log10(vd),'Color','red');
%plot(log10(abs(eid)),'Color','blue');
%plot(log10(abs(eid)./vd),'Color','green');
%%plot(t,px,'Color','cyan')
%legend('dual','rhs','solution');%,'interpolation');

%figure;
%plot(log10(abs(eid)./vd),'Color','black');

%% Dual prec
chidp = eidp./vdp; chiD = chidp; niter = size(chiD,1);
%[indm2,pol] = findPicard2(log10(abs(chiD)), 3, 1);
%n = size(pol,1);
%t = 1:.05:niter; tt = zeros(n,20*(niter-1)+1);
%for j=1:n
%   tt(j,:) = t.^(n-j);
%end
%px = pol'*tt;

%figure;
%hold on;
%plot(log10(vdp),'Color','red');
%plot(log10(abs(eidp)),'Color','blue');
%plot(log10(abs(eidp)./vdp),'Color','green');
%%plot(t,px,'Color','cyan')
%legend('dual prec','rhs','solution');%,'interpolation');

figure;
hold on;
plot(log10(abs(eid)./vd),'Color','blue');
plot(log10(abs(eidp)./vdp),'Color','black');
legend('no prec','prec');

uptot = Vp*chib;
up = Vp(:,1:30)*chib(1:30); tp = S1*up - Kgj*Ikjj*Krj'*ur;
upp = Vpp(:,1:30)*chibp(1:30); tpp = S1*upp - Kgj*Ikjj*Krj'*ur;
%chib2 = inv(Vpp'*Vpp)*Vpp'*uptot;

tdtot = Vd*chid;
td = Vd(:,1:30)*chid(1:30); ud = D1*td + Ikgg*Kgj*Il*Kjr*ur;
tdp = Vdp(:,1:30)*chidp(1:30); udp = D1*tdp + Ikgg*Kgj*Il*Kjr*ur;
%chid2 = inv(Vdp'*Vdp)*Vdp'*tdtot;

%figure;
%hold on;
%plot(log10(abs(chib)),'Color','blue');
%plot(log10(abs(chib2)),'Color','cyan');
%plot(log10(abs(chibp)),'Color','black');
%legend('no prec','no prec on the energy basis','prec');
%
%figure;
%hold on;
%plot(log10(abs(chid)),'Color','blue');
%plot(log10(abs(chid2)),'Color','cyan');
%plot(log10(abs(chidp)),'Color','black');
%legend('no prec','no prec on the energy basis','prec');

%figure;
%hold on;
%plot(up,'Color','red');
%plot(upp,'Color','magenta');
%plot(um,'Color','blue');
%
%figure;
%hold on;
%plot(tp,'Color','red');
%plot(tpp,'Color','magenta');
%plot(tm,'Color','blue');
%
%figure;
%hold on;
%plot(ud,'Color','red');
%plot(udp,'Color','magenta');
%plot(um,'Color','blue');
%
%figure;
%hold on;
%plot(tp,'Color','blue');
%plot(tpp,'Color','cyan');
%plot(td,'Color','red');
%plot(tdp,'Color','magenta');
%plot(tm,'Color','green');
%legend('SPP', 'SPP+prec', 'SPD', 'SPD+prec', 'ref');

%figure;
%hold on;
%plot(log10(svd(Korth)),'Color','red');
%plot(log10(svd(D1^(1/2)*Korth)),'Color','magenta');
%plot(log10(svd(Korth2)),'Color','blue');
%plot(log10(svd(S2^(1/2)*Korth2)),'Color','cyan');
%legend('Korth primal','Korth precond','Korth dual','Korth precond');
