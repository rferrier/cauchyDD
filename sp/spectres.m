% 30/03/2017
% Calculs de spectres

close all;
clear all;

addpath(genpath('./tools'))

% Parameters
E       = 70000;  % MPa : Young modulus
nu      = 0.3;    % Poisson ratio

mat = [2, E, nu, 1, 1];
%mat = [0, E, nu];

% Import the mesh
[ nodes,elements,ntoelem,boundary,order ] = readmesh( 'meshes/platet.msh' );
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
[K,C,nbloq] = Krig2 (nodes,elements,mat,order,boundary,[]);

Kinter = K(1:2*nnodes, 1:2*nnodes);
M      = mass_mat(nodes, elements);
[node2b4, b2node4] = mapBound(4, boundaryp1, nnodes);
[node2b3, b2node3] = mapBound(3, boundaryp1, nnodes);
[node2b1, b2node1] = mapBound(1, boundaryp1, nnodes);
[node2b2, b2node2] = mapBound(2, boundaryp1, nnodes);
index = [2*b2node2-1;2*b2node2];
%%
%K1 = Kinter;
%K1(:,[2*b2node1-1;2*b2node1;2*b2node3-1;2*b2node3;2*b2node4-1;2*b2node4]) = [];
%K1([2*b2node1-1;2*b2node1;2*b2node3-1;2*b2node3;2*b2node4-1;2*b2node4],:) = [];
%K2 = Kinter;
%K2(:,[2*b2node1-1;2*b2node1;2*b2node3-1;2*b2node3]) = [];
%K2([2*b2node1-1;2*b2node1;2*b2node3-1;2*b2node3],:) = [];
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

%hold on;
%plot(log10(svd(K1k)),'Color','red');
%plot(log10(svd(Stot)),'Color','magenta');
%plot(log10(svd(K2k)),'Color','blue');
%plot(log10(svd(Dtot)),'Color','cyan');
%plot(log10(svd(D1*K1k)),'Color','green');
%plot(log10(svd(S2*K2k)),'Color','black');
%plot(log10(svd(Korth)),'Color','yellow');
%legend('Sd-Sn','Sd-Sn cancel','Dn-Dd','Dn-Dd cancel','Dd(Sd-Sn)','Sn(Dn-Dd)','Korth')

hold on;
plot(log10(svd(K1k)),'Color','red');
plot(log10(svd(K2k)),'Color','blue');
plot(log10(svd(D1*K1k)),'Color','green');
plot(log10(svd(S2*K2k)),'Color','black');
legend('Sd-Sn','Dn-Dd','Dd(Sd-Sn)','Sn(Dn-Dd)')
