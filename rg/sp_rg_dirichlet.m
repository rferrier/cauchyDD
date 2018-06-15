% 22/05/2018
% Détection de fissure 1D plane par Cauchy puis écart à la réciprocité domaine carré, condition de Dirichlet

close all;
clear all;

% Parameters
E           = 210000; % MPa : Young modulus
nu          = 0.3;    % Poisson ratio
fscalar     = 250;    % N.mm-1 : Loading on the plate
mat         = [0, E, nu];
br          = .0;      % Noise level

[ nodesg, elementsg, ntoelemg, boundaryg, order, physical ] = readmesh( 'meshes/rg_sp_squared/plate_c_squared3_noref.msh' );
nnodesg = size(nodesg,1);

%% Direct problem
dirichlet  = [8,1,0 ; 8,2,0];

neumann1   = [1,2,-fscalar];
neumann2   = [2,1,fscalar ; 7,1,fscalar ; 4,1,-fscalar; 9,1,-fscalar];
neumann3   = [1,1,-fscalar ; 1,2,-fscalar ; 4,1,-fscalar ; 4,2,-fscalar; 9,1,-fscalar ; 9,2,-fscalar];
neumann4   = [1,1,fscalar ; 1,2,-fscalar ; 2,1,fscalar ; 2,2,-fscalar ; 7,1,fscalar ; 7,2,-fscalar];

%neumann1   = [8,2,fscalar ; 1,2,-fscalar];
%neumann2   = [2,1,fscalar ; 7,1,fscalar ; 4,1,-fscalar; 9,1,-fscalar];
%neumann3   = [8,1,fscalar ; 8,2,fscalar ; 2,1,fscalar ; 2,2,fscalar ; 7,1,fscalar ; 7,2,fscalar ; ...
%              1,1,-fscalar ; 1,2,-fscalar ; 4,1,-fscalar ; 4,2,-fscalar; 9,1,-fscalar ; 9,2,-fscalar];
%neumann4   = [8,1,-fscalar ; 8,2,fscalar ; 2,1,fscalar ; 2,2,-fscalar ; 7,1,fscalar ; 7,2,-fscalar ; ...
%              1,1,fscalar ; 1,2,-fscalar ; 4,1,-fscalar ; 4,2,fscalar ; 9,1,-fscalar ; 9,2,fscalar];

[K,C,nbloq,node2c,c2node] = Krig2 (nodesg,elementsg,mat,order,boundaryg,dirichlet);
Kinter = K( 1:2*nnodesg, 1:2*nnodesg );

Xmax = max(nodesg(:,1)); Xmin = min(nodesg(:,1)); Xmoy = (Xmax+Xmin)/2;
Ymax = max(nodesg(:,2)); Ymin = min(nodesg(:,2)); Ymoy = (Ymax+Ymin)/2;
Lx = Xmax-Xmin; Ly = Ymax-Ymin;

f1  = loading(nbloq,nodesg,boundaryg,neumann1);
f2  = loading(nbloq,nodesg,boundaryg,neumann2);
f3  = loading(nbloq,nodesg,boundaryg,neumann3);
f4  = loading(nbloq,nodesg,boundaryg,neumann4);

uin = K\[f1,f2,f3,f4];
u1 = uin(1:2*nnodesg,1); u2 = uin(1:2*nnodesg,2);
u3 = uin(1:2*nnodesg,3); u4 = uin(1:2*nnodesg,4);
f1 = Kinter*u1; f2 = Kinter*u2; f3 = Kinter*u3; f4 = Kinter*u4;

ui = reshape(u1,2,[])';  ux = ui(:,1);  uy = ui(:,2);

u1ref = u1; u2ref = u2; u3ref = u3; u4ref = u4;
f1ref = f1; f2ref = f2; f3ref = f3; f4ref = f4;

sigma = stress(u1,E,nu,nodesg,elementsg,order,1,ntoelemg);
plotGMSH({u1,'U1';u2,'U2';u3,'U3';u4,'U4'}, elementsg, nodesg, 'output/reference');

%% Split the mesh

indexu = find( physical == 2); indexd = find( physical == 1 );
elementsu = elementsg(indexu,:); elementsd = elementsg(indexd,:);

nodesu_glo = unique(elementsu(:)); % Global numbering of the nodes of the upper part.
[~, nodesg_up] = ismember (1:nnodesg, nodesu_glo ); % Local numbering of the global nodes
nnodesu = size(nodesu_glo,1);
nodesd_glo = unique(elementsd(:)); % Global numbering of the nodes of the lower part.
[~, nodesg_do] = ismember (1:nnodesg, nodesd_glo ); % Local numbering of the global nodes
nnodesd = size(nodesd_glo,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the reference position of the crack
x1 = nodesg(5,1); y1 = nodesg(5,2);
x2 = nodesg(6,1); y2 = nodesg(6,2);

xmin = min(nodesg(:,1)); xmax = max(nodesg(:,1));
ymin = min(nodesg(:,2)); ymax = max(nodesg(:,2));

% 4 intersections of the line with the boundaries
t1 = (xmin-x1)/(x2-x1); t2 = (xmax-x1)/(x2-x1);
t3 = (ymin-y1)/(y2-y1); t4 = (ymax-y1)/(y2-y1);
xy11 = [xmin;y1+(y2-y1)*t1];
xy22 = [xmax;y1+(y2-y1)*t2];
xy33 = [x1+(x2-x1)*t3;ymin];
xy44 = [x1+(x2-x1)*t4;ymax];
xy1234 = [xy11,xy22,xy33,xy44];
% limit to those that are inside the square
elim1 = find(xy1234(1,:)>1); elim2 = find(xy1234(2,:)>1);
elim3 = find(xy1234(1,:)<0); elim4 = find(xy1234(2,:)<0);
total = setdiff( [1,2,3,4] , union(union(elim1,elim2),union(elim3,elim4)) );
xyt = xy1234(:,total); % Normally, this one should be of size 2
xy1 = xyt(:,1); xy2 = xyt(:,2); xy1r = xy1; xy2r = xy2;

%% Recover the reference
step = (xy2r-xy1r)/100;
n = [-step(2);step(1)]; n = n/norm(n); % Normal
nodes3r = [ xy1r(1):step(1):xy2r(1) ; xy1r(2):step(2):xy2r(2) ];
nodes3r = nodes3r'; nnodes3r = size(nodes3r,1);
curvr = sqrt( (nodes3r(:,1)-nodes3r(1,1)).^2 + (nodes3r(:,2)-nodes3r(1,2)).^2 );
%% Check if need to reverse the nodes (for visu purposes)
%if norm(nodes3r(end,:)-nodes3b(1,:)) < norm(nodes3r(1,:)-nodes3b(1,:))
%   nodes3r = nodes3r(end:-1:1,:);
%end
nodes3s = nodes3r + 1e-3*ones(size(nodes3r,1),1)*n';
nodes3i = nodes3r - 1e-3*ones(size(nodes3r,1),1)*n';
urs = passMesh2D( nodesg, elementsg, nodes3s, [], [u1,u2,u3,u4] );
uri = passMesh2D( nodesg, elementsg, nodes3i, [], [u1,u2,u3,u4] );
urg = -(uri-urs);  % Vectorial gap
urg([1,2,end-1,end],:) = 0; % Overwrite the strange stuff that comes from the fact that we get out of the domain

rn1 = urg(1:2:end-1,1)*n(1) + urg(2:2:end,1)*n(2); % Reference normal gaps
rn2 = urg(1:2:end-1,2)*n(1) + urg(2:2:end,2)*n(2);
rn3 = urg(1:2:end-1,3)*n(1) + urg(2:2:end,3)*n(2);
rn4 = urg(1:2:end-1,4)*n(1) + urg(2:2:end,4)*n(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First step: resolution of the Cauchy problem

% Use the upper mesh
nnodes    = nnodesu;
nodes     = nodesg(nodesu_glo,:);
elements  = nodesg_up(elementsg);
toremove1 = find(elements(:,1)==0); toremove2 = find(elements(:,2)==0); toremove3 = find(elements(:,3)==0);
toremove = union(toremove1,union(toremove2,toremove3));
elements(toremove,:) = [];

boundary = [ boundaryg(:,1) , nodesg_up(boundaryg(:,[2,3])) ];
toremove1 = find(boundary(:,3)==0); toremove2 = find(boundary(:,2)==0);
toremove = union(toremove1,toremove2);
boundary(toremove,:) = [];

% List of dofs (same as nodes, but with the 2*nnodes-1;2*nnodes)
dofu_glo = [2*nodesu_glo-1,2*nodesu_glo];
dofu_loc = [(1:2:2*nnodes-1)',(2:2:2*nnodes)'];

% Recover the fields
u1u = zeros(2*nnodes,1); u2u = zeros(2*nnodes,1);
u3u = zeros(2*nnodes,1); u4u = zeros(2*nnodes,1);
u1u(dofu_loc) = u1(dofu_glo); u2u(dofu_loc) = u2(dofu_glo);
u3u(dofu_loc) = u3(dofu_glo); u4u(dofu_loc) = u4(dofu_glo);

[Kinter,~,~,~,~] = Krig2 (nodes,elements,mat,order,boundary,[]);
f1u = Kinter*u1u; f2u = Kinter*u2u; f3u = Kinter*u3u; f4u = Kinter*u4u;

%% Debug plot
%figure;
%patch('Faces',elements,'Vertices',nodes,'FaceAlpha',0);
%nodesdefo = nodes + 1e2*reshape(u2u,2,[])';
%patch('Faces',elements,'Vertices',nodesdefo,'FaceVertexCData',u1u(2:2:2*nnodes),'FaceColor','interp');

%% Build the requested operators
% First, remove the corner points in boundary
xmax = max(nodes(:,1)); xmin = min(nodes(:,1));
ymax = max(nodes(:,2)); ymin = min(nodes(:,2));
no1  = findNode(xmin, ymin, nodes, 1e-5);
no2  = findNode(xmax, ymin, nodes, 1e-5);
no3  = findNode(xmax, ymax, nodes, 1e-5);
no4  = findNode(xmin, ymax, nodes, 1e-5);
boundaryp1 = suppressBound( boundary, no1, 3 );
boundaryp1 = suppressBound( boundaryp1, no2, 3 );
boundaryp1 = suppressBound( boundaryp1, no1, 9 );
boundaryp1 = suppressBound( boundaryp1, no4, 9 );
boundaryp1 = suppressBound( boundaryp1, no2, 7 );
boundaryp1 = suppressBound( boundaryp1, no3, 7 );

% Operators
dirichlet1d = [ 7,1,0;7,2,0 ; 8,1,0;8,2,0 ; 9,1,0;9,2,0 ];
[K1d,C1d,nbloq1d,node2c1d,c2node1d] = Krig2 (nodes,elements,mat,order,boundary,dirichlet1d);
dirichlet2d = [ 8,1,0;8,2,0 ];
[K2d,C2d,nbloq2d,node2c2d,c2node2d] = Krig2 (nodes,elements,mat,order,boundary,dirichlet2d);

% Restriction on 3
dirichletdummy = [ 3,1,0;3,2,0 ];
[~,Cu,~,~,~] = Krig2 (nodes,elements,mat,order,boundaryp1,dirichletdummy);
% Restriction on 7 & 9
dirichletdummy = [ 7,1,0;7,2,0 ; 9,1,0;9,2,0 ];
[~,Ck,~,~,~] = Krig2 (nodes,elements,mat,order,boundaryp1,dirichletdummy);

% Dirichlet and Neumann Rhs
%fD = [ zeros(2*nnodes,4) ; C1d'*C1d*C1d'*[u1u,u2u,u3u,u4u] ]; % Dirichlet rhs on the redondant boundary (for K1d)
%fN = [ Ck*Ck'*[f1u,f2u,f3u,f4u] ; C2d'*C2d*C2d'*[u1u,u2u,u3u,u4u] ]; % Neumann rhs on the redondant boundary
fD = [ zeros(2*nnodes,2) ; C1d'*C1d*C1d'*[u1u,u3u] ];
fN = [ Ck*Ck'*[f1u,f3u] ; C2d'*C2d*C2d'*[u1u,u3u] ];
%fD = [ zeros(2*nnodes,2) ; C1d'*C1d*C1d'*[u1u,u2u] ];
%fN = [ Ck*Ck'*[f1u,f2u] ; C2d'*C2d*C2d'*[u1u,u2u] ];
%fD = [ zeros(2*nnodes,1) ; C1d'*C1d*C1d'*u4u ];
%fN = [ Ck*Ck'*f2u ; C2d'*C2d*C2d'*u4u ];

clear Ck; clear C2d; clear C1d;
% Call the SPD function
tic
[uc,ind] = spd_mrhs(2*nnodes, nbloq1d, nbloq2d, K1d, K2d, fD, fN, Cu, 20, 0, 0, Cu*Cu'*[f1u,f3u]);%,f4u]);
fc = Kinter*uc;
toc

fc1 = fc(:,1); fc3 = fc(:,2);
uc1 = uc(:,1); uc3 = uc(:,2);

clear K1d; clear K2d;

%% Debug plot
%figure;
%patch('Faces',elements,'Vertices',nodes,'FaceAlpha',0);
%nodesdefo = nodes + 1e2*reshape(uc1,2,[])';
%patch('Faces',elements,'Vertices',nodesdefo,'FaceVertexCData',uc1(2:2:2*nnodes),'FaceColor','interp');

% Comparison
[~, b2node3] = mapBound(3, boundaryp1, nnodes);
u1ref = u1u([2*b2node3-1, 2*b2node3]); u1id = uc1([2*b2node3-1, 2*b2node3]);
error1u = norm(u1ref-u1id,'fro')/norm(u1ref,'fro');
f1ref = f1u([2*b2node3-1, 2*b2node3]); f1id = fc1([2*b2node3-1, 2*b2node3]);
error1f = norm(f1ref-f1id,'fro')/norm(f1ref,'fro');

% And on the other side
[~, b2node8] = mapBound(8, boundaryp1, nnodes);
f1ref8 = f1u([2*b2node8-1, 2*b2node8]); f1id8 = fc1([2*b2node8-1, 2*b2node8]);
error1f8 = norm(f1ref8-f1id8,'fro')/norm(f1ref8,'fro');

try
figure;
hold on;
plot(u1ref(:,1),'Color', 'red');
plot(u1id(:,1),'Color', 'blue');
legend('U : reference','U : solution');
end

try
figure;
hold on;
plot(f1ref(:,2),'Color', 'red');
plot(f1id(:,2),'Color', 'blue');
legend('F : reference','F : solution');
end

try
figure;
hold on;
plot(f1ref8(:,2),'Color', 'red');
plot(f1id8(:,2),'Color', 'blue');
legend('F : referenceon the top','F : solution on the top');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Identify the crack on the large domain

% Use the global mesh
nnodes   = nnodesg;
nodes    = nodesg;
elements = elementsg;
boundary = boundaryg;

% Transpose F (U is known : no need to re-transpose it)
zenodes = nodesu_glo(b2node8); %nodes at which the f should be replaced
toreplace = [2*zenodes-1, 2*zenodes]; % Dofs in the upper numerotation

f1N = f1; f3N = f3;
f1N(toreplace) = fc1([2*b2node8-1, 2*b2node8]);
f3N(toreplace) = fc3([2*b2node8-1, 2*b2node8]);

%plotGMSH({f1-f1N,'F1';f3-f3N,'F3'}, elementsg, nodesg, 'output/fields'); % Debug plot

%% Manage the boundary stuff
boundary(find(boundary(:,1)==3),:) = []; % Remove the interior bounds
boundary(find(boundary(:,1)==5),:) = [];
boundary(find(boundary(:,1)==6),:) = [];

nboun1 = size(boundary,1); nelem1 = size(elements,1);
extnorm1 = zeros( nboun1, 2 );
frr = zeros(2*nboun1,2);
for i=1:nboun1
   % Volumic element
   no1 = boundary(i,2); no2 = boundary(i,3); % only with 2 nodes even if order > 1
   cand1 = rem( find(elements==no1),nelem1 ); % find gives line + column*size
   cand2 = rem( find(elements==no2),nelem1 );
   elt = intersect(cand1, cand2); % If everything went well, there is only one
  
   % Exterior normal
   no3 = setdiff( elements( elt, 1:3 ), [no1,no2]);
   x1 = nodes(no1,1); y1 = nodes(no1,2);
   x2 = nodes(no2,1); y2 = nodes(no2,2);
   x3 = nodes(no3,1); y3 = nodes(no3,2);
   extnorm1(i,:) = [ y1-y2 , -(x1-x2) ];
   extnorm1(i,:) = extnorm1(i,:)/norm(extnorm1(i,:));
   if (x3-x2)*(y1-y2) - (y3-y2)*(x1-x2) > 0 % Check that the normal is exterior
      extnorm1(i,:) = -extnorm1(i,:);
   end
   
   % Value of F
   len = sqrt( (x1-x2)^2 + (y1-y2)^2 );
   frr([2*i-1,2*i],1) = 1/(2*len)*(f1N([2*no1-1,2*no1]) + f1N([2*no2-1,2*no2])); % Pass f on the element
   frr([2*i-1,2*i],2) = 1/(2*len)*(f3N([2*no1-1,2*no1]) + f3N([2*no2-1,2*no2]));
end

% Call the function
tic
[normal,Cte1,Cte2,ug,x] = rg_poly_crack_2d( nodes, extnorm1, boundary, 1, mat, [u1,u3], frr, 3, 1, 0, 2 );
toc
Q = [ normal(2), normal(1) ; - normal(1), normal(2) ];

%%%%

%% Compute the computed position of the crack
Vp1 = [2;-Cte1]; Vp2 = [-2;-Cte1];
Vm1 = [2;-Cte2]; Vm2 = [-2;-Cte2];
vp1 = Q*Vp1; vm1 = Q*Vm1; vp2 = Q*Vp2; vm2 = Q*Vm2;

xmin = min(nodesg(:,1)); xmax = max(nodesg(:,1));
ymin = min(nodesg(:,2)); ymax = max(nodesg(:,2));

%% First one
x1 = vp1(1); y1 = vp1(2);
x2 = vp2(1); y2 = vp2(2);
% 4 intersections of the line with the boundaries
t1 = (xmin-x1)/(x2-x1); t2 = (xmax-x1)/(x2-x1);
t3 = (ymin-y1)/(y2-y1); t4 = (ymax-y1)/(y2-y1);
xy11 = [xmin;y1+(y2-y1)*t1];
xy22 = [xmax;y1+(y2-y1)*t2];
xy33 = [x1+(x2-x1)*t3;ymin];
xy44 = [x1+(x2-x1)*t4;ymax];
xy1234 = [xy11,xy22,xy33,xy44];
% limit to those that are inside the square
elim1 = find(xy1234(1,:)>1); elim2 = find(xy1234(2,:)>1);
elim3 = find(xy1234(1,:)<0); elim4 = find(xy1234(2,:)<0);
total = setdiff( [1,2,3,4] , union(union(elim1,elim2),union(elim3,elim4)) );
xyt = xy1234(:,total); % Normally, this one should be of size 2
if max(size(total)) < 2
   warning('crack is outside the domain');
   xyp1 = [0;0]; xyp2 = [0;0];
else
   xyp1 = xyt(:,1); xyp2 = xyt(:,2);
end

%% Second one
x1 = vm1(1); y1 = vm1(2);
x2 = vm2(1); y2 = vm2(2);
% 4 intersections of the line with the boundaries
t1 = (xmin-x1)/(x2-x1); t2 = (xmax-x1)/(x2-x1);
t3 = (ymin-y1)/(y2-y1); t4 = (ymax-y1)/(y2-y1);
xy11 = [xmin;y1+(y2-y1)*t1];
xy22 = [xmax;y1+(y2-y1)*t2];
xy33 = [x1+(x2-x1)*t3;ymin];
xy44 = [x1+(x2-x1)*t4;ymax];
xy1234 = [xy11,xy22,xy33,xy44];
% limit to those that are inside the square
elim1 = find(xy1234(1,:)>1); elim2 = find(xy1234(2,:)>1);
elim3 = find(xy1234(1,:)<0); elim4 = find(xy1234(2,:)<0);
total = setdiff( [1,2,3,4] , union(union(elim1,elim2),union(elim3,elim4)) );
xyt = xy1234(:,total); % Normally, this one should be of size 2
if max(size(total)) < 2
   warning('crack is outside the domain');
   xym1 = [0;0]; xym2 = [0;0];
else
   xym1 = xyt(:,1); xym2 = xyt(:,2);
end

% Vizualize the crack's line
x6 = nodes(6,1); y6 = nodes(6,2); x5 = nodes(5,1); y5 = nodes(5,2); 
try
figure; hold on;
x1m = xym1(1); y1m = xym1(2);
x2m = xym2(1); y2m = xym2(2);
x1p = xyp1(1); y1p = xyp1(2);
x2p = xyp2(1); y2p = xyp2(2);
x1r = xy1r(1); y1r = xy1r(2);
x2r = xy2r(1); y2r = xy2r(2);
plot( [xmin,xmax], [ymin,ymin], 'Color', 'black');
plot( [xmax,xmax], [ymin,ymax], 'Color', 'black');
plot( [xmax,xmin], [ymax,ymax], 'Color', 'black');
plot( [xmin,xmin], [ymax,ymin], 'Color', 'black');
plot( [x1r,x2r], [y1r,y2r], 'Color', 'black', 'LineWidth', 3 );
plot( [x5,x6], [y5,y6] ,'Color', 'magenta', 'LineWidth',5);
plot( [x1m,x2m], [y1m,y2m], 'Color', 'green', 'LineWidth', 3 );
plot( [x1p,x2p], [y1p,y2p], 'Color', 'red', 'LineWidth', 3 );
axis equal;
end

%% Outputs plots
%% Plot the normal
%try
%figure
%hold on;
%%ret = patch('Faces',elements2(:,1:3),'Vertices',nodes2,'FaceAlpha',0);
%x6 = nodes(6,1); y6 = nodes(6,2); x5 = nodes(5,1); y5 = nodes(5,2); 
%Xmin = min(nodes(:,1)); Xmax = max(nodes(:,1));
%Ymin = min(nodes(:,2)); Ymax = max(nodes(:,2));

%xc = .5*(x6+x5);
%yc = .5*(y6+y5);

%% Rectangle
%plot( [Xmin,Xmax], [Ymin,Ymin], 'Color', 'black' );
%plot( [Xmax,Xmax], [Ymin,Ymax], 'Color', 'black' );
%plot( [Xmax,Xmin], [Ymax,Ymax], 'Color', 'black' );
%plot( [Xmin,Xmin], [Ymax,Ymin], 'Color', 'black' );

%x2 = xc + .5*normal(1); y2 = yc + .5*normal(2);
%plot( [xc,x2], [yc,y2] ,'Color', 'red', 'LineWidth',3);
%plot( [x5,x6], [y5,y6] ,'Color', 'magenta', 'LineWidth',3);
%axis('equal');
%end

%% Plot the line
%try
%figure
%hold on;

%Vp1 = [2;-Cte1]; Vp2 = [-2;-Cte1];
%Vm1 = [2;-Cte2]; Vm2 = [-2;-Cte2];
%vp1 = Q*Vp1; vm1 = Q*Vm1; vp2 = Q*Vp2; vm2 = Q*Vm2;

%% Rectangle
%plot( [Xmin,Xmax], [Ymin,Ymin], 'Color', 'black' );
%plot( [Xmax,Xmax], [Ymin,Ymax], 'Color', 'black' );
%plot( [Xmax,Xmin], [Ymax,Ymax], 'Color', 'black' );
%plot( [Xmin,Xmin], [Ymax,Ymin], 'Color', 'black' );

%plot( [vp1(1), vp2(1)], [vp1(2), vp2(2)] ,'Color', 'red', 'LineWidth',3);
%plot( [vm1(1), vm2(1)], [vm1(2), vm2(2)] ,'Color', 'green', 'LineWidth',3);
%plot( [x5,x6], [y5,y6] ,'Color', 'magenta', 'LineWidth',4);
%axis('equal');
%end

% Plot the gap
try
figure
hold on;
plot(curvr,rn1,'Color','red');
plot(x-x(1),ug*sign(normal'*n),'Color','blue');
legend('Reference','Solution');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Identify the crack on the small domain
% Use the lower mesh
nnodes    = nnodesd;
nodes     = nodesg(nodesd_glo,:);
elements  = nodesg_do(elementsg);
toremove1 = find(elements(:,1)==0); toremove2 = find(elements(:,2)==0); toremove3 = find(elements(:,3)==0);
toremove = union(toremove1,union(toremove2,toremove3));
elements(toremove,:) = [];

boundary = [ boundaryg(:,1) , nodesg_do(boundaryg(:,[2,3])) ];
toremove1 = find(boundary(:,3)==0); toremove2 = find(boundary(:,2)==0);
toremove = union(toremove1,toremove2);
boundary(toremove,:) = [];

% List of dofs (same as nodes, but with the 2*nnodes-1;2*nnodes)
dofd_glo = [2*nodesd_glo-1,2*nodesd_glo];
dofd_loc = [(1:2:2*nnodes-1)',(2:2:2*nnodes)'];

zenodes = nodesu_glo(b2node3); %nodes at which the f should be replaced
toreplace = [2*zenodes-1, 2*zenodes]; % Dofs in the upper numerotation

f1N = f1; f3N = f3; u1D = u1; u3D = u3;
f1N(toreplace) = fc1([2*b2node3-1, 2*b2node3]);
f3N(toreplace) = fc3([2*b2node3-1, 2*b2node3]);
u1D(toreplace) = uc1([2*b2node3-1, 2*b2node3]);
u3N(toreplace) = uc3([2*b2node3-1, 2*b2node3]);

% Reference fields
u1r = zeros(2*nnodes,1); u3r = zeros(2*nnodes,1);
f1r = zeros(2*nnodes,1); f3r = zeros(2*nnodes,1);
u1r(dofd_loc) = u1(dofd_glo); u3r(dofd_loc) = u3(dofd_glo);
[Kinter,~,~,~,~] = Krig2 (nodes,elements,mat,order,boundary,[]);
f1r = Kinter*u1r; f3r = Kinter*u3r;

% Recover the fields
u1d = zeros(2*nnodes,1); u3d = zeros(2*nnodes,1);
f1d = zeros(2*nnodes,1); f3d = zeros(2*nnodes,1);
u1d(dofd_loc) = u1D(dofd_glo); u3d(dofd_loc) = u3D(dofd_glo);
f1d(dofd_loc) = -f1N(dofd_glo); f3d(dofd_loc) = -f3N(dofd_glo); % The - is because f(1->2) = -f(2->1)

%% Debug : replace by the reference fields
%u1d = u1r; u3d = u3r; f1d = f1r; f3d = f3r;

%% Debug plot
%mesh2GMSH( nodes, elements, boundary, 'output/downmesh' );
%plotGMSH({u1d,'U1';u3d,'U3';f1d,'F1';f3d,'F3'}, elements, nodes, 'output/fields');

%% Manage the boundary stuff
boundary(find(boundary(:,1)==5),:) = []; % Remove the crack's bounds
boundary(find(boundary(:,1)==6),:) = [];

nboun1 = size(boundary,1); nelem1 = size(elements,1);
extnorm1 = zeros( nboun1, 2 );
frr = zeros(2*nboun1,2);
for i=1:nboun1
   % Volumic element
   no1 = boundary(i,2); no2 = boundary(i,3); % only with 2 nodes even if order > 1
   cand1 = rem( find(elements==no1),nelem1 ); % find gives line + column*size
   cand2 = rem( find(elements==no2),nelem1 );
   elt = intersect(cand1, cand2); % If everything went well, there is only one
  
   % Exterior normal
   no3 = setdiff( elements( elt, 1:3 ), [no1,no2]);
   x1 = nodes(no1,1); y1 = nodes(no1,2);
   x2 = nodes(no2,1); y2 = nodes(no2,2);
   x3 = nodes(no3,1); y3 = nodes(no3,2);
   extnorm1(i,:) = [ y1-y2 , -(x1-x2) ];
   extnorm1(i,:) = extnorm1(i,:)/norm(extnorm1(i,:));
   if (x3-x2)*(y1-y2) - (y3-y2)*(x1-x2) > 0 % Check that the normal is exterior
      extnorm1(i,:) = -extnorm1(i,:);
   end
   
   % Value of F
   len = sqrt( (x1-x2)^2 + (y1-y2)^2 );
   frr([2*i-1,2*i],1) = 1/(2*len)*(f1d([2*no1-1,2*no1]) + f1d([2*no2-1,2*no2])); % Pass f on the element
   frr([2*i-1,2*i],2) = 1/(2*len)*(f3d([2*no1-1,2*no1]) + f3d([2*no2-1,2*no2]));
end

% Call the function
tic
[normal,Cte1,Cte2,ug,x] = rg_poly_crack_2d( nodes, extnorm1, boundary, 1, mat, [u1d,u3d], frr, 1, 1, 0, 1 );
toc
Q = [ normal(2), normal(1) ; - normal(1), normal(2) ];

%%%%
%% Compute the computed position of the crack
Vp1 = [2;-Cte1]; Vp2 = [-2;-Cte1];
Vm1 = [2;-Cte2]; Vm2 = [-2;-Cte2];
vp1 = Q*Vp1; vm1 = Q*Vm1; vp2 = Q*Vp2; vm2 = Q*Vm2;

xmin = min(nodesg(:,1)); xmax = max(nodesg(:,1));
ymin = min(nodesg(:,2)); ymax = max(nodesg(:,2));

%% First one
x1 = vp1(1); y1 = vp1(2);
x2 = vp2(1); y2 = vp2(2);
% 4 intersections of the line with the boundaries
t1 = (xmin-x1)/(x2-x1); t2 = (xmax-x1)/(x2-x1);
t3 = (ymin-y1)/(y2-y1); t4 = (ymax-y1)/(y2-y1);
xy11 = [xmin;y1+(y2-y1)*t1];
xy22 = [xmax;y1+(y2-y1)*t2];
xy33 = [x1+(x2-x1)*t3;ymin];
xy44 = [x1+(x2-x1)*t4;ymax];
xy1234 = [xy11,xy22,xy33,xy44];
% limit to those that are inside the square
elim1 = find(xy1234(1,:)>1); elim2 = find(xy1234(2,:)>1);
elim3 = find(xy1234(1,:)<0); elim4 = find(xy1234(2,:)<0);
total = setdiff( [1,2,3,4] , union(union(elim1,elim2),union(elim3,elim4)) );
xyt = xy1234(:,total); % Normally, this one should be of size 2
if max(size(total)) < 2
   warning('crack is outside the domain');
   xyp1 = [0;0]; xyp2 = [0;0];
else
   xyp1 = xyt(:,1); xyp2 = xyt(:,2);
end

%% Second one
x1 = vm1(1); y1 = vm1(2);
x2 = vm2(1); y2 = vm2(2);
% 4 intersections of the line with the boundaries
t1 = (xmin-x1)/(x2-x1); t2 = (xmax-x1)/(x2-x1);
t3 = (ymin-y1)/(y2-y1); t4 = (ymax-y1)/(y2-y1);
xy11 = [xmin;y1+(y2-y1)*t1];
xy22 = [xmax;y1+(y2-y1)*t2];
xy33 = [x1+(x2-x1)*t3;ymin];
xy44 = [x1+(x2-x1)*t4;ymax];
xy1234 = [xy11,xy22,xy33,xy44];
% limit to those that are inside the square
elim1 = find(xy1234(1,:)>1); elim2 = find(xy1234(2,:)>1);
elim3 = find(xy1234(1,:)<0); elim4 = find(xy1234(2,:)<0);
total = setdiff( [1,2,3,4] , union(union(elim1,elim2),union(elim3,elim4)) );
xyt = xy1234(:,total); % Normally, this one should be of size 2
if max(size(total)) < 2
   warning('crack is outside the domain');
   xym1 = [0;0]; xym2 = [0;0];
else
   xym1 = xyt(:,1); xym2 = xyt(:,2);
end

% Vizualize the crack's line
try
figure; hold on;
x1m = xym1(1); y1m = xym1(2);
x2m = xym2(1); y2m = xym2(2);
x1p = xyp1(1); y1p = xyp1(2);
x2p = xyp2(1); y2p = xyp2(2);
x1r = xy1r(1); y1r = xy1r(2);
x2r = xy2r(1); y2r = xy2r(2);
plot( [xmin,xmax], [ymin,ymin], 'Color', 'black');
plot( [xmax,xmax], [ymin,ymax], 'Color', 'black');
plot( [xmax,xmin], [ymax,ymax], 'Color', 'black');
plot( [xmin,xmin], [ymax,ymin], 'Color', 'black');
plot( [x1r,x2r], [y1r,y2r], 'Color', 'black', 'LineWidth', 3 );
plot( [x5,x6], [y5,y6] ,'Color', 'magenta', 'LineWidth',5);
plot( [x1m,x2m], [y1m,y2m], 'Color', 'green', 'LineWidth', 3 );
plot( [x1p,x2p], [y1p,y2p], 'Color', 'red', 'LineWidth', 3 );
axis equal;
end

%% Outputs plots
%% Plot the normal
%try
%figure
%hold on;

%plot( [Xmin,Xmax], [Ymin,Ymin], 'Color', 'black' );
%plot( [Xmax,Xmax], [Ymin,Ymax], 'Color', 'black' );
%plot( [Xmax,Xmin], [Ymax,Ymax], 'Color', 'black' );
%plot( [Xmin,Xmin], [Ymax,Ymin], 'Color', 'black' );

%x2 = xc + .5*normal(1); y2 = yc + .5*normal(2);
%plot( [xc,x2], [yc,y2] ,'Color', 'red', 'LineWidth',3);
%plot( [x5,x6], [y5,y6] ,'Color', 'magenta', 'LineWidth',3);
%axis('equal');
%end

%% Plot the line
%try
%figure
%hold on;

%Vp1 = [2;-Cte1]; Vp2 = [-2;-Cte1];
%Vm1 = [2;-Cte2]; Vm2 = [-2;-Cte2];
%vp1 = Q*Vp1; vm1 = Q*Vm1; vp2 = Q*Vp2; vm2 = Q*Vm2;

%% Rectangle
%plot( [Xmin,Xmax], [Ymin,Ymin], 'Color', 'black' );
%plot( [Xmax,Xmax], [Ymin,Ymax], 'Color', 'black' );
%plot( [Xmax,Xmin], [Ymax,Ymax], 'Color', 'black' );
%plot( [Xmin,Xmin], [Ymax,Ymin], 'Color', 'black' );

%plot( [vp1(1), vp2(1)], [vp1(2), vp2(2)] ,'Color', 'red', 'LineWidth',3);
%plot( [vm1(1), vm2(1)], [vm1(2), vm2(2)] ,'Color', 'green', 'LineWidth',3);
%plot( [x5,x6], [y5,y6] ,'Color', 'magenta', 'LineWidth',4);
%axis('equal');
%end

% Plot the gap
try
figure
hold on;
plot(curvr,rn1,'Color','red');
plot(x-x(1),ug*sign(normal'*n),'Color','blue');
legend('Reference','Solution');
end
