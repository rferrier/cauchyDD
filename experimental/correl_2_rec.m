% Transformations de données de corrélation
clear all;
close all;

E     = 70000; Ey = E;
nu    = 0.3;
mat   = [0,E,nu];

% load("../MASTER/Master-longi_ter/img_0001-0002-Mesh.mat");
%load("../MASTER/Master-longi_sec/img_0001-0002-Mesh.mat");
load("../MASTER/Master-trans_sec/img_0010-0011-Mesh.mat");
%load("../MASTER/Master-longi/img_0001-0002-Mesh.mat");

nodes_c = Mesh.Znode;
eleme_c = Mesh.TRI;


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

[K,C,nbloq] = Krig2 (nodes,elements,mat,order,boundary,[]);

% Rigid modes
R = rigidModes(nodes);
P = eye(2*nnodes) - R*((R'*R)\R'); % Orthogonal projector to modes

% Deformation modes
uxx = zeros(2*nnodes,1); uxx(1:2:2*nnodes-1) = nodes(:,1); uxx = P*uxx;
uyy = zeros(2*nnodes,1); uyy(2:2:2*nnodes)   = nodes(:,2); uyy = P*uyy;
% Exx = uxx*((uxx'*uxx)\uxx'); % Projector onto uxx
% Eyy = uyy*((uyy'*uyy)\uyy');
xmax = max(nodes(:,1)); xmin = min(nodes(:,1)); xmoy = .5*(xmax+xmin); Lx = xmax-xmin;
ymax = max(nodes(:,2)); ymin = min(nodes(:,2)); ymoy = .5*(ymax+ymin); Ly = ymax-ymin;
a = .5*Lx/2; b = .5*Ly/2;

% Restricted deformation modes
uxr = zeros(2*nnodes,1); uyr = zeros(2*nnodes,1);
for i=1:nnodes
   rad = (nodes(i,1)-xmoy)^2/a^2 + (nodes(i,2)-ymoy)^2/b^2;
   if rad < 10
      uxr(2*i-1) = nodes(i,1);
      uyr(2*i)   = nodes(i,2);
   end
end
uxr = P*uxr; uyr = P*uyr;

u_c = zeros(2*nnodes,39);
Exx = zeros(1,39); Eyy = zeros(1,39);
for i=1:31 % 39
   ind = i + 10;
   no = num2str( ind );
   if ind<10
      no = strcat("0",no);
   end
    load( strcat("../MASTER/Master-trans_sec/img_0010-00",no,"-Mesh.mat") );
%    load( strcat("../MASTER/Master-longi_ter/img_0001-00",no,"-Mesh.mat") );
%   load( strcat("../MASTER/Master-longi_sec/img_0001-00",no,"-Mesh.mat") );
   %load( strcat("../MASTER/Master-longi/img_0001-00",no,"-Mesh.mat") );
   % Read u
   u_c(:,i) = P*U; % Remove rigid modes
   u = u_c(:,i);

%    Exx(i) = ((uxx'*uxx)\uxx') * u;
%    Eyy(i) = ((uyy'*uyy)\uyy') * u;
   Exx(i) = ((uxr'*uxr)\uxr') * u;
   Eyy(i) = ((uyr'*uyr)\uyr') * u;
   
   %sigma  = stress(u,Ey,nu,nodes,elements,order,1,ntoelem);
   %f = K*u;
end

% Exx : 22-27
% Eyy : 12-22
t = 1:39;
ixx = 9:17;%27:33;%9:17;%27:33;
iyy = 9:17;%27:33;%9:17;%27:33;
pxx = polyfit(t(ixx),Exx(ixx),1); kxx = pxx(1);
pyy = polyfit(t(iyy),Eyy(iyy),1); kyy = pyy(1);

figure;
hold on;
plot(Exx,'Color','blue');
plot( t(ixx), pxx(2) + pxx(1)*t(ixx),'Color','blue' );
plot(Eyy,'Color','red');
plot( t(iyy), pyy(2) + pyy(1)*t(ixx),'Color','red' );
legend('Exx','','Eyy','');

uchoice = 17;%22
ux22 = uxx*Exx(uchoice); uy22 = uyy*Eyy(uchoice); %22
sigma  = stress(u_c(:,uchoice),Ey,nu,nodes,elements,order,1,ntoelem);
sx = sigma(1:3:end-2); sy = sigma(2:3:end-1); sxy = sigma(3:3:end);

% Export u
%plotGMSH({u,'Vect_U';f,'Vect_F';sigma,'stress'},elements, nodes, 'output/correlation');
plotGMSH({u_c(:,uchoice),'Vect_U';sx,'sx';sy,'sy';sxy,'sxy';...
          ux22,'Uxx';uy22,'Uyy'},elements, nodes, 'output/correlation');
