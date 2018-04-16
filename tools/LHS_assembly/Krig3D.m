function [K,C,ntot,node2c,c2node] = Krig3D (nodes, elem, mat, order,...
    boundary, bc, varargin)
 % Computes the global stiffness matrix
 % input : nodes    : list of nodes : [x1,y1;x2,y2;.....]
 %         elem     : list of elements
 %         mat      : Material : 0=isotropic, 1=orthotropic
 %         order    : order of the elements
 %         bc       : boundary conditions
 %         varargin : if 1 : plane deformations else plane constrains
 %         varargin2: multiplicator for C (better condition number)
 
 % mat = [0, E, nu] if isotropic
 % mat = [1, Ex, Ey, nuxy, Gxy] TODO : add rotation of the axis
 
 % output : K    : global stiffness matrix
 %          ntot : number of nodes on the uy=0 border

 % Only for debug purposes :
 %Kin = zeros(2*size(nodes,1), 2*size(nodes,1));
 
 % First build the model's stiffness matrix
 mul = 1;
 if numel(varargin)>1
     mul = cell2mat(varargin(2));
 end
 
 if mat(1) == 0
    E = mat(2); nu = mat(3);
    S = 1/E*[1,-nu,-nu,0,0,0 ; -nu,1,-nu,0,0,0 ; -nu,-nu,1,0,0,0 ;...
             0,0,0,2*(1+nu),0,0 ; 0,0,0,0,2*(1+nu),0 ; 0,0,0,0,0,2*(1+nu)];
    Sm1 = inv(S);
 elseif mat(1) == 1
    error('Orthotropic material not implemented in 3D (altrough its easier than in 2D...)')
 else
    error('The required material behaviour is not implemented')
 end
 
% Kin = sparse(3*size(nodes,1), 3*size(nodes,1)); % 'cause I love when things go fast
 ndof = 3*size(nodes(elem(1,:),:),1); % /!\ Assume that all the elements have the same ndof
 indI = zeros( 1, size(elem,1)*ndof^2 ); 
 indJ = zeros( 1, size(elem,1)*ndof^2 );
 Koef = zeros( 1, size(elem,1)*ndof^2 );
 for i=1:size(elem,1)
     Xloc1 = nodes(elem(i,:),:);    % Extract and adapt coords
     nnodes = size(Xloc1,1);
     Xloc = zeros(3*nnodes,1);
     
     map = zeros(1,3*nnodes);
     for j=1:nnodes
         map(1,[3*j-2,3*j-1,3*j]) = [3*elem(i,j)-2, 3*elem(i,j)-1, 3*elem(i,j)];
         Xloc([3*j-2,3*j-1,3*j],1) = [Xloc1(j,1);Xloc1(j,2);Xloc1(j,3)];
     end

     Ke = stifmat3D(Xloc,order,Sm1,0);

     % Build K
     map1 = map'*ones(1,size(map,2)); map1 = map1(:);
     map2 = ones(size(map,2),1)*map;  map2 = map2(:);
%     Kin(map,map) = Kin(map,map) + Ke;
     indI( (i-1)*ndof^2 + (1:ndof^2) ) = map1;
     indJ( (i-1)*ndof^2 + (1:ndof^2) ) = map2;
     Koef( (i-1)*ndof^2 + (1:ndof^2) ) = Ke(:);
 end

 Kin = sparse( indI, indJ, Koef, 3*size(nodes,1), 3*size(nodes,1) );
 
%  % Boundary conditions with Lagrange multiplicators
%  d_lines = [];
%  for i = 1:size(bc,1)
%      d_lines(i,1) = bc(i,1);
%  end

 % First compute the number of elements with Dirichlet bc
% inodes = zeros(size(nodes,1),4); % Vector that stores witch coords are imposed
% ntot = 0;
% for i=1:size(boundary,1)
%     entity = boundary(i,1);
%     for j=1:size(bc,1)
%         if bc(j,1) == entity
%             for k=2:size(boundary,2)
%
%                 if bc(j,2) == 1     % Mark this node as x
%                     inodes(boundary(i,k),1) = 1;
%                 elseif bc(j,2) == 2 % Mark this node as y
%                     inodes(boundary(i,k),2) = 1;
%                 elseif bc(j,2) == 3 % Mark this node as z
%                     inodes(boundary(i,k),3) = 1;
%                 else
%                     error('Unable to read Dirichlet BC axis')
%                 end
%                 
%             end
%         end
%     end
% end
% 
% % Number of dof with BCs (an other loop because corner)
% for i=1:size(nodes,1)
%     if inodes(i,1) ~= 0
%         ntot = ntot+1;
%     end
%     if inodes(i,2) ~= 0
%         ntot = ntot+1;
%     end
%     if inodes(i,3) ~= 0
%         ntot = ntot+1;
%     end
% end
%
% node2c = zeros( 1, 3*size(nodes,1) ); % This list makes lines of C match with ddls
% c2node = zeros( 1, ntot );
%
% for j=1:size(bc,1)
%     if bc(j,1) == 0 % Regularisation/x or y
%         ntot = ntot+1;
%     end
% end
% 
% C = zeros(3*size(nodes,1),ntot);
%
% % Second time : populate C (I don't like growing matrixes)
% j=1;
% for i=1:size(nodes,1)
%     if inodes(i,1) ~= 0
%         C(3*i-2,j) = inodes(i,1);
%         node2c(1,3*i-2) = j;
%         c2node(1,j) = 3*i-2;
%         j=j+1;
%     end
%     if inodes(i,2) ~= 0
%         C(3*i-1,j) = inodes(i,2);
%         node2c(1,3*i-1) = j;
%         c2node(1,j) = 3*i-1;
%         j=j+1;
%     end
%     if inodes(i,3) ~= 0
%         C(3*i,j) = inodes(i,3);
%         node2c(1,3*i) = j;
%         c2node(1,j) = 3*i;
%         j=j+1;
%     end
% end
% 
% % Regularisation on x : sum(ux) = 0
% for k=1:size(bc,1)
%     if bc(k,1) == 0 && bc(k,2) == 1 % Regularisation/x
%         for i=1:size(nodes,1)
%             C(3*i-2,j) = C(3*i-2,j) + 1;
%         end
%         j=j+1;
%     end
%     if bc(k,1) == 0 && bc(k,2) == 2 % Regularisation/y
%         for i=1:size(nodes,1)
%             C(3*i-1,j) = C(3*i-1,j) + 1;
%         end
%         j=j+1;
%     end
%     if bc(k,1) == 0 && bc(k,2) == 3 % Regularisation/z
%         for i=1:size(nodes,1)
%             C(3*i,j) = C(3*i,j) + 1;
%         end
%         j=j+1;
%     end
%     if bc(k,1) == 0 && bc(k,2) == 4 % Regularisation/thetax
%         % First, find the ~=barycenter of the solid
%         z0 = sum( nodes(:,3) ) / size(nodes,3);
%         y0 = sum( nodes(:,2) ) / size(nodes,2);
%         for i=1:size(nodes,1)
%             z = nodes(i,3);
%             y = nodes(i,2);
%
%             C(3*i,j) = C(3*i,j) + (y-y0);
%             C(3*i-1,j) = C(3*i-1,j) - (z-z0);
%         end
%         j=j+1;
%     end
%     if bc(k,1) == 0 && bc(k,2) == 5 % Regularisation/thetay
%         % First, find the ~=barycenter of the solid
%         z0 = sum( nodes(:,3) ) / size(nodes,3);
%         x0 = sum( nodes(:,1) ) / size(nodes,1);
%         for i=1:size(nodes,1)
%             z = nodes(i,3);
%             x = nodes(i,1);
%
%             C(3*i,j) = C(3*i,j) + (x-x0);
%             C(3*i-2,j) = C(3*i-2,j) - (z-z0);
%         end
%         j=j+1;
%     end
%      if bc(k,1) == 0 && bc(k,2) == 6 % Regularisation/thetaz
%         % First, find the ~=barycenter of the solid
%         x0 = sum( nodes(:,1) ) / size(nodes,1);
%         y0 = sum( nodes(:,2) ) / size(nodes,2);
%         for i=1:size(nodes,1)
%             x = nodes(i,1);
%             y = nodes(i,2);
%
%             C(3*i-2,j) = C(3*i,j) + (y-y0);
%             C(3*i-1,j) = C(3*i-1,j) - (x-x0);
%         end
%         j=j+1;
%     end
% end
 [ C, node2c, c2node, ntot ] = Cbound ( nodes, bc, boundary );
 K = [Kin,mul*C;mul*C',zeros(ntot)];
 
end