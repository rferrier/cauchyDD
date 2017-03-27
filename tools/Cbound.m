function [ C, node2c, c2node, ntot ] = Cbound ( nodes, bc, boundary )
% This function computes the trace matrix C such that C'*u = u|bound

 % First compute the number of elements with Dirichlet bc
 inodes = zeros(size(nodes,1),4); % Vector that stores witch coords are imposed
 ntot = 0;
 for i=1:size(boundary,1)
     entity = boundary(i,1);
     for j=1:size(bc,1)
         if bc(j,1) == entity
             for k=2:size(boundary,2)

                 if bc(j,2) == 1     % Mark this node as x
                     inodes(boundary(i,k),1) = 1;
                 elseif bc(j,2) == 2 % Mark this node as y
                     inodes(boundary(i,k),2) = 1;
                 elseif bc(j,2) == 3 % Mark this node as z
                     inodes(boundary(i,k),3) = 1;
                 else
                     error('Unable to read Dirichlet BC axis')
                 end
                 
             end
         end
     end
 end
 
 % Number of dof with BCs (an other loop because corner)
 for i=1:size(nodes,1)
     if inodes(i,1) ~= 0
         ntot = ntot+1;
     end
     if inodes(i,2) ~= 0
         ntot = ntot+1;
     end
     if inodes(i,3) ~= 0
         ntot = ntot+1;
     end
 end

 node2c = zeros( 1, 3*size(nodes,1) ); % This list makes lines of C match with ddls
 c2node = zeros( 1, ntot );

 for j=1:size(bc,1)
     if bc(j,1) == 0 % Regularisation/x or y
         ntot = ntot+1;
     end
 end
 
 C = zeros(3*size(nodes,1),ntot);

 % Second time : populate C (I don't like growing matrixes)
 j=1;
 for i=1:size(nodes,1)
     if inodes(i,1) ~= 0
         C(3*i-2,j) = inodes(i,1);
         node2c(1,3*i-2) = j;
         c2node(1,j) = 3*i-2;
         j=j+1;
     end
     if inodes(i,2) ~= 0
         C(3*i-1,j) = inodes(i,2);
         node2c(1,3*i-1) = j;
         c2node(1,j) = 3*i-1;
         j=j+1;
     end
     if inodes(i,3) ~= 0
         C(3*i,j) = inodes(i,3);
         node2c(1,3*i) = j;
         c2node(1,j) = 3*i;
         j=j+1;
     end
 end
 
 % Regularisation on x : sum(ux) = 0
 for k=1:size(bc,1)
     if bc(k,1) == 0 && bc(k,2) == 1 % Regularisation/x
         for i=1:size(nodes,1)
             C(3*i-2,j) = C(3*i-2,j) + 1;
         end
         j=j+1;
     end
     if bc(k,1) == 0 && bc(k,2) == 2 % Regularisation/y
         for i=1:size(nodes,1)
             C(3*i-1,j) = C(3*i-1,j) + 1;
         end
         j=j+1;
     end
     if bc(k,1) == 0 && bc(k,2) == 3 % Regularisation/z
         for i=1:size(nodes,1)
             C(3*i,j) = C(3*i,j) + 1;
         end
         j=j+1;
     end
     if bc(k,1) == 0 && bc(k,2) == 4 % Regularisation/thetax
         % First, find the ~=barycenter of the solid
         z0 = sum( nodes(:,3) ) / size(nodes,3);
         y0 = sum( nodes(:,2) ) / size(nodes,2);
         for i=1:size(nodes,1)
             z = nodes(i,3);
             y = nodes(i,2);

             C(3*i,j) = C(3*i,j) + (y-y0);
             C(3*i-1,j) = C(3*i-1,j) - (z-z0);
         end
         j=j+1;
     end
     if bc(k,1) == 0 && bc(k,2) == 5 % Regularisation/thetay
         % First, find the ~=barycenter of the solid
         z0 = sum( nodes(:,3) ) / size(nodes,3);
         x0 = sum( nodes(:,1) ) / size(nodes,1);
         for i=1:size(nodes,1)
             z = nodes(i,3);
             x = nodes(i,1);

             C(3*i,j) = C(3*i,j) + (x-x0);
             C(3*i-2,j) = C(3*i-2,j) - (z-z0);
         end
         j=j+1;
     end
      if bc(k,1) == 0 && bc(k,2) == 6 % Regularisation/thetaz
         % First, find the ~=barycenter of the solid
         x0 = sum( nodes(:,1) ) / size(nodes,1);
         y0 = sum( nodes(:,2) ) / size(nodes,2);
         for i=1:size(nodes,1)
             x = nodes(i,1);
             y = nodes(i,2);

             C(3*i-2,j) = C(3*i,j) + (y-y0);
             C(3*i-1,j) = C(3*i-1,j) - (x-x0);
         end
         j=j+1;
     end
 end

end
