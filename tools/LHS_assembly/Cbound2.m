function [ C, ntot ] = Cbound2 ( nodes, bc, boundary )
% This function computes the trace matrix C such that C'*u = u|bound
% Does not work for rigid body motions

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

 cdofs = zeros(ntot,1); %cindex = zeros(ntot,1); cvalues = zeros(ntot,1);
 cindex = 1:ntot; cvalues = 1;
 j=1;
 for i=1:size(nodes,1)
     if inodes(i,1) ~= 0
         cdofs(j) = 3*i-2;
%         cindex(j) = j;
%         cvalues(j) = inodes(i,1);
         j=j+1;
     end
     if inodes(i,2) ~= 0
         cdofs(j) = 3*i-1;
%         cindex(j) = j;
%         cvalues(j) = inodes(i,2);
         j=j+1;
     end
     if inodes(i,3) ~= 0
         cdofs(j) = 3*i;
%         cindex(j) = j;
%         cvalues(j) = inodes(i,3);
         j=j+1;
     end
 end

 C = sparse( cdofs, cindex, cvalues, 3*size(nodes,1), ntot  );
 
end
