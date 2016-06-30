function [ K2, f2 ] = suppressLine( K1, f1, ent, boundary, nnodes )
 % This function suppresses the lines and columns of a matrix, that
 % represent dof on the entuity.
 
 map = [];
 j = 1;
 inodes = zeros(nnodes, 1);
 for k=1:size(ent,2)
     entity = ent(1,k);
     for i=1:size(boundary,1)
         if boundary(i,1) == entity
             if inodes(boundary(i,2)) == 0 % Add this one
                map(1,[2*j-1,2*j]) = [2*boundary(i,2)-1,2*boundary(i,2)];
                inodes(boundary(i,2)) = 1;
                j = j+1;
             end
             if inodes(boundary(i,3)) == 0
                map(1,[2*j-1,2*j]) = [2*boundary(i,3)-1,2*boundary(i,3)];
                inodes(boundary(i,3)) = 1;
                j = j+1;
             end
         end
     end
 end
 
 K2 = K1;
 f2  = f1;
 K2(map,:) = [];
 K2(:,map) = [];
 f2(map,:)  = [];

end

