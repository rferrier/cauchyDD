function K = penalise(  K1, ent, boundary, nnodes, param )
 % This function penalises a matrix, that ont dof on the entity.
 
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

 penamatr = param * eye(size(map,2));
 K = K1;
 K(map,map) = K(map,map) + penamatr;

end

