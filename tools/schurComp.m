function [ S, b, map ] = schurComp( K, f, entity, boundary, nnodes )
 % This function computes the Schur complement and the Schur Rhs :
 
 % K = [Kii,Kib;Kbi,Kbb]
 % f = [fi;fb]
 map = [];
 j = 1;
 inodes = zeros(nnodes, 1);
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

 Kbb = K(map,map);
 fb  = f(map,1);
 
 Kii = K;
 fi  = f;
 Kii(map,:) = [];
 Kii(:,map) = [];
 fi(map,:)  = [];
 
 Kbi = K(map,:);
 Kbi(:,map) = [];
 
 Kib = K(:,map);
 Kib(map,:) = [];
 
 S = full(Kbb - Kbi*inv(Kii)*Kib);
 b = fb - Kbi*inv(Kii)*fi;
end