function ret = plotOnBound( u, ent, boundary, nnodes, axis )
 % This functions plots the field u on the entity

 map = [];
 inodes = zeros(nnodes);
 j = 1;
 
 for i=1:size(boundary,1)
     if boundary(i,1) == ent
         if inodes(boundary(i,2), 1) == 0
            inodes(boundary(i,2), 1) = 1;
            map(1,[2*j-1,2*j]) = [2*boundary(i,2)-1,2*boundary(i,2)];
            j = j+1;
         end
         if inodes(boundary(i,3), 1) == 0
            inodes(boundary(i,3), 1) = 1;
            map(1,[2*j-1,2*j]) = [2*boundary(i,3)-1,2*boundary(i,3)];
            j = j+1;
         end
     end
 end
 
 urb = reshape(u(map,1),2,[])';
 ret = plot(urb(:,axis));
 figure;
end

