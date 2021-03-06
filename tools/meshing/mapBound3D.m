function [ node2b, b2node ] = mapBound3D( entity, boundary, nnodes, varargin )
 % This function builds the equivalence table between a boundary and the
 % global nodal notation
 % if entity = 0 : maps all the boundaries

 % I GUESS THIS ONE COULD BE MERGED WITH mapBound (2D)
 
 % Excluded nodes
 exNo = [];
 if numel( varargin ) > 0
     exNo = cell2mat( varargin(1) );
 end

 node2b = zeros(nnodes,1);
 b2node = zeros(nnodes,1);
 j = 1;
 inodes = zeros(nnodes,1);
 for i=1:size(boundary,1)
     if boundary(i,1) == entity || entity == 0
         ok = ones(size(boundary,2)); % Is node i2 in exNo ?
         for i2 = 2:size(boundary,2)
            for k=1:size(exNo)
               if exNo(k) == boundary(i,i2)
                  ok1(i2-1) = 0;
                  break
               end
            end
            
            if inodes( boundary(i,i2) ) == 0 && ok(i2-1) == 1
               node2b(boundary(i,i2),1) = j;
               b2node(j) = boundary(i,i2);
               inodes( boundary(i,i2) ) = 1;  % Don't re-return this node again
               j = j+1;
            end
         end
%         ok1 = 1;  % Is node 1 in exNo ?
%         ok2 = 1;
%         ok3 = 1;
%         for k=1:size(exNo)
%            if exNo(k) == boundary(i,2)
%                ok1 = 0;
%                break
%            end
%         end
%         for k=1:size(exNo)
%            if exNo(k) == boundary(i,3)
%                ok2 = 0;
%                break
%            end
%         end
%         for k=1:size(exNo)
%            if exNo(k) == boundary(i,4)
%                ok3 = 0;
%                break
%            end
%         end
%         
%         if inodes( boundary(i,2) ) == 0 && ok1 == 1
%            node2b(boundary(i,2),1) = j;
%            b2node(j) = boundary(i,2);
%            inodes( boundary(i,2) ) = 1;
%            j = j+1;
%         end
%         if inodes( boundary(i,3) ) == 0 && ok2 == 1
%            node2b(boundary(i,3),1) = j;
%            b2node(j) = boundary(i,3);
%            inodes( boundary(i,3) ) = 1;
%            j = j+1;
%         end
%         if inodes( boundary(i,4) ) == 0 && ok3 == 1
%            node2b(boundary(i,4),1) = j;
%            b2node(j) = boundary(i,4);
%            inodes( boundary(i,4) ) = 1;
%            j = j+1;
%         end
     end
 end
 
 % reshape b2node
 b2node(j:nnodes) = [];

end
