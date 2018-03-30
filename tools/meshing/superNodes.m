function [ b1to2, b2to1 ] = superNodes( nodes1, nodes2, epsilon )
 % This function builds the tables of coincident nodes between mesh1 and
 % mesh2

 b1to2 = zeros( size(nodes1,1), 1 );
 b2to1 = zeros( size(nodes2,1), 1 );
 
 for i=1:size(nodes1,1)
     for j=1:size(nodes2,1)
         if nodes1(i,1) <= nodes2(j,1)+epsilon &&...
                 nodes1(i,1) >= nodes2(j,1)-epsilon &&...
                 nodes1(i,2) <= nodes2(j,2)+epsilon &&...
                 nodes1(i,2) >= nodes2(j,2)-epsilon
             b1to2(i,1) = j;
             b2to1(j,1) = i;
         end
     end
 end

end

