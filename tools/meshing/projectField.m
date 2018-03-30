function u2 = projectField( u1, nodes1, nodes2, epsilon )
 % This function takes the values of u1 on nodes1 and puts it on nodes2
 % for every coincident node
 
 u2 = zeros( 2*size(nodes1,1), 1 );
 
 for i=1:size(nodes1,1)
     for j=1:size(nodes2,1)
         if nodes1(i,1) <= nodes2(j,1)+epsilon &&...
                 nodes1(i,1) >= nodes2(j,1)-epsilon &&...
                 nodes1(i,2) <= nodes2(j,2)+epsilon &&...
                 nodes1(i,2) >= nodes2(j,2)-epsilon
             u2(2*j-1) = u1(2*i-1);
             u2(2*j) = u1(2*i);
         end
     end
 end
end

