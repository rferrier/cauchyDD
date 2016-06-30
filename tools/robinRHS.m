function [K, f] = robinRHS( ntot, nodes, boundary, u, k, ent )
 % This function builds the right hand side for Robin problem and adds a
 % matrix K to add to the stiffness.

 f = zeros(2*size(nodes,1)+ntot,1);
 K = sparse(2*size(nodes,1)+ntot, 2*size(nodes,1)+ntot);
 
 % Populate f and K
 for i=1:size(boundary,1)
     entity = boundary(i,1);
     if ent == entity
        node1 = boundary(i,2);
        node2 = boundary(i,3);
        map = [2*node1-1,2*node1,2*node2-1,2*node2];

        f(map,1) = k * u(map,1);
        K(map,map) = k*eye(4);
     end
 end
 
end

