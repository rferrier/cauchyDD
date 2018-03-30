function M = mass_mat(nodes, elements)
 % This function computes the mass matrix in 2D
 
 nnodes = size(nodes,1);
 M = sparse(2*nnodes, 2*nnodes);
 
 for i=1:size(elements)
     no1 = elements(i,1); no2 = elements(i,2); no3 = elements(i,3);
     x1 = nodes(no1,1); x2 = nodes(no2,1); x3 = nodes(no3,1);
     y1 = nodes(no1,2); y2 = nodes(no2,2); y3 = nodes(no3,2);
     S = abs( .5*((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)) );% element area
     map = [2*no1-1, 2*no1, 2*no2-1, 2*no2, 2*no3-1, 2*no3];
     
     Be = 1/3 * [1,0,1,0,1,0;0,1,0,1,0,1];
     M(map,map) = M(map,map) + S*(Be)'*Be;
 end
end

