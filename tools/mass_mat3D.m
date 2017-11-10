function M = mass_mat3D(nodes, elements)
 % This function computes the mass matrix in 3D
 
 nnodes = size(nodes,1);
 M = sparse(3*nnodes, 3*nnodes);
 
 for i=1:size(elements)
     no1 = elements(i,1); no2 = elements(i,2);
     no3 = elements(i,3); no4 = elements(i,4);
     x1 = nodes(no1,1); x2 = nodes(no2,1); x3 = nodes(no3,1); x4 = nodes(no4,1);
     y1 = nodes(no1,2); y2 = nodes(no2,2); y3 = nodes(no3,2); y4 = nodes(no4,2);
     z1 = nodes(no1,3); z2 = nodes(no2,3); z3 = nodes(no3,3); z4 = nodes(no4,3);
     
    V = ( (y2-y1).*(z3-z1)-(y3-y1).*(z2-z1) ).*(x4-x1) + ...
        (-(x2-x1).*(z3-z1)+(x3-x1).*(z2-z1) ).*(y4-y1) + ...
        ( (x2-x1).*(y3-y1)-(x3-x1).*(y2-y1) ).*(z4-z1);
          
     map = [ 3*no1-2, 3*no1-1, 3*no1, 3*no2-2, 3*no2-1, 3*no2, ...
             3*no3-2, 3*no3-1, 3*no3, 3*no4-2, 3*no4-1, 3*no4 ];
     
     Be = 1/4 * [ 1,0,0,1,0,0,1,0,0,1,0,0; ...
                  0,1,0,0,1,0,0,1,0,0,1,0; ...
                  0,0,1,0,0,1,0,0,1,0,0,1 ];
     M(map,map) = M(map,map) + V*(Be)'*Be;
 end
end