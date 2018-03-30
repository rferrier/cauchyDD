function [M] = bMass_mat (nodes, boundary, index)
% This function computes the boundary mass matrix
% such that u'Mv = int_{Gamma}(u*v) 
% provided u and v are interpolated with shape functions

% /!\ Problem at corner nodes

 M = zeros(2*size(nodes,1));

 for i=1:size(index,1)
    bound = index(i);
    ibound = zeros(size(boundary,1), 1);
    
    % Look for an element in boundary
    for j=1:size(boundary)
       if bound == boundary(j,1) && ibound(j) == 0
          no1 = boundary(j,2);
          no2 = boundary(j,3);
          map = [2*no1-1, 2*no1, 2*no2-1, 2*no2];
          x1 = nodes(no1,1); y1 = nodes(no1,2);
          x2 = nodes(no2,1); y2 = nodes(no2,2);
          leng = sqrt( (x2-x1)^2 + (y2-y1)^2 );
          M(map,map) = leng/3 * [1,0,.5,0;0,1,0,.5;.5,0,1,0;0,.5,0,1];
       end
    end
    
 end

end