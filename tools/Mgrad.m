function M = Mgrad( nb, nodes, boundary, entity )
%Mgrad computes the mass matrix of the gradient operator

 M = zeros( nb );
 for i=1:size(boundary)
     if boundary(i,1) == entity
         no1 = boundary(i,2); no2 = boundary(i,3);
         x1 = nodes(no1,1); y1 = nodes(no1,2);
         x2 = nodes(no2,1); y2 = nodes(no2,2);
         map = [2*no1-1, 2*no1, 2*no2-1, 2*no2];
         len = sqrt((x1-x2)^2 + (y1-y2)^2);
         M(map,map) = M(map,map) + 1/(2*len)*[ 1,0,-1,0; 0,1,0,-1;...
                                               -1,0,1,0; 0,-1,0,1];
     end
 end

end

