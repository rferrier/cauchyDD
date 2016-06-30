function Mu = norm_bound(u, nodes, boundary, entity)
 % This function computes a norm on a boundary
 % M is such as u'Mu is a mesure of the magnitude of u
 
 Mu = zeros(size(u,1),1);
 
 for i=1:size(boundary)
     if boundary(i,1) == entity
         no1 = boundary(i,2); no2 = boundary(i,3);
         x1 = nodes(no1,1); y1 = nodes(no1,2);
         x2 = nodes(no2,1); y2 = nodes(no2,2);
         len = sqrt((x1-x2)^2 + (y1-y2)^2);
         map = [2*no1-1, 2*no1, 2*no2-1, 2*no2];
         Mu(map,1) = Mu(map,1) + len/2*[ 1,0,1/3,0; 0,1,0,1/3;...
                                  1/3,0,1,0; 0,1/3,0,1] * u(map,1);
     end
 end
end

