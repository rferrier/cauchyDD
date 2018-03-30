function f = pressureLoad( ntot, nodes, boundary, pressure, ent )
 % This function returns the loading on the structure
 % pressure is an array of polynomial values
 % pressure = [a,b;c,d] means p = a+bx+cy+dxy

 f = zeros(2*size(nodes,1)+ntot,1);
 
 % User put a scalar as imput
 if size(pressure,1) == 1
     pressure = [pressure, 0 ; 0, 0];
 end
 
 % Populate f : add a vertical loading for each boundary element
 for i=1:size(boundary,1)
     entity = boundary(i,1);
     for j=1:size(ent,1) % This elt needs a Neumann BC
         if ent(j)==entity
            node1 = boundary(i,2);
            node2 = boundary(i,3);
            onetoto = nodes(node2,:)' - nodes(node1,:)'; %Vector from n1 to n2
            leng = norm(onetoto);
            theta = [onetoto(2);-onetoto(1)];  % cos(theta); sin(theta)
            Xg = 0.5*(nodes(node2,:)' + nodes(node1,:)'); % Gauss point
            pr = pressure(1,1) + pressure(1,2)*Xg(1) +...
                pressure(2,1)*Xg(2) + pressure(2,2)*Xg(1)*Xg(2); % Value at Gauss point
            mapx = [2*node1-1,2*node2-1];
            mapy = [2*node1,2*node2];
            
            f(mapx,1) = f(mapx,1) + leng * pr * theta(1) * [1/2;1/2] ;
            f(mapy,1) = f(mapy,1) + leng * pr * theta(2) * [1/2;1/2] ;
         end
     end
 end
 
end
