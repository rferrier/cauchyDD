function f = pressureLoad3D( ntot, nodes, boundary, pressure, ent )
 % This function returns the loading on the structure
 % pressure is an array of polynomial values
 % pressure = [a,b,c,d,e] means p = a+bx+cy+dz+ex^2 (not very general...)

 f = zeros(3*size(nodes,1)+ntot,1);
 
 % User put a scalar as imput
 if size(pressure,1) == 1
     pressure = [pressure, 0, 0, 0, 0];
 end
 
 % Populate f : add a vertical loading for each boundary element
 for i=1:size(boundary,1)
     entity = boundary(i,1);
     for j=1:size(ent,1) % This elt needs a Neumann BC
         if ent(j)==entity
            node1 = boundary(i,2); node2 = boundary(i,3); node3 = boundary(i,4);
            x1 = nodes( boundary(i,2), 1) ;
            y1 = nodes( boundary(i,2), 2 );
            z1 = nodes( boundary(i,2), 3 );
            x2 = nodes( boundary(i,3), 1 );
            y2 = nodes( boundary(i,3), 2 );
            z2 = nodes( boundary(i,3), 3 );
            x3 = nodes( boundary(i,4), 1 );
            y3 = nodes( boundary(i,4), 2 );
            z3 = nodes( boundary(i,4), 3 );
            
            % Prepare vectors for surface computation
            a = [ x2-x1 ; y2-y1 ; z2-z1 ];
            b = [ x3-x1 ; y3-y1 ; z3-z1 ];
            n = [ a(2)*b(3)-a(3)*b(2) ;...
                  a(3)*b(1)-a(1)*b(3) ;...
                  a(1)*b(2)-a(2)*b(1) ];
            S = 1/2*norm(n); % Aera of the element
            n = n/norm(n);
            
            Xg = 1/3*([x1;y1;z1] + [x2;y2;z2] + [x3;y3;z3]); % Gauss point
            pr = pressure(1) + pressure(2)*Xg(1) + pressure(3)*Xg(2) + ...
                 pressure(4)*Xg(3) + pressure(5)*Xg(1)^2; % Value at Gauss point
            mapx = [3*node1-2,3*node2-2,3*node3-2];
            mapy = [3*node1-1,3*node2-1,3*node3-1];
            mapz = [3*node1,3*node2,3*node3];
            
            f(mapx,1) = f(mapx,1) + S * pr * n(1) * [1/3;1/3;1/3] ;
            f(mapy,1) = f(mapy,1) + S * pr * n(2) * [1/3;1/3;1/3] ;
            f(mapz,1) = f(mapz,1) + S * pr * n(3) * [1/3;1/3;1/3] ;
         end
     end
 end
 
end
