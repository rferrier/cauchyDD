function f = loading3D( ntot, nodes, boundary, neumann )
 % This function returns the loading on the structure

 f = zeros(3*size(nodes,1)+ntot,1);
 
 % Populate f : add a vertical loading for each boundary element
 for i=1:size(boundary,1)
     entity = boundary(i,1);
     for j=1:size(neumann,1) % This elt needs a Neumann BC
         if neumann(j,1)==entity
            if neumann(j,2) == 1
               map = [3*boundary(i,2)-2, 3*boundary(i,3)-2, 3*boundary(i,4)-2];
            elseif neumann(j,2) == 2
               map = [3*boundary(i,2)-1, 3*boundary(i,3)-1, 3*boundary(i,4)-1];
            elseif neumann(j,2) == 3
               map = [3*boundary(i,2), 3*boundary(i,3), 3*boundary(i,4)];
            else
               error('Unable to read the axis of the surfacic load')
            end
   
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
   
            S = 1/2*norm( [ a(2)*b(3)-a(3)*b(2) ;...
                        a(3)*b(1)-a(1)*b(3) ;...
                        a(1)*b(2)-a(2)*b(1) ] ); % Aera of the element
            f(map,1) = f(map,1) + neumann(j,3)*S/3*[1;1;1];
         end
     end
 end
 
end