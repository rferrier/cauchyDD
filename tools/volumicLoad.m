function f = volumicLoad( ntot, nodes, elems, axis, fscalar )
 % This function computes the volumic RHS associated to a mesh and a value
 % of volumic loading.
 
 % imput : ntot : number of constrained dofs (to match the dimension)
 %         nodes : list of nodes of the mesh
 %         elems : list of elements of the mesh
 %         axis  : if axis = 1 : on x, if axis = 2, on y
 %         fscalar : array of polynomial coeffs :
 %         [a,b;c,d] means a + bx + cy + dxy
 %         note : the integration is exact only if the terms under the
 %         bottom-left diagonal are 0 : ex :
 %         [a,b,c;d,e,0;g,0,0]; is exact
 %         [a,b;c,d]; is inexact, but [a,b,0;c,d,0;0,0,0]; is exact
 %         The maximum order is 8

 f = zeros(2*size(nodes,1) + ntot, 1);
 Ng = size(fscalar,1);
 
 if Ng == 1 % In order to speed up the most simple case
     for i=1:size(elems,1)
         if axis == 1
            map = [2*elems(i,1)-1, 2*elems(i,2)-1, 2*elems(i,3)-1];
         elseif axis == 2
            map = [2*elems(i,1), 2*elems(i,2), 2*elems(i,3)];
         else
            error('Unable to read the axis of the volumic load')
         end

         x1 = nodes(elems(i,1),1);
         y1 = nodes(elems(i,1),2);
         x2 = nodes(elems(i,2),1);
         y2 = nodes(elems(i,2),2);
         x3 = nodes(elems(i,3),1);
         y3 = nodes(elems(i,3),2);

         S = abs( .5*((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)) ); % Aera of the element
         f(map,1) = f(map,1) + fscalar*S/3*[1;1;1];
     end
 else
     for i=1:size(elems,1)
         if axis == 1
            map = [2*elems(i,1)-1, 2*elems(i,2)-1, 2*elems(i,3)-1];
         elseif axis == 2
            map = [2*elems(i,1), 2*elems(i,2), 2*elems(i,3)];
         else
            error('Unable to read the axis of the volumic load')
         end

         x1 = nodes(elems(i,1),1);
         y1 = nodes(elems(i,1),2);
         x2 = nodes(elems(i,2),1);
         y2 = nodes(elems(i,2),2);
         x3 = nodes(elems(i,3),1);
         y3 = nodes(elems(i,3),2);

         [ Xg, Wg ] = gaussPt( Ng ); % Compute the Gauss points
         ng         = size(Wg,1);    % Effective number of Gauss points
         ftot       = [0;0;0];       % Total value
         for j=1:ng
             xl = Xg(j,1);
             yl = Xg(j,2);
             sf = shapefunc( [xl;yl],0,1 );
             x = sf*[x1;x2;x3];  % Global coordinates of Gauss points
             y = sf*[y1;y2;y3];
             xy = [1,   x,     x^2,     x^3,     x^4,     x^5,     x^6,     x^7,     x^8
                   y,   x*y,   x^2*y,   x^3*y,   x^4*y,   x^5*y,   x^6*y,   x^7*y,   x^8*y
                   y^2, x*y^2, x^2*y^2, x^3*y^2, x^4*y^2, x^5*y^2, x^6*y^2, x^7*y^2, x^8*y^2
                   y^3, x*y^3, x^2*y^3, x^3*y^3, x^4*y^3, x^5*y^3, x^6*y^3, x^7*y^3, x^8*y^3
                   y^4, x*y^4, x^2*y^4, x^3*y^4, x^4*y^4, x^5*y^4, x^6*y^4, x^7*y^4, x^8*y^4
                   y^5, x*y^5, x^2*y^5, x^3*y^5, x^4*y^5, x^5*y^5, x^6*y^5, x^7*y^5, x^8*y^5
                   y^6, x*y^6, x^2*y^6, x^3*y^6, x^4*y^6, x^5*y^6, x^6*y^6, x^7*y^6, x^8*y^6
                   y^7, x*y^7, x^2*y^7, x^3*y^7, x^4*y^7, x^5*y^7, x^6*y^7, x^7*y^7, x^8*y^7
                   y^8, x*y^8, x^2*y^8, x^3*y^8, x^4*y^8, x^5*y^8, x^6*y^8, x^7*y^8, x^8*y^8];
             fgi = fscalar .* xy( 1:Ng, 1:Ng );
             ftot = ftot + Wg(j)*sum(sum(fgi)) * [sf(1) ; sf(2) ; sf(3)];
         end

         S = abs( .5*((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)) ); % Aera of the element
         f(map,1) = f(map,1) + ftot*S;
     end
 end
end

