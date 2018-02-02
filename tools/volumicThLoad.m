function f = volumicThLoad( ntot, nodes, elems, fscalar )
 % This function computes the volumic RHS associated to a mesh and a value
 % of volumic loading.
 
 % imput : ntot : number of constrained dofs (to match the dimension)
 %         nodes : list of nodes of the mesh
 %         elems : list of elements of the mesh
 %         axis  : if axis = 1 : on x, if axis = 2, on y
 %         fscalar : scalar

 f = zeros(size(nodes,1) + ntot, 1);
 Ng = size(fscalar,1);
 
 if Ng == 1 % In order to speed up the most simple case
     for i=1:size(elems,1)
         map = [elems(i,1), elems(i,2), elems(i,3)];


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
    warning('Gauss number unimplemented');
 end
end
