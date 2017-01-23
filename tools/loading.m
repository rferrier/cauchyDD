function f = loading( ntot, nodes, boundary, neumann )
 % This function returns the loading on the structure

 f = zeros(2*size(nodes,1)+ntot,1);
 
 % Determine the order
 if size(boundary,2) == 3
    order = 1;
 elseif size(boundary,2) == 4
    order = 2;
    %s3 = sqrt(3); % 'cause time is money (or not)
 else
    warning('boundary elements have not the expected number of nodes.');
 end
 
 % Populate f : add a vertical loading for each boundary element
 for i=1:size(boundary,1)
     entity = boundary(i,1);
     for j=1:size(neumann,1) % This elt needs a Neumann BC
         if neumann(j,1)==entity
            node1 = boundary(i,2);
            node2 = boundary(i,3);
            onetoto = nodes(node2,:)' - nodes(node1,:)'; %Vector from n1 to n2
            leng = norm(onetoto);
            if order==2
               node3 = boundary(i,4);
               map = [2*node1,2*node2,2*node3];
            else % Standard (good for order 1)
               map = [2*node1,2*node2];
            end
            
            if neumann(j,2) == 1          % X
                map = map-1;
            elseif neumann(j,2) == 2      % Y
                %map = [2*node1,2*node2];
            else
                error('Unable to read axis of Neumann BC')
            end
            
            if size(neumann,2) == 3
               if order == 2
                  f(map,1) = f(map,1) + 1/6*leng*neumann(j,3) * [1;1;4];% Assume segments are regular...
                             %([1+s3;1-s3;4] + [1-s3;1+s3;4]);
               else
                  f(map,1) = f(map,1) + leng * neumann(j,3) * [1/2;1/2];
               end
            else  % Size = 5 for now (1 GaussPoint) [~,~,a,b,c], f = a+bx+cy
                  % In ths case, order 2 mesh is not supported
               if order >= 2
                  warning('order 2 mesh unsupported for this loading : precision will be inferior.')
               end
               onepto = nodes(node2,:)' + nodes(node1,:)';
               f(map,1) = f(map,1) + leng * (neumann(j,3) +...
                          neumann(j,4)*onepto(1)/2 +...
                          neumann(j,5)*onepto(2)/2) * [1/2;1/2];
            end
         end
     end
 end
 
end