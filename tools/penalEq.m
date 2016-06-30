function K = penalEq( K1, K2, param, nodes1, nodes2, epsilon, varargin )
 % This function merges 2 meshes by penalisation
 [i1,j1] = size(K1);
 [i2,j2] = size(K2);
 
 K = [K1,zeros(i1,j2) ; zeros(i2,j1),K2];
 
 % Get the nodes in front of other nodes
 for i=1:size(nodes1,1)
     for j=1:size(nodes2,1)
         if nodes1(i,1) <= nodes2(j,1)+epsilon &&...
                 nodes1(i,1) >= nodes2(j,1)-epsilon &&...
                 nodes1(i,2) <= nodes2(j,2)+epsilon &&...
                 nodes1(i,2) >= nodes2(j,2)-epsilon
             map = [2*i-1,2*i,2*j-1+i1,2*j+i1];
             K(map,map) = K(map,map) +...
                 param*[1,0,-1,0;0,1,0,-1;-1,0,1,0;0,-1,0,1];
         end
     end
 end
 
end

