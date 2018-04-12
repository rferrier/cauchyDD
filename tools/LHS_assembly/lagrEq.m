function [K,nb] = lagrEq( K1, K2, nodes1, nodes2, epsilon, c2node1, c2node2 )
 % This function merges 2 meshes with Lagrange multiplicators
 [i1,j1] = size(K1);
 [i2,j2] = size(K2);
 
 K = [K1,zeros(i1,j2) ; zeros(i2,j1),K2];
 C = [];
 
 % Get the nodes in front of other nodes
 for i=1:size(nodes1,1)
     for j=1:size(nodes2,1)
         if nodes1(i,1) <= nodes2(j,1)+epsilon &&...
                 nodes1(i,1) >= nodes2(j,1)-epsilon &&...
                 nodes1(i,2) <= nodes2(j,2)+epsilon &&...
                 nodes1(i,2) >= nodes2(j,2)-epsilon
             
             % Test if dof isn't already imposed at 0
             if size(find(c2node1==2*i-1),2)<1 &&...
                     size(find(c2node2==2*j-1),2)<1
                Ci = zeros(1,i1+i2);
                Ci(1,2*i-1) = 1;
                Ci(1,2*j-1+i1) = -1;
                C = [C;Ci];
             else
%                  [2*i-1,2*j-1]
%                  [size(find(c2node1==2*i-1)),size(find(c2node2==2*j-1))]
             end
             if size(find(c2node1==2*i),2)<1 &&...
                     size(find(c2node2==2*j),2)<1
                Ci = zeros(1,i1+i2);
                Ci(1,2*i) = 1;
                Ci(1,2*j+i1) = -1;
                C = [C;Ci];
             else
%                  [2*i,2*j]
%                  [size(find(c2node1==2*i)),size(find(c2node2==2*j))]
                 
                 
             end
         end
     end
 end
 
 K = [K,C';C,zeros(size(C,1))];
 nb = size(C,1);
end

