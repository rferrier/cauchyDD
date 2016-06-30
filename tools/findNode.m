function index = findNode( x, y, nodes, epsilon )
 % This function finds the asked node
 index = 0;
 
 for i=1:size(nodes)
     if nodes(i,1) < x+epsilon && nodes(i,1) > x-epsilon &&...
             nodes(i,2) < y+epsilon && nodes(i,2) > y-epsilon
         index = i;
         break
     end
 end
 
 if index == 0
     error('umpossible to find the node')
 end

end

