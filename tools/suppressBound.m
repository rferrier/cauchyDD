function bound2 = suppressBound( bound1, exNo, entity )
% This function removes all the boundary elements that belongs to an entity
% and that contain a node in exNo

 bound2 = bound1;

 for i=1:size(exNo,1)
     j = 1;
     while j <= size(bound2,1) % bound2 changes size
         if bound2(j,1) == entity &&...
                ( bound2(j,2) == exNo(i,1) || bound2(j,3) == exNo(i,1) )
             bound2(j,:) = [];
             j = j-1;
         end
         j = j+1;
     end
 end

end

