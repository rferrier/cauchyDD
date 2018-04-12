function fdir2 = assembleDirichlet( fdir1 )
% This function assemblies the Didichlet Rhs contained in fdir1 as columns.

 fdir2 = zeros(size(fdir1,1),1); % zero
 
 for i=1:size(fdir1,1)
     for j=1:size(fdir1,2)
         if abs(fdir1(i,j)) > 1e-6
             fdir2(i) = fdir1(i,j);
         end
     end
 end

end

