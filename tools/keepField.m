function u2 = keepField( u1, entity, boundary, varargin )
 % This function copies u1, but only keeps the values on the entity
 % varargin represents the coordinate (by default, both coords are kept)

 % component to keep
 co = 0;
 if numel(varargin)>0
     co = cell2mat(varargin(1));
 end
 
 u2 = u1-u1;
 
 if co==1
     for i=1:size(boundary,1)
         if boundary(i,1) == entity
             u2(2*boundary(i,2)-1) = u1(2*boundary(i,2)-1);
             u2(2*boundary(i,3)-1) = u1(2*boundary(i,3)-1);
         end
     end
 elseif co==2
     for i=1:size(boundary,1)
         if boundary(i,1) == entity
             u2(2*boundary(i,2)) = u1(2*boundary(i,2));
             u2(2*boundary(i,3)) = u1(2*boundary(i,3));
         end
     end
 else % co==0
     for i=1:size(boundary,1)
         if boundary(i,1) == entity
             u2(2*boundary(i,2)-1) = u1(2*boundary(i,2)-1);
             u2(2*boundary(i,2)) = u1(2*boundary(i,2));
             u2(2*boundary(i,3)-1) = u1(2*boundary(i,3)-1);
             u2(2*boundary(i,3)) = u1(2*boundary(i,3));
         end
     end
 end
end

