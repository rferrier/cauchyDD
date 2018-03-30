function bound2 = suppressBound( bound1, exNo, entity, varargin )
% This function removes all the boundary elements that belongs to an entity
% and that contains a node in exNo
% if exNo is the string 'extreme', the nodes with the extreme abscissae are
% eliminated, and varargin is the tolerance

 bound2 = bound1;

 nodes = [];
 if numel(varargin) > 0
    nodes = varargin{1};
 end
 epsilon = 1e-9;
 if numel(varargin) > 1
    epsilon = varargin{2};
 end
 
 if strcmp( exNo, 'extreme' )
    if size(nodes,1) == 0
       error('Please specify nodes for this option');
    end
    
    % Find the extrema
    dim = size(nodes,2)-1; % in n-1 D /!\ it only works for horizontal plane
    extrem = zeros( dim, 2 );
    for i=1:dim
       extrem(i,1) = max(nodes(:,i))-epsilon;
       extrem(i,2) = min(nodes(:,i))+epsilon;
    end

    j = 1;
    while j <= size(bound2,1) % bound2 changes size
       if bound2(j,1) == entity
          killhim = 0;
          for i=1:dim  % loop over the dimensions and nodes of the element to see if it must be removed
             for k=1:size(bound2,2)-1
                no = bound2(j,k+1);  % Remember first index in bound is not a node
                if nodes(no,i) >= extrem(i,1) || nodes(no,i) <= extrem(i,2)
                   killhim = 1;
                end % TODO : place break to speed up the stuff
             end
          end
          
          if killhim == 1
            bound2(j,:) = [];
             j = j-1;
          end
       end
       j = j+1;
    end
    
 else
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
 
end