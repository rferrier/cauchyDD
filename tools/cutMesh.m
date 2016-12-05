function [nodes1, elements1, boundary1, map1...
          nodes2, elements2, boundary2, map2, newbound] =...
                                     cutMesh (nodes, elements, boundary, line)
% This function creates 2 sub-meshes from a mesh and a line
% nodes    : nodes of the mesh
% elements : elements of the mesh
% boundary : boundary elements of the mesh
% line     : separation line defined by a segment [ xa ; ya ; xb ; yb ]

 nelem  = size(elements,1);
 nnodes = size(nodes,1);
 %gotow = zeros( nelem , 1 ); % This list stores the subdomain
 nodes1 = nodes; nodes2 = nodes; boundary1 = boundary; boundary2 = boundary;

 % A loop over the elements to determine to witch subdomain each element 
 % belongs.
 xa = line(1); ya = line(2); xb = line(3); yb = line(4);
 xba = xb-xa; yba = yb-ya;   % small time gain
 j = 1; k = 1; % indices of the elements in 1 and 2
 for i=1:nelem
    n1 = elements(i,1); n2 = elements(i,2); n3 = elements(i,3);
    
    x1 = nodes(n1,1); y1 = nodes(n1,2);
    x2 = nodes(n2,1); y2 = nodes(n2,2);
    x3 = nodes(n3,1); y3 = nodes(n3,2);
    
    % ~ ugly estimate of the aera : if>0, goto 1 else, goto 2.
    aera = ( x1+x2+x3 - 3*xa )*yba - ( y1+y2+y3 - 3*ya )*xba;
    
    if aera >= 0
       elements1(j,:) = elements(i,:);
       j = j+1;
    else
       elements2(k,:) = elements(i,:);
       k = k+1;
    end
 end

 % Now, kill floating nodes
 j = 1; k = 1; % indices of the suppressed nodes
 map1 = zeros(nnodes,1);  map2 = zeros(nnodes,1);
 for i=1:nnodes
    map1(i) = i-j+1; % associates new indices to old ones
    map2(i) = i-k+1;
    if isempty( find (elements1==i) )
       killed1(j) = i;
       j = j+1;
    end
    if isempty( find (elements2==i) )
       killed2(k) = i;
       k = k+1;
    end
 end
 
 % Ensure the last index
 map1(end) = nnodes-j+1;
 map2(end) = nnodes-k+1;

 nodes1( killed1, : ) = [];
 nodes2( killed2, : ) = [];
 
 newbound = setdiff( (1:nnodes)', union(killed1,killed2) ) ; % The boundary
 
 % Re-index the nodes in the elements
 elements1 = map1(elements1);
 elements2 = map2(elements2);
 
 % Suppress floating boundary elements
 for i=1:size(boundary,1)
    if ~isempty( find( killed1 == boundary(i,2) ) ) ||...
                        ~isempty( find( killed1 == boundary(i,3) ) )
       tosuppr1(end+1) = i;
    end
    if ~isempty( find( killed2 == boundary(i,2) ) ) ||...
                        ~isempty( find( killed2 == boundary(i,3) ) )
       tosuppr2(end+1) = i;
    end
 end

 boundary1(tosuppr1,:) = [];  boundary2(tosuppr2,:) = [];
 
 boundary1(:,[2,3]) = map1( boundary1(:,[2,3]) );
 boundary2(:,[2,3]) = map2( boundary2(:,[2,3]) );
 
end
