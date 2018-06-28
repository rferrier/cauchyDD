function u2 = passMesh3Db (nodes1, elements1, nodes2, elements2, u1, varargin)
% This function moves the field u1 from mesh 1 to mesh 2 by FE interpolation
% Optimized if the new mesh size is smaller
% /!\ Only works for linear fields
% nodes1    : nodes of the old mesh
% elements1 : elements of the olf mesh
% nodes2    : nodes of the new mesh
% elements2 : elements of the new mesh
% u1        : field on the old mesh
% varargin  : option to pass only values on the boundaries

% Formula from : http://steve.hollasch.net/cgindex/geometry/ptintet.html
 
 s2 = size(u1,2); % Keep in mind u1 can be a concatenation of many fields
 nnodes2 = size(nodes2,1);
 nelem1 = size(elements1,1);
 u2 = zeros( 3*nnodes2, s2 );
 
 if numel(varargin) > 0
    warning('Extra argument not recognized');
 end
 
 for i = 1:nelem1 % Find the tethraedron for this node
    x = nodes2(:,1); y = nodes2(:,2); z = nodes2(:,3);
    
    x1 = nodes1(elements1(i,1),1); y1 = nodes1(elements1(i,1),2); z1 = nodes1(elements1(i,1),3);
    x2 = nodes1(elements1(i,2),1); y2 = nodes1(elements1(i,2),2); z2 = nodes1(elements1(i,2),3);
    x3 = nodes1(elements1(i,3),1); y3 = nodes1(elements1(i,3),2); z3 = nodes1(elements1(i,3),3);
    x4 = nodes1(elements1(i,4),1); y4 = nodes1(elements1(i,4),2); z4 = nodes1(elements1(i,4),3);
       
    D1 = ( (y2-y).*(z3-z)-(y3-y).*(z2-z) ).*(x4-x) + ...
         (-(x2-x).*(z3-z)+(x3-x).*(z2-z) ).*(y4-y) + ...
         ( (x2-x).*(y3-y)-(x3-x).*(y2-y) ).*(z4-z);
    D2 = ( (y-y1).*(z3-z1)-(y3-y1).*(z-z1) ).*(x4-x1) + ...
         (-(x-x1).*(z3-z1)+(x3-x1).*(z-z1) ).*(y4-y1) + ...
         ( (x-x1).*(y3-y1)-(x3-x1).*(y-y1) ).*(z4-z1);
    D3 = ( (y2-y1).*(z-z1)-(y-y1).*(z2-z1) ).*(x4-x1) + ...
         (-(x2-x1).*(z-z1)+(x-x1).*(z2-z1) ).*(y4-y1) + ...
         ( (x2-x1).*(y-y1)-(x-x1).*(y2-y1) ).*(z4-z1);
    D4 = ( (y2-y1).*(z3-z1)-(y3-y1).*(z2-z1) ).*(x-x1) + ...
         (-(x2-x1).*(z3-z1)+(x3-x1).*(z2-z1) ).*(y-y1) + ...
         ( (x2-x1).*(y3-y1)-(x3-x1).*(y2-y1) ).*(z-z1);
    
    D0 = D1+D2+D3+D4;
       
    j1 = find(D0.*D1 >= 0); j2 = find(D0.*D2 >= 0);
    j3 = find(D0.*D3 >= 0); j4 = find(D0.*D4 >= 0);
    jj = intersect( intersect( intersect(j1,j2) , j3) , j4);

    % Loop over all the nodes that are inside the element
    if size(jj) >= 1
       for k=1:max(size(jj))
         j = jj(k);

         elt = elements1(i,:);
         a = D1(j)/D0(j); b = D2(j)/D0(j); c = D3(j)/D0(j); d = D4(j)/D0(j);

         old = [3*elt(1)-2,3*elt(1)-1,3*elt(1)
               3*elt(2)-2,3*elt(2)-1,3*elt(2)
               3*elt(3)-2,3*elt(3)-1,3*elt(3)
               3*elt(4)-2,3*elt(4)-1,3*elt(4)];
         new = [3*j-2,3*j-1,3*j];
         u2(new,:) = a*u1(old(1,:),:) + b*u1(old(2,:),:) + c*u1(old(3,:),:) + d*u1(old(4,:),:);
      end
    end
 end
 
end
