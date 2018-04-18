function u2 = passMesh2D3d (nodes1, elements1, nodes2, elements2, u1, varargin)
% This function moves the field u1 from mesh 1 to mesh 2 by FE interpolation
% Fields are on 2D meshes, but the space is 3D
% /!\ Only works for linear fields
% nodes1    : nodes of the old mesh
% elements1 : elements of the old mesh
% nodes2    : nodes of the new mesh
% elements2 : elements of the new mesh
% u1        : field on the old mesh
% varargin  : option to pass only values on the boundaries
 
 s2 = size(u1,2); % Keep in mind u1 can be a concatenation of many fields
 nnodes2 = size(nodes2,1);
 nelem1 = size(elements1,1);
 u2 = zeros( 3*nnodes2, s2 );
 
 if numel(varargin) > 0
    onlybound = cell2mat(varargin(1));
 else
    onlybound = 0;
 end
 
 if onlybound == 1
    boundary = cell2mat(varargin(1));
    [ ~ , list ] = mapBound( 0, boundary, nnodes2 ); list = list';
 else
    list = 1:nnodes2;
 end

 % Parfor is only for ornament on Octave
 for i = list % Find the triangle for this node
    x = nodes2(i,1); y = nodes2(i,2); z = nodes2(i,3);

    x1 = nodes1(elements1(:,1),1); y1 = nodes1(elements1(:,1),2); z1 = nodes1(elements1(:,1),3);
    x2 = nodes1(elements1(:,2),1); y2 = nodes1(elements1(:,2),2); z2 = nodes1(elements1(:,2),3);
    x3 = nodes1(elements1(:,3),1); y3 = nodes1(elements1(:,3),2); z3 = nodes1(elements1(:,3),3);

    % Compute the full vectroial product
    D1x = (y2-y).*(z3-z) - (y3-y).*(z2-z);
    D1y = (z2-z).*(x3-x) - (z3-z).*(x2-x);
    D1z = (x2-x).*(y3-y) - (x3-x).*(y2-y);
    
    D2x = (y3-y).*(z1-z) - (y1-y).*(z3-z);
    D2y = (z3-z).*(x1-x) - (z1-z).*(x3-x);
    D2z = (x3-x).*(y1-y) - (x1-x).*(y3-y);
    
    D3x = (y1-y).*(z2-z) - (y2-y).*(z1-z);
    D3y = (z1-z).*(x2-x) - (z2-z).*(x1-x);
    D3z = (x1-x).*(y2-y) - (x2-x).*(y1-y);
    
    % D0 is no longer the sum (x may be out of the plane)
    D0x = (y2-y1).*(z3-z1) - (y3-y1).*(z2-z1);
    D0y = (z2-z1).*(x3-x1) - (z3-z1).*(x2-x1);
    D0z = (x2-x1).*(y3-y1) - (x3-x1).*(y2-y1);
    
    % Compute the norms (the actual aera (X2))
    D1 = D1x.*D1x + D1y.*D1y + D1z.*D1z;
    D2 = D2x.*D2x + D2y.*D2y + D2z.*D2z;
    D3 = D3x.*D3x + D3y.*D3y + D3z.*D3z;
    D0 = D0x.*D0x + D0y.*D0y + D0z.*D0z;

    j1 = find(D0x.*D1x + D0y.*D1y + D0z.*D1z >= 0);
    j2 = find(D0x.*D2x + D0y.*D2y + D0z.*D2z >= 0);
    j3 = find(D0x.*D3x + D0y.*D3y + D0z.*D3z >= 0);
    j  = intersect( intersect(j1,j2) , j3);
    
%    j, D1(j), D2(j), D3(j), D0(j), bug;

    % A node can be at a frontier, in that case, choose a random element (the first one)
    if size(j) >= 1
       j = j(1);
       
       elt = elements1(j,:);
       a = D1(j)/D0(j); b = D2(j)/D0(j); d = D3(j)/D0(j);

       old = [3*elt(1)-2,3*elt(1)-1,3*elt(1)
              3*elt(2)-2,3*elt(2)-1,3*elt(2)
              3*elt(3)-2,3*elt(3)-1,3*elt(3)];
       new = [3*i-2,3*i-1,3*i];

       u2(new,:) = a*u1(old(1,:),:) + b*u1(old(2,:),:) + d*u1(old(3,:),:);
    end
 end
 
end
