function u2 = passMesh2D (nodes1, elements1, nodes2, elements2, u1, varargin)
% This function moves the field u1 from mesh 1 to mesh 2 by FE interpolation
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
 u2 = zeros( 2*nnodes2, s2 );
 
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
    x = nodes2(i,1); y = nodes2(i,2);

    x1 = nodes1(elements1(:,1),1); y1 = nodes1(elements1(:,1),2);
    x2 = nodes1(elements1(:,2),1); y2 = nodes1(elements1(:,2),2);
    x3 = nodes1(elements1(:,3),1); y3 = nodes1(elements1(:,3),2);
    
    D1 = (x2-x).*(y3-y) - (x3-x).*(y2-y);
    D2 = (x3-x).*(y1-y) - (x1-x).*(y3-y);
    D3 = (x1-x).*(y2-y) - (x2-x).*(y1-y);
    D0 = D1+D2+D3;

    j1 = find(D0.*D1 >= 0);
    j2 = find(D0.*D2 >= 0);
    j3 = find(D0.*D3 >= 0);
    j = intersect( intersect(j1,j2) , j3);

    % A node can be at a frontier, in that case, choose a random element (the first one)
    if size(j) >= 1
       j = j(1);
       
       elt = elements1(j,:);
       a = D1(j)/D0(j); b = D2(j)/D0(j); d = D3(j)/D0(j);

       old = [2*elt(1)-1,2*elt(1)
              2*elt(2)-1,2*elt(2)
              2*elt(3)-1,2*elt(3)];
       new = [2*i-1,2*i];

       u2(new,:) = a*u1(old(1,:),:) + b*u1(old(2,:),:) + d*u1(old(3,:),:);
    end
 end
 
end
