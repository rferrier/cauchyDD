function u2 = passMesh2D3d (nodes1, elements1, nodes2, elements2, u1, varargin)
% This function moves the field u1 from mesh 1 to mesh 2 by FE interpolation
% Fields are on 2D meshes, but the space is 3D
% /!\ Only works for linear fields, and mesh 1 should be on a plane
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

 % Find the plane equation from the first element
 x1 = nodes1(elements1(1,1),1); y1 = nodes1(elements1(1,1),2); z1 = nodes1(elements1(1,1),3);
 x2 = nodes1(elements1(1,2),1); y2 = nodes1(elements1(1,2),2); z2 = nodes1(elements1(1,2),3);
 x3 = nodes1(elements1(1,3),1); y3 = nodes1(elements1(1,3),2); z3 = nodes1(elements1(1,3),3);

 mat = [ x1,y1,z1,1 ; ...
          x2,y2,z2,1 ;
          x3,y3,z3,1 ;
          1,1,1,1 ];

 if cond(mat) > 1e4 % In case of bad conditionning (close to x-y+z-1=0) change the extra condition
    mat = [ x1,y1,z1,1 ; ...
             x2,y2,z2,1 ;
             x3,y3,z3,1 ;
             1,-1,1,-1 ];
 end

 abcd = mat \ [0;0;0;1];

 % Parfor is only for ornament on Octave
 for i = list % Find the triangle for this node
    x = nodes2(i,1); y = nodes2(i,2); z = nodes2(i,3);

    x1 = nodes1(elements1(:,1),1); y1 = nodes1(elements1(:,1),2); z1 = nodes1(elements1(:,1),3);
    x2 = nodes1(elements1(:,2),1); y2 = nodes1(elements1(:,2),2); z2 = nodes1(elements1(:,2),3);
    x3 = nodes1(elements1(:,3),1); y3 = nodes1(elements1(:,3),2); z3 = nodes1(elements1(:,3),3);

    % Projection point
    alpha = -( abcd(1)*x + abcd(2)*y + abcd(3)*z + abcd(4) ) /...
             (abcd(1)^2 + abcd(2)^2 + abcd(3)^2);
    xn = x+alpha*abcd(1);
    yn = y+alpha*abcd(2);
    zn = z+alpha*abcd(3);

    % Compute the weights
    D1x = (y2-yn).*(z3-zn) - (y3-yn).*(z2-zn);
    D1y = (z2-zn).*(x3-xn) - (z3-zn).*(x2-xn);
    D1z = (x2-xn).*(y3-yn) - (x3-xn).*(y2-yn);
    D2x = (y3-yn).*(z1-zn) - (y1-yn).*(z3-zn);
    D2y = (z3-zn).*(x1-xn) - (z1-zn).*(x3-xn);
    D2z = (x3-xn).*(y1-yn) - (x1-xn).*(y3-yn);
    D3x = (y1-yn).*(z2-zn) - (y2-yn).*(z1-zn);
    D3y = (z1-zn).*(x2-xn) - (z2-zn).*(x1-xn);
    D3z = (x1-xn).*(y2-yn) - (x2-xn).*(y1-yn);
    D0x = (y1-y3).*(z2-z3) - (y2-y3).*(z1-z3);
    D0y = (z1-z3).*(x2-x3) - (z2-z3).*(x1-x3);
    D0z = (x1-x3).*(y2-y3) - (x2-x3).*(y1-y3);

    D1 = sqrt(D1x.*D1x + D1y.*D1y + D1z.*D1z);
    D2 = sqrt(D2x.*D2x + D2y.*D2y + D2z.*D2z);
    D3 = sqrt(D3x.*D3x + D3y.*D3y + D3z.*D3z);
    D0 = D1+D2+D3;

    % And find the right one
    j1 = find(D0x.*D1x + D0y.*D1y + D0z.*D1z >= 0);
    j2 = find(D0x.*D2x + D0y.*D2y + D0z.*D2z >= 0);
    j3 = find(D0x.*D3x + D0y.*D3y + D0z.*D3z >= 0);
    j  = intersect( intersect(j1,j2) , j3);
    
    %j, D0(j), D1(j), D2(j), D3(j), bug;

    % A node can be at a frontier, in that case, choose a random element (the first one)
    if size(j) >= 1
       j = j(1);
       elt = elements1(j,:);

    else % Extrapolate
 
       xm = (x1+x2+x3)/3; ym = (y1+y2+y3)/3; zm = (z1+z2+z3)/3; % Coordiantes of the barycenters of the elements
       [~,j] = min( (x-xm).^2 + (y-ym).^2 + (z-zm).^2 ); % Find the closest one
       elt = elements1(j,:);

    end

    a = D1(j)/D0(j); b = D2(j)/D0(j); d = D3(j)/D0(j);
    old = [3*elt(1)-2,3*elt(1)-1,3*elt(1)
           3*elt(2)-2,3*elt(2)-1,3*elt(2)
           3*elt(3)-2,3*elt(3)-1,3*elt(3)];
    new = [3*i-2,3*i-1,3*i];

    u2(new,:) = a*u1(old(1,:),:) + b*u1(old(2,:),:) + d*u1(old(3,:),:);

 end
 
end
