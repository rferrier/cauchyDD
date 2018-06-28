function u2 = passMesh3DF (nodes1, elements1, nodes2, elements2, u1, varargin)
% This function moves the field u1 from mesh 1 to mesh 2 by FE interpolation
% In this function, u1 is a vector of nodal forces
% nodes1    : nodes of the old mesh
% elements1 : elements of the olf mesh
% nodes2    : nodes of the new mesh
% elements2 : elements of the new mesh
% u1        : field on the old mesh
% varargin  : option to pass only values on the boundaries

 s2 = size(u1,2); % Keep in mind u1 can be a concatenation of many fields
 nnodes1 = size(nodes1,1);
 nnodes2 = size(nodes2,1);
 nelem1 = size(elements1,1);
 nelem2 = size(elements2,1);
 u2 = zeros( 3*nnodes2, s2 );
 
 if numel(varargin) > 0
    warning('Extra argument is not recognized');
 end
 
 %% First task : pass the forces from the nodes1 to the elements1
 ntoelem = zeros(nnodes1,1);
% for i=1:nelem1
%    no1 = elements1(i,1); no2 = elements1(i,2); no3 = elements1(i,3); no4 = elements1(i,4);
%    x1 = nodes1(no1,1); y1 = nodes1(no1,2); z1 = nodes1(no1,3);
%    x2 = nodes1(no2,1); y2 = nodes1(no2,2); z2 = nodes1(no2,3);
%    x3 = nodes1(no3,1); y3 = nodes1(no3,2); z3 = nodes1(no3,3);
%    x4 = nodes1(no4,1); y4 = nodes1(no4,2); z4 = nodes1(no4,3);
%       
%    a = [ x2-x1 ; y2-y1 ; z2-z1 ];
%    b = [ x3-x1 ; y3-y1 ; z3-z1 ];
%    c = [ x4-x1 ; y4-y1 ; z4-z1 ];
%    V = (1/6)*det( [a,b,c] );

%    ntoelem(no1) = ntoelem(no1) + V;
%    ntoelem(no2) = ntoelem(no2) + V;
%    ntoelem(no3) = ntoelem(no3) + V;
%    ntoelem(no4) = ntoelem(no4) + V;
% end

 no1 = elements1(:,1); no2 = elements1(:,2); no3 = elements1(:,3); no4 = elements1(:,4);
 x1 = nodes1(no1,1); y1 = nodes1(no1,2); z1 = nodes1(no1,3);
 x2 = nodes1(no2,1); y2 = nodes1(no2,2); z2 = nodes1(no2,3);
 x3 = nodes1(no3,1); y3 = nodes1(no3,2); z3 = nodes1(no3,3);
 x4 = nodes1(no4,1); y4 = nodes1(no4,2); z4 = nodes1(no4,3);
       
 a = [ x2-x1 , y2-y1 , z2-z1 ];
 b = [ x3-x1 , y3-y1 , z3-z1 ];
 c = [ x4-x1 , y4-y1 , z4-z1 ];
 V = (1/6) * ( a(:,1).*(b(:,2).*c(:,3)-b(:,3).*c(:,2)) -...
               a(:,2).*(b(:,1).*c(:,3)-c(:,1).*b(:,3)) +...
               a(:,3).*(b(:,1).*c(:,2)-c(:,1).*b(:,2)) );
       
 ntoelem(no1) = ntoelem(no1) + V;
 ntoelem(no2) = ntoelem(no2) + V;
 ntoelem(no3) = ntoelem(no3) + V;
 ntoelem(no4) = ntoelem(no4) + V;

 ntoelemm1 = zeros(nnodes1,1);
 ntoelemm1(find(ntoelem)) = 1./ntoelem(find(ntoelem)); % Avoid the 0

 % Duplicate ntoelemm1
 ntoelemm = zeros(3*nnodes1,1);
 ntoelemm(1:3:3*nnodes1-2) = ntoelemm1;
 ntoelemm(2:3:3*nnodes1-1) = ntoelemm1;
 ntoelemm(3:3:3*nnodes1)   = ntoelemm1;

 fR = u1.*(ntoelemm*ones(1,s2));
 fe1 = zeros(3*size(elements1,1),s2);

% for i=1:nelem1      
%    no1 = elements1(i,1); no2 = elements1(i,2); no3 = elements1(i,3); no4 = elements1(i,4);
%    fe(3*i-2,:) = (fR(3*no1-2,:) + fR(3*no2-2,:) + fR(3*no3-2,:) + fR(3*no4-2,:) );
%    fe(3*i-1,:) = (fR(3*no1-1,:) + fR(3*no2-1,:) + fR(3*no3-1,:) + fR(3*no4-1,:));
%    fe(3*i,:)   = (fR(3*no1,:)   + fR(3*no2,:)   + fR(3*no3,:)   + fR(3*no4,:));
% end
 no1 = elements1(:,1); no2 = elements1(:,2); no3 = elements1(:,3); no4 = elements1(:,4);
 fe1(3*(1:nelem1)-2,:) = (fR(3*no1-2,:) + fR(3*no2-2,:) + fR(3*no3-2,:) + fR(3*no4-2,:) );
 fe1(3*(1:nelem1)-1,:) = (fR(3*no1-1,:) + fR(3*no2-1,:) + fR(3*no3-1,:) + fR(3*no4-1,:));
 fe1(3*(1:nelem1),:)   = (fR(3*no1,:)   + fR(3*no2,:)   + fR(3*no3,:)   + fR(3*no4,:));

 %% Second task : pass the field element-to-element
 % For this, we do a 1 Gauss-point integration : we identify for each barycenter 
 % of the mesh 1 the element in witch it is and give all the force to it.

 fe2 = zeros(3*size(elements2,1),s2);

 for i = 1:nelem2 % Find the tethraedrons that are inside this element
    x1 = nodes1(elements1(:,1),1); y1 = nodes1(elements1(:,1),2); z1 = nodes1(elements1(:,1),3);
    x2 = nodes1(elements1(:,2),1); y2 = nodes1(elements1(:,2),2); z2 = nodes1(elements1(:,2),3);
    x3 = nodes1(elements1(:,3),1); y3 = nodes1(elements1(:,3),2); z3 = nodes1(elements1(:,3),3);
    x4 = nodes1(elements1(:,4),1); y4 = nodes1(elements1(:,4),2); z4 = nodes1(elements1(:,4),3);

    x = (x1+x2+x3+x4)/4; y = (y1+y2+y3+y4)/4; z = (z1+z2+z3+z4)/4; % Barycenter
    
    x1 = nodes2(elements2(i,1),1); y1 = nodes2(elements2(i,1),2); z1 = nodes2(elements2(i,1),3);
    x2 = nodes2(elements2(i,2),1); y2 = nodes2(elements2(i,2),2); z2 = nodes2(elements2(i,2),3);
    x3 = nodes2(elements2(i,3),1); y3 = nodes2(elements2(i,3),2); z3 = nodes2(elements2(i,3),3);
    x4 = nodes2(elements2(i,4),1); y4 = nodes2(elements2(i,4),2); z4 = nodes2(elements2(i,4),3);
       
    %xc = (x1+x2+x3+x4)/4; yc = (y1+y2+y3+y4)/4; zc = (z1+z2+z3+z4)/4;

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
    if size(jj) >= 1
       for j=1:max(size(jj))
          old = [3*jj(j)-2,3*jj(j)-1,3*jj(j)];
          new = [3*i-2,3*i-1,3*i];
          fe2(new,:) = fe2(new,:) + fe1(old,:)*V(jj(j)); % Add it
       end
    end
 end

 %% Third and last task : pass the forces from the elements 2 to the nodes 2
 for i=1:nelem2
   no1 = elements2(i,1); no2 = elements2(i,2);
   no3 = elements2(i,3); no4 = elements2(i,4);

%   x1 = nodes2(no1,1); x2 = nodes2(no2,1); x3 = nodes2(no3,1); x4 = nodes2(no4,1);
%   y1 = nodes2(no1,2); y2 = nodes2(no2,2); y3 = nodes2(no3,2); y4 = nodes2(no4,2);
%   z1 = nodes2(no1,3); z2 = nodes2(no2,3); z3 = nodes2(no3,3); z4 = nodes2(no4,3);

%   a = [ x2-x1 ; y2-y1 ; z2-z1 ];
%   b = [ x3-x1 ; y3-y1 ; z3-z1 ];
%   c = [ x4-x1 ; y4-y1 ; z4-z1 ];
%   V = (1/6) * ( a(1)*(b(2)*c(3)-b(3)*c(2)) - a(2)*(b(1)*c(3)-c(1)*b(3)) + a(3)*(b(1)*c(2)-c(1)*b(2)) );

   ind1 = [ 3*no1-2, 3*no1-1, 3*no1 ];
   ind2 = [ 3*no2-2, 3*no2-1, 3*no2 ];
   ind3 = [ 3*no3-2, 3*no3-1, 3*no3 ];
   ind4 = [ 3*no4-2, 3*no4-1, 3*no4 ];
   inde = [ 3*i-2, 3*i-1, 3*i ];

   u2(ind1,:) = u2(ind1,:) + .25*fe2(inde,:);
   u2(ind2,:) = u2(ind2,:) + .25*fe2(inde,:);
   u2(ind3,:) = u2(ind3,:) + .25*fe2(inde,:);
   u2(ind4,:) = u2(ind4,:) + .25*fe2(inde,:);
 end 

% no1 = elements2(:,1); no2 = elements2(:,2);
% no3 = elements2(:,3); no4 = elements2(:,4);

% x1 = nodes2(no1,1); x2 = nodes2(no2,1); x3 = nodes2(no3,1); x4 = nodes2(no4,1);
% y1 = nodes2(no1,2); y2 = nodes2(no2,2); y3 = nodes2(no3,2); y4 = nodes2(no4,2);
% z1 = nodes2(no1,3); z2 = nodes2(no2,3); z3 = nodes2(no3,3); z4 = nodes2(no4,3);

% a = [ x2-x1 , y2-y1 , z2-z1 ];
% b = [ x3-x1 , y3-y1 , z3-z1 ];
% c = [ x4-x1 , y4-y1 , z4-z1 ];
% V = (1/6) * ( a(:,1).*(b(:,2).*c(:,3)-b(:,3).*c(:,2)) -...
%               a(:,2).*(b(:,1).*c(:,3)-c(:,1).*b(:,3)) +...
%               a(:,3).*(b(:,1).*c(:,2)-c(:,1).*b(:,2)) );

% % Dispatch the V
% Vv = zeros(3*nelem2,1);
% Vv(1:3:3*nelem2-2) = V;
% Vv(2:3:3*nelem2-1) = V;
% Vv(3:3:3*nelem2)   = V;

% ind1 = [ 3*no1-2, 3*no1-1, 3*no1 ];
% ind2 = [ 3*no2-2, 3*no2-1, 3*no2 ];
% ind3 = [ 3*no3-2, 3*no3-1, 3*no3 ];
% ind4 = [ 3*no4-2, 3*no4-1, 3*no4 ];
% inde = [ 3*(1:nelem2)-2, 3*(1:nelem2)-1, 3*(1:nelem2) ];

% u2(ind1,:) = u2(ind1,:) + .25*Vv.*fe2(inde,:);
% u2(ind2,:) = u2(ind2,:) + .25*Vv.*fe2(inde,:);
% u2(ind3,:) = u2(ind3,:) + .25*Vv.*fe2(inde,:);
% u2(ind4,:) = u2(ind4,:) + .25*Vv.*fe2(inde,:);

end
