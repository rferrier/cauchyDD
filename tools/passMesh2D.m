function u2 = passMesh2D (nodes1, elements1, nodes2, elements2, u1, varargin)
% This function moves the field u1 from mesh 1 to mesh 2 by FE interpolation
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
 u2 = zeros( 2*nnodes2, s2 );
 
 if numel(varargin) > 0
    boundary = cell2mat(varargin(1));
    [ ~ , list ] = mapBound3D( 0, boundary, nnodes2 ); list = list';
 else
    list = 1:nnodes2;
 end
 
 % Parfor is only for ornament on Octave
 parfor i = list % Find the tethraedron for this node
    x = nodes2(i,1); y = nodes2(i,2);
    for j = 1:nelem1
       elt = elements1(j,:);
       x1 = nodes1(elt(1),1); y1 = nodes1(elt(1),2);
       x2 = nodes1(elt(2),1); y2 = nodes1(elt(2),2);
       x3 = nodes1(elt(3),1); y3 = nodes1(elt(3),2);

       D0 = det([x1,y1,1;x2,y2,1;x3,y3,1]);
       D1 = det([x,y,1;x2,y2,1;x3,y3,1]);
       D2 = det([x1,y1,1;x,y,1;x3,y3,1]);
       D3 = det([x1,y1,1;x2,y2,1;x,y,1]);
       
       if D0*D1>=0 && D0*D2>=0 && D0*D3>=0
          % Find the interpolation coefficients
          X = [x-x3;y-y3];
          A = [x1-x3,x2-x3; y1-y3,y2-y3];
          C = A\X; a=C(1); b=C(2); d = 1-a-b;

          old = [2*elt(1)-1,2*elt(1)
                 2*elt(2)-1,2*elt(2)
                 2*elt(3)-1,2*elt(3)];
          new = [2*i-1,2*i];
%          old = [3*elt(1)-2,3*elt(1)-1,3*elt(1) % /!\ if the interpolated field is 3D
%                 3*elt(2)-2,3*elt(2)-1,3*elt(2)
%                 3*elt(3)-2,3*elt(3)-1,3*elt(3)
%                 3*elt(4)-2,3*elt(4)-1,3*elt(4)];
%          new = [3*i-2,3*i-1,3*i];
          %a=.25; b=.25; c=.25; Rough debug
          u2(new,:) = a*u1(old(1,:),:) + b*u1(old(2,:),:) + d*u1(old(3,:),:);
          break; % 'Cause I don't want to die
       end
    end
 end
 
end
