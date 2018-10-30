function nodes2 = morph2D( nodes, elements, boundnodes, geometry )
 % This function morphs the given mesh onto the requested geometry
 % Mesh cam be in a 2D or 3D space. If the latter is true, it must be plane
 % @input  : nodes      : nodes of the mesh
 %           elements   : elements of the mesh
 %           boundnodes : boundary nodes of the mesh
 %           geometry   : points giving the polygonal geometry
 % @output : nodes2     : new coordinates of the nodes (connectivity is unchanged)

 nnodes = size(nodes,1);
 nnboun = size(boundnodes,1);
 ng     = size(geometry,1);

 if size(nodes,2) == 2
    nodes = [nodes,zeros(nnodes,1)];
 end

 if size(geometry,2) == 2
    geometry = [geometry,zeros(ng,1)];
 end

 %% Translations
 t1 = nodes(1,:)';    nodes = nodes-t1';
 t2 = geometry(1,:)'; geometry = geometry-t2';

 %% Rotate the mesh in the Oxy plane
 theta = 0; phi = 0; %Initialize the angles
 for i=1:100 % Newton iterations
    vec = [ cos(theta)*sin(phi) ; sin(theta) ; cos(theta)*cos(phi) ];
    res = nodes*vec; % Rem : we use all the nodes for a better condition number

    if norm(res) < 1e-12, break; end

    dft = nodes * [ -sin(theta)*sin(phi) ; cos(theta) ; -sin(theta)*cos(phi) ]; % Derivative / theta
    dfp = nodes * [ cos(theta)*cos(phi) ; 0 ; -cos(theta)*sin(phi) ]; % Derivative / theta
    H   = [dft,dfp]'*[dft,dfp]; % Tangent matrix

    dt    = H\([dft,dfp]'*res); % Actualize
    theta = theta-dt(1);
    phi   = phi-dt(2);
 end
 R = [ cos(phi), 0, -sin(phi) ; ...
       -sin(theta)*sin(phi), cos(theta), -sin(theta)*cos(phi) ; ...
       cos(theta)*sin(phi), sin(theta), cos(theta)*cos(phi) ];
 nodes = transpose(R*nodes'); % Apply the rotation
 %nodes(:,3) = 0; % Clean the last row

%norm(res)
%norm(nodes(:,3))
%theta
%phi

 %% Rotate the geometry
 theta = pi/4; phi = pi/4; %Initialize the angles
 for i=1:100 % Newton iterations
    vec = [ cos(theta)*sin(phi) ; sin(theta) ; cos(theta)*cos(phi) ];
    res = geometry*vec; % Rem : we use all the nodes for a better condition number

    if norm(res) < 1e-12, break; end

    dft = geometry * [ -sin(theta)*sin(phi) ; cos(theta) ; -sin(theta)*cos(phi) ]; % Derivative / theta
    dfp = geometry * [ cos(theta)*cos(phi) ; 0 ; -cos(theta)*sin(phi) ]; % Derivative / theta
    H   = [dft,dfp]'*[dft,dfp]; % Tangent matrix

    dt    = H\([dft,dfp]'*res); % Actualize
    theta = theta-dt(1);
    phi   = phi-dt(2);
 end
% R = [ cos(phi), -sin(phi), 0 ; ...
%       cos(theta)*sin(phi), cos(theta)*cos(phi), -sin(theta) ; ...
%       sin(theta)*sin(phi), sin(theta)*cos(phi), cos(theta) ];
 R = [ cos(phi), 0, -sin(phi) ; ...
       -sin(theta)*sin(phi), cos(theta), -sin(theta)*cos(phi) ; ...
       cos(theta)*sin(phi), sin(theta), cos(theta)*cos(phi) ];
 geometry = transpose(R*geometry'); % Apply the rotation
 geometry(:,3) = 0; % Clean the last row

%norm(geometry(:,3))
%theta
%phi

 %% Compute the curvilign abscissae
 lcurv = zeros(nnboun,1);
 lcurg = zeros(ng,1);

 cum = 0;
 for i=2:nnboun
    xp1 = nodes(boundnodes(i-1),1); yp1 = nodes(boundnodes(i-1),2); zp1 = nodes(boundnodes(i-1),3);
    x1  = nodes(boundnodes(i),1);   y1  = nodes(boundnodes(i),2);   z1  = nodes(boundnodes(i),3);
    cum = cum + sqrt( (x1-xp1)^2 + (y1-yp1)^2 + (z1-zp1)^2 );
    lcurv(i) = cum;
 end
 % And the last point
 xp1 = nodes(boundnodes(end-1),1); yp1 = nodes(boundnodes(end-1),2); zp1 = nodes(boundnodes(end-1),3);
 x1  = nodes(boundnodes(end),1);   y1  = nodes(boundnodes(end),2);   z1  = nodes(boundnodes(end),3);
 cum = cum + sqrt( (x1-xp1)^2 + (y1-yp1)^2 + (z1-zp1)^2 );
 lcurv = lcurv/cum; % Normalize

 cum = 0;
 for i=2:ng
    xp1 = geometry(i-1,1); yp1 = geometry(i-1,2); zp1 = geometry(i-1,3);
    x1  = geometry(i,1);   y1  = geometry(i,2);   z1  = geometry(i,3);
    cum = cum + sqrt( (x1-xp1)^2 + (y1-yp1)^2 + (z1-zp1)^2 );
    lcurg(i) = cum;
 end
 % And the last point
 xp1 = geometry(end-1,1); yp1 = geometry(end-1,2); zp1 = geometry(end-1,3);
 x1  = geometry(end,1);   y1  = geometry(end,2);   z1  = geometry(end,3);
 cum = cum + sqrt( (x1-xp1)^2 + (y1-yp1)^2 + (z1-zp1)^2 );
 lcurg = lcurg/cum;

 %% Attribute any bounday node to a part of the geometry
 nodes1 = zeros(nnboun,3); % New abscissae on the boundary
 % First, choose the points
 isattributed = zeros(nnboun);
 avaliable = 1:nnboun; index = 0;

 for i=1:ng
    lc = lcurg(i);
    % Look for the closest inattributed
    [~,zebest] = min(abs(lcurv(avaliable)-lc));
    lasti = index; index = avaliable(zebest);
    nodes1( index, : ) = geometry(i,:); % place the node
    avaliable(zebest) = []; % Eliminate it
    % Place all the others
    if i > 1
       lcc = ( lcurv(lasti+1:index-1) - lcurg(i-1) ) / ( lcurg(i) - lcurg(i-1) ); % Re-normalize
       nodes1(lasti+1:index-1,:) = lcc*geometry(i,:) + (1-lcc)*geometry(i-1,:);
    end
    if i == ng
       lcc = ( lcurv(index+1:end) - 1 ) / ( lcurg(i) - 1 );
       nodes1(index+1:end,:) = lcc*geometry(i,:) + (1-lcc)*geometry(1,:);
    end
 end
 nodes11 = zeros(nnodes,3);
 nodes11(boundnodes,:) = nodes1; %

 %% Morph
 mat = [0,1,0];
 K = Krig2(nodes,elements,mat,1,[],[]);

 % The matrix for the Dirichlet equations
 C = zeros(2*nnboun,2*nnodes);
 for i=1:nnboun
    C ( 2*i-1, 2*boundnodes(i)-1 ) = 1;
    C ( 2*i,2*boundnodes(i) )     = 1;
 end

 udir = zeros(2*nnodes,1);
 udir(1:2:2*nnodes-1) = nodes11(:,1)-nodes(:,1); udir(2:2:2*nnodes) = nodes11(:,2)-nodes(:,2);
 Kk = [ K, C' ; C, zeros(2*nnboun) ];
 fk = [ zeros(2*nnodes,1) ; C*udir ];

 du = Kk\fk;

 nodes2 = nodes;
 nodes2(:,1) = nodes2(:,1) + du(1:2:2*nnodes-1);
 nodes2(:,2) = nodes2(:,2) + du(2:2:2*nnodes);

 %% Re-rotate
 nodes2 = transpose(R'*nodes2');
 nodes2 = nodes2 + (t2' - nodes2(1,:));

end
