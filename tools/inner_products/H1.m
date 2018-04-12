function K = H1 (nodes, elem, order)
 % Computes the surfacic H1 norm matrix (in 3D)
 % input : nodes    : list of nodes : [x1,y1;x2,y2;.....]
 %         elem     : list of elements
 %         order    : order of the elements
 
 % output : K       : inner product matrix
 
 K = sparse(3*size(nodes,1), 3*size(nodes,1));
 
 if order > 1
    warning('H1 norm function is not ready for T>3 elements')
 end
 
 for i=1:size(elem,1)
     Xloc1 = nodes(elem(i,:),:);    % Extract and adapt coords
     nnodes = size(Xloc1,1);
     Xloc = zeros(3*nnodes,1);

     map = zeros(1,3*nnodes);
     for j=1:nnodes
         map(1,[3*j-2,3*j-1,3*j]) = [3*elem(i,j)-2, 3*elem(i,j)-1, 3*elem(i,j)];
         Xloc([3*j-2,3*j-1,3*j],1) = [Xloc1(j,1);Xloc1(j,2);Xloc1(j,3)];
     end
     
     x11 = Xloc(1,1); x21 = Xloc(4,1); x31 = Xloc(7,1); % nodal coordinates
     x12 = Xloc(2,1); x22 = Xloc(5,1); x32 = Xloc(8,1);
     x13 = Xloc(3,1); x23 = Xloc(6,1); x33 = Xloc(9,1);
     
     % Build a basis of the element
     t  = [ x21-x11 ; x22-x12 ; x23-x13 ]; t  = t/norm(t);
     v1 = [ x31-x11 ; x32-x12 ; x33-x13 ]; v1 = v1/norm(v1);
     n  = [ t(2)*v1(3) - t(3)*v1(2)
            -t(1)*v1(3) + t(3)*v1(1)
            t(1)*v1(2) - t(2)*v1(1)  ]; 
     n  = n/norm(n); % normal
     v  = [ n(2)*t(3) - n(3)*t(2)
            -n(1)*t(3) + n(3)*t(1)
            n(1)*t(2) - n(2)*t(1)  ];
     v = v/norm(v);
            
     Q = [t,n,v];    % Change basis
     Xlq = [ Q',zeros(3,6) ; zeros(3),Q',zeros(3) ; zeros(3,6),Q' ]*Xloc;
    
     x11 = Xlq(1,1); x21 = Xlq(3,1); x31 = Xlq(5,1);  % nodal coordinates
     x12 = Xlq(2,1); x22 = Xlq(4,1); x32 = Xlq(6,1);
    
     S = abs( .5*((x21-x11)*(x32-x12)-(x31-x11)*(x22-x12)) );% element area
     Be=[x22-x32,0,0,x32-x12,0,0,x12-x22,0,0;
         x31-x21,0,0,x11-x31,0,0,x21-x11,0,0;
         0,x22-x32,0,0,x32-x12,0,0,x12-x22,0;
         0,x31-x21,0,0,x11-x31,0,0,x21-x11,0;
         0,0,x22-x32,0,0,x32-x12,0,0,x12-x22;
         0,0,x31-x21,0,0,x11-x31,0,0,x21-x11  ]/(2*S);
    
     Ke=S*Be'*Be;
     
     % Build K
     K(map,map) = K(map,map) + Ke;
 end
end