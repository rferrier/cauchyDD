function K = stifmat3D (Xloc,order,A,format)
 % Computes the stiffness matrix of the element
 % input : Xloc     : coord of the points ot the element
 %         E        : Young Modulus
 %         nu       : Fish ratio (because I'm funny)
 %         order    : order of the element
 %         A        : Stiffness of the model (\epsilon A = \sigma)
 %         format   : if 0, Ku = f
 %                    if 1, Ku = sigma
 
 % output : K     : elementary stiffness matrix
 
 if order == 1
    x11 = Xloc(1,1); x21 = Xloc(4,1); x31 = Xloc(7,1); x41 = Xloc(10,1);  % nodal coordinates (I bet you didn't guess)
    x12 = Xloc(2,1); x22 = Xloc(5,1); x32 = Xloc(8,1); x42 = Xloc(11,1);
    x13 = Xloc(3,1); x23 = Xloc(6,1); x33 = Xloc(9,1); x43 = Xloc(12,1);
    
    a = [ x21-x11 ; x22-x12 ; x23-x13 ];
    b = [ x31-x11 ; x32-x12 ; x33-x13 ];
    c = [ x41-x11 ; x42-x12 ; x43-x13 ];
    V = 1/6*det( [a,b,c] );  % Volume of the tetraedron
    
    % Invert a 4*4 matrix to find grad(u)=A*epsilon 
    Am1 = [x11, x12, x13, 1
           x21, x22, x23, 1
           x31, x32, x33, 1
           x41, x42, x43, 1];

    Ap1 = inv(Am1); % inv is faster than comatrix stuff

    Al = Ap1(1:3,:); % We don't need the constant
    Onmat = [1,0,0,0,0,0,0,0,0
             0,0,0,0,1,0,0,0,0
             0,0,0,0,0,0,0,0,1
             0,1,0,1,0,0,0,0,0
             0,0,1,0,0,0,1,0,0
             0,0,0,0,0,1,0,1,0];% Builds epsilon
             
    Permu = [1,0,0,0,0,0,0,0,0,0,0,0
             0,0,0,1,0,0,0,0,0,0,0,0
             0,0,0,0,0,0,1,0,0,0,0,0
             0,0,0,0,0,0,0,0,0,1,0,0
             0,1,0,0,0,0,0,0,0,0,0,0
             0,0,0,0,1,0,0,0,0,0,0,0
             0,0,0,0,0,0,0,1,0,0,0,0
             0,0,0,0,0,0,0,0,0,0,1,0
             0,0,1,0,0,0,0,0,0,0,0,0
             0,0,0,0,0,1,0,0,0,0,0,0
             0,0,0,0,0,0,0,0,1,0,0,0
             0,0,0,0,0,0,0,0,0,0,0,1]; % Permutation operator (u1x,u2x... instead of
                                        %                       u1x,u1y...)
             
    Be = Onmat*[Al, zeros(3,8)
                zeros(3,4), Al, zeros(3,4)
                zeros(3,8), Al] * Permu;
                
    % This one is 10/20% slower
%    Xloc1 = [Xloc(1:3:10) , Xloc(2:3:11) , Xloc(3:3:12)];
%
%    DN=[-1 1 0 0 ; -1 0 1 0 ; -1 0 0 1 ]'; % derivative of shape functions
%    J=Xloc1'*DN;                               % jacobian matrix
%    detJ=det(J);                               % jacobian
%    invJ=inv(J);                               % inverse of jacobian matrix
%
%    GN=DN*invJ;                                % gradient of shape functions
%    Be=[GN(1,1) 0 0 GN(2,1) 0 0 GN(3,1) 0 0 GN(4,1) 0 0 ;
%        0 GN(1,2) 0 0 GN(2,2) 0 0 GN(3,2) 0 0 GN(4,2) 0 ;
%        0 0 GN(1,3) 0 0 GN(2,3) 0 0 GN(3,3) 0 0 GN(4,3) ;
%        GN(1,2) GN(1,1) 0 GN(2,2) GN(2,1) 0 GN(3,2) GN(3,1) 0 GN(4,2) GN(4,1) 0 ;
%        GN(1,3) 0 GN(1,1) GN(2,3) 0 GN(2,1) GN(3,3) 0 GN(3,1) GN(4,3) 0 GN(4,1) ;
%        0 GN(1,3) GN(1,2) 0 GN(2,3) GN(2,2) 0 GN(3,3) GN(3,2) 0 GN(4,3) GN(4,2) ] ;

    if format == 0
        K=abs(V)*Be'*A*Be;%Be'*A*Be*detJ/6;%1/6 is the Gauss weight
    elseif format == 1
        K=A*Be;
    else
        error('unknown format asked')
    end

 elseif order == 2
    % Yay, I did it by myself ! (thank you, Aster)
    s5 = sqrt(5);
    coe1 = (5-s5)/20; coe2 = (5+3*s5)/20;

    a_gauss = [ coe1 coe1 coe1
                coe1 coe1 coe2
                coe1 coe2 coe1
                coe2 coe1 coe1 ];             % Gauss abscissae
    w_gauss=[1/24 1/24 1/24 1/24];                 % Gauss weights

    Xloc1 = [Xloc(1:3:28) , Xloc(2:3:29) , Xloc(3:3:30)];  % Transform the data
    
    if format == 0
        K=zeros(30);
    elseif format == 1
        K=zeros(6,30);
    end
    
    for g=1:4                                    % loop over Gauss points
      a=a_gauss(g,:);                            % coordinates for gauss point

      DN=[-3+4*(a(1)+a(2)+a(3)) 4*a(1)-1 0 0 4*(1-2*a(1)-a(2)-a(3))... % derivative of shape functions
          4*a(2) -4*a(2) -4*a(3) 0 4*a(3) ;
          -3+4*(a(1)+a(2)+a(3)) 0 4*a(2)-1 0  -4*a(1)...         % GMSH numerotation
          4*a(1) 4*(1-a(1)-2*a(2)-a(3)) -4*a(3) 4*a(3) 0;
          -3+4*(a(1)+a(2)+a(3)) 0 0 4*a(3)-1 -4*a(1) ...
          0 -4*a(2) 4*(1-a(1)-a(2)-2*a(3)) 4*a(2) 4*a(1) ]';
      J=Xloc1'*DN;                               % jacobian matrix
      detJ=det(J);                               % jacobian
      invJ=inv(J);                               % inverse of jacobian matrix

      GN=DN*invJ;                                % gradient of shape functions
      Be=[GN(1,1) 0 0 GN(2,1) 0 0 GN(3,1) 0 0 GN(4,1) 0 0 GN(5,1) 0 0 ...
          GN(6,1) 0 0 GN(7,1) 0 0 GN(8,1) 0 0 GN(9,1) 0 0 GN(10,1) 0 0 ;
          0 GN(1,2) 0 0 GN(2,2) 0 0 GN(3,2) 0 0 GN(4,2) 0 0 GN(5,2) 0 ...
          0 GN(6,2) 0 0 GN(7,2) 0 0 GN(8,2) 0 0 GN(9,2) 0 0 GN(10,2) 0 ;
          0 0 GN(1,3) 0 0 GN(2,3) 0 0 GN(3,3) 0 0 GN(4,3) 0 0 GN(5,3) ...
          0 0 GN(6,3) 0 0 GN(7,3) 0 0 GN(8,3) 0 0 GN(9,3) 0 0 GN(10,3) ;
          GN(1,2) GN(1,1) 0 GN(2,2) GN(2,1) 0 GN(3,2) GN(3,1) 0 ...
          GN(4,2) GN(4,1) 0 GN(5,2) GN(5,1) 0 GN(6,2) GN(6,1) 0 ...
          GN(7,2) GN(7,1) 0 GN(8,2) GN(8,1) 0 GN(9,2) GN(9,1) 0 ...
          GN(10,2) GN(10,1) 0 ;
          GN(1,3) 0 GN(1,1) GN(2,3) 0 GN(2,1) GN(3,3) 0 GN(3,1) ...
          GN(4,3) 0 GN(4,1) GN(5,3) 0 GN(5,1) GN(6,3) 0 GN(6,1) ...
          GN(7,3) 0 GN(7,1) GN(8,3) 0 GN(8,1) GN(9,3) 0 GN(9,1) ...
          GN(10,3) 0 GN(10,1) ;
          0 GN(1,3) GN(1,2) 0 GN(2,3) GN(2,2) 0 GN(3,3) GN(3,2) ...
          0 GN(4,3) GN(4,2) 0 GN(5,3) GN(5,2) 0 GN(6,3) GN(6,2) ...
          0 GN(7,3) GN(7,2) 0 GN(8,3) GN(8,2) 0 GN(9,3) GN(9,2) ...
          0 GN(10,3) GN(10,2) ;
          ]; % I don't love you, Be
      
       if format == 0
          K=K+Be'*A*Be*detJ*w_gauss(g);     % contribution to stiffness matrix
       elseif format == 1
          K=K+A*Be/4;   
          % Rem : this function returns 1/4 of the sum of sigma at each 
          % Gauss point (there are 4 Gauss points).
       else
          error('unknown format asked')
       end
    end
    
 else
     error('Order of element not implemented');
 end

end