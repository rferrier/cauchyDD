function K = stifmat (Xloc,order,A,format)
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
%    a = [1/3;1/3]; % Coordinates of the gauss point
%    Belem = shapefunc(a,1,1);
%
%    BJ1 = [Belem(1,1),0,Belem(1,2),0,Belem(1,3),0 ;...
%           0,Belem(1,1),0,Belem(1,2),0,Belem(1,3)] ; % derivative/a1
%    BJ2 = [Belem(2,1),0,Belem(2,2),0,Belem(2,3),0 ;...
%           0,Belem(2,1),0,Belem(2,2),0,Belem(2,3)] ; % derivative/a2
%
%    J = [BJ1*Xloc,BJ2*Xloc];
%    Ja = det(J);  % Jacobian
%    Jm1 = inv(J);
%    B1 = Jm1*Belem;
%    
%    B = [B1(1,1),0,B1(1,2),0,B1(1,3),0 ;...
%         0,B1(2,1),0,B1(2,2),0,B1(2,3) ;
%         B1(1,1)/2,B1(2,1)/2,B1(1,3)/2,...
%         B1(2,3)/2,B1(1,2)/2,B1(2,2)/2] ;
%    
%    K = B'*A*B*Ja;
    % TODO : fix because it's KK
    
    % Other method, dirtier but works
    % Adaptated from a codebase from ?? (someone I don't know)
    x11 = Xloc(1,1); x21 = Xloc(3,1); x31 = Xloc(5,1);  % nodal coordinates
    x12 = Xloc(2,1); x22 = Xloc(4,1); x32 = Xloc(6,1);
    
    S = abs( .5*((x21-x11)*(x32-x12)-(x31-x11)*(x22-x12)) );% element area
    Be=[x22-x32,0,x32-x12,0,x12-x22,0;
        0,x31-x21,0,x11-x31,0,x21-x11;
        x31-x21,x22-x32,x11-x31, ...
        x32-x12,x21-x11,x12-x22]/(2*S);
%        (x31-x21)/2,(x22-x32)/2,(x11-x31)/2, ...
%        (x32-x12)/2,(x21-x11)/2,(x12-x22)/2]/(2*S);
    
    if format == 0
        K=S*Be'*A*Be;
    elseif format == 1
        K=A*Be;
    else
        error('unknown format asked')
    end

 elseif order == 2
    % Again adapted from the same code
    
    a_gauss=1/6*[4 1 1; 1 4 1; 1 1 4];           % Gauss abscissae
    w_gauss=[1/6 1/6 1/6];                       % Gauss weights

    Xloc1 = [Xloc(1:2:11) , Xloc(2:2:12)];       % Transform the data
    
    if format == 0
        K=zeros(12);
    elseif format == 1
        K=zeros(3,12);
    end
    
    for g=1:3,                                   % loop over Gauss points
      a=a_gauss(g,:);                            % coordinates for gauss point
      DN=[4*a(1)-1 0 -4*a(3)+1 4*a(2)...         % derivative of shape functions
          -4*a(2) 4*(a(3)-a(1));                 % w.r.t. a_1, a_2 and a3
          0 4*a(2)-1 -4*a(3)+1 4*a(1) ...
          4*(a(3)-a(2)) -4*a(1)]';
      J=Xloc1'*DN;                               % jacobian matrix
      detJ=J(1,1)*J(2,2)-J(1,2)*J(2,1);          % jacobian
      invJ=1/detJ*[ J(2,2) -J(1,2); ...          % inverse jacobian matrix
                   -J(2,1)  J(1,1)];
      GN=DN*invJ;                                % gradient of shape functions
      Be=[GN(1,1) 0 GN(2,1) 0 GN(3,1) 0 ...
          GN(4,1) 0 GN(5,1) 0 GN(6,1) 0;
          0 GN(1,2) 0 GN(2,2) 0 GN(3,2)...
          0 GN(4,2) 0 GN(5,2) 0 GN(6,2);
          GN(1,2) GN(1,1) GN(2,2) GN(2,1) ...
          GN(3,2) GN(3,1) GN(4,2) GN(4,1) ...
          GN(5,2) GN(5,1) GN(6,2) GN(6,1)];
      
       if format == 0
          K=K+Be'*A*Be*detJ*w_gauss(g);        % contribution to stiffness matrix
       elseif format == 1
          K=K+A*Be/3;   
          % Rem : this function returns 1/3 of the sum of sigma at each 
          % Gauss point (there are 3 Gauss points).
       else
          error('unknown format asked')
       end
    end
    
 else
     error('Order of element not implemented');
 end

end
