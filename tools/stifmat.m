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
    a = [1/3;1/3]; % Coordinates of the gauss point
    Belem = shapefunc(a,1,1);

    BJ1 = [Belem(1,1),0,Belem(1,2),0,Belem(1,3),0 ;...
           0,Belem(1,1),0,Belem(1,2),0,Belem(1,3)] ; % derivative/a1
    BJ2 = [Belem(2,1),0,Belem(2,2),0,Belem(2,3),0 ;...
           0,Belem(2,1),0,Belem(2,2),0,Belem(2,3)] ; % derivative/a2

    J = [BJ1*Xloc,BJ2*Xloc];
    Ja = det(J);  % Jacobian
    Jm1 = inv(J);
    B1 = Jm1*Belem;
    
    B = [B1(1,1),0,B1(1,2),0,B1(1,3),0 ;...
         0,B1(2,1),0,B1(2,2),0,B1(2,3) ;
         B1(1,1)/2,B1(2,1)/2,B1(1,3)/2,...
         B1(2,3)/2,B1(1,2)/2,B1(2,2)/2] ;
    
    K = B'*A*B*Ja;
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
    
    if format == 0
        K=S*Be'*A*Be;
    elseif format == 1
        K=A*Be;
    else
        error('unknown format asked')
    end

 else
     error('Order of element not implemented');
 end

end

