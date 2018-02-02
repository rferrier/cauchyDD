function K = thermat (Xloc,order,A,format)
 % Computes the conductivity matrix of the element
 % input : Xloc     : coord of the points ot the element
 %         order    : order of the element
 %         A        : Stiffness of the model (\epsilon A = \sigma)
 %         format   : if 0, Ku = f
 %                    if 1, Ku = sigma
 
 % output : K     : elementary stiffness matrix
 
 if order == 1
    % Adaptated from a codebase from ?? (someone I don't know)
    x11 = Xloc(1,1); x21 = Xloc(3,1); x31 = Xloc(5,1);  % nodal coordinates
    x12 = Xloc(2,1); x22 = Xloc(4,1); x32 = Xloc(6,1);
    
    S = abs( .5*((x21-x11)*(x32-x12)-(x31-x11)*(x22-x12)) );% element area
    Be=[x22-x32,x32-x12,x12-x22;
        x31-x21,x11-x31,x21-x11]/(2*S);
    
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
