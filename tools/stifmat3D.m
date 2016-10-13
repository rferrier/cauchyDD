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
    % Adaptated from a codebase from ?? (someone I don't know)
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
    % Comatrix stuff (inv is faster)
%    Ap1 = 1/det(Am1) * [ det(Am1(2:4,2:4)),     -det(Am1(2:4,[1,3,4])),...
%                               det(Am1(2:4,[1,2,4])),     -det(Am1(2:4,1:3))
%                         -det(Am1([1,3,4],2:4)), det(Am1([1,3,4],[1,3,4])),...
%                               -det(Am1([1,3,4],[1,2,4])), det(Am1([1,3,4],1:3))
%                         det(Am1([1,2,4],2:4)), -det(Am1([1,2,4],[1,3,4])),...
%                               det(Am1([1,2,4],[1,2,4])), -det(Am1([1,2,4],1:3))
%                         -det(Am1(1:3,2:4)),     det(Am1(1:3,[1,3,4])),...
%                                -det(Am1(1:3,[1,2,4])),     det(Am1(1:3,1:3)) ]'

    Ap1 = inv(Am1);

    Al = Ap1(1:3,:); % We don't need the constant
    Onmat = [1,0,0,0,0,0,0,0,0
             0,0,0,0,1,0,0,0,0
             0,0,0,0,0,0,0,0,1
             0,.5,0,.5,0,0,0,0,0
             0,0,.5,0,0,0,.5,0,0
             0,0,0,0,0,.5,0,.5,0];% Builds epsilon
             
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

    if format == 0
        K=V*Be'*A*Be;
    elseif format == 1
        K=A*Be;
    else
        error('unknown format asked')
    end

 else
     error('Order of element not implemented');
 end

end

