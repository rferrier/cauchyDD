function B = shapefunc( a,gradient,order )
 % Computes the shape function
 % input : a = |a_1 : coordinates of the point in the ref. base
 %             |a_2
 %         gradient : order of gradient needed
 %         order : order of the element
 
 % output : B : matrix of shape functions
 
 if order == 1
     a1 = a(1,1);
     a2 = a(2,1);
     if gradient == 0 % Only give the shape functions
         B = [1-a1-a2,a1,a2];
     elseif gradient == 1 % symetric part of the 1st gradient of shape
                          %     functions (epsilon)
         B = [-1,1,0;-1,0,1];
     else
         error('Second gradient not implemented for shape functions');
     end
 else
     error('Order of element not implemented');
 end
     
end

